import scipy as SP
import scipy.stats as ST
import os
import cPickle
import statsmodels
from statsmodels.stats.multitest import fdrcorrection
from collections import Counter
from io_tools import read_suppression_values
from common import *

SP.random.seed(42)


def read_complexes(tsq_genes, filename, i_gene, i_cplx, min_alleles=2):
    ifh = file(filename, 'r')
    gene_i = {}
    for i in range(len(tsq_genes)):
        g = tsq_genes[i]
        if g not in gene_i: gene_i[g] = []
        gene_i[g].append(i)
        
    res = {}
    for l in ifh:
        d = l.strip().split("\t")
        if len(d) < max(i_gene, i_cplx) + 1:
            continue
        gene, cplx = d[i_gene], d[i_cplx]
        if gene not in gene_i: continue
        if cplx not in res: res[cplx] = []
        res[cplx].extend(gene_i[gene])
    return {cplx:res[cplx] for cplx in res if len(res[cplx]) >= min_alleles}


def read_complexes_from_table(tsq_genes, filename, cutoff=0.5, min_alleles=2):
    ifh = file(filename, 'r')
    gene_i = {}
    for i in range(len(tsq_genes)):
        g = tsq_genes[i]
        if g not in gene_i: gene_i[g] = []
        gene_i[g].append(i)
        
    complexes = ifh.next().strip().split("\t")[2:]
    res = {cplx:[] for cplx in complexes}
    
    for l in ifh:
        d = l.strip().split("\t")
        gene = d[0]
        if gene not in gene_i: continue
        for i in range(len(complexes)):
            if float(d[i+2]) > cutoff:
                res[complexes[i]].extend(gene_i[gene])
    return {cplx:res[cplx] for cplx in res if len(res[cplx]) >= min_alleles}
    

def _cor(x,y):
    J = SP.isnan(x) | SP.isnan(y)
    return SP.corrcoef(x[~J], y[~J])[0,1]


def _get_random_cors(gene_counts, x, genes, n_random, multigene):
    frac_shared = 1.*sum(map(lambda x:x*(x-1)/2, gene_counts))/(len(x)*(len(x)-1)/2)

    cors = []
    for i in range(int(n_random*frac_shared)): # pairs of alleles from same gene
        g = SP.random.choice(multigene) # genes with at least two alleles
        j1, j2 = SP.random.choice(SP.where(x[genes == g])[0], 2, replace=False) # pick two alleles for a gene
        cors.append(_cor(x[j1], x[j2]))
        
    for i in range(int(n_random*(1.-frac_shared))): # any pair of alleles
        j1, j2 = SP.random.choice(range(len(x)), 2, replace=False)
        cors.append(_cor(x[j1], x[j2]))
    return cors


def _get_random_means(x, n_alleles, n_random):
    means = []
    for i in range(n_random):
        I_rnd = SP.random.choice(range(len(x)), n_alleles, replace=False)
        means.append(SP.nanmedian(x[I_rnd]))
    return SP.array(means)
        

def _get_cors(x):
    c = []
    for i in range(len(x)-1):
        for j in range(i+1, len(x)):
            c.append(_cor(x[i], x[j]))
    return c


def summarize_genesets(genesets, scores, meta, strains, max_cors=150, n_random=1000, do_print=True, do_plot=True, min_plot_m=0.2, min_plot_r=0.3, min_plot_size=4, max_plot_heatmap=41):
    stats = {}
    all_gene_counts = Counter(meta[:, 1])
    multigene = [g for g in all_gene_counts if all_gene_counts[g] > 1]
    
    set_keys = SP.array(genesets.keys())
    J = SP.argsort([len(genesets[k]) for k in set_keys]) # indexes of gene sets sorted by size
    if do_print: print "M(R)\tM(Rrnd)\tp_diff\tM(S)\tM(Srnd)\tp_low\tp_high\tAlleles\tGeneset"

    for geneset in set_keys[J]:
        I, n_alleles = genesets[geneset], len(genesets[geneset])
        if n_alleles == 1: continue # a single allele not a gene set
            
        # 1. random correlations and means [former respecting gene sharing of alleles]
        cors_rnd = _get_random_cors(Counter(meta[I, 1]).values(), scores, meta[:,1], n_random, multigene)
        means_rnd = _get_random_means(scores, n_alleles, n_random*2)
            
        # 2. real correlations and mean
        cors = [0] 
        if n_alleles < max_cors: # skip some _huge_ categories, like "cytoplasmic"
            cors = _get_cors(scores[I])
            
        # 3. calculate summary stats
        cm, cm_rnd, m, m_rnd = SP.nanmean(cors), SP.nanmean(cors_rnd), SP.nanmedian(scores[I]), SP.nanmedian(means_rnd)
        if SP.isnan(cm): cm = cm_rnd
        pc_diff = ST.ttest_ind_from_stats(cm, ST.nanstd(cors), len(cors), cm_rnd, ST.nanstd(cors_rnd), n_random, False)[1]
        pm_low, pm_high = SP.mean(m > means_rnd), SP.mean(m < means_rnd)
        stats[geneset] = [cm, cm_rnd, pc_diff, m, m_rnd, pm_low, pm_high, n_alleles]
        if do_print: print "%.2f\t%.2f\t%.0e\t%.2f\t%.2f\t%.0e\t%.0e\t%d"%tuple(stats[geneset]), geneset

        if not do_plot:
            continue
            
        # 4. plot 
        if (n_alleles >= min_plot_size) and ((cm-cm_rnd >= min_plot_r) or (abs(m-m_rnd) >= min_plot_m)):
            PL.figure(None, [9,3])
            PL.subplot(131)
            J = [j for j in I]
            if len(J) > max_plot_heatmap:
                J = SP.random.choice(I, max_plot_heatmap, replace=False)
            PL.imshow(scores[J], interpolation="none", vmin=-1, vmax=1)
            PL.yticks(range(len(J)), meta[J,0])
            PL.xticks(range(len(strains)), strains, rotation=90)
            PL.title(geneset)
            
            PL.subplot(132)            
            if len(cors) > 1:
                PL.hist(cors, normed=True, range=(-1,1), bins=20, alpha=0.7, color='b')
            PL.axvline(cm, linestyle='dashed', color='b', linewidth=1)
            if len(cors_rnd) > 0:
                PL.hist(cors_rnd, normed=True, range=(-1,1), bins=20, alpha=0.1, color='k')
            PL.axvline(cm_rnd, linestyle='dashed', color='k', linewidth=1)
            PL.title("%s\nN=%d G=%d p=%.0e\nm(R)=%.2f m(R0)=%.2f"%(geneset, len(I), len(SP.unique(meta[I,1])), pc_diff, cm, cm_rnd))
            PL.xlabel("Pearson's R")
            
            PL.subplot(133)
            PL.hist(means_rnd, normed=True, range=(-1,1), bins=20, alpha=0.1, color='k')
            PL.axvline(m, linestyle='dashed', color='b', linewidth=1)
            PL.axvline(m_rnd, linestyle='dashed', color='k', linewidth=1)
            PL.xlabel("Mean suppression")
            PL.title("p_low=%.0e p_high=%.0e"%(pm_low, pm_high))
    return stats            

    
def reduce_to_gene(scores, meta, datasets):
    orfs = SP.unique(meta[:,1])
    tsqs_used = {}
    dataset_genes = {}
    gene_score, gene_meta = SP.zeros([len(orfs), scores.shape[1]]), []
    for o, orf in enumerate(orfs):
        I = SP.where(meta[:,1] == orf)[0]
        meds = SP.nanmedian(scores[I],axis=1)
        J = SP.argsort(meds)
        gene_score[o] = scores[I[J[-1]]]
        gene_meta.append(meta[I[J[-1]]])
        tsqs_used[I[J[-1]]] = o
    for k in datasets:
        dataset_genes[k] = {}
        for geneset in datasets[k]:
            dataset_genes[k][geneset] = [tsqs_used[i] for i in datasets[k][geneset] if i in tsqs_used]
    return gene_score, SP.array(gene_meta), dataset_genes
    

def print_paper_complex_enrichments(fdr=0.2, effect=0.2, max_sd_lim=1, max_sd=0.5, do_print=False, rerun=False):
    scores, meta, strains = read_suppression_values()
    I = ((meta[:,6] == "0") & (meta[:,8] == "0") & (~(SP.isnan(scores[:,:,0,0]).all(axis=0))))# TS, not translocated
    strains = list(strains)
    scores[strains.index("Y14278"), (meta[:,7] == "1")] = SP.nan # blank out the duplicated chrII for one strain
    for i in range(len(scores)):
        for j in range(len(scores[i])):
            if scores[i,j,0,1] < max_sd_lim and scores[i,j,1,1] > max_sd: scores[i,j,:,1] = SP.nan
                
    scores, meta = scores[:,I,0,1].T, meta[I]
    tsqs = list(meta[:,0])
    tsq_genes = list(meta[:,1])
    meta = SP.array([meta[:,2], meta[:,1], meta[:,0]]).T

    if (not rerun) and os.path.exists("precomputed_complex_stats.pickle"):
        (datasets, all_stats, dataset_genes, gene_stats) = cPickle.load(file("precomputed_complex_stats.pickle", 'rb'))
    else:
        datasets = {'EBI': read_complexes(tsq_genes, "%s/paper/meta/180518_Yeast_complexes_EBI.tab"%DIR_DATA, i_gene=2, i_cplx=1),
        'GO slim': read_complexes(tsq_genes, "%s/paper/meta/go_slim_mapping.tab"%DIR_DATA, i_gene=0, i_cplx=4),
        'KEGG': read_complexes(tsq_genes, "%s/paper/meta/190401_Yeast_KEGGpathways.tab"%DIR_DATA, i_gene=0, i_cplx=3),
        'FunCats': read_complexes(tsq_genes, "%s/paper/meta/160408_Yeast_functional_categories.tab"%DIR_DATA, i_gene=0, i_cplx=2),
        'CoLoc': read_complexes_from_table(tsq_genes, "%s/paper/meta/191008_Yeast_localization_CYCLoPs_WT1_LOCscore.tab"%DIR_DATA)}

        all_stats = {k: summarize_genesets(datasets[k], scores, meta, strains, max_cors=150, n_random=1000, do_print=do_print, do_plot=False) for k in datasets}
        gene_scores, gene_meta, dataset_genes = reduce_to_gene(scores, meta, datasets)
        gene_stats = {k: summarize_genesets(dataset_genes[k], gene_scores, gene_meta, strains, max_cors=150, n_random=1000, do_print=do_print, do_plot=False) for k in datasets}
        cPickle.dump((datasets,all_stats, dataset_genes, gene_stats), file('precomputed_complex_stats.pickle', 'wb'), -1)

    for k in ['EBI','CoLoc', 'KEGG', 'FunCats', 'GO slim']:
        mu = SP.array(gene_stats[k].values())[:,0] # c, crnd, m, mrnd, n  # correlation of allele values across set
        stats = SP.array(gene_stats[k].values())[:,2] # c, crnd, m, mrnd, n  # correlation of allele values across set
        I = ~SP.isnan(stats)
        pcorr = statsmodels.stats.multitest.fdrcorrection(stats[I])[1]
        N = 1.*len(pcorr)
        J = (pcorr < fdr) & (abs(mu[I]) > effect)
        print "%s: correlation of values across gene set\t%.2f (%d of %d at FDR=%.2f)"%(k, sum(J)/N, sum(J), N, fdr)


def _old_print_paper_complex_enrichments_old(fdr=0.2):
    scores, meta, strains = read_suppression_values()
    I = ((meta[:,3] == "0") & (meta[:,5] == "0") & (~(SP.isnan(scores[:,:,0]).all(axis=1))))# TS, not translocated
    scores[(meta[:,6] == "1"), strains.index("Y14278")] = SP.nan # blank out the duplicated chrII for one strain
    
    scores, meta, strains = scores[I,0:10,0], meta[I], strains[0:10]
    tsqs = list(meta[:,2])
    all_gene_counts = Counter(meta[:, 1])
    multigene = [g for g in all_gene_counts if all_gene_counts[g] > 1]
        
    if os.path.exists("complex_stats.pickle"):
        (datasets, all_stats, dataset_genes, gene_stats) = cPickle.load(file("complex_stats.pickle", 'rb'))
    else:
        datasets = {'EBI': read_complexes(tsq_genes, "%s/paper/meta/180518_Yeast_complexes_EBI.tab"%DATA_DIR, i_gene=2, i_cplx=1),
        'GO slim': read_complexes(tsq_genes, "%s/paper/meta/go_slim_mapping.tab"%DATA_DIR, i_gene=0, i_cplx=4),
        'KEGG': read_complexes(tsq_genes, "%s/paper/meta/190401_Yeast_KEGGpathways.tab"%DATA_DIR, i_gene=0, i_cplx=3),
        'FunCats': read_complexes(tsq_genes, "%s/paper/meta/160408_Yeast_functional_categories.tab"%DATA_DIR, i_gene=0, i_cplx=2),
        'CoLoc': read_complexes_from_table(tsq_genes, "%s/paper/meta/191008_Yeast_localization_CYCLoPs_WT1_LOCscore.tab"%DATA_DIR)}

        all_stats = {k: summarize_genesets(datasets[k], scores, meta, strains, max_cors=150, n_random=1000, do_print=False, do_plot=False) for k in datasets}
        gene_scores, gene_meta, dataset_genes = reduce_to_gene(scores, meta, datasets)
        gene_stats = {k: summarize_genesets(dataset_genes[k], gene_scores, gene_meta, strains, max_cors=150, n_random=1000, do_print=False, do_plot=False) for k in datasets}

    for k in ['EBI','CoLoc', 'KEGG', 'FunCats', 'GO slim']:
        stats = SP.array(gene_stats[k].values())[:,2] # c, crnd, m, mrnd, n  # correlation of allele values across set
        I = ~SP.isnan(stats)
        pcorr = statsmodels.stats.multitest.fdrcorrection(stats[I])[1]
        N = 1.*len(pcorr)
        J = pcorr < fdr
        print "%s: correlation of values across gene set\t%.2f (%d of %d at FDR=%.2f)"%(k, sum(J)/N, sum(J), N, fdr)
