import scipy as SP
import scipy.stats as ST
import pylab as PL
from common import *
from io_tools import read_suppression_values, read_fitness_values, read_followup_phenotypes, read_wild_strain_name

def print_strain_stats(nan_sick=False):
    total_strains_colonysize = len(file("%s/paper/tables/combined_colonysizes.tab"%DIR_DATA, 'r').readlines()) - 1
    total_strains_fitness = len(file("%s/paper/tables/combined_fitnesses.tab"%DIR_DATA, 'r').readlines()) - 1
    total_strains_suppression = len(file("%s/paper/tables/combined_suppression.tab"%DIR_DATA, 'r').readlines()) - 1
    print "Strains: %d with [some] fitnesses: %d with [some] suppression: %d"%(total_strains_colonysize, total_strains_fitness, total_strains_suppression)
    
    supp, meta_all, strains = read_suppression_values(keep_translocated=True)
    meta = SP.array(meta_all[:,6:], bool)
    print "Strains with suppression values: %d for %d genes (%d bad or irreproducible)."%(total_strains_suppression, len(set(meta_all[:,1])), total_strains_colonysize-total_strains_suppression)
    print "chrXVI translocation: %d, chrVIII: %d, non-TS: %d"%(sum(meta[:,0]), sum(meta[:,1]), sum(meta[:,2]))
    I = ~meta[:,2]
    print "TS strains: %d for %d genes"%(sum(I), len(set(meta_all[I,1])))
    I = ~(meta[:,2] | meta[:,0])
    print "TS strains, no chrXVI: %d for %d genes"%(sum(I), len(set(meta_all[I,1])))
    

def _get_i(orf):
    n = 0
    i = 3
    while i < len(orf) and orf[i] in '0123456789':
        n = 10*n + int(orf[i])
        i += 1
    return n


def plot_allele_supp_genomic_xvi(cutoff=0.75, print_gene_means=False):
    supp, meta, strains = read_suppression_values(keep_translocated=True)

    # 1. set up observation matrix
    Ic = [o[1] == "P" for o in meta[:,1]] # chr XVI
    Its = meta[:,8] == "0"
    gene_i = []
    for orf in meta[Ic&Its,1]:
        if orf[2] == "R": gene_i.append(_get_i(orf) + 2000)
        else: gene_i.append(1000 - _get_i(orf))
    Isort = SP.argsort(gene_i)
    Itrans = (meta[:,6] == "1")[Ic&Its][Isort]
    y = supp[:,Ic&Its,0,1].T[Isort]

    # 2. plot
    PL.figure(None, [3,1.5])
    pi = 1
    N = len(y)
    x = SP.arange(N)
    titles = "No translocation: %d of %d", "With translocation %d of %d"
    colors = 'br'
    ymi, yma = -0.75, 2.5
    for I in [[0,1,2,3,5,7,8], [4,6,9]]:
        a, b = x, SP.nanmean(y[:,I],axis=1)
        #print I, SP.array(10*SP.nanmean(y[Itrans][:,I], axis=1),int)
        if print_gene_means:
            n1, n2 = 0,0
            m = meta[Ic&Its][Isort][Itrans]
            for orf in sorted(set(m[:,1])):
                Iorf = m[:,1] == orf
                n1 += sum(~SP.isnan(b[Itrans][Iorf]))
                n2 += sum(b[Itrans][Iorf] < 0.5)
                print "%s: %d m=%.2f low=%d"%(orf, sum(~SP.isnan(b[Itrans][Iorf])), SP.nanmean(b[Itrans][Iorf]), sum(b[Itrans][Iorf] < 0.5))
            print n1, n2
        PL.plot(a, b, colors[4 in I]+".", markersize=6, alpha=0.7)
        i1, i2 = min(SP.where(Itrans)[0]), max(SP.where(Itrans)[0])
        PL.fill([i1,i2,i2,i1],[ymi,ymi,yma,yma], color='k', alpha=0.05)
        PL.xlim(0,N); PL.ylim(ymi, yma)
        PL.yticks([0,1,2])
        PL.ylabel("Average\nsuppression")
        PL.xticks([])
        if pi == 3: PL.xlabel("TS alleles for queries along chrXVI")
        PL.plot([0,N], [0,0], 'k--', alpha=1, linewidth=1)
        PL.plot([0,N], [cutoff, cutoff], 'k-', alpha=1, linewidth=1)
        n_supp, total = sum(SP.nanmean(y[:,I],axis=1)[Itrans] > cutoff), sum(Itrans) - sum(SP.isnan(b))
        n_lenient = sum(SP.nanmean(y[:,I],axis=1)[Itrans] > cutoff-0.25)
        print titles[4 in I]%(n_supp, total), "of measured alleles suppressed on chrXVI across strains on average, %d leniently"%(n_lenient)
        PL.title(titles[pi-2]%(n_supp, total))
    PL.tight_layout()
    PL.savefig("%s/paper/plots/2B_chr_xvi_translocation.svg"%DIR_DATA)
        

def print_paper_max_suppression_profile_corrs():
    supp, meta, strains = read_suppression_values(keep_translocated=False)
    I = meta[:,8] == "0" # TS strains only
    
    print "Max and Average correlation to others:"                  
    for s in [5,1,3,7,8,0,2,4,6,9]:
        cors = []
        for i in range(len(strains)):
            if i == s: continue
            x,y = supp[i,I,0,1], supp[s,I,0,1]            
            J = ~(SP.isnan(x) | SP.isnan(y))
            cors.append(SP.corrcoef(x[J],y[J])[0,1])
        print "\t%s: %.3f\t%.3f"%(strains[s], SP.nanmax(cors), SP.mean(cors))


def plot_trees():
    strains = SP.array(["Y%d"%i for i in range(14273, 14283)])
    PL.figure(None, [4,3.5])

    # 1. Create SNP distance matrix and dendrogram
    d = _calc_dists(strains)
    dendrog = dendrogram(average(d), labels=strains, no_plot=True)
    # Reorder a bit
    xsg, ysg = [[a for a in x] for x in dendrog['icoord']], [[b for b in y] for y in dendrog['dcoord']]
    labelxg, labeltxtg = SP.arange(5,100,10), strains[dendrog['leaves']]
    xs, ys, Ireorder, flipped, targeted = _replot_order(xsg, ysg, labelxg, [4,5], verbose=False)
    # Plot
    PL.subplot(211)
    _plot_dendro(xs, ys, labelxg, labeltxtg[Ireorder], flipped, targeted, colors='ggg')
    PL.yticks([0,25000,50000,75000], ['0k','25k','50k','75k'])
    PL.ylabel("Number of sites"); PL.title("Genotype tree")

    # 2. Create suppression distance matrix and dendrogram
    supp, meta, strains = read_suppression_values(keep_translocated=False) # supp: S x K x M/V x 26/34
    I = SP.ones(supp.shape[1], bool)
    I[(meta[:,1] == "CDC28") | (meta[:,8] == "1") | (meta[:,6] == "1")] = False # filter out CDC28 and non-ts
    J = SP.array([2,0,6,9,8,4,7,3,5,1])
    pd = _pdistnan(supp[:,I][J,:,0,1])
    dendrop = dendrogram(average(pd), labels=strains[J], no_plot=True)
    # Reorder a bit
    xsp, ysp = [[a for a in x] for x in dendrop['icoord']], [[b for b in y] for y in dendrop['dcoord']]
    labelxp, labeltxtp = SP.arange(5,100,10), strains[J[dendrop['leaves']]]
    x, y, Ireorder, flipped, targeted = _replot_order(xsp, ysp, labelxp, [6,1], verbose=False) # was [2,5]
    # Plot
    PL.subplot(212)
    _plot_dendro(x, 1.-SP.array(y), labelxp, labeltxtp[Ireorder], flipped, targeted, colors='bbb')
    PL.yticks([0.4,0.6,0.8,1.0],["0.6","0.4","0.2","0"])
    PL.ylabel("1-R(suppression)"); PL.title("Phenotype tree")
    PL.tight_layout()

    PL.savefig("%s/paper/plots/2C_geno_pheno_tree.svg"%DIR_DATA)



import scipy.cluster.hierarchy as hierarchy
from scipy.cluster.hierarchy import average, dendrogram
from scipy.spatial.distance import pdist, squareform

def _read_sites(strain, do_rm=False):
    import os
    sitefile = "sites_%s.tab"%strain
    if not os.path.exists(sitefile):
        cmd = "gunzip -c %s/seq/afs/af_calls/%s_all_samples.joint_varcall.freebayes.annotated.SNP.afs.vcf.gz | cut -f1,2 | grep -v \\# > %s"%(DIR_DATA, strain, sitefile)
        os.system(cmd)

    ifh = file(sitefile, 'r')
    sites = {l: True for l in ifh if l[0:3] == "chr"}
    if do_rm:
        os.system("rm %s"%sitefile)
    return sites


def _calc_dists(strains):
    sites = {s: set(_read_sites(s)) for s in strains}
    d = SP.zeros([len(sites), len(sites)])
    for i in range(len(sites) - 1):
        for j in range(i+1, len(sites)):
            #d[i,j] = d[j,i] = 1. - 1.*len(sites[strains[i]] & sites[strains[j]])/len(sites[strains[i]] | sites[strains[j]])
            d[i,j] = d[j,i] = len(sites[strains[i]] | sites[strains[j]]) - len(sites[strains[i]] & sites[strains[j]])
    return squareform(d)


EPS = 1e-10

def _flip(to_flip, xs, ys, i_parent, i_xrange, verbose=False):
    xmin, xmax = i_xrange[to_flip]
    center = 0.5*(xmin+xmax) # center for flipping
    if verbose: print xmin, xmax, center
    flipped = []

    # 1. Downwards: children's coordinates are flipped around the center
    for i, x in enumerate(xs):
        if (x[0] >= xmin) and (x[-1] <= xmax): # in range
            if verbose:
                print "flipping", x, "to",
            for j in range(len(x)):
                x[j] = xmin + (xmax - x[j]) # flip around center
            if verbose: print x
            flipped.append(i)

    # 2. Upwards: update parent coordinates
    new_root = 0.5*(xs[to_flip][1] + xs[to_flip][2]) # current root with updated coordinates
    old_root = xmin + (xmax - new_root) # pre-flipping
    parent = i_parent[to_flip]
    while parent is not None:
        if verbose: print parent, old_root, new_root
        parent_root = 0.5*(xs[parent][1] + xs[parent][2])
        if old_root > parent_root: # to the right
            xs[parent][2] = xs[parent][3] = new_root
        else:
            xs[parent][0] = xs[parent][1] = new_root

        old_root = parent_root
        new_root = 0.5*(xs[parent][1] + xs[parent][2])
        parent = i_parent[parent]

    return xs, ys, flipped


def _construct_dag(xs, ys):
    N = len(xs)
    i_parent, i_xrange = [None]*N, [[1e10,-1e10] for i in range(N)]

    for i in range(N):
        for j in range(N):
            if (abs(ys[j][0] - ys[i][1]) < EPS) or (abs(ys[j][3] - ys[i][1]) < EPS):
                i_parent[i] = j

    for i in range(N): # for each node, update all its parents to include it in its xrange
        cur_xrange = (xs[i][0], xs[i][3])
        j = i
        while j is not None:
            i_xrange[j][0] = min(i_xrange[j][0], cur_xrange[0])
            i_xrange[j][1] = max(i_xrange[j][1], cur_xrange[1])
            cur_xrange = i_xrange[j]
            j = i_parent[j]

    return i_parent, i_xrange


def _update_order(I, leaf_x, x_span, verbose=False):
    i1 = SP.where(leaf_x >= x_span[0])[0][0]
    i2 = None
    if x_span[1] >= max(leaf_x): i2 = len(leaf_x)
    else:  i2 = SP.where(leaf_x > x_span[1])[0][0]
    if verbose: print leaf_x, x_span, i1, i2, I,
    Inew = I[0:i1] + I[i1:i2][::-1] + I[i2:]
    if verbose: print Inew
    return Inew


def _replot_order(xs, ys, leaf_x, flip_list, verbose=False):
    # re-order leaves, flipping all joins ordered from bottom
    y_order = SP.argsort([y[1] for y in ys])
    i_parent, i_xrange = _construct_dag(xs, ys)
    Iorder = range(len(leaf_x))

    if verbose:
        print y_order
        for i in y_order: print "\t", i, y_order[i], ys[y_order[i]][1:3]
        print i_parent
        for i in range(len(xs)):
            j = i_parent[i]
            if j is not None:
                print "\t", i, int(0.5*xs[i][1]+0.5*xs[i][2]), ys[i][1], "=>", j, int(0.5*xs[j][1]+0.5*xs[j][2]), ys[j][1]
        print i_xrange

    all_flipped = set([])
    for i in y_order: # go through tree bottom to top
        if i in flip_list: # if want to flip
            xs, ys, flipped = _flip(i, xs, ys, i_parent, i_xrange, verbose)
            Iorder = _update_order(Iorder, leaf_x, i_xrange[i], verbose)
            all_flipped = all_flipped | set(flipped)

    return xs, ys, Iorder, list(all_flipped), flip_list


def _plot_dendro(xs, ys, labelx, labeltxt, flipped=None, targeted=None, colors='gbr'):
    for i in range(len(xs)):
        x,y = xs[i], ys[i]
        color = colors[0]
        if (flipped is not None) and (i in flipped): color = colors[1]
        if (targeted is not None) and (i in targeted): color = colors[2]
        for j in range(3):
            PL.plot(x[j:j+2], y[j:j+2], color+"-")
    PL.xticks(labelx, labeltxt, rotation=90)
    

def _pdistnan(d):
    N = d.shape[0]
    pds = []
    for i in range(N-1):
        for j in range(i+1, N):
            x,y = d[i], d[j]
            In = ~(SP.isnan(x) | SP.isnan(y))
            x, y = x[In], y[In]
            pds.append(1-SP.corrcoef(x,y)[0,1])
    return pds


def _read_replicate_fitnesses(strain, temp=34):
    ifh = file("%s/paper/tables/fitness_%s.tab"%(DIR_DATA, strain), 'r')
    header = ifh.next().strip().split("\t")
    i1, i2 = header.index("R1_%d_mean"%temp), header.index("R2_%d_mean"%temp)
    rep_means = []
    for l in ifh:
        d = l.strip().replace("*","").split("\t")
        rep_means.append([d[i1],d[i2]])
    x,y = SP.array(rep_means, float).T
    I = SP.isnan(x) | SP.isnan(y) | SP.isinf(x) | SP.isinf(y)
    return x, y, SP.corrcoef(x[~I],y[~I])[0,1]
    


def _cor(x,y):
    I = SP.isnan(x) | SP.isnan(y) | SP.isinf(x) | SP.isinf(y)
    return SP.corrcoef(x[~I], y[~I])[0,1]


def plot_fitnesses():
    data, meta, strains = read_fitness_values()
    strain_names = read_wild_strain_name()
    J = meta[:,8] == "1" # not temperature sensitive strains
    S = len(strains)
    P = 5
    PL.figure(None, [2.4*5, 2.0*P])
    pi = 1
    
    for s,strain in enumerate(strains):
        if s == 5:
            PL.tight_layout()
            PL.savefig("%s/paper/plots/FigS1_FitnessOverview1_210305.png"%DIR_DATA, dpi=300)
            PL.savefig("%s/paper/plots/FigS1_FitnessOverview1_210305.svg"%DIR_DATA)
            PL.savefig("%s/paper/plots/FigS1_FitnessOverview1_210305.pdf"%DIR_DATA, dpi=300)
            
            P = 6
            PL.figure(None, [2.4*5, 2.0*P])
            pi = 1
        x0,x1 = data[s,:,0,0], data[s,:,0,1] # means at two temperatures
        for temp in [34]:#[26,34]:
            PL.subplot(P,5,pi); pi += 1
            a,b,r = _read_replicate_fitnesses(strain, temp)
            PL.plot(a, b, ".", alpha=0.2)
            PL.plot([7.5,11.5],[7.5,11.5], 'k--', linewidth=1)
            PL.xticks(SP.arange(8,12,1)); PL.yticks(SP.arange(8,12,1))
            PL.xlabel("%dC: Replicate 1"%temp); PL.ylabel("%dC: Replicate 2"%temp)
            PL.title("R=%.2f"%r)
            PL.xlim(7.5, 11.5); PL.ylim(7.5, 11.5)

            a, b = data[s,:,0,int(temp==34)], data[10,:,0,int(temp==34)]
            r = _cor(a,b)
            PL.subplot(P,5,pi); pi += 1
            PL.plot(a, b, "m.", alpha=0.2)
            PL.plot([7.5,11.5],[7.5,11.5], 'k--', linewidth=1)
            PL.xticks(SP.arange(8,12,1)); PL.yticks(SP.arange(8,12,1))
            PL.xlabel("%dC: %s"%(temp, strain_names[strain])); PL.ylabel("%dC: Control"%temp)
            PL.title("R=%.2f"%r)
            PL.xlim(7.5, 11.5); PL.ylim(7.5, 11.5)
        PL.subplot(P,5,pi); pi += 1
        PL.hist(x0, range=(7.5,11.5), bins=20, alpha=0.2) # Strains x TSQ x 26/34 x M/Sd
        PL.hist(x1, range=(7.5,11.5), bins=20, alpha=0.2)
        PL.xticks(SP.arange(8,12,1))
        PL.ylabel("# alleles"); PL.legend(["26 C", "34 C"], loc="upper left"); PL.xlabel("Growth = log2(colony size)")
        PL.subplot(P,5,pi); pi += 1
        PL.plot(x0[~J], x1[~J], ".", markersize=6, alpha=0.2)
        PL.plot(x0[J], x1[J], "k.", markersize=6, alpha=0.2)
        PL.legend(['TS','Not TS'], numpoints=1, loc='upper left')
        PL.plot([6,12],[6,12], 'k--', alpha=0.4); PL.xlim(7.5,11.5); PL.ylim(7.5,11.5); PL.xticks(SP.arange(8,12,1)); PL.yticks(SP.arange(8,12,1))
        PL.xlabel("Growth at 26 C"); PL.ylabel("Growth at 34 C"); PL.title(strain_names[strain])
        PL.subplot(P,5,pi); pi += 1
        PL.hist(x0-x1, range=(-1,3), bins=16, alpha=0.5); PL.xticks(SP.arange(-1,4,1))
        PL.axvline(0, linestyle='dashed', color='k', linewidth=1)
        PL.xlabel("Growth defect (26 C - 34 C)"); #PL.title(strain)
    PL.tight_layout()
    PL.savefig("%s/paper/plots/FigS1_FitnessOverview2_210305.png"%DIR_DATA, dpi=300)
    PL.savefig("%s/paper/plots/FigS1_FitnessOverview2_210305.svg"%DIR_DATA)
    PL.savefig("%s/paper/plots/FigS1_FitnessOverview2_210305.pdf"%DIR_DATA, dpi=300)


def print_allele_summaries(s_cutoff=0.75, z_cutoff=4.5, verbose=False):
    supp, meta, strains = read_suppression_values(keep_translocated=False)
    Iok = (meta[:,8] == "0") # TS
    Inots = (meta[:,8] == "1") # Not TS
    supp_nots = supp[:,Inots]
    supp, meta = supp[:,Iok], meta[Iok]
    
    hits, tested = {}, {}
    hits['s_34'] = supp[:,:,0,1] > s_cutoff
    tested['s_34'] = ~SP.isnan(supp[:,:,0,1])
    hits['s_26'] = supp[:,:,0,0] > s_cutoff
    tested['s_26'] = ~SP.isnan(supp[:,:,0,0])
    hits['z_34'] = (supp[:,:,0,1] > s_cutoff) & (supp[:,:,0,1]/supp[:,:,1,1] > z_cutoff)
    tested['z_34'] = ~SP.isnan(supp[:,:,0,1]/supp[:,:,1,1])
    hits['z_26'] = (supp[:,:,0,0] > s_cutoff) & (supp[:,:,0,0]/supp[:,:,1,0] > z_cutoff)
    tested['z_26'] = ~SP.isnan(supp[:,:,0,0]/supp[:,:,1,0])
    hits['0_34'] = supp_nots[:,:,0,1] > s_cutoff
    hits['1_34'] = (supp_nots[:,:,0,1] > s_cutoff) & (supp_nots[:,:,0,1]/supp_nots[:,:,1,1] > z_cutoff)
    tested['0_34'] = ~SP.isnan(supp_nots[:,:,0,1])
    tested['1_34'] = ~SP.isnan(supp_nots[:,:,0,1]/supp_nots[:,:,1,1])
    
    if verbose:
        print "s > %.2f z > %.1f"%(s_cutoff, z_cutoff)
        for k in sorted(hits):
            n_hits, n_tested = SP.sum(hits[k]), SP.sum(tested[k])
            n_strain, n_allele = SP.sum(hits[k], axis=1), SP.sum(hits[k], axis=0)
            print "%s\t%d/%d  \t%d\t%d/%d"%(k, n_hits, n_tested, SP.median(n_strain), SP.sum(n_allele>0), SP.sum(~SP.isnan(n_allele)))
    else:
        k = 'z_34'
        p = ST.norm.cdf(-z_cutoff)
        n_hits, n_tested = SP.sum(hits[k]), SP.sum(tested[k])
        n_strain, n_allele = SP.sum(hits[k], axis=1), SP.sum(hits[k], axis=0)
        print "Suppression at s > %.2f, z > %.1f (Bonferroni p after correcting for %d tests=%.2e; nominal %.2e)"%(s_cutoff, z_cutoff, n_tested, p*n_tested, p)
        print "Alleles:"
        print "\tAll tests: %d/%d (%.3f)"%(n_hits, n_tested, 1.*n_hits/n_tested)
        print "\tAny strain: %d/%d (%.2f)"%(SP.sum(n_allele>0), SP.sum(~SP.isnan(n_allele)), 1.*SP.sum(n_allele>0)/SP.sum(~SP.isnan(n_allele)))
        print "\tMedian per strain: %d (%.2f)"%(SP.median(n_strain), 1.*SP.median(n_strain)/SP.sum(~SP.isnan(n_allele)))
        genehits, genetested = _reduce_hits_to_gene(hits[k], tested[k], meta)
        n_hits, n_tested = SP.sum(genehits), SP.sum(genetested)
        n_strain, n_allele = SP.sum(genehits, axis=1), SP.sum(genehits, axis=0)
        print "Genes:"
        print "\tAll tests: %d/%d (%.3f)"%(n_hits, n_tested, 1.*n_hits/n_tested)
        print "\tAny strain: %d/%d (%.2f)"%(SP.sum(n_allele>0), SP.sum(~SP.isnan(n_allele)), 1.*SP.sum(n_allele>0)/SP.sum(~SP.isnan(n_allele)))
        print "\tMedian per strain: %d"%(SP.median(n_strain))


def _reduce_hits_to_gene(hits, tested, meta):
    orfs = meta[:,2]
    genehits = SP.zeros([hits.shape[0], len(SP.unique(orfs))], bool)
    genetested = SP.zeros([hits.shape[0], len(SP.unique(orfs))], bool)
    for o,orf in enumerate(SP.unique(orfs)):
        Iorf = (orfs == orf)
        genehits[:,o] = (hits[:,Iorf].sum(axis=1) > 0)
        genetested[:,o] = (tested[:,Iorf].sum(axis=1) > 0)
    return genehits, genetested

        
def print_paper_s288c_only(s_cutoff=0.75, z_cutoff=4.5, s288c_cutoff=0.5, min_fraction=0.51, do_plot_count=False, do_plot_scatter=False, do_plot_hist=True):
    supp, meta, strains = read_suppression_values(keep_translocated=False)
    I = (meta[:,8] == "0") # TS
    supp, meta = supp[:,I], meta[I] # apply filter
    s = supp[:,:,0,1]
    z = s/supp[:,:,1,1]
    if do_plot_count:
        PL.figure(None, [3,1])
        xs = SP.arange(0.3,2,0.01)
        PL.plot(xs, ([(s > x).any(axis=0).sum() for x in xs]))
        PL.plot(xs, ([(s*(z>z_cutoff) > x).any(axis=0).sum() for x in xs]))
        pval = ST.norm.cdf(-z_cutoff)*s.shape[1]
        fdr = pval/(s*(z>z_cutoff) > x).any(axis=0).sum()
        PL.text(0.77, 400, "p < %.2f\nFDR=%.2f"%(pval, fdr))
        PL.xlim(0.3,1)
        PL.yticks([0,200,400,600])
        PL.ylabel("Alleles")
        PL.xlabel("Suppression cutoff")
        PL.axvline(0.75, linestyle='dashed', linewidth=1)
        
    if do_plot_scatter:
        PL.figure(None, [4,3.5])
        for i in range(10): PL.plot(s[i],abs(z[i]), ".", alpha=0.1)
        PL.xlim(-2,2); PL.ylim(0,20)
        
    Js = SP.where((s > s_cutoff).any(axis=0))[0]
    Jz = SP.where(((s > s_cutoff) & (z > z_cutoff)).any(axis=0))[0]
    if do_plot_hist:
        PL.figure(None, [3,1.6])
        
    for J in (Js, Jz):
        n_ref = 0
        alln = []
        for j in J:
            x = s[:,j]
            x = x[~SP.isnan(x)]
            alln.append(1.*sum(x>s288c_cutoff))
            if sum(x > s288c_cutoff) >= min_fraction*len(strains): # *len(x) to have a relative number
                n_ref += 1
        if (J is Jz):
            print "Of %d suppressed alleles, %d (%d%%) have at least %d wild strains suppressing at %.2f"%(len(J), n_ref, 100*n_ref//len(J), int(min_fraction*len(strains))+1, s288c_cutoff)
            if not do_plot_hist: continue
            PL.hist(alln, range=[0.5,10.5], bins=10, alpha=0.7, )
            PL.xlim(0.5,10.5)
            #PL.axvline(5.5,linewidth=1, linestyle='dashed')
            PL.xticks(range(1,11))
            PL.yticks([0,10,20,30])
            PL.xlabel("Suppressing strains (s > 0.50)"); PL.ylabel("Alleles")
            PL.text(5.9,23,"%d/%d in at\nleast 6 strains"%(sum(SP.array(alln) > 5), len(alln)))
    PL.tight_layout()
    PL.savefig("%s/paper/plots/2D_suppression_strain_count.png"%(DIR_DATA), dpi=300)
    PL.savefig("%s/paper/plots/2D_suppression_strain_count.svg"%(DIR_DATA), dpi=300)
        
    
def print_paper_gene_concordance(s_cutoff=0.75, z_cutoff=4.5):
    supp, meta, strains = read_suppression_values(keep_translocated=False)
    I = (meta[:,8] == "0") & (meta[:,2] != "CDC28") # TS, not CDC28
    supp, meta = supp[:,I], meta[I] # apply filter
    
    from collections import Counter
    orf_count = Counter(meta[:,1])
    any_vals, any_seconds, supp_seconds = [], [], []
    strongest_vals_all = [] # overview across all values
    strongest_vals = [] # null model - picking a second value is the same as first; just pick strongest per strain
    second_vals = [] # alternative - picking a second value is skewed; picking strongest for strain suppressing across other alleles more biased
    for orf in orf_count:
        for s in range(len(strains)):
            x = supp[s,meta[:,1] == orf,0,1] # mean suppression at 34 for all strains and TSQs of this ORF
            x = sorted(x[~SP.isnan(x)])[::-1]
            any_vals.extend(x)
            if len(x) > 0:
                strongest_vals_all.append(x[0]) # random selection for any strain is strongest allele
            if len(x) < 2: continue # only multi-allele orfs after filtering
                
            for i in range(len(x)):
                if x[i] > s_cutoff:
                    if i == 0:
                        supp_seconds.append(x[1])
                    else:
                        supp_seconds.append(x[0])
            strongest_vals.append(x[0]) # random selection for any strain is strongest allele
            if x[0] > s_cutoff:
                second_vals.append(x[1])
            else:
                any_seconds.append(x[1])
                
    PL.figure(None, [7.5,4])
    xmi,xma = -1,2
    PL.hist(any_vals, alpha=0.75, color='k', normed=True, range=(xmi,xma), bins=30, histtype="step")
    PL.hist(strongest_vals_all, alpha=0.95, color='c', normed=True, range=(xmi,xma), bins=30, histtype="step")
    PL.hist(strongest_vals, alpha=0.95, color='b', normed=True, range=(xmi,xma), bins=30, histtype="step")
    PL.hist(supp_seconds, alpha=0.95, color='r', normed=True, range=(xmi,xma), bins=30, histtype="step")
    print "Any allele+strain: %.2f\nStrongest strain per allele: %.2f\nStrongest strain per allele (>1 allele): %.2f\nGiven suppression of one allele by a strain, strongest other allele of same gene for same strain : %.2f"%(SP.nanmean(any_vals), SP.nanmean(strongest_vals_all), SP.nanmean(strongest_vals), SP.nanmean(supp_seconds))
    PL.legend(["Any allele per wild strain (%d)"%len(any_vals),
               "Strongest allele per gene per wild strain (%d)"%len(strongest_vals_all),
               "Strongest allele per gene wild strain (>1 allele) (%d)"%len(strongest_vals),
               "Strongest other allele per gene per wild strain for suppressing allele (%d)"%len(supp_seconds)], ncol=1)
#    for x in SP.nanmean(any_vals), SP.nanmean(strongest_vals), SP.nanmean(supp_seconds):
#        PL.axvline(x, linestyle='dashed', linewidth=1)
#        PL.text(x+0.02-0.36*(x<0.5),1.2,"%.2f"%x)
    PL.axvline(0.75, linestyle='dashed', linewidth=1, color='k')
    PL.xlabel("Suppression"); PL.ylabel("Density")
    PL.yticks([0,0.5,1,1.5]); PL.ylim(0,1.75)
    PL.savefig("%s/paper/plots/2A_suppression_dist.png"%(DIR_DATA), dpi=300)
    PL.savefig("%s/paper/plots/2A_suppression_dist.svg"%(DIR_DATA))


def print_reference_dependence(s_cutoff=0.75, z_cutoff=4.5):
    supp, meta, strains = read_suppression_values(keep_translocated=False)
    strain_names = read_wild_strain_name()
    I = (meta[:,8] == "0") # TS
    z = (supp[:,I,0,1]/supp[:,I,1,1]).T
    x = supp[:,I,0,1].T*(z > z_cutoff)
    supp, meta = supp[:,I,0,1], meta[I] # apply filter
    nums = []
    nums.append((x > s_cutoff).sum(axis=0))
    N = (x > s_cutoff).any(axis=1).sum()
    print "Ref.   R=%.2f N=%d med=%d Per strain:"%(1, N, SP.nanmedian(nums[-1])), nums[-1]
    ystrains = ["Ref (%d)"%N]
    for i in [3,5,1,4,7,8,0,2,9,6]:
        x = supp - supp[i]
        nums.append((x > s_cutoff).sum(axis=1)*1.)
        nums[-1][i] = SP.nan
        n1,n2 = list(nums[0]), list(nums[-1])
        n1 = n1[0:i] + n1[i+1:] # skip focal strain - all 0
        n2 = n2[0:i] + n2[i+1:]
        N = (x > s_cutoff).any(axis=0).sum()
        print "%s R=%.2f N=%d med=%d Per strain:"%(strains[i], SP.corrcoef(n1,n2)[0,1], N, SP.nanmedian(nums[-1])), nums[-1]
        ystrains.append("%s (%d)"%(strain_names[strains[i]], N))
    nums = SP.array(nums)
    PL.imshow(nums)
    PL.xticks(range(len(strains)), [strain_names[s] for s in strains], rotation=90)
    PL.yticks(range(len(strains)+1), ystrains, rotation=0)
    PL.colorbar()
    PL.tight_layout()
    PL.savefig("%s/paper/plots/AS3_suppression_ref_dependence.png"%(DIR_DATA), dpi=300)
    PL.savefig("%s/paper/plots/AS3_suppression_ref_dependence.svg"%(DIR_DATA))
    


def _reduce_to_gene(scores, meta):
    orfs = SP.unique(meta[:,1])
    tsqs_used = {}
    gene_score, gene_meta = SP.zeros([len(orfs), scores.shape[0]]), []
    for o, orf in enumerate(orfs):
        I = SP.where(meta[:,1] == orf)[0]
        meds = SP.nanmedian(scores[:,I],axis=0)
        J = SP.argsort(meds)
        gene_score[o] = scores[:,I[J[-1]]]
        gene_meta.append(meta[I[J[-1]]])
        tsqs_used[I[J[-1]]] = o
    return gene_score, SP.array(gene_meta)


def _get_supp_colormap():
    yellow = [255, 205, 2]
    white = [255,255,255]
    cyan = [53, 153, 204]
    colors = [list(SP.linspace(yellow[i], white[i], 30)) + list(SP.linspace(white[i],cyan[i],30)) for i in range(3)]
    from matplotlib.colors import ListedColormap
    return ListedColormap(SP.array(colors).T/255.)


def plot_paper_allele_suppression_examples(genes=("GAB1","NSE4"), outfilename="suppression_example.svg", vmi=-1.5, vma=1.5):
    supp, meta, strains = read_suppression_values(keep_translocated=False)
    I = []
    for gene in genes:
        I0 = (meta[:,2] == gene) & (meta[:,6] == "0") & (meta[:,8] == "0") # TS, not aneuploid
        I.extend(SP.where(I0)[0])
    x,names = supp[:,I,0,1], meta[I,0]
    I = [1,0,2,6,5,4,3]

    PL.figure(None,[3.5,3.5])
    PL.imshow(x[:,I], cmap=_get_supp_colormap(), interpolation='nearest', vmin=vmi, vmax=vma)
    PL.xticks(range(len(I)), names[I], rotation=90)
    PL.yticks(range(10), strains)
    PL.colorbar()
    PL.savefig("%s/paper/plots/%s"%(DIR_DATA, outfilename))


def plot_paper_gene_suppression_examples(genes, outfilename="2A_gene_suppression_example.svg", vmi=-1.5, vma=1.5):
    supp, meta, strains = read_suppression_values(keep_translocated=False)
    I = ((meta[:,6] == "0") & (meta[:,8] == "0")) # TS, not translocated
    supp[list(strains).index("Y14278"), (meta[:,7] == "1"), :,:] = SP.nan # blank out the duplicated chrII for one strain
    supp, meta = supp[:,I,0,1], meta[I]
    gene_supp, gene_meta = _reduce_to_gene(supp,meta)
    
    I = []
    for gene in genes:
        I0 = (gene_meta[:,1] == gene) # TS, not aneuploid
        I.extend(SP.where(I0)[0])
        
    PL.figure(None,[3.5,3.5])
    PL.imshow(gene_supp[I].T, interpolation='nearest', cmap=_get_supp_colormap(), vmin=vmi, vmax=vma)
    PL.xticks(range(len(I)), gene_meta[I,2], rotation=90)
    PL.yticks(range(10), strains)
    PL.colorbar()
    PL.savefig("%s/paper/plots/%s"%(DIR_DATA, outfilename))


def print_replicate_count_impact(s_cutoff=0.75, z_cutoff=4.5):
    scores, meta, strains = read_suppression_values(keep_translocated=False)
    strains = SP.array(strains)
    I = (meta[:,8] == "0") # not translocated, TS
    s = scores[:,I,0,1]
    sd = scores[:,I,1,1]
    z = s/sd
    I1 = [6,8,9]
    I2 = [0,2,4]
    print "One rep:", strains[I1], "Two reps:", strains[I2]

    for i,I in enumerate([I1,I2]):
        print "%d reps: %.1f alleles s>%.2f on average, %.1f also z>%.1f"%(i+1, SP.nanmean((s[I]>s_cutoff).sum(axis=1)), s_cutoff, SP.nanmean(((s[I]>0.75)&(z[I]>z_cutoff)).sum(axis=1)), z_cutoff)
    v1 = SP.nanstd(s[I1], axis=0)
    v2 = SP.nanstd(s[I2], axis=0)
    I = SP.isnan(v1) | SP.isnan(v2)
    print "%d/%d means of two reps more variable than one"%(sum(v2[~I]>v1[~I]), sum(~I))



def get_wild_strain_names():
    ifh = file("%s/meta/strains.txt"%DIR_DATA, 'r')
    strains = {}
    for i in range(10):
        strains[ifh.next().split("\t")[2]] = i
    return strains

# chr1	63992	0	0	0	0	YAL041W	H385Y	DELETERIOUS	0.04	2.78	27	27
def read_sift(filename="meta/SIFT-predictions-SGRP2-cerevisiae.txt"): #
    num_deleterious = {}
    orf_scores = {}
    strains = get_wild_strain_names()
    strain_i = {}
    ifh = file("%s/%s"%(DIR_DATA, filename), 'r')
    for l in ifh:
        if l[0:6] == "#CHROM":
            h = l.split("\t")
            for s in strains:
                if s in h:
                    strain_i[strains[s]] = h.index(s) - 2
        elif l[0] != "#":
            d = l.strip().split()
            orf = d[-7]
            if orf not in num_deleterious: 
                num_deleterious[orf] = SP.zeros(10)*SP.nan
                orf_scores[orf] = [[] for i in range(10)]
                #for i in strain_i: num_deleterious[orf][i] = 0
            is_del = (d[-5] == "DELETERIOUS") # is it deleterious?
            for i in strain_i:
                if (d[strain_i[i]] != ".") and SP.isnan(num_deleterious[orf][i]):
                    num_deleterious[orf][i] = 0
                if d[strain_i[i]] == "1":
                    num_deleterious[orf][i] += is_del
                    if d[-4] != ".":
                        orf_scores[orf][i].append(float(d[-4]))
    return num_deleterious, orf_scores


def calc_suppression_stratified_sift():
    ndel, osco = read_sift()
    supp, meta, strains = read_suppression_values(keep_translocated=False)
    dsupp = supp[:,:,0,1].T
    z = (supp[:,:,0,1]/supp[:,:,1,1]).T
    from collections import Counter
    counts = Counter(meta[:,1])
    suppressed_scores, suppressed_delcount, restricted_scores, restricted_delcount = [],[],[],[]

    for i in range(len(meta)): # for each strain
        orf = meta[i,1]
        if counts[orf] > 1: continue # look at 1 strain only for now
        if orf not in ndel: continue # must have SIFT values
        
        for j in range(10):
            if SP.isnan(ndel[orf][j]): continue
            if dsupp[i,j] > 0.75 and z[i,j] > 4.5: # TSQ i can be suppressed by strain j:
                suppressed_scores.extend(osco[orf][j])
                suppressed_delcount.append(ndel[orf][j])
            else: # TSQ i cannot be suppressed: 
                restricted_delcount.append(ndel[orf][j])
                restricted_scores.extend(osco[orf][j])
    return suppressed_scores, suppressed_delcount, restricted_scores, restricted_delcount


def plot_sift():
    PL.figure(None, [8,2])
    PL.subplot(121)
    supp_s, supp_n, rest_s, rest_n = map(SP.array, calc_suppression_stratified_sift())
    PL.hist(supp_n, log=True, range=(-0.5,10.5), bins=11, alpha=0.6, normed=True)
    PL.hist(rest_n, log=True, range=(-0.5,10.5), bins=11, alpha=0.6, normed=True)
    PL.xlim(-1,11)
    PL.legend(["Suppressed", "Other"], ncol=2); PL.ylabel("Number of genes"); PL.xlabel("Deleterious mutations")
    print "Suppressed genes: %.3f deleterious mutations on average; restricted:%.3f"%(SP.nanmean(supp_n), SP.nanmean(rest_n))    
    PL.subplot(122)
    PL.hist(supp_s, log=True, range=(0,1), bins=20, alpha=0.6, normed=True)
    PL.hist(rest_s, log=True, range=(0,1), bins=20, alpha=0.6, normed=True)
    PL.xlabel("SIFT score"); PL.ylabel("Number of genes")
    PL.legend(["Suppressed", "Other"], ncol=2)


def plot_followup_suppression(verbose=False, filename="S2E_suppression_followup.png"):
    x, y1, refy1, y2, refy2, z, col = read_followup_phenotypes()
    PL.figure(None, [2.5,2.2])
    vals = []
    for i in range(len(x)):
        color = "rb"[col[i]]
        marker = "."
        if refy2[i] > -1: marker = "x" # not TS - count ratio of control strain not skewed
        vals.append([x[i], y2[i]-refy2[i]])
        PL.plot(vals[i][0], vals[i][1], color + marker, markersize=12*(1.-2**refy2[i]), alpha=0.7)
        if x[i] > 0.75 and y2[i]-refy2[i] < 0 and color == "r" and marker == ".": 
            if verbose:
                print strains[i], x[i],2**y2[i], 2**refy2[i]
    xv, xmi, xma, ymi, yma = 0.75, -1, 2.5, -3,7
    PL.plot([xv,xv], [ymi, yma], 'k--', alpha=0.7, linewidth=1.2)
    PL.plot([0,0], [ymi, yma], 'k', alpha=0.5, linewidth=0.5)
    PL.plot([xmi, xma], [0,0], 'k', alpha=0.5, linewidth=0.5)
    PL.xlim(xmi,xma); PL.ylim(ymi, yma)
    PL.xticks([-1,0,1,2])
    PL.yticks([-2,0,2,4,6])
    PL.xlabel("Screen"); PL.ylabel("Followup colony count"); PL.title("Suppression")
    PL.savefig("%s/paper/plots/%s"%(DIR_DATA, filename), dpi=300)
    PL.savefig("%s/paper/plots/%s"%(DIR_DATA, filename.replace(".png",".svg")))


def write_suppression_values_followup():
    supp, meta, strains = read_suppression_values(keep_translocated=False)
    ifh = file("%s/paper/input_data/followup_strains.tab"%DIR_DATA, 'r')
    ofh = file("%s/paper/tables/followup_strain_suppression.tab"%DIR_DATA, 'w')
    
    for l in ifh:
        ofh.write(l.strip())
        tsq, _, _, strain = l.strip().split("\t")
        if strain == "DMA1":
            ofh.write("\t#N/A\n")
            continue
        s = int(strain[1:]) - 14273
        i = list(meta[:,0]).index(tsq)
        ofh.write(("\t%.3f\n"%(supp[s,i,0,1])).replace("nan", "#N/A")) # 0,1 = mean of 34C.
    ofh.close()
