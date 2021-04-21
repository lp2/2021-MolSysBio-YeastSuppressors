import scipy as SP
import pylab as PL
import glob
import os
import cPickle
from common import *
from io_tools import *
from sqtl.io_tools import read_vcf


DIR_DATA = "/Users/lp2/data/projects/wildsupp"
DIR_AFS = "%s/seq/Joint-call/afs"%DIR_DATA
DIR_META = "%s/paper/meta"%DIR_DATA
DIR_TABLE = "%s/paper/tables"%DIR_DATA
DIR_PLOT = "%s/paper/plots"%DIR_DATA
DIR_QTL = "%s/paper/tables/qtls"%DIR_DATA
DIR_QTLPLOT = "%s/paper/plots/qtls"%DIR_DATA
CHRMS = "I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI".split(" ")


def _sort_chrms(chrms):
    ROMAN = "I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI Mito MT".split(" ")
    if "XII" in chrms: return [c for c in ROMAN if c in chrms]
    if "chrXII" in chrms: return ["chr%s"%c for c in ROMAN if "chr%s"%c in chrms]
    return sorted(chrms)


def plot_example_qtl_map(tsq="TSQ200", strain="Y14281", rep="2"):
    import glob
    vcf = read_vcf(glob.glob("%s/*%s*.vcf.gz"%(DIR_AFS, strain))[0])
    sample_34 = "%s_%s_34.%s"%(tsq,strain,rep)
    sample_26 = "%s_%s_26.%s"%(tsq,strain,rep)
    i34 = list(vcf['chrII']['samples']).index(sample_34)
    i26 = list(vcf['chrII']['samples']).index(sample_26)

    skip_edges = 5
    right_end = 0
    colors = 'br'
    PL.figure(None, [183/25.4,1])

    for chrm in _sort_chrms(vcf): # for chromosomes in ascending order
        for i in [i34,i26]:
            I = SP.ones(len(vcf[chrm]['PS'][i]), bool)
            I[0:skip_edges] = False
            I[-skip_edges:] = False
            col = colors[i == i34]
            line = PL.plot(SP.array(vcf[chrm]['L'])[I] + right_end, 1.-SP.array(vcf[chrm]['PM'][i])[I], alpha=0.4, lw=2, color=col)[0]  # offset x-coordinate by the end of last chromosome

        if right_end > 0: PL.plot([right_end, right_end], [0,1], 'k--', lw=0.2, alpha=0.8) # plot separators between chromosomes
        new_right = right_end + max(vcf[chrm]['L']) # add length of chromosome
        PL.text(right_end + 0.5*(new_right - right_end), -0.17, str(chrm[3:]), horizontalalignment='center')
        right_end = new_right # update rightmost end

    PL.plot([0,right_end], [0.5,0.5], 'k', alpha=0.8, linewidth=1) # 0.5 frequency guide
    PL.xlim(0,right_end) # Decorators
    PL.xticks([]); PL.yticks([0,0.5,1], [0,0.5,1])
    PL.ylabel("Wild allele\nfrequency")
    PL.savefig("%s/paper/plots/3A_example_qtl.svg"%DIR_DATA)


def _read_qtls(strain, tsq):
    return _read_file_qtls("%s/qtls/%s_%s_34_qtls.tab"%(DIR_TABLE, tsq, strain))

    
def _read_file_qtls(filename):
    ifh = file(filename, "r")
    res = []
    for l in ifh:
        if l[0] == "#": continue
        qtlset, chrm, peak,start,c_start,c_end, end, length, af_peak, sd_peak, numsd_peak, genes = l.strip().split("\t")
        res.append((chrm, int(peak), float(af_peak), float(sd_peak), float(numsd_peak), int(start), int(end), int(c_start), int(c_end), genes))
    return res



def print_paper_qtl_stats():
    genes = read_genes_gff(DIR_META)
    genelocs = read_gene_locs(DIR_META)
    tsq_names = read_tsq_names(DIR_META)
    good_qtlmaps = read_good_qtl_maps(DIR_TABLE)
    qtls = {}
    qtl_afs = {}
    for (strain, tsq) in good_qtlmaps:
        if good_qtlmaps[(strain, tsq)] and tsq != "SN851":
            qtls[(strain, tsq)] = _read_qtls(strain, tsq)
            qtl_afs[(strain,tsq)] = SP.array([q[2] for q in qtls[(strain,tsq)]])
            
    tot_strong, tot_weak = 0, 0
    cutoff = 0.2
    for k in qtl_afs: 
        tot_strong += sum(abs(qtl_afs[k]) >= cutoff)
        tot_weak += sum(abs(qtl_afs[k]) < cutoff)
    print "%d good crosses; %d strong QTLs [%.2f average]; %d weak QTLs [%.2f average]"%(len(qtls), tot_strong, 1.*tot_strong/len(qtls), tot_weak, 1.*tot_weak/len(qtls))
    tot_strong_wild, tot_weak_wild = 0, 0
    for k in qtl_afs: 
        tot_strong_wild += sum(qtl_afs[k] < -cutoff)
        tot_weak_wild += sum((qtl_afs[k] < 0) & (qtl_afs[k] > -cutoff))
    print "Of strong QTLs, %d towards wild [%.2f of all]; of weak ones %d towards wild [%.2f of all]"%(tot_strong_wild, 1.*tot_strong_wild/tot_strong, tot_weak_wild, 1.*tot_weak_wild/tot_weak)    
    


def plot_qtl_summary_statistics(strong_qtl=0.2):
    good_qtlmaps = read_good_qtl_maps(DIR_TABLE)
    afs, good_afs = [], []
    num_crosses, num_good_crosses = 0,0

    for qtl_file in glob.glob("%s/paper/tables/qtls/*.tab"%DIR_DATA):
        qtls = _read_file_qtls(qtl_file)
        allele, strain = qtl_file.split("/")[-1].split("_")[0:2]
        is_good = good_qtlmaps[(strain,allele)]
        for qtl in qtls:
            afs.append(qtl[2])
            if is_good: good_afs.append(qtl[2])
        num_crosses += 1
        num_good_crosses += is_good

    lims = SP.arange(0.12,0.55,0.01)
    afs, good_afs = SP.array(afs), SP.array(good_afs)
    N1, N2 = len(afs), len(good_afs)
    PL.figure(None, [60/25.4,40/25.4])
    npos = SP.array([1.*sum(afs>l) for l in lims])
    nneg = SP.array([1.*sum(afs<-l) for l in lims])
    PL.plot(lims, nneg/num_crosses, 'y')
    PL.plot(lims, npos/num_crosses, 'b')
    PL.legend(["Wild", "Reference"])
    PL.axvline(strong_qtl, linewidth=1, color='k', linestyle='dashed')
    PL.xlim(0.1,0.4); PL.ylim(0,3)
    PL.yticks([0,0.5,1,1.5,2,2.0,2.5])
    PL.xticks(SP.arange(0.1,0.5,0.1))
    PL.xlabel("Allele frequency change"); PL.ylabel("Suppressor loci\nper cross")
    PL.savefig("%s/paper/plots/3B_qtl_count.svg"%DIR_DATA)


def plot_singlerepqtl_reproducibility(vcfs, sample_qtls, good_qtlmaps, tsq_names, strong_qtl=0.2, do_print_r=False, outfilename=None):
    PL.figure(None, [8,3.5])
    PL.subplot(121)
    (allx, ally, good) = _plot_singlerepqtl_afs(vcfs, sample_qtls, good_qtlmaps, strong_qtl, do_print_r)
    PL.subplot(122)
    _plot_qtl_sign_concordance(allx, allx,ally,good)
    PL.tight_layout()
    if outfilename is not None:
        PL.savefig(outfilename, dpi=300)
    return allx,ally,good

        
def _plot_qtl_sign_concordance(qtlx, allx, ally, good, x_start=0.1, step=0.01): 
    x = SP.arange(x_start-step, 0.4, step)
    I = SP.array(good, bool)
    allfrac, goodfrac = [], []
    for thresh in x:
        allfrac.append(((allx*ally > 0) & (abs(qtlx) >= thresh)).sum()*1./((abs(qtlx) >= thresh).sum()))
        goodfrac.append(((allx*ally > 0) & (abs(qtlx) >= thresh))[I].sum()*1./((abs(qtlx) >= thresh)[I].sum()))
    PL.plot(x, goodfrac, 'b', alpha=0.7)
    #PL.plot(x, allfrac, 'r', alpha=0.7)
    #PL.xlabel("QTL allele frequency change")
    #PL.ylabel("Fraction correct sign\nin other replicate")
    PL.xlim(min(x),max(x)-step)
    PL.ylim(0,1.02)
    PL.axhline(1, linestyle='dashed', linewidth=1, color='k')
    #PL.legend(["Both replicates good", "One replicate bad"], loc='lower right')


def _plot_paper_qtl_replicate_afs(vcfs, sample_qtls, tsq_names, strong_qtl=0.2, do_print_r=False):
    combx, allx, ally, good = [], [], [], []
    legend_lines, legend_texts = [None,None], [None,None]
    total, correct_sign = [0,0],[0,0]
    
    for strain in vcfs:
        for sample in sample_qtls[strain]:
            allele, strain = sample.split("_")[0:2]
            for qtl in sample_qtls[strain][sample]:
                chrm, peak, combdelta = qtl[0], int(qtl[1]), float(qtl[2])
                d = vcfs[strain][chrm]
                i = list(d['L']).index(peak)
                deltas = []
                col = 'b'
                
                for rep in [1,2]:
                    s34 = list(d['samples']).index("%s.%d"%(sample,rep))
                    s26 = list(d['samples']).index("%s.%d"%(sample.replace("_34","_26"),rep))
                    delta = d['PM'][s34][i] - d['PM'][s26][i]
                    delta = min(max(delta, -0.4), 0.4)
                    deltas.append(delta)
                    if not _read_cross_sequence_rep_good(strain, allele, str(rep)): col = 'r'

                l = PL.plot([deltas[0]], [deltas[1]], col+".", markersize=12, alpha=0.5)[0]
                if col == 'b':
                    total[0] += 1
                    correct_sign[0] += (deltas[0]*deltas[1] > 0)
                    legend_lines[0], legend_texts[0] = l, "Both replicates good"
                else:
                    total[1] += 1
                    correct_sign[1] += (deltas[0]*deltas[1] > 0)
                    legend_lines[1], legend_texts[1] = l, "One replicate bad"
                combx.append(combdelta)
                allx.append(deltas[0])
                ally.append(deltas[1])
                good.append(col=='b')
                
    xmi,xma = -0.45, 0.45
    PL.plot([xmi,xma], [xmi,xma], 'k--', linewidth=1)
    PL.plot([-strong_qtl,-strong_qtl], [xmi,xma], 'k--', linewidth=1, alpha=0.2)
    PL.plot([xmi,xma],[-strong_qtl,-strong_qtl], 'k--', linewidth=1, alpha=0.2)
    PL.plot([strong_qtl,strong_qtl], [xmi,xma], 'k--', linewidth=1, alpha=0.2)
    PL.plot([xmi,xma], [strong_qtl,strong_qtl], 'k--', linewidth=1, alpha=0.2)
    PL.xlim(xmi,xma); PL.ylim(xmi,xma)
    #PL.xlabel("Allele frequency change (Rep 1)")
    #PL.ylabel("Allele frequency change (Rep 2)")
    PL.legend(legend_lines, legend_texts, loc="upper left", numpoints=1)

    if do_print_r:
        x,y, Igood = SP.array(allx), SP.array(ally), SP.array(good)
        for lim in [0.3,0.25,0.2,0.18,0.16,0.14, 0.13]:
            I = (abs(x) < lim)
            f1 = 1.*sum(x[I]*y[I] > 0)/sum(I)
            f2 = 1.*sum(x[I&Igood]*y[I&Igood] > 0)/sum(I&Igood)
            print "Max signal=%.2f: R=%.2f Good only: %.2f Fraction consistent sign: %.2f Good only: %.2f"%(lim, SP.corrcoef(x[I], y[I])[0,1], SP.corrcoef(x[I&Igood], y[I&Igood])[0,1], f1, f2)
    print "Good total QTLs=%d, correct sign in other replicate=%d (%d%%)"%(total[0], correct_sign[0], 100*correct_sign[0]/total[0])
    print "Total QTLs with one replicate bad=%d, correct sign in other replicate=%d (%d%%)"%(total[1], correct_sign[1], 100*correct_sign[1]/total[1])
    return map(SP.array, [combx, allx, ally, good])



def _plot_singlerepqtl_afs(vcfs, sample_qtls, tsq_names, strong_qtl=0.2, do_print_r=False):
    allx, ally, good = [], [], []
    legend_lines, legend_texts = [None,None], [None,None]
    total, correct_sign = [0,0],[0,0]
    
    for strain in vcfs:
        for sample in sample_qtls[strain]:
            sample_norep = sample[0:-2]
            allele, strain = sample.split("_")[0:2]
            for qtl in sample_qtls[strain][sample]:
                chrm, peak = qtl[0], int(qtl[1])
                d = vcfs[strain][chrm]
                i = list(d['L']).index(peak)
                deltas = []
                col = 'b'
                
                for rep in [int(sample[-1]),3-int(sample[-1])]: # 1, 2 or 2,1
                    s34 = list(d['samples']).index("%s.%d"%(sample_norep,rep))
                    s26 = list(d['samples']).index("%s.%d"%(sample_norep.replace("_34","_26"),rep))
                    delta = d['PM'][s34][i] - d['PM'][s26][i]
                    delta = min(max(delta, -0.4), 0.4)
                    deltas.append(delta)
                    if not _read_cross_sequence_rep_good(strain, allele, str(rep)): col = 'r'   # if one replicate not good, color red
                    if not _read_cross_sequence_rep_good(strain, allele, str(3-rep)): col = 'r'

                l = PL.plot([deltas[0]], [deltas[1]], col+".", markersize=12, alpha=0.5)[0]
                if col == 'b':
                    total[0] += 1
                    correct_sign[0] += (deltas[0]*deltas[1] > 0)
                    legend_lines[0], legend_texts[0] = l, "Both replicates good"
                else:
                    total[1] += 1
                    correct_sign[1] += (deltas[0]*deltas[1] > 0)
                    legend_lines[1], legend_texts[1] = l, "One replicate bad"
                allx.append(deltas[0])
                ally.append(deltas[1])
                good.append(col=='b')
    xmi,xma = -0.45, 0.45
    PL.plot([xmi,xma], [xmi,xma], 'k--', linewidth=1)
    PL.plot([-strong_qtl,-strong_qtl], [xmi,xma], 'k--', linewidth=1, alpha=0.2)
    PL.plot([xmi,xma],[-strong_qtl,-strong_qtl], 'k--', linewidth=1, alpha=0.2)
    PL.plot([strong_qtl,strong_qtl], [xmi,xma], 'k--', linewidth=1, alpha=0.2)
    PL.plot([xmi,xma], [strong_qtl,strong_qtl], 'k--', linewidth=1, alpha=0.2)
    PL.xlim(xmi,xma)
    PL.ylim(xmi,xma)
    #PL.xlabel("Allele frequency change (discovery)")
    #PL.ylabel("Allele frequency change (replication)")
    PL.legend(legend_lines, legend_texts, loc="upper left", numpoints=1)

    if do_print_r:
        x,y, Igood = SP.array(allx), SP.array(ally), SP.array(good)
        for lim in [0.3,0.25,0.2,0.18,0.16,0.14, 0.13]:
            I = (abs(x) < lim)
            f1 = 1.*sum(x[I]*y[I] > 0)/sum(I)
            f2 = 1.*sum(x[I&Igood]*y[I&Igood] > 0)/sum(I&Igood)
            print "Max signal=%.2f: R=%.2f Good only: %.2f Fraction consistent sign: %.2f Good only: %.2f"%(lim, SP.corrcoef(x[I], y[I])[0,1], SP.corrcoef(x[I&Igood], y[I&Igood])[0,1], f1, f2)
    print "Good total QTLs=%d, correct sign in other replicate=%d (%d%%)"%(total[0], correct_sign[0], 100*correct_sign[0]/(1e-6+total[0]))
    print "Total QTLs with one replicate bad=%d, correct sign in other replicate=%d (%d%%)"%(total[1], correct_sign[1], 100*correct_sign[1]/(total[1]+1e-6))
    return map(SP.array, [allx, ally, good])


def _plot_paper_qtl_other_query_afs(vcfs, sample_qtls, tsq_names, strong_qtl):
    multiallele_strain = [("Y14281", ["TSQ200", "TSQ1254"]),
                          ("Y14277", ["TSQ658", "TSQ811"]), 
                          ("Y14275", ["TSQ200", "TSQ1254"]), 
                          ("Y14274", ["TSQ521", "TSQ534"])]
    
    allx, ally, good = [], [], []
    cols = 'rbgm'
    c = -1
    total, correct_sign = 0,0
    legend_texts, legend_lines = [],[]
    
    for strain, alleles in multiallele_strain:
        c += 1
        legend_done = False
        for a,allele in enumerate(alleles):
            sample = "%s_%s"%(allele, strain)
            other_sample = "%s_%s"%(alleles[1-a], strain)
            for qtl in sample_qtls[strain][sample + "_34"]:
                chrm, peak = qtl[0], int(qtl[1])
                d = vcfs[strain][chrm]
                i = SP.argmin(abs(d['L']-peak))
                m1, m2 = _get_sample_diff(d, sample, i), _get_sample_diff(d, other_sample, i)
                l = PL.plot([m1], [m2], cols[c]+".", markersize=12, alpha=0.5)[0]
                total += 1
                correct_sign += (m1*m2 > 0)
                if not legend_done:
                    legend_lines.append(l)
                    legend_texts.append("%s %s"%(strain, tsq_names[allele]))
                    legend_done = True
    #PL.xlabel("Allele frequency change (query 1)")
    #PL.ylabel("Allele frequency change (query 2)")
    xmi,xma = -0.45, 0.45
    PL.plot([xmi,xma], [xmi,xma], 'k--', linewidth=1)
    PL.xlim(xmi,xma)
    PL.ylim(xmi,xma)

    PL.legend(legend_lines, legend_texts, loc="upper left", numpoints=1)
    print "Total QTLs with multiple queries=%d, correct sign in other query=%d (%d%%)"%(total, correct_sign, 100*correct_sign/total)

    
def _plot_paper_qtl_other_strain_afs(vcfs, sample_qtls, tsq_names, strong_qtl):
    multistrain_allele = [("TSQ200", ["Y14275","Y14281"]),
                          ("TSQ1254", ["Y14275","Y14281"]), 
                          ("TSQ2884x", ["Y14279","Y14281"])]

    allx, ally, good = [], [], []
    cols = 'rbgm'
    c = -1
    legend_texts, legend_lines = [],[]
    total, correct_sign = 0,0
    
    for allele, strains in multistrain_allele:
        c += 1
        legend_done = False
        for s,strain in enumerate(strains):
            sample = "%s_%s"%(allele, strain)
            other_sample = "%s_%s"%(allele, strains[1-s])
            for qtl in sample_qtls[strain][sample + "_34"]:
                chrm, peak = qtl[0], int(qtl[1])
                d = vcfs[strain][chrm]
                m1 = _get_sample_diff(d, sample, SP.argmin(abs(d['L']-peak)))
                d = vcfs[strains[1-s]][chrm]
                m2 = _get_sample_diff(d, other_sample, SP.argmin(abs(d['L']-peak)))
                l = PL.plot([m1], [m2], cols[c]+".", markersize=12, alpha=0.5)[0]
                total += 1
                correct_sign += (m1*m2 > 0)
                if not legend_done:
                    legend_lines.append(l)
                    legend_texts.append("%s %s"%(tsq_names[allele], allele))
                    legend_done = True
    xmi,xma = -0.45, 0.45
    PL.plot([xmi,xma], [xmi,xma], 'k--', linewidth=1)
    PL.xlim(xmi,xma)
    PL.ylim(xmi,xma)
    #PL.xlabel("Allele frequency change (strain 1)")
    #PL.ylabel("Allele frequency change (strain 2)")
    PL.legend(legend_lines, legend_texts, loc="upper left", numpoints=1)
    print "Total QTLs with multiple strains per query=%d, correct sign in other wild strain=%d (%d%%)"%(total, correct_sign, 100*correct_sign/total)


def read_af_vcfs(force_reread=False):
    import cPickle
    file_pickle = "%s/paper/input_data/final_vcfs.pickle"%DIR_DATA
    
    if force_reread or (not os.path.exists(file_pickle)):
        vcfs = {vcf_file.split("/")[-1].split("_")[0]:read_vcf(vcf_file) for vcf_file in glob.glob("%s/*.vcf.gz"%DIR_AFS)}
        cPickle.dump(vcfs, file(file_pickle, 'wb'), -1)
        
    return cPickle.load(file(file_pickle, 'rb'))
    

def paper_call_qtls(af_lenient=0.08, af_stringent=0.12, strong_qtl=0.2, verbose=False):
    vcfs = read_af_vcfs()
    genes = read_genes_gff(DIR_META)
    genelocs = read_gene_locs(DIR_META)
    tsq_names = read_tsq_names(DIR_META)
    good_qtlmaps = read_good_qtl_maps(DIR_TABLE)
    sample_qtls, samplerep_qtls = {}, {}
    ofh = file("%s/paper/tables/TableS5_QTLs_20210305.xls"%DIR_DATA, 'w')
    ofh.write("#Set\tGene\tStrain\tChrm\tPeak\tStart\tCentre_start\tCentre_end\tEnd\tLength\tAF_peak\tSD_peak\tnumSD_peak\tCentre_genes\n") # write header

    for strain in ["Y%d"%i for i in range(14273,14283)]: # cleanest to re-call all QTLs here, as single rep QTLs not stored
        sample_reps, sample_ctrl = read_repmap("%s/%s_repmap.xls"%(DIR_AFS, strain)) # read contrast and replicate structure
        sample_qtls[strain] = call_qtls(strain, vcfs[strain], sample_reps, sample_ctrl, genes, good_qtlmaps, "%s/qtls"%DIR_TABLE, af_stringent=af_stringent, af_lenient=af_lenient, combined_qtl_file=ofh, verbose=verbose)
    ofh.close()


def write_simple_afs():
    vcfs = read_af_vcfs()
    for strain in vcfs:
        v = vcfs[strain]    
        samples = v['chrI']['samples']
        ofh = file("%s/paper/tables/seq_%s.xls"%(DIR_DATA, strain), 'w')
        ofh.write("#Chromosome\tLocation\tRef_allele\tAlt_allele")
        for s in samples:
            ofh.write("\t%s_ref_count\t%s_alt_count\t%s_posterior_mean\t%s_posterior_SD"%(s,s,s,s))
        ofh.write("\n")
        for chrm in ['chr' + c for c in "I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI".split(" ")]:
            assert (SP.array(samples) == SP.array(v[chrm]['samples'])).all() # make sure order correct        
            for i in range(len(v[chrm]['L'])): # for each site
                ofh.write("%s\t%s\t%s\t%s"%(chrm, v[chrm]['L'][i], v[chrm]['seq'][i][0], v[chrm]['seq'][i][1])) # output data on it
                for j in range(len(samples)):
                    ofh.write("\t%d\t%d"%(v[chrm]['D'][j,i,0], v[chrm]['D'][j,i,1])) # counts
                    ofh.write("\t%.3f\t%.3f"%(v[chrm]['PM'][j,i], v[chrm]['PS'][j,i])) # posterior mean/SD
                ofh.write("\n")
        ofh.close()
    os.system("tar cvfz %s/paper/tables/Data_D4.tar.gz %s/paper/tables/seq*.xls"%(DIR_DATA, DIR_DATA))


        
def plot_paper_qtl_reproducibility(af_lenient=0.08, af_stringent=0.12, strong_qtl=0.2):
    vcfs = read_af_vcfs()
    genes = read_genes_gff(DIR_META)
    genelocs = read_gene_locs(DIR_META)
    tsq_names = read_tsq_names(DIR_META)
    good_qtlmaps = read_good_qtl_maps(DIR_TABLE)
    sample_qtls, samplerep_qtls = {}, {}

    for strain in ["Y%d"%i for i in range(14273,14283)]: # cleanest to re-call all QTLs here, as single rep QTLs not stored
        sample_reps, sample_ctrl = read_repmap("%s/%s_repmap.xls"%(DIR_AFS, strain)) # read contrast and replicate structure
        sample_qtls[strain] = call_qtls(strain, vcfs[strain], sample_reps, sample_ctrl, genes, good_qtlmaps, "%s/qtls"%DIR_TABLE, af_stringent=af_stringent, af_lenient=af_lenient)
        samplerep_qtls[strain] = call_singlerep_qtls(strain, vcfs[strain], sample_reps, sample_ctrl, genes, good_qtlmaps, "%s/rep_qtls"%DIR_TABLE, af_stringent=0.1, af_lenient=af_lenient)
        
    PL.figure(None, [150/25.,45/25.])
    PL.subplot(131)
    (combx, allx, ally, good) = _plot_paper_qtl_replicate_afs(vcfs, sample_qtls, good_qtlmaps, strong_qtl, False)
    PL.subplot(132)
    _plot_paper_qtl_other_query_afs(vcfs, sample_qtls, tsq_names, strong_qtl)
    PL.subplot(133)
    _plot_paper_qtl_other_strain_afs(vcfs, sample_qtls, tsq_names, strong_qtl)
    PL.tight_layout()
    PL.savefig("%s/paper/plots/3CD_qtl_reproducibility.svg"%DIR_DATA)


STR_VCF_HEADER_SMOOTH_METADATA = """##smoothing=sQTLmapper v%s
##sQTLmapper_inference_rounds=%d
##sQTLmapper_recombination_rate=%.2e (events/bp)
##sQTLmapper_recombination_rate_cutoff=%.2f
##sQTLmapper_nearby_SNP_cutoff=%d bp (min distance between nearby sites)
##sQTLmapper_coverage_cap=%.1f of the median coverage
"""
STR_VCF_HEADER_SMOOTH_FORMAT = """##FORMAT=<ID=ML,Number=1,Type=Float,Description="Observed reference allele frequency (maximum likelihood estimate)">
##FORMAT=<ID=PM,Number=1,Type=Float,Description="Inferred reference allele frequency (posterior mean)">
##FORMAT=<ID=PS,Number=1,Type=Float,Description="Standard deviation of inferred reference allele frequency">
##FORMAT=<ID=FF,Number=1,Type=Integer,Description="Site deemed of good quality?">
"""


def _read_vcf_header(vcf):
    dp4_genotype, samples = False, None
    format_index = 0
    l = vcf.next()

    while l[0] == "#":
        if l[0:6] == "#CHROM":
            samples = l.strip().split("\t")[9:]
            return samples, dp4_genotype            
        dp4_genotype = dp4_genotype or l.count("##INFO=<ID=DP4") > 0 # check info field for DP4 values
        l = vcf.next()
            
    assert False, "No line that starts with #CHROM to indicate end of VCF header"

    
def _read_vcf_header_old(vcf):
    full_genotype, dp4_genotype, samples = False, False, None
    format_index = 0
    ref_index, alt_index = None, None
    l = vcf.next()

    while l[0] == "#":
        if l[0:6] == "#CHROM":
            samples = l.strip().split("\t")[9:]
            return samples, ((ref_index is not None) and (alt_index is not None)), dp4_genotype, (ref_index, alt_index)

        if l[0:9] == "##FORMAT=": # check output format info for ref/alt counts
            if l[10:15] == "ID=RO": # this index gives reference count
                ref_index = format_index
            elif l[10:15] == "ID=AO": # this index gives alt count
                alt_index = format_index
            format_index += 1
            
        dp4_genotype = dp4_genotype or l.count("##INFO=<ID=DP4") > 0 # check info field for DP4 values
        
        l = vcf.next()
            
    assert False, "No line that starts with #CHROM to indicate end of VCF header"


def _extract_ref_alt(filter_str):
    assert filter_str.count("DP4") > 0, "DP4 field does not exist in genotype information" + filter_str
    for tok in filter_str.split(";"):
        if tok.count("=") == 0: continue
        key, val = tok.split("=")
        if key == "DP4":
            counts = map(int, val.split(","))
            return counts[0] + counts[1], counts[2] + counts[3]

        
""" Read site allele frequency data from vcf file
@param file_name vcf file; can be gzipped and have .gz extension
@return map of chromosome name => {
    'L': list of loci for segregating sites,
    'seq': reference and alternative allele at these sites,
    'D': reference and alternative allele counts at these sites,
    'ML': maximum likelihood allele frequency estimate at these sites [reference allele count / (reference allele count + alternative allele count)],
    'PM': posterior mean estimate of allele frequency. Empty list if not yet calculated
    'PS': posterior mean standard deviation. Empty list if not yet calculated
    'FF': if site is considered an outlier [ML allele frequency estimate discordant from the neighbours, and estimated local mean.]
    }
@require exists(file_name)
"""
def read_vcf_old(file_name, quality_cutoff=1, index_adjust=0):
    ifh = file(file_name, 'r')
    if file_name[-3:] == ".gz": # if zipped, use gzip file open
        import gzip
        ifh = gzip.open(file_name, 'r')

    samples, full_genotype, dp4_genotype, (ref_index, alt_index) = _read_vcf_header(ifh)
    ref_index, alt_index = ref_index + index_adjust, alt_index+index_adjust
    #print file_name, samples, full_genotype, dp4_genotype, (ref_index, alt_index)
    res = {} # set up structure to return    
    
    for l in ifh:
        assert (full_genotype or dp4_genotype), "Either full genotype should be given in output, or DPI4 allele counts in INFO field"
        d = l.strip().split("\t") # read vcf line, split into constituents
        chrm, loc, ref, alt, qual, geno_filter, geno_desc, genotype = d[0], d[1], d[3], d[4], float(d[5]), d[7], d[8], map(lambda x: x.split(":"), d[9:]) # read out interesting fields, especially genotype one.
        if qual < quality_cutoff: continue
        if chrm not in res:
            res[chrm] = {'L':[], 'seq':[], 'D':[[] for i in range(len(samples))], 'ML':[[] for i in range(len(samples))], 'PM':[[] for i in range(len(samples))], 'PS':[[] for i in range(len(samples))], 'FF':[[] for i in range(len(samples))], 'samples':samples} # initialize chromosome input if needed
        res[chrm]['L'].append(int(loc)) # populate location
        res[chrm]['seq'].append([ref, alt]) # sequence
        for j in range(len(samples)):
            if genotype[j][0] == ".": # zero coverage in this sample
                res[chrm]['D'][j].append([0,0])
            elif full_genotype:
                res[chrm]['D'][j].append([int(genotype[j][ref_index]), sum(map(int, genotype[j][alt_index].split(",")))]) # allele counts (sum alternative alleles)
            elif dp4_genotype:
                res[chrm]['D'][j].append(_extract_ref_alt(geno_filter))
            if geno_desc.count("ML:PM:PS:FF") > 0: # if file genotype info already decorated with posterior estimates
                res[chrm]['ML'][j].append(float(genotype[j][-4])) # populate as well
                res[chrm]['PM'][j].append(float(genotype[j][-3]))
                res[chrm]['PS'][j].append(float(genotype[j][-2]))
                res[chrm]['FF'][j].append(bool(genotype[j][-1]))

    for chrm in res:
        for k in res[chrm]: res[chrm][k] = SP.array(res[chrm][k])
            
    return res


""" Read site allele frequency data from vcf file
@param file_name vcf file; can be gzipped and have .gz extension
@return map of chromosome name => {
    'L': list of loci for segregating sites,
    'seq': reference and alternative allele at these sites,
    'D': reference and alternative allele counts at these sites,
    'ML': maximum likelihood allele frequency estimate at these sites [reference allele count / (reference allele count + alternative allele count)],
    'PM': posterior mean estimate of allele frequency. Empty list if not yet calculated
    'PS': posterior mean standard deviation. Empty list if not yet calculated
    'FF': if site is considered an outlier [ML allele frequency estimate discordant from the neighbours, and estimated local mean.]
    }
@require exists(file_name)
"""
def read_vcf(file_name, quality_cutoff=1):
    ifh = file(file_name, 'r')
    if file_name[-3:] == ".gz": # if zipped, use gzip file open
        import gzip
        ifh = gzip.open(file_name, 'r')

    samples, dp4_genotype = _read_vcf_header(ifh)
    res = {} # set up structure to return    
    
    for l in ifh:
        d = l.strip().split("\t") # read vcf line, split into constituents
        chrm, loc, ref, alt, qual, geno_filter, geno_desc, genotype = d[0], d[1], d[3], d[4], float(d[5]), d[7], d[8], map(lambda x: x.split(":"), d[9:]) # read out interesting fields, especially genotype one.
        if qual < quality_cutoff: continue
        if chrm not in res:
            res[chrm] = {'L':[], 'seq':[], 'D':[[] for i in range(len(samples))], 'ML':[[] for i in range(len(samples))], 'PM':[[] for i in range(len(samples))], 'PS':[[] for i in range(len(samples))], 'FF':[[] for i in range(len(samples))], 'samples':samples} # initialize chromosome input if needed
        res[chrm]['L'].append(int(loc)) # populate location
        res[chrm]['seq'].append([ref, alt]) # sequence
        
        for j in range(len(samples)): # for each sample
            if genotype[j][0] == ".": # zero coverage in this sample
                res[chrm]['D'][j].append([0,0])
            elif geno_desc.count("AO") > 0:
                ref_index, alt_index = geno_desc.split(":").index("RO"), geno_desc.split(":").index("AO")
                res[chrm]['D'][j].append([int(genotype[j][ref_index]), sum(map(int, genotype[j][alt_index].split(",")))]) # allele counts (sum alternative alleles)
            elif dp4_genotype:
                res[chrm]['D'][j].append(_extract_ref_alt(geno_filter))
                
            if geno_desc.count("ML:PM:PS:FF") > 0: # if file genotype info already decorated with posterior estimates
                gi = geno_desc.split(":").index("ML")
                res[chrm]['ML'][j].append(float(genotype[j][gi])) # populate as well
                res[chrm]['PM'][j].append(float(genotype[j][gi+1]))
                res[chrm]['PS'][j].append(float(genotype[j][gi+2]))
                res[chrm]['FF'][j].append(bool(genotype[j][gi+3]))
    for chrm in res:
        for k in res[chrm]: res[chrm][k] = SP.array(res[chrm][k])
            
    return res


""" Helper function to open files for reading and writing, and set up necessary directory structures
@param vcf_filename input vcf filename
@param opts structure of options provided as input to overall program
@effects creates opts.outdir if specified, and does not exist
@return ifh input file handle [either standard or gzip file]
@return ofh output file handle [either standard or gzip file]
"""
def _open_posterior_files(vcf_filename, opts):
    outdir = "/".join(vcf_filename.split("/")[0:-1]) # outdir is output directory without the trailing /
    if len(outdir) == 0: outdir = "." # if not specified for input file name, use current directory
    outfilename = vcf_filename.split("/")[-1].replace(".vcf",".afs.vcf") # outfilename ignores directory
    
    if 'outdir' in dir(opts) and opts.outdir is not None: # if outdir specified elsewhere, use it
        outdir = opts.outdir
        if not os.path.exists(outdir): 
            LOG.info("Output directory %s does not exist; creating"%outdir)
            os.system("mkdir -p %s"%outdir)
            
    if vcf_filename[-3:] == ".gz": # if zipped file, both input and output will be opened as zipped
        import gzip
        return gzip.open(vcf_filename, 'r'), gzip.open("%s/%s"%(outdir, outfilename), 'w')
    else:
        return file(vcf_filename, 'r'), file("%s/%s"%(outdir, outfilename), 'w')


""" Write inferred allele frequencies
@param out_file output file name
@param header str(file header information)
@param means [init_mean, posterior_mean, posterior_var, posterior_stats, bad, loc, coverage]
@param counts sample->chrm->{D,L,seq}
"""
def write_posterior_vcf(vcf_filename, means, opts, qual_cutoff=1):
    ifh, ofh = _open_posterior_files(vcf_filename, opts)
    site_i = {chrm:0 for chrm in means} # current index per chromosome
    header_meta_done, header_format_done, prev_meta = False, False, "[blank to start with]"
    for l in ifh:
        if l[0] == "#": # comment line 
            if (not header_format_done) and (l.lower()[0:8] != "##format") and (prev_meta.lower()[0:8] == "##format"):
                ofh.write(STR_VCF_HEADER_SMOOTH_FORMAT)
                header_format_done = True
            ofh.write(l)
            if (not header_meta_done) and l.lower().count("fileformat") > 0: # put sQTL mapper information after file format
                ofh.write(STR_VCF_HEADER_SMOOTH_METADATA%(SQTL_VERSION, opts.rounds, opts.rec_rate, opts.rec_cutoff, opts.snp_cutoff, opts.max_num_median_coverage))
                header_meta_done = True
            prev_meta = l
        else: # data line
            d = l.strip("\n").split("\t")
            c,qual = d[0], float(d[5]) # chromsosome
            if qual < qual_cutoff: continue 
            if c in means: # output all existing vcf data, updating the 10th and 11th columns with posteriors
                i = site_i[c] # keep track of how many sites already output for this chromsosome; take next one from the list
                ofh.write("\t".join(d[0:9]) + ":ML:PM:PS:FF")
                for s in range(len(means[c][0])): # for each sample
                    ofh.write("\t" + d[9+s] + ":%.3f:%.3f:%.3f:%d"%(means[c][0][s][i], means[c][1][s][i], means[c][2][s][i], means[c][3][s][i]))
                ofh.write("\n")
                site_i[c] += 1
    ofh.close()


def call_qtls(strain, vcf, sample_reps, sample_ctrl, genes, good_qtlmaps, out_dir=".", af_lenient=0.1, sd_lenient=3, af_stringent=0.15, sd_stringent=5, length_cutoff=1000, peak_cutoff=0.03, combined_qtl_file=None, verbose=False):
    sample_qtls = {}
    for sample in sample_ctrl: # for each sample with a control; sample = [allele]_[strain]_[temp]
        if (sample.count("_34") == 0) or (sample.count("SN851") > 0):  # permissive temperature, or a control
            continue
        
        allele = sample.split("_")[0] # good quality cross?
        if (strain, allele) not in good_qtlmaps or not good_qtlmaps[(strain, allele)]:
            if verbose: print (strain, allele), "skipped for QTL mapping, as not a good cross"
            continue
            
        sample_qtls[sample] = [] # restrictive temperature, good cross - step through chromosomes and map QTLs
        for chrm in vcf:
            if chrm == "S288c_plasmid": continue # we do not map on plasmid
                
            delta_mu, delta_sd = _get_samplectrl_data(sample, sample_ctrl[sample], vcf[chrm], sample_reps) # calculate allele frequency change
            chrm_qtls = _call_chrm_qtls(SP.array(vcf[chrm]['L']), delta_mu, delta_sd, chrm, af_lenient, sd_lenient, af_stringent, sd_stringent, length_cutoff, peak_cutoff, genes[chrm])
            genelocs = read_gene_locs(DIR_META)
            query_loc = get_query_loc(read_tsq_names(DIR_META)[allele], DIR_META) # filter out query and SGA linkage
            for loc in [genelocs[g] for g in ["MAT", "CAN1", "LYP1"]] + [query_loc]:
                if (loc is not None) and (loc[0] == chrm): # if query on same chromosome
                    chrm_qtls = [qtl for qtl in chrm_qtls if (loc[1] < int(qtl[5])) or (loc[1] > int(qtl[6]))] # only leave in QTLs not in linkage with a selection marker
            sample_qtls[sample].extend(chrm_qtls)

        out_filename = "%s/%s_qtls.tab"%(out_dir, sample)
        outfile_header = STR_HEADER_CALL_REGIONS%(SQTL_VERSION, out_filename, "%s.annotated.afs.vcf.gz"%strain, sample, sample_ctrl[sample], af_lenient, sd_lenient, af_stringent, sd_stringent, length_cutoff)
        write_qtls(out_filename, sample, sample_qtls[sample], outfile_header)
        if combined_qtl_file is not None:
            append_qtls(combined_qtl_file, sample, sample_qtls[sample])            
        
    return sample_qtls


def call_singlerep_qtls(strain, vcf, sample_reps, sample_ctrl, genes, good_qtlmaps, out_dir=".", af_lenient=0.1, sd_lenient=3, af_stringent=0.15, sd_stringent=5, length_cutoff=1000, peak_cutoff=0.03):
    sample_qtls = {}

    # Compile samples - all non-controls at 34 degrees
    samples = []
    for s in sample_reps:
        for r in sample_reps[s]:
            if r.count("34") > 0 and r.count("SN851") == 0:
                samples.append(r)
                
    for sample in samples: # for each non-control sample at 34; sample = [allele]_[strain]_34.[Rep]
        ctrl_sample = sample.replace("34","26")
        allele = sample.split("_")[0] # good quality cross?
        
        if (strain, allele) not in good_qtlmaps:
            print (strain, allele), "skipped for QTL mapping, as not a good cross"
            continue
            
        sample_qtls[sample] = [] # restrictive temperature, good cross - step through chromosomes and map QTLs
        for chrm in vcf:
            if chrm == "S288c_plasmid": continue # we do not map on plasmid
                
            delta_mu, delta_sd = _get_samplectrl_data(sample, ctrl_sample, vcf[chrm], {sample:[sample], ctrl_sample:[ctrl_sample]}) # calculate allele frequency change
            chrm_qtls = _call_chrm_qtls(SP.array(vcf[chrm]['L']), delta_mu, delta_sd, chrm, af_lenient, sd_lenient, af_stringent, sd_stringent, length_cutoff, peak_cutoff, genes[chrm])
            
            query_loc = get_query_loc(read_tsq_names(DIR_META)[allele], DIR_META) # filter out query linkage
            if (query_loc is not None) and (query_loc[0] == chrm): # if query on same chromosome
                chrm_qtls = [qtl for qtl in chrm_qtls if (abs(int(qtl[1]) - query_loc[1]) > 25000)] # only leave in linkages that are at least 25kb away from query
            sample_qtls[sample].extend(chrm_qtls)

        out_filename = "%s/%s_qtls.tab"%(out_dir, sample)
        outfile_header = STR_HEADER_CALL_REGIONS%(SQTL_VERSION, out_filename, "%s.annotated.afs.vcf.gz"%strain, sample, ctrl_sample, af_lenient, sd_lenient, af_stringent, sd_stringent, length_cutoff)
        #write_qtls(out_filename, sample, sample_qtls[sample], outfile_header)
        
    return sample_qtls


def _get_samplectrl_data(test, ctrl, data, samplereps):
    mus, sds = [[],[]], [[],[]]
    for s,sample in enumerate([test, ctrl]):
        for i,sam in enumerate(data['samples']):
            if sam in samplereps[sample]:
                mus[s].append(data['PM'][i])
                sds[s].append(data['PS'][i])
    mus, precs = map(SP.array, mus), map(lambda(x): SP.array(x)**(-2), sds)
    mu = SP.array([(precs[m]*mus[m]).sum(axis=0)/precs[m].sum(axis=0) for m in [0,1]])
    var_comb = SP.array([prec.sum(axis=0)**(-1) for prec in precs])
    delta_mu, delta_sd = mu[0]-mu[1], (var_comb.sum(axis=0))**0.5
    return delta_mu, delta_sd


    
""" Workhorse function for QTL calling at given location
@param locs length-L array of segregating sites
@param delta length-L array of allele frequency change at the sites
@param sd length-L array of standard deviations of AF change at the sites
@param chrm chromosome name
@param af_lenient minimum allele frequency change to be considered interesting [0.1]
@param sd_lenient minimum allele frequency change in SD units to be considered interesting [3]
@param af_stringent minimum allele frequency change to be considered very interesting [0.15]
@param sd_stringent minimum allele frequency change in SD units to be considered interesting [5]
@param length_cutoff minimum span of region with lenient allele frequency change to be considered a QTL [1000]
@param peak_cutoff allele frequency change within which a site is considered to be in the peak region [0.03]
@param chrm_genes map of midpoint location => (gene name, gene length) for genes on the chromosome
@return list of qtls [chrm, peak loc, af difference at peak, sd of af difference at peak, af difference in SD units at peak, lenient start, lenient end, peak start, peak end, [list of genes within 10kb]]
"""
def _call_chrm_qtls(locs, delta, sd, chrm, af_lenient=0.1, sd_lenient=3, af_stringent=0.15, sd_stringent=5, length_cutoff=1000, peak_cutoff=0.03, chrm_genes=None):
    qtls = []
    lenient = ((abs(delta) > af_lenient) & (abs(delta)/sd > sd_lenient)) | SP.isnan(delta) # and whether the differences satisfy the stringent
    stringent = ((abs(delta) > af_stringent) & (abs(delta)/sd > sd_stringent)) | SP.isnan(delta) # and lenient criteria for large change.
        
    for i in range(len(locs)): # scan loci left to right
        if lenient[i] and (not SP.isnan(delta[i])): # if change at the locus leniently significant by all measures
            is_qtl = False # not a QTL yet, but start scanning
            start, end, peak = i,i,i # start, end, strongest signal of peak region, region we're confident in having the gene in
            while (end < len(locs) - 1) and lenient[end]: # while can extend peak to right, do so
                end += 1
                if (not SP.isnan(delta[end])) and stringent[end]:
                    is_qtl = True # actual QTL if the stretch of leniently significant loci has at least one properly significant change
                    if abs(delta[end]) > abs(delta[peak]):
                        peak = end # store location of strongest peak

            if is_qtl and (locs[end] - locs[start] > length_cutoff): # if strong enough change exists, and the region is long enough,
                Icentre = SP.where((locs >= locs[start]) & (locs <= locs[end]) & (abs(delta[peak]) <= abs(delta) + peak_cutoff))[0] # all loci with change within peak_cutoff
                genes = [""]
                if chrm_genes is not None:
                    genes = [chrm_genes[mid][0] for mid in chrm_genes if abs(mid - locs[peak]) < 10000]                    
                qtls.append((chrm, locs[peak], delta[peak], sd[peak], abs(delta[peak])/(sd[peak] + 1e-6), locs[start], locs[end], locs[Icentre.min()], locs[Icentre.max()], genes)) # store QTL
            elif is_qtl:
                print "Skipping QTL %s %d (af change=%.2f) as too short (%d)"%(chrm, locs[peak], delta[peak], locs[end] - locs[start])
                
            lenient[start:end+1] = stringent[start:end+1] = False # remove the stretch just observed from future consideration
    return qtls

    

def plot_combined_heatmap(bin_window=5000, sd_cutoff=2, min_allele_signal=0.1, max_shown_signal=0.4, control_only=False, sort_gene=False, good_only=False, filename=None):
    vcfs = read_af_vcfs()
    good_qtlmaps = read_good_qtl_maps(DIR_TABLE)
    # read all genome-wide differences
    samples, diffs, locs = _read_genomewide_differences(vcfs, control_only, sd_cutoff)
    chrms = ["chr%s"%c for c in "I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI".split(" ")]
    genelocs = read_gene_locs(DIR_META)
    for i in range(len(samples)):
        allele = samples[i].split("_")[1]
        query_loc = get_query_loc(read_tsq_names(DIR_META)[allele], DIR_META) 
        for chrm,loc in [genelocs[g] for g in ["MAT", "CAN1", "LYP1"]] + [query_loc]:
            j = chrms.index(chrm)
            I = SP.where(abs(SP.array(locs[i][j]) - loc) < 25000)[0]
            diffs[i][j][I] = SP.nan
    
    # bin in chunks
    data = None
    ends, mids = [0],[]
    for i in range(len(locs[0])):
        d = _get_binned(diffs, locs, i, bin_window, max_entry=max_shown_signal)
        if data is None: data = d
        else: data = SP.concatenate((data, d), axis=1)
        ends.append(data.shape[1])
        mids.append(0.5*(ends[-1] + ends[-2]))
    I = SP.nanmax(abs(data),axis=1) >= min_allele_signal # retain only regions with some signal
    data = data[I]
    samples = SP.array(samples)[I]
    
    if good_only:
        I = [i for i in range(len(samples)) if good_qtlmaps[tuple(samples[i].replace("1864", "2885x").split("_"))]] ## udpate here in read_genomewide_differences
        data = data[I]
        samples = samples[I]
        
    # plot
    tsq_names = read_tsq_names(DIR_META)
    tsq_names['SN851'] = "Control"
    strains = SP.array([s.split("_")[0] for s in samples])
    genes = SP.array([tsq_names[s.split("_")[1]] for s in samples])
    if sort_gene:
        I = SP.argsort(genes)
        genes, data, samples, strains = genes[I], data[I], samples[I], strains[I]
    PL.figure(None, [12,0.3*len(samples)+0.5])
    from matplotlib.colors import ListedColormap
    yellow, white, cyan = [255, 205, 2], [255,255,255], [53,153,204]
    #cmap = ListedColormap(SP.array([list(SP.linspace(yellow[i], white[i], 30)) + list(SP.linspace(white[i],cyan[i],30)) for i in range(3)]).T/255.)    
    cmap = ListedColormap(SP.array([list(SP.linspace(cyan[i], white[i], 30)) + list(SP.linspace(white[i],yellow[i],30)) for i in range(3)]).T/255.)    
    PL.imshow(data, interpolation='nearest', cmap=cmap, aspect='auto', vmin=-max_shown_signal, vmax=max_shown_signal)
    for e in ends[1:len(ends)-1]: PL.axvline(e, color='k', linewidth=1)
    for i in range(len(strains) - 1):
        if (sort_gene and (genes[i] != genes[i+1])) or ((not sort_gene) and (strains[i] != strains[i+1])):
            PL.axhline(i+0.5, linewidth=0.5, color='k')
    if control_only:
        samples = [s.split("_")[0] for s in samples]
    else:
        samples = ["%s %s"%(s.replace("_", " "), tsq_names[s.split("_")[1]]) for s in samples]

    PL.yticks(range(len(samples)), samples)
    PL.xticks(mids, CHRMS)
    PL.xlim(0, ends[-1])
    PL.colorbar()
    if control_only: PL.title("Allele frequency change in control cross")
    else: PL.title("Allele frequency change in cross")
    if filename is not None:
        PL.savefig(filename, dpi=300)
        PL.savefig(filename.replace("png","svg"))
    else: PL.show()



def _read_genomewide_differences(vcfs, control_only=False, sd_cutoff=2, rep="all"):
    filename = "genomewide_differences_%s.pickle"%(["tsqs", "controls"][control_only])
    
    diffs, locs = [], []
    samples = []
    if os.path.exists(filename):
        return cPickle.load(file(filename, 'rb'))
    
    for strain in vcfs:
        samplereps, samplectrl = read_repmap("%s/%s_repmap.xls"%(DIR_AFS, strain))
        for sample in samplereps:
            allele, strain, temp = sample.split("_")
            if temp != "34": # only look for restrictive temperature samples
                continue
            if control_only:
                if allele != "SN851": continue
            else:
                if allele == "SN851": continue
                
            diffs.append([])
            locs.append([])
            samples.append("%s_%s"%(strain, allele))
            test, ctrl = "%s_%s_34"%(allele, strain), "%s_%s_26"%(allele, strain)
            for chrm in _sort_chrms(vcfs[strain]): # for chromosomes in ascending order
                data = vcfs[strain][chrm]
                delta_mu, delta_sd = _get_samplectrl_data(test, ctrl, data, samplereps)
                delta_mu[SP.where(abs(delta_mu)/delta_sd < sd_cutoff)[0]] = 0
                diffs[-1].append(delta_mu)
                locs[-1].append(data['L'])
    cPickle.dump((samples, diffs, locs), file(filename, 'wb'), -1)
    return (samples, diffs, locs)


def _get_binned(data, locs, c, bin_window, max_entry=0.4):
    d = [x[c] for x in data]
    l = [x[c] for x in locs]
    maxl = max(map(SP.nanmax, l))
    n_bins = maxl//bin_window + 1
    res = SP.zeros([len(data), n_bins,2])*SP.nan # keep track of largest allele frequency changes in both direction (0:increase, 1:decrease)
#    res[:,:,1] = 10 # want smallest value, initialize to large
    
    for i in range(len(d)):
        for j in range(len(l[i])):
            k = int(l[i][j]/bin_window)
            res[i,k,0] = SP.nanmax([res[i,k,0], d[i][j]]) # largest increase
            res[i,k,1] = SP.nanmin([res[i,k,1], d[i][j]]) # largest decrease
            
    to_return = SP.zeros(res.shape[0:2]) # retain one with strongest deviation
    for i in range(res.shape[0]):
        for k in range(res.shape[1]):
            if abs(res[i,k,0]) > abs(res[i,k,1]):
                to_return[i,k] = min(res[i,k,0], max_entry)
            else:
                to_return[i,k] = max(res[i,k,1], -max_entry)
    return to_return


def _read_cross_sequence_rep_good(strain, allele, rep):
    ifh = file("%s/qtl_sequencing_qc.tab"%DIR_TABLE)
    ok26, ok34 = True, True
    for l in ifh:
        d = l.strip().split("\t")
        if (d[0] == strain) and (d[1] == allele) and (d[4] == rep):
            if d[3] == "34": ok34 = ok34 and (d[-1] == "TRUE") and (d[-2] == "TRUE") and (d[-3] == "TRUE")
            if d[3] == "26": ok26 = ok26 and (d[-1] == "TRUE") and (d[-2] == "TRUE") and (d[-3] == "TRUE")
    return ok26 and ok34


def _get_sample_diff(data, sample, i, lim_delta=0.4):
    mus, precs = {}, {}
    
    for temp in [26,34]:
        m,p = [],[]
        for rep in [1,2]:
            s_i = list(data['samples']).index('%s_%s.%d'%(sample, temp, rep))
            m.append(data['PM'][s_i][i])
            p.append((data['PS'][s_i][i])**(-2))
        m, p = SP.array(m), SP.array(p)
        mus[temp] = (m*p).sum()/(p.sum())
        precs[temp] = p.sum()
    d = mus[34] - mus[26]
    return min(lim_delta, max(-lim_delta, d))


def call_afs(do_run=True):
    # 1. Call allele frequencies for each VCF file
    cmd_sqtlm = "python /Users/lp2/prog/projects/sqtl/streamlined/sqtlm.py calc_afs -a 0.000005 -d %s %s" # recombination rate of 60 events per 12 MB; output to DIR_AFS
    DIR_VCFS = "%s/seq/afs/annotated_vcfs"%DIR_DATA
    DIR_AFS = "%s/seq/afs/af_calls"%DIR_DATA

    for vcf_file in glob.glob(DIR_VCFS + "/*.vcf.gz"):
        cmd = cmd_sqtlm%(DIR_AFS, vcf_file)
        print cmd
        if do_run:
            os.system(cmd)

    # 2. Combine VCF files as supplementary dataset D4
    if do_run:
        cmd_combine = "tar cvfz DataS4_AlleleFrequencyVCFs_201109.tar.gz Y*.gz"
        cmd = cmd_combine%(DIR_AFS, DIR_PAPER)
        print cmd
        os.system(cmd) # full VCFs
        write_simple_afs() # reduced tables



def _read_vcfs(picklefile="%s/paper/input_data/final_vcfs.pickle"%DIR_DATA):
    import cPickle    
    if not os.path.exists(picklefile):
        from sqtl.io_tools import read_vcf
        vcfs = {vcf_file.split("/")[-1].split("_")[0]:read_vcf(vcf_file) for vcf_file in glob.glob("%s/*.vcf.gz"%DIR_AFS)}
        cPickle.dump(vcfs, file(picklefile, 'wb'), -1)
    return cPickle.load(file(picklefile, 'rb'))


def _calc_mean_coverage(vcf, s):
    n, tot = 0, 0
    for chrm in vcf.keys():
        n += len(vcf[chrm]['D'][s])
        tot += vcf[chrm]['D'][s].sum()
    return 1.*tot/n


def _calc_frac_extreme_af(vcf, s, limit):
    afs = []
    for chrm in vcf:
        a = vcf[chrm]['PM'][s]
        afs.extend(a-SP.nanmean(a))
    frac_ok = SP.mean(abs(afs-SP.nanmean(afs)) < limit)
    return 1.- frac_ok


def calc_sequencing_qc_metrics(vcfs=None, do_print=False, outfilename="%s/qtl_sequencing_qc.tab"%DIR_TABLE, 
                min_marker_linkage = 0.666, # require at least 0.666 linkage ...
                marker_linkage_window=20000, # ... in 20kb window around selection marker
                min_coverage=7, # and 7x sequencing coverage
                extreme_af_limit=0.1, max_extreme_af_fraction=0.2):
    
    if vcfs is None: 
        vcfs = _read_vcfs()    
        
    tsq_names = read_tsq_names(DIR_META)
    tsq_names ["SN851"] = "Control"
    genelocs = read_gene_locs(DIR_META)

    ofh = file(outfilename, 'w')
    ofh.write("#Strain\tTS_allele\tTS_gene\tTemperature\tReplicate\tCoverage\tMAT_linkage\tCAN1_linkage\tLYP1_linkage\tTS_allele_linkage\tExtreme_AFS\tCoverage_ok?\tLinkage_ok?\tManual_ploidy_ok?\tAF_variation_ok?\n")

    for strain in sorted(vcfs):
        if do_print:
            print
            print "*", strain, "*"
            print "Sample\t\t\tGene\tCoverage MAT\tCAN\tLYP\tTSQ     Extr_AFs Coverage Linkage Manual_ploidy AF_variation"
        v = vcfs[strain]
        samples = v['chrII']['samples']

        for s in SP.argsort(samples):
            sam = samples[s] # e.g. TSQ200_Y14278_26.1
            tsq, strain, temp = sam.split("_")
            gene = tsq_names[tsq]
            ofh.write("%s\t%s\t%s\t%s\t%s\t"%(strain, tsq, gene, temp.split(".")[0], temp.split(".")[1]))
            cov = _calc_mean_coverage(v, s)
            is_low_coverage = cov < min_coverage
            if do_print: print "%s\t%s\t%d\t"%(sam, gene, cov),#, [" ", "(C)"][is_low_coverage]),
            ofh.write("%d\t"%(cov))

            for g in ["MAT", "CAN1", "LYP1", gene]: # markers that should have linkage
                if g == "Control":
                    ofh.write("N/A\t")
                    if do_print: print "N/A\t",
                    continue

                chrm, pos = genelocs[g]
                I, w = [], marker_linkage_window
                while len(I) < 6: # use at least 6 SNPs to estimate selection strength
                    I = SP.where(abs(v[chrm]['L'] - pos) < w)[0]
                    w = 1.5*w

                strongest_linkage = v[chrm]['PM'][s,I].max()
                if g == "MAT": strongest_linkage = 1-v[chrm]['PM'][s,I].min() # selected towards non-reference
                is_low_linkage = (strongest_linkage < min_marker_linkage)

                if do_print: print "%.2f%s\t"%(strongest_linkage, " *"[is_low_linkage]),
                ofh.write("%.2f\t"%(strongest_linkage))
                
            is_ploidy_bad = "%s_%s"%(tsq, strain) in ["TSQ860_Y14273", "TSQ2266_Y14276"]
            extreme_af = _calc_frac_extreme_af(v, s, extreme_af_limit)
            is_too_af_variable = extreme_af > max_extreme_af_fraction
            if do_print: print " %.2f%s "%(extreme_af," *"[is_too_af_variable]),
            ofh.write("%.2f\t"%(extreme_af))

            if do_print:
                print ["    ","(C) "][is_low_coverage and (gene != "Control")],
                print ["    ","(L) "][is_low_linkage and (gene != "Control")],
                print ["    ","(P) "][is_ploidy_bad and (gene != "Control")],
                print ["    ","(A) "][is_too_af_variable and (gene != "Control")]
            ofh.write("%s\t%s\t%s\t%s\n"%(["TRUE", "FALSE"][is_low_coverage], ["TRUE", "FALSE"][is_low_linkage], ["TRUE", "FALSE"][is_ploidy_bad], ["TRUE", "FALSE"][is_too_af_variable]))
    ofh.close()
