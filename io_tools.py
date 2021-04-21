import scipy as SP
import scipy.stats as ST
import pylab as PL
import glob
import os
from common import *


def _read_data(filename): # read raw colony size data
    ifh = file(filename, 'r')
    data = []
    for l in ifh:
        if l[0] in "(#": continue # skip comments
        d = l.strip().split("\t")
        row, col, val = int(d[0]), int(d[1]), float(d[2].replace("NA","nan"))
        data.append((row, col, val))
    return SP.array(data)


"""
Read all raw pixel size data from Bryan's or Meredith's screen
@return map of wild strain -> temperature -> plate -> 32x48 colony sizes
"""
def read_screen_colonysize_data(replicate="R1"):
    experimenter = ["bry","mer"][replicate == "R2"]
    temps = [26,34]
    strains = list(set([x.split("_")[-4] for x in glob.glob("%s/paper/input_data/mer/*/*.JPG.dat"%DIR_DATA)]))
    if experimenter == "bry": strains = list(set([x.split("_")[-1].split(".")[0] for x in glob.glob("%s/screen_data/bry/*/*.JPG.dat"%DIR_DATA)]))
    result = {strain: {temp: SP.zeros([6,32,48],float) for temp in temps} for strain in strains}

    for strain in strains:
        for temp in temps:
            for plate in range(1,7):
                pattern = "%s/screen_data/mer/T%d/laver_%s_TSQ_%d_%d.JPG.dat*"%(DIR_DATA, temp, strain, plate, temp)
                if experimenter == "bry": pattern = "%s/screen_data/bry/T%d/JVL_*_00%d_T%d_%s.JPG.dat"%(DIR_DATA, temp, plate, temp, strain)
                if len(glob.glob(pattern)) == 0:  print "No files: %s"%pattern
                else:  result[strain][temp][plate-1] = _read_data(glob.glob(pattern)[0])[:,2].reshape([32,48])
                    
    return {k.replace("BY4741", "YOR202W"):result[k] for k in result} # make sure control names aligned


def _compute_transloc_meta(orf):
    is_xvi = (orf[1:3] == "PL") and (int(orf[3:6]) >= 92) # chrXVI left arm, beyond YPL092 are included in translocation
    is_viii = orf[1] == "H" # chromosome VIII
    return "01"[is_xvi], "01"[is_viii]

    
""" Write a combined data file from all gitter output
@param data: map of replicate(experimenter)->wild strain->temperature->D colony sizes
@param meta: list of D TSQ array strain IDs
@return None
@effect creates output file of raw data
"""
def write_combined_colonysize_data(data,meta, max_sd=0.5, outfilename="%s/paper/tables/combined_colonysizes.tab"%DIR_DATA, paperfilename="%s/paper/tables/TableS1_RawColonySize_20210305.xls"%DIR_DATA, filters=(("chrIV", 0, 120000, "R2"), ("chrV", 0, 230000, "R1"),  ("chrV", 0, 230000, "R2"), ("chrXV", 710000, 765000, "R1"))):
    new_tsq_id = {"TSQ2353": "TSQ2884x",
               "TSQ1864": "TSQ2885x",
               "TSQ1877": "TSQ2886x", 
               "TSQ1879": "TSQ2887x"}
    bad_alleles = _get_bad_alleles()
    tsqnames = _read_tsq_orfnames()
    genelocs = read_genes_loc_gff()
    
    ofh = file(outfilename, "w")
    ofh.write("#TSQ_strain\tGene\tORF\tChrm\tStart\tEnd\tis_XVI_transloc?\tis_VIII?")
    for setname in sorted(data): # output header: replicate(experimenter) - wild strain - temperature
        for strain in sorted(data[setname]):
            for temp in sorted(data[setname][strain]):
                for sta in ["mean","SD"]:
                    ofh.write("\t%s-%s-T%d-%s"%(setname,strain,temp, sta))
    ofh.write("\n")
    
    bads_skipped, linkage_skipped = 0, 0
    nans = 0
    for i in range(len(meta)): # populate data
        tsq = meta[i]
        if tsq in new_tsq_id: tsq = new_tsq_id[tsq]
        if tsq in bad_alleles: 
            bads_skipped += 1
            continue # skip bad TSQs
            
        orf, name = tsqnames[tsq]
        if len(name.strip()) < 2: name = orf
            
        to_filter = {setname: False for setname in data}
        chrm, start, end = "None", "NaN", "NaN"
        
        if orf not in ["Control", "Empty"]:
            chrm, (start,end) = "chr" + "I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI".split(" ")[ord(orf[1]) - ord('A')], genelocs[orf]
            loc = 0.5*(float(start) + float(end)) # gene loc
            
            # check which wild strains to linkage filter out
            for (f_chrm, f_start, f_end, f_set) in filters: # filter range
                if (f_chrm == chrm) and (f_start <= loc) and (f_end >= loc): # if orf matches filter
                    to_filter[f_set] = True
                    
        if SP.all(to_filter.values()): # if both datasets filtered, do not output allele
            linkage_skipped += 10
            continue
            
        is_transloc, is_chrviii = _compute_transloc_meta(orf)
        ofh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%(tsq, orf, name, chrm, str(start), str(end), is_transloc, is_chrviii)) # TSQ array strain ID, ORF, gene name, chromosome, and midpoint
        for setname in sorted(data):
            for strain in sorted(data[setname]):
                for temp in sorted(data[setname][strain]):
                    m,sd = data[setname][strain][temp][i]
                    if to_filter[setname]:
                        ofh.write("\tnan\tnan")
                        linkage_skipped += 1
                    elif sd < max_sd:
                        ofh.write("\t%.2f\t%.4f"%(m,sd))
                    else:
                        ofh.write("\tnan\tnan")
                        nans += 1
        ofh.write("\n")
    LOG.debug("Skipped %d alleles designated bad, %d linkage-filtered, and wrote %d nans"%(bads_skipped, linkage_skipped, nans))
    ofh.close()
    os.system("cp %s %s"%(outfilename, paperfilename))
    

def _get_bad_alleles():
    bad_alleles = {s.strip():True for s in file("%s/bad_TSQ_strain_list.txt"%DIR_META, 'r')}
    for l in file("%s/190311_tsq-strain-status.txt"%DIR_META, 'r'):
        s, _, gene, _, _, status = l.strip().split("\t")
        if status == "Bad" or gene == "BRR2": bad_alleles[s] = True # updated list; also ignore all BRR2, as these are likely bad as well. 
    bad_alleles["EMPTY"] = True # skip empty strain
    return bad_alleles


def _read_tsq_orfnames(dir_meta=DIR_META, tsq_filename="TSQ_strains.txt"):
    res = {'SN851':('Control','Control'), 'EMPTY':('Empty','Empty')}
    
    for l in file("%s/%s"%(dir_meta, tsq_filename), 'r'):
        d = l.strip().split("\t")
        res[d[5]] = (d[4], d[2])
    return res


def read_genes_loc_gff(dir_meta=DIR_META, filename="saccharomyces_cerevisiae_R64-2-1_20150113.gff"):
    res = {'Control':('nan','nan'), 'Empty':('nan','nan')}

    for l in file("%s/%s"%(dir_meta, filename), 'r'): # step through input
        d = l.strip().split("\t")
        if len(d) < 8 or d[2] != "gene": continue # only look at genes

        chrm, start, end, annot, gene = d[0], d[3], d[4], d[8], None # parse information from GFF
        if chrm not in res: res[chrm] = {} # initialize chromosome info if needed

        for token in annot.split(";"):  # annotation tokens hold the common name of the gene in form "Name=YBR002;gene=IRA2;ENSEMBL_ID=ENSG1230001;..."
            name, val = token.split("=") # parse 'key=value' pair
            if (gene is None) and (name == "Name"): # as default, take YBR002 type form; this guarantees a value for each
                gene = val
#            elif name == "gene": # Pick gene name if present
#                gene = val
        mid = 0.5*(int(start) + int(end)) # midpoint of the gene
        length = int(end) - int(start)
        res[gene] = (int(start), int(end))
    return res


""" Read colony sizes after postprocessing
@param filters [chrm, start, end, replicate] to not consider """
def read_combined_colonysize_data():
    chrms = "I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI".split()    
    data, meta = {}, []
    
    # open file, read header, initialize observations
    ifh = file("%s/paper/tables/combined_colonysizes.tab"%DIR_DATA, 'r')
    header = ifh.next().strip().split("\t")
    sets = [h.split("-") for h in header[8:]]
    for s,st,t,stat in sets: # for each replicate-strain-temperature-[mu/sd], init
        if s not in data: data[s] = {}
        if st not in data[s]: data[s][st] = {}
        if int(t[1:]) not in data[s][st]: data[s][st][int(t[1:])] = [[],[]]
    
    # fill in data
    for l in ifh:
        d = l.strip().split("\t")
        tsq, orf, name, chrm, start, end, is_xvi, is_viii = d[0:8]
                        
        # store data according to filters
        meta.append(d[0:8])
        for i, (s,st,t,stat) in enumerate(sets): # replicate, strain, temperature
            #if s not in to_filter:
            data[s][st][int(t[1:])][stat == "SD"].append(float(d[i+8]))
            #else: data[s][st][int(t[1:])][stat == "SD"].append(SP.nan)

    # cast into arrays and return
    for s in data:
        for st in data[s]:
            for t in data[s][st]: 
                #print s, st, t, map(len, data[s][st][t])
                data[s][st][t] = SP.array(data[s][st][t], float)
                
    return data, meta


""" Write aligned, combined data. Do not output bad strains
@param strain strain to output
@param meta Kx6 array of strain metadata (gene/orf/strain/chrm/start/end)}
@param x KxTxE array of aligned values
@param fitness 2xKxT array of mean/standard deviation estimates for median phenotypes
@param supp 2xKxT array of suppression values (differences to YOR202W, aligned) """
def write_combined_strain_data(strain, meta, x, fitness, supp):
    ofh = file("%s/paper/tables/fitness_%s.tab"%(DIR_DATA,strain), 'w')
    ofh.write("#TSQ\tGene\tORF\tChrm\tStart\tEnd\tis_XVI_transloc?\tis_VIII?\tis_non-TS?")
    for t in range(2):
        for rep in "12":
            for s in ["mean", "SD"]:
                ofh.write("\tR%s_%d_%s"%(rep, [26,34][t], s))
    ofh.write("\tMean_26\tSD_26\tMean_34\tSD_34")
    ofh.write("\tMean_Supp_26\tSD_Supp_26\tMean_Supp_34\tSD_Supp_34\n")

    mus, sds = fitness[0], fitness[1]
    for i in range(x.shape[0]):
        # 1. filter
        if SP.isnan(mus[i,1]):
            LOG.debug("Skipping %s as no mean fitness estimate at 34"%(" ".join(meta[i])))
            continue

        # 2. write metadata
        ofh.write("%s"%("\t".join(meta[i])))

        # 3. write values
        for t in range(2): # two temperatures
            for j in range(2): # two experiments
                for v in range(2): # mean and sd:
                    ofh.write("\t%.4f"%(x[i,t,j,v]))
        ofh.write("\t%.4f\t%.4f\t%.4f\t%.4f"%(mus[i,0], sds[i,0], mus[i,1], sds[i,1]))
        for t in range(2): # two temperatures
            for v in range(2): # mean and sd:
                ofh.write("\t%.4f"%(supp[v,i,t]))
                if str(meta[i][6]) == "1": # chrXVI translocation region - 77,79,82 affected
                    if strain in ("Y14277", "Y14279", "Y14282"):
                        ofh.write("*")
                if str(meta[i][7]) == "1": # chrVIII aneuploid - Y14278 affected
                    if strain == "Y14278":
                        ofh.write("*")
        ofh.write("\n")
        
    ofh.close()

    
""" Write aligned, combined data. Do not output bad strains
@param meta Kx6 array of strain metadata (gene/orf/strain/chrm/start/end)}
@param fitness {strain:2xKxT array of fitness values} """
def write_combined_fitness_data(fitness, meta):
    ofh = file("%s/paper/tables/combined_fitnesses.tab"%DIR_DATA, 'w')
    ofh.write("#TSQ\tGene\tORF\tChrm\tStart\tEnd\tis_XVI_transloc?\tis_VIII?\tis_non-TS?")
    strains = sorted(fitness.keys())
    
    for strain in strains:
        for temp in [26,34]:
            for val in ["Mean", "SD"]:
                ofh.write("\t%s_%d_%s"%(strain, temp, val))
    ofh.write("\n")
    
    for i in range(len(meta)):
        x = SP.array([fitness[strain][:,i,:] for strain in strains]) # Sx2x2 values [strains x mean/sd x temps]
        
        # 1. filter for missing or uncertain data
        if SP.isnan(x[:,0]).all(): # if no values to output => skip
            LOG.debug("Skipping %s as no fitness estimate at 34 in any strain"%(" ".join(meta[i])))
            continue

        # 2. write metadata
        ofh.write("%s"%("\t".join(meta[i])))
        
        # 3. write values
        for s in range(len(strains)):
            for t in range(2): # temps 26, 34
                for v in range(2): # mean/SD
                    ofh.write("\t%.4f"%x[s,v,t])
        ofh.write("\n")
                
    ofh.close()


def read_fitness_values(n_meta_col=9, keep_translocated=True):
    ifh = file("%s/paper/tables/combined_fitnesses.tab"%DIR_DATA, 'r')
    header = ifh.next()
    strains = [h.split("_")[0] for h in header.strip().split("\t")[n_meta_col:] if h.count("_26_SD") > 0]
    meta, data = [], [[] for i in range(len(strains))]
    for l in ifh:
        d = l.strip().split("\t")
        meta.append(d[0:n_meta_col])
        for i in range(len(strains)):
            s = SP.zeros(4)*SP.nan
            for j in range(4):
                v = d[n_meta_col+4*i+j]
                if (v.count("*") == 0) or keep_translocated: # this combination from 
                    s[j] = float(v.replace("*",""))
            data[i].append([[s[0],s[2]],[s[1],s[3]]]) # First means, then SDs
    return SP.array(data), SP.array(meta), SP.array(strains)


""" Write aligned, combined data. Do not output bad strains
@param meta Kx6 array of strain metadata (gene/orf/strain/chrm/start/end)}
@param supp {strain:2xKxT array of suppression values (difference to YOR202W, aligned)} """
def write_combined_suppression_data(supp, meta, paperfilename="%s/paper/tables/TableS2_Suppression_20210305.xls"%DIR_DATA):
    outfilename = "%s/paper/tables/combined_suppression.tab"%DIR_DATA
    ofh = file(outfilename, 'w')
    ofh.write("#TSQ\tGene\tORF\tChrm\tStart\tEnd\tis_XVI_transloc?\tis_VIII?\tis_non-TS?")
    strains = ["Y%d"%s for s in range(14273, 14283)]
    
    for strain in strains:
        for temp in [26,34]:
            for val in ["Mean", "SD"]:
                ofh.write("\t%s_%d_%s"%(strain, temp, val))
    ofh.write("\n")
    
    for i in range(len(meta)):
        x = SP.array([supp[strain][:,i,:] for strain in strains]) # Sx2x2 values [strains x mean/sd x temps]

        # 1. filter for missing data
        if SP.isnan(x[:,0,1]).all(): # if no values to output at 34 => skip
            LOG.debug("Skipping %s as no suppression estimates in any strain"%(" ".join(meta[i])))
            continue

        # 2. write metadata
        ofh.write("%s"%("\t".join(meta[i])))

        # 3. write values
        for s in range(len(strains)):
            for t in range(2): # temps 26, 34
                for v in range(2): # mean/SD
                    ofh.write("\t%.4f"%x[s,v,t])
                    if str(meta[i][6]) == "1": # chrXVI translocation region - 77,79,82 affected
                        if s in (14277-14273, 14279-14273, 14282-14273):
                            ofh.write("*")
                    elif str(meta[i][7]) == "1": # chrVIII aneuploid - Y14278 affected
                        if s == (14278-14273):
                            ofh.write("*")
        ofh.write("\n")
                
    ofh.close()
    os.system("cp %s %s"%(outfilename, paperfilename))

    
def read_suppression_values(n_meta_col=9, keep_translocated=False):
    ifh = file("%s/paper/tables/combined_suppression.tab"%DIR_DATA, 'r')
    header = ifh.next()
    strains = [h.split("_")[0] for h in header.strip().split("\t")[n_meta_col:] if h.count("_26_SD") > 0]
    meta, data = [], [[] for i in range(len(strains))]
    for l in ifh:
        d = l.strip().split("\t")
        meta.append(d[0:n_meta_col])
        #x = SP.array(d[n_meta_col:], float)
        for i in range(len(strains)):
            s = SP.zeros(4)*SP.nan
            for j in range(4):
                v = d[n_meta_col+4*i+j]
                if (v.count("*") == 0) or keep_translocated: # this combination from 
                    s[j] = float(v.replace("*",""))
            #s = x[4*i:4*i+4] # previous
            data[i].append([[s[0],s[2]],[s[1],s[3]]]) # First means, then SDs
    return SP.array(data), SP.array(meta), SP.array(strains)


def _read_followup(filename="%s/paper/tables/TableS3_RandomSpo.tab"%DIR_DATA):
    ifh = file(filename, 'r')
    header = ifh.next().strip().split("\t")
    res = {}
    for l in ifh:
        d = l.replace("#N/A", "NaN").strip().split("\t")
        tsq, strain = d[1], d[4]
        if tsq == "SN851": continue
        # SuppScore	TotalArea26C	TotalArea34C	AreaSD26C	AreaSD34C	MeanArea26C	MeanArea34C	MedianArea26C	MedianArea34C	Count26C	Count34C	TotalAreaRatio	MedianAreaRatio	CountRatio
        supp = SP.array(d[6:-1],float)
        confirmed = (d[-1] == "v")
        if tsq not in res: res[tsq] = {}
        res[tsq][strain] = (supp, confirmed)
    return res, header[6:-1]


def read_followup_phenotypes():
    d,h = _read_followup()
    strains = []
    x,y1,y2,z, refy1, refy2, col = [], [], [], [], [], [],[]
    for tsq in d:
        for strain in d[tsq]:
            if strain == "DMA1": continue
            if "DMA1" not in d[tsq]: 
                print tsq
                continue
            x.append(d[tsq][strain][0][0]) # Suppression from screen
            y1.append(SP.log2(d[tsq][strain][0][-3])) # Total Area Ratio
            refy1.append(SP.log2(d[tsq]["DMA1"][0][-3]))
            y2.append(SP.log2(d[tsq][strain][0][-1])) # Count Ratio
            refy2.append(SP.log2(d[tsq]["DMA1"][0][-1]))
            z.append(SP.log2(d[tsq][strain][0][-4]+1)) # Count at *34*
            col.append(d[tsq][strain][1]) # confirmed or not
            strains.append([tsq, strain])
    return x, y1, refy1, y2, refy2, z, col


def read_wild_strain_name():
    ifh = file("%s/paper/meta/strains.txt"%DIR_DATA, 'r')
    names = {"YOR202W":"S288C"}
    for l in ifh:
        d = l.strip().split("\t")
        names[d[0]] = d[2]
    return names


def read_tsq_names(dir_meta=DIR_META,tsq_filename="TSQ_strains.txt"):
    res = {}
    for l in file("%s/%s"%(dir_meta, tsq_filename), 'r'):
        d = l.strip().split("\t")
        res[d[5]] = d[2]
    return res


def read_gene_locs(dir_meta=DIR_META, filename="gene_loc.tab"):
    res = {}
    for l in file("%s/%s"%(dir_meta, filename), 'r'):
        gene, loc = l.strip().split("\t")
        kb, chrm = loc.split(" ")
        if "chr" not in chrm: chrm = "chr" + chrm
        res[gene] = (chrm, 1000*int(kb[0:-2]))
    return res


def read_genes_gff(dir_meta=DIR_META, filename="saccharomyces_cerevisiae_R64-2-1_20150113.gff"):
    res = {}

    for l in file("%s/%s"%(dir_meta, filename), 'r'): # step through input
        d = l.strip().split("\t")
        if len(d) < 8 or d[2] != "gene": continue # only look at genes

        chrm, start, end, annot, gene = d[0], d[3], d[4], d[8], None # parse information from GFF
        if chrm not in res: res[chrm] = {} # initialize chromosome info if needed

        for token in annot.split(";"):  # annotation tokens hold the common name of the gene in form "Name=YBR002;gene=IRA2;ENSEMBL_ID=ENSG1230001;..."
            name, val = token.split("=") # parse 'key=value' pair
            if (gene is None) and (name == "Name"): # as default, take YBR002 type form; this guarantees a value for each
                gene = val
            elif name == "gene": # Pick gene name if present
                gene = val
        mid = 0.5*(int(start) + int(end)) # midpoint of the gene
        length = int(end) - int(start)
        res[chrm][mid] = (gene, length)
    return res


def read_good_qtl_maps(dir_table=DIR_TABLE, good_filename="qtl_sequencing_qc.tab"):
    pair_ok = {}
    ifh = file("%s/%s"%(dir_table, good_filename), 'r')
    h = ifh.next()
    for l in ifh:
        d = l.strip().split("\t")
        strain, allele, gene = d[0:3]
        pair = (strain, allele)
        temp, rep = d[3:5]
        ti = temp.count("34")
        coverage_ok, linkage_ok, ploidy_ok, afs_ok = (d[-4] == "TRUE"), (d[-3] == "TRUE"), (d[-2] == "TRUE"), (d[-1] == "TRUE")
        if pair not in pair_ok: pair_ok[pair] = [False, False]
        pair_ok[pair][ti] = pair_ok[pair][ti] or (ploidy_ok and linkage_ok and coverage_ok and afs_ok)
    return {p: pair_ok[p][0] and pair_ok[p][1] for p in pair_ok}


def read_repmap(filename):
    samplereps = {}
    samplectrl = {}
    for l in file(filename, 'r'):
        if l[0] == "#": continue # skip header
        rep, sample, ctrl = l.strip().split("\t")
        samplereps[sample] = samplereps.get(sample, []) + [rep]
        samplectrl[sample] = ctrl
    return samplereps, samplectrl


def get_query_loc(gene, dir_meta=DIR_META, filename="gene_loc.tab"):
    ifh = file("%s/%s"%(dir_meta, filename), 'r')
    for l in ifh:
        d = l.strip().split("\t")
        loc, chrm = d[1].split()
        if d[0] == gene: return ("chr"+chrm, 1000*int(loc.replace("kb","")))
    return None


def write_qtls(out_file, experiment, qtls, header):
    ofh = file(out_file, "w")
    ofh.write(header)
    if len(qtls) == 0:
        ofh.write("# NO QTLS FOUND\n")
        ofh.close()
        return

    ofh.write("#Set\tChrm\tPeak\tStart\tCentre_start\tCentre_end\tEnd\tLength\tAF_peak\tSD_peak\tnumSD_peak\tCentre_genes\n") # write header
    for q, (chrm, peak, d, s, sds, start, end, c_start, c_end, genes) in enumerate(qtls): # for each QTL
        ofh.write("%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.1f\t%s\n"%("%s-QTL-%d"%(experiment, q+1), chrm, peak, start, c_start, c_end, end, end-start, d, s, sds, ",".join(genes))) # output it
    ofh.close()


def append_qtls(ofh, experiment, qtls):
    tsq_names = read_tsq_names(DIR_META)
    strain_names = read_strain_names(DIR_META)
    tsq, strain = experiment.split("_")[0:2]
    for q, (chrm, peak, d, s, sds, start, end, c_start, c_end, genes) in enumerate(qtls): # for each QTL
        ofh.write("%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\t%.1f\t%s\n"%("%s-QTL-%d"%(experiment, q+1), tsq_names[tsq], strain_names[strain], chrm, peak, start, c_start, c_end, end, end-start, d, s, sds, ",".join(genes))) # output it


def read_strain_names(dir_meta=DIR_META):
    names = {}
    for l in file("%s/strains.txt"%dir_meta, 'r'):
        d = l.strip().split("\t")
        if d[0][0] == "Y":
            names[d[0]] = d[2]
    return names
