import scipy as SP
import numpy as NP
import scipy.stats as ST
import glob
from common import *
from io_tools import read_screen_colonysize_data, write_combined_colonysize_data, read_combined_colonysize_data, write_combined_strain_data, write_combined_fitness_data, write_combined_suppression_data, read_fitness_values


""" Combine technical replicates for colony size estimates
@param max_sd: maximum standard deviation between technical replicates that is output; NaNs output if larger
"""
def create_combined_colonysize_data(max_sd):
    # 1. read metadata to align screens
    platelocs = _create_platelocs()

    # 2. Create data matrix
    data = {}
    for setname in ["R1", "R2"]: # for each screen
        data[setname] = combine_locations(combine_techreps(read_screen_colonysize_data(setname), platelocs))

    # 3. output
    write_combined_colonysize_data(data, sorted(platelocs), max_sd)
    LOG.debug("Done creating and writing two-experiment comparison")


# Return TSQ => [(plate, row, col)]
def _create_platelocs():
    meta = []
    platelocs = {}
    for f in glob.glob("%s/paper/meta/plate?.tab"%DIR_DATA):
        ifh = file(f, 'r')
        p = int(f.split("/")[-1][5:6])
        for i in range(5): ifh.next()
        for l in ifh:
            c, r, tsq = l.strip("\n").split("\t")
            c, r = int(c), int(r)
            if tsq not in platelocs: platelocs[tsq] = []
            platelocs[tsq].append(((p-1,r-1,c-1)))
    return platelocs


""" Calculate one value for each array strain
@param data: strain->temp->Px2Rx2C float array of colony sizes in pixels
@param meta length-D array of TSQ strain IDs
@param platelocs map of TSQ strain ID -> list of (plate, row, col) triplets where strain is found
@param var_explained_cutoff: amount of variance explained by one colony to be considered an outlier
@return d' : strain->temp->length-D array; d'[s][t][p,r,c] = median(d[s][t][p, 2r:2r:2, 2c:2c+2])"""
def combine_techreps(data, platelocs, var_explained_cutoff=0.9):
    res = {}
    
    for strain in data:
        res[strain] = {}
        for temp in data[strain]:
            plates,rows,cols = data[strain][temp].shape
            d = SP.zeros([plates,rows/2,cols/2,2])
            for p in range(plates):
                for r in range(rows/2):
                    for c in range(cols/2):
                        x = list(data[strain][temp][p,(2*r):(2*r+2),(2*c):(2*c+2)].reshape(4))
                        d[p,r,c] = _lg1p_median_filter_techrep(x, var_explained_cutoff=var_explained_cutoff)
            res[strain][temp] = [[d[p,r,c] for (p,r,c) in platelocs[ts]] for ts in sorted(platelocs)]
    return res


""" Filter outliers, return median of technical replicates
@param x list of replicate values
@param var_explained_cutoff total amount of variance explained by one colony to be considered an outlier
@return log2-scale median of non-outlier colonies, and standard deviation of posterior mean estimate [sqrt(var/n)] in log2-scale
"""
def _lg1p_median_filter_techrep(x, var_explained_cutoff=0.9, min_size=8, prior=1.):
    to_remove = None
    v = SP.var(x)
    for i in range(len(x)):
        if SP.var(x[0:i] + x[i+1:]) < (1-var_explained_cutoff)*v:
            to_remove = i
    if to_remove is not None:
        x[to_remove] = SP.nan
    if min_size is not None:
        x[x < min_size] = SP.nan
    return SP.log2(ST.nanmedian(x)+prior), ST.nanstd(SP.log2(SP.array(x,float)+prior))/(sum(~SP.isnan(x))**0.5)

    
""" Combine locations, if a TSQ is present multiple times in an array
@param data: map of strain->temperature->TSQ->[[log2 scale colony size; sd] x locations on plates]
@return map of strain->temperature->TSQ->[median log2 scale colony size, SD estimate]
"""
def combine_locations(data, var_prior=1e-4):
    d = {}
    for strain in data: # for each strain
        d[strain] = {}
        for t in data[strain]: # for each temperature
            d[strain][t] = SP.zeros([len(data[strain][t]),2])
            for i in range(len(data[strain][t])): # for each TSQ
                x = SP.array(data[strain][t][i])
                d[strain][t][i][0] = ST.nanmedian(x[:,0], axis=0) # median across positions of techreps
                # for combined variance of mean, treat them as technical replicates, and reduce variance thanks to N reps as sqrt(var/N)
                d[strain][t][i][1] = ((ST.nanmean(x[:,1]**2, axis=0) + var_prior)/(sum(~SP.isnan(x[:,1]))))**0.5 # mean variance
    return d



""" Fit a line y=ax that optimizes the _distance_ of points to the line (rather than simply MSE) """
def _fit_ts_values(x,y):
    x2, xy, y2 = SP.dot(x,x), SP.dot(x,y), SP.dot(y,y)
    a0,a1,a2 = xy, y2-x2, -xy
    s = (a1*a1-4*a0*a2)**0.5 # solve the quadratic
    return 0.5*(-a1-s)/a2


""" Calculate fitness difference between two strains
@param f: fitness values of strain [mean/SD x strains x temperatures]
@param f0: fitness values of control [mean/SD x strains x temperatures]
@param max_bulk_diff: maximum fitness difference between f and f0 to be considered in alignment
@param ts_cutoff: minimum fitness difference between permissive and restrictive temperatures in f0 to be considered temperature sensitive
"""
def calc_fitness_diff(f, f0, max_bulk_diff=2, ts_cutoff=0.2):
    diff = SP.zeros([2, f.shape[1], 2])*SP.nan # Output suppression: [mean, SD] x strains x [26, 34]
    Its = f0[0,:,0] - f0[0,:,1] > ts_cutoff # TS if (control strain at 26) - (control strain at 34) > (TS cutoff)
    
    for t in range(2): # for 26 and 34 degrees:
        x, y = f[0,:,t], f0[0,:,t] # means
        xv, yv = f[1,:,t]**2, f0[1,:,t]**2 # variances

        # 0. identify reasonable strains
        I = ~(SP.isnan(x) | SP.isnan(y)) & (abs(x-y) < max_bulk_diff) # filter for not NaN and small bulk differences
        if t == 0: # if permissive temperature, 
            I = I & (y > 10.2) & (x > 10) # ignore small colony sizes that skew results
            
        # 1. align overall - medians between non-TS strains; slopes of reasonable strains
        mx, my = ST.nanmedian(x[I&(~Its)]), ST.nanmedian(y[I&(~Its)]) # medians of reasonable controls - non-TS and not a large difference in fitting the bulk
        r = ST.linregress(x[I],y[I])[2]
        slope = _fit_ts_values(x[I]-mx, y[I]-my)
        diff[0,:,t] = (x - mx)*slope - (y - my) # suppression mean value is growth in strain - growth in control
        diff[1,:,t] = (yv + xv*slope*slope - 2*slope*r*((yv*xv)**0.5))**0.5 # standard deviation is the sqrt of variance of the mean.
        # var((x-mx)*slope - (y-my)) = var(x-mx)*slope*slope + var(y-my) - 2*slope*cov(x-mx, y-my)
        # cov(x-mx, y-my) = r*sd(x-mx)*sd(y-my); leading to the experssion above.
    return diff


def _update_meta_is_ts(meta, ctrl_fitness, ts_cutoff):
    for i in range(len(meta)):
        is_nonts = ctrl_fitness[0,i,0] - ctrl_fitness[0,i,1] < ts_cutoff
        meta[i].append("01"[is_nonts]) # mean ctrl fitness at 26 is not much more than mean ctrl fitness at 34
    return meta


""" For each query wild strain, calculate TS allele cross phenotypes and suppression effects
@param None
@return None
@effects creates a per-strain output file, combined fitness and suppression output files
"""
def create_suppression_values(ts_cutoff=0.2, do_rep_plot=False, do_rep_mean_plot=False):
    data, meta = read_combined_colonysize_data()
    temps = [26,34] # temperatures used for combining
    expts = ["R1","R2"] # experimenters
    aligned, fitness, suppression = {}, {}, {} # outputs
    
    for strain in ["YOR202W"] + ["Y%d"%s for s in range(14273, 14283)]:
        aligned[strain] = calc_repaligned_vals(strain, data, meta, expts, temps, do_plot=do_rep_plot) # align replicates (linear transformation)
        fitness[strain] = calc_repaveraged_vals(aligned[strain]) # calculate fitness values based on aligned replicates
        suppression[strain] = calc_fitness_diff(fitness[strain], fitness["YOR202W"])
        if strain == "YOR202W": meta = _update_meta_is_ts(meta, fitness[strain], ts_cutoff=ts_cutoff)
        write_combined_strain_data(strain, meta, aligned[strain], fitness[strain], suppression[strain]) # output per strain
        LOG.debug("Wrote strain fitness data for %s"%strain)

    if do_rep_mean_plot: _plot_rep_vs_mean(aligned, fitness)
        
    write_combined_fitness_data(fitness, meta)
    write_combined_suppression_data(suppression, meta)

        
""" Align values from different screens, and calculate the median as the final phenotype
@param strain strain to analyse
@param data map experimenter->strain->temperature->[value]
@param expts list of experimenters
@param temps list of temperatures
@return aligned_vals - experimental data aligned to to_align
"""
def calc_repaligned_vals(strain, data, meta, expts, temps, prior=SP.log2(16+1.), do_plot=False):
    aligned_vals = SP.zeros([len(meta), len(temps), len(expts), 2])*SP.nan # values aligned to most dense screen for both temperatures; mu/SD
    
    # align second rep to 1st
    for i,rep in enumerate(["R1", "R2"]):
        if strain not in data[rep]: 
            LOG.debug("Skipping %s %s - no data"%(strain, rep))
            continue # this strain not measured by this experimentor
        for t in range(2): # for both temperatures
            x = SP.array(data[rep][strain][temps[t]]).T
            x[x[:,0] < prior-1] = SP.nan  # "-1" to keep consistency with a previous version
            # the 2**x in next line is to keep code consistency with previous linear versions
            if i == 0: # Rep1 - copy over
                aligned_vals[:,t,i,0] = SP.log2(2**x[:,0] + 2**prior - 2) # both v and prior have a "+1"; this only to keep consistency with previous code
                aligned_vals[:,t,i,1] = x[:,1] # SD
            else:      # Rep2 - align to Rep1
                xval = SP.log2(2**x[:,0]+2**prior-2) # -2 as above
                ref = aligned_vals[:,t,0,0] # Rep1 xval
                aligned_vals[:,t,i] = align(xval, ref, xsd=x[:,1], sdref=aligned_vals[:,t,0,1], Inonts=aligned_vals[:,0,0,0]-aligned_vals[:,1,0,0] < 0.2)
                
    if not do_plot:
        return aligned_vals # the rest is plotting
        
    for t in range(2):
        PL.subplot(1,2,t+1)
        PL.plot(aligned_vals[:,t,0,0], aligned_vals[:,t,1,0], ".", alpha=0.2)
        #I = aligned_vals[:,0,0,0] - aligned_vals[:,1,0,0] < 0.25
        #PL.plot(aligned_vals[I,t,0,0], aligned_vals[I,t,1,0], "k.", markersize=3)        
        PL.plot([7,12],[7,12], 'k--')
        PL.xlim(7,12); PL.ylim(7,12)
        
    return aligned_vals


def _cor(x,y):
    I = SP.isnan(x) | SP.isnan(y) | SP.isinf(x) | SP.isinf(y)
    return SP.corrcoef(x[~I], y[~I])[0,1]

"""
Calculate final fitness values given aligned measurements
@param aligned_vals: KxTxEx2 array of values [K strains, T temperatures, E experimentors, [mean/sd]]
@param control_size: value to normalize control fitnesses to
@param sd_to_filter: minimum SD of values to consider for jackknife test of outliers (no need to finetune non-varying strains more) [0.2]
@param frac_var_explained: fraction of variance explained by a single measurement for it to be called an outlier and filtered out [0.9]
@return final KxT array of mean fitness values
@return sds KxT array of SD of mean fitness values
@return filtered KxTxE boolean array of whether a reading was filtered out
"""
def calc_repaveraged_vals(aligned_vals, control_size=11.):
    K,T = aligned_vals.shape[0:2]
    permissive_med = ST.nanmedian(aligned_vals[:,0,0,0]) # replicate 1, 26 degrees, mean
    if not SP.isnan(permissive_med):
        aligned_vals[:,:,:,0] -= (permissive_med - control_size) # adjust values so that controls at 26 are exactly the same size for all BRY strains and temperatures if present
        
    mus = SP.zeros([K,T])*SP.nan
    sds = SP.zeros([K,T])*SP.nan

    for t in range(T): # for each temperature
        I = SP.where((~SP.isnan(aligned_vals[:,t,:,0])).sum(axis=0) > 0)[0] # at least some mean measurements present at this temperature
        for k in range(K): # for each strain
            v = SP.array([aligned_vals[k,t,i] for i in I]) # take values present
            if len(v) == 1:
                mus[k,t] = v[0,0]
                sds[k,t] = v[0,1]
            else:
                mus[k,t] = NP.nanmedian(v[:,0]) # Record median and next, uncertainty
                sds[k,t] = SP.nanmedian(v[:,1]) # uncertainty is for replicate 1 [always measured; 4 strains have rep2 unmeasured]
                delta = abs(v[0,0]-v[1,0]) # difference between aligned values
                # if all values and sds present - makes sense to compute combined SD across biological nad technical replicates
                if (~SP.isnan(aligned_vals[k,t])).all():
                    v1, v2, vall = v[0,1]**2, v[1,1]**2, 2*(delta/2.)**2 # rep1 mean ~ N(mu, var_tech1 + var_all); rep2 mean ~ N(mu, var_Tech2 + var_all)
                    a1, a2 = 1./(vall + v1), 1./(vall + v2)
                    vpost = 1./(a1 + a2)
                    mpost = (a1*v[0,0] + a2*v[1,0])/(a1+a2)
                    mus[k,t] = mpost
                    sds[k,t] = vpost**0.5
                    
    return SP.array([mus, sds])


def _plot_rep_vs_mean(aligned, fitness):        
    import pylab as PL
    PL.figure(None, [15,6])
    pi = 1
    
    for rep in range(2):
        for strain in sorted(aligned):
            if strain in ['Y14282','Y14281','Y14278','Y14279']: continue # one rep
            PL.subplot(3,7,pi); pi += 1
            PL.plot(fitness[strain][0,:,1], aligned[strain][:,1,rep], ".", alpha=0.2) # [:,1] gives values at 34; fitness[0] is mean
            PL.plot([4,11],[4,11], 'k--', linewidth=1)
            PL.xlim(7,11); PL.ylim(7.5,11.5)
            PL.xticks(range(7,12)); PL.yticks(range(7,12))
            PL.ylabel("Rep %d"%(rep+1))
            if rep == 0: PL.title(strain)
            else: PL.xlabel("Mean 34C")
    for strain in sorted(aligned):
            if strain in ['Y14282','Y14281','Y14278','Y14279']: continue # one rep
            PL.subplot(3,7,pi); pi += 1
            x = fitness[strain][1,:,1]
            x[x>0.5] = 0.5
            PL.hist(x, bins=10, range=(0,0.5)) # [:,1] gives values at 34; fitness[0] is mean
            PL.ylim(0,350); PL.yticks([0,100,200,300])
    PL.tight_layout()

    

""" Align vector x to reference ref, matching mean and slope
@param xval vector to align
@param ref vector to align to
@param lim max difference between aligned and reference data to be considered in finetuning iterations
@param minsize minimum size of colony to consider (empties huge variance, small ones also large variance)
@param rounds number of refinement rounds
@param debug whether to plot the aligned values
@return x0 offset and scaled to match ref as well as possible in mean squared error sense
"""
def align(xval, ref, lim=0.75, minsize=SP.log2(8+1), rounds=3, Inonts=None, debug=False, xsd=None, sdref=None):    
    I0 = SP.isnan(xval) | SP.isnan(ref) | (xval < minsize) | (ref < minsize) | (sdref > 0.5) # filter on ok values - finite, not nan, well measured
    if xsd is not None: I0 = I0 | (xsd > 0.5)
    x, y = xval[~I0], ref[~I0] # discard rest for fit
    do_plot = False
    if do_plot:
        import pylab as PL
        PL.figure(None, [12,2.5])
        PL.subplot(141)
        PL.plot(x,y,".", alpha=0.2)
        PL.plot(x[Inonts[~I0]],y[Inonts[~I0]],"k.", markersize=2)
        PL.plot([7,12],[7,12],'k--')
        PL.xlim(7,12); PL.ylim(7,12)
    mx, my = 0, 0
    if Inonts is not None: # if given set medians of non-ts to be equal
        nx, ny = len(x)*0.1, len(y)*0.1
        I = abs(xval-ref) < lim
        mx, my = ST.nanmedian(xval[I & Inonts & (~I0)]), ST.nanmedian(ref[I & Inonts & (~I0)])
        mx, my = (sorted(xval[I & (~I0)]))[int(nx)], sorted(ref[I & (~I0)])[int(ny)]
        
    x -= mx
    y -= my # set y median same as x
    slope = _fit_ts_values(x,y)
    intercept = 0
    #slope, intercept = ST.linregress(x, y)[0:2] # get first alignment
    if do_plot:
        PL.subplot(142)
        PL.plot(x,y,".", alpha=0.2)
        PL.plot(x[Inonts[~I0]],y[Inonts[~I0]],"k.", markersize=2)
        PL.plot([-4,1],[-4,1],'k--')
        PL.xlim(-4,1); PL.ylim(-4,1)
#    slope, intercept = 1, 0 # first alignment is identity
    I = abs(x-y) < lim
    for i in []: #range(rounds):
        slope, intercept, rv, pv, se = ST.linregress(x[I], y[I]) # re-align
        I = abs(slope*x + intercept - y) < lim # retain only nicely fitting ones
    if do_plot:
        PL.subplot(143)
        PL.plot(x*slope+intercept,y,".", alpha=0.2)
        PL.plot(slope*x[Inonts[~I0]]+intercept,y[Inonts[~I0]],"k.", markersize=2)
        PL.plot([-4,1],[-4,1],'k--')
        PL.xlim(-4,1); PL.ylim(-4,1)
    
    xval = (xval-mx)*slope
    xval += (intercept + my)
    if do_plot:
        PL.subplot(144)
        PL.plot(xval,ref,".", alpha=0.2)
        PL.plot(xval[Inonts],ref[Inonts],"k.", markersize=2)
        PL.plot([7,12],[7,12],'k--')
        PL.xlim(7,12); PL.ylim(7,12)
    
    if xsd is None: # if no SD, return mean
        return xval
    result = SP.zeros([len(xval),2]) # otherwise construct array of mean and variance
    result[:,0] = xval
    result[:,1] = xsd*slope # var(ax) = a^2var(x)
    return result



def plot_detailed_ref_ts(ts_cutoff=0.2):
    fitness, _, strains = read_fitness_values()
    ref_fitness = fitness[list(strains).index("YOR202W"),:,0,:] # 0 = mean; 1 = SD
    import pylab as PL
    PL.figure(None, [3,1.5])
    PL.hist(ref_fitness[:,0] - ref_fitness[:,1], range=(-0.5,2.5), bins=20, alpha=0.7)
    PL.axvline(ts_cutoff, color='r', linewidth=1)   
    PL.xlabel("Control fitness defect"); PL.ylabel("Number of alleles")
