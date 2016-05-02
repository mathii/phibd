# Estimating ibd and relatedness in pseudo-haploid samples

from __future__ import division, print_function
import sys, argparse, pyEigenstrat, itertools
from multiprocessing import Pool
import phibd_hmm, phibd_interpret
import numpy as np

################################################################################
#Default human-centric!

AUTOSOMES=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
           "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21"]

################################################################################

def parse_options():
    """
    arguments -h for help
    """
    parser=argparse.ArgumentParser()
    parser.add_argument('-e', '--eigenstrat', type=str, default="", help=
                        "Root for eigenstrat files - i.e {root.snp, root.geno, root.ind}")
    parser.add_argument('-i', '--individuals', type=str, default="", help=
                        "List of individuals to include (default all)")
    parser.add_argument('-p', '--pairs', type=str, default="", help=
                        "List of pairs to test (default individuals*individuals)")
    parser.add_argument('-m', '--min_chunk', type=float, default=10.0, help=
                        "Filtered results remove chunks less than this size")
    parser.add_argument('-n', '--ncore', type=int, default=1, help=
                        "number of cores to parallelize IBD computation")
    parser.add_argument('-t', '--max_iters', type=int, default=15, help=
                        "Maximum number of training iterations")
    parser.add_argument('-l', '--tolerance', type=float, default=0.001, help=
                        "Heterozygosity tolerance in training algorithm")
    parser.add_argument('-a,', '--auto', type=str, default=",".join(AUTOSOMES), help=
                        "Comma separated list of chromosomes to treat as autosomes")
    parser.add_argument('-x,', '--chrx', type=str, default="23".join(AUTOSOMES), help=
                        "Chromosome to treat as X")
    parser.add_argument('-y,', '--chry', type=str, default="24".join(AUTOSOMES), help=
                        "Chromosome to treat as Y")

    options=parser.parse_args()
    options.auto=options.auto.split(",")
    options.chromosomes=options.auto+[options.chrx, options.chry]

    return options

################################################################################

def make_pairs(data, options):
    """
    If we passed a list of pairs, use all those pairs
    Otherwise if we passed a list of individuals, use all pairs containing those individuals
    Otherwise use all pairs.
    """
    if options.pairs:
        return [x[:-1].split() for x in open(options.pairs, "r").readlines()]
    elif options.individuals:
        raise Exception("Individuals not implemented")
    else:
        return [x for x in itertools.combinations(data.ind["IND"], 2)]
    
################################################################################

def get_job(data, pair, options):
    """
    Get the match/mismatch and position matrix for chromosome each pair
    """
    p0=pair[0]
    i0=np.where(data.ind["IND"]==p0)[0][0]
    g0=data.geno()[:,i0]
    p1=pair[1]
    i1=np.where(data.ind["IND"]==p1)[0][0]
    g1=data.geno()[:,i1]
    
    job={"pair":pair, "options":options}
    for chrom in  options.chromosomes:
        include = data.snp["CHR"]==chrom
        g0c=g0[include]
        g1c=g1[include]
        nonmissing=np.logical_and(g0c!=9, g1c!=9)
        states=(g0c==g1c)[nonmissing]
        pos=data.snp["POS"][include][nonmissing]
        job["chr"+chrom]={"states":states, "pos":pos}
    
    return job
    
################################################################################

def make_jobs(data, options):
    """
    Create all the jobs
    """
    pairs=make_pairs(data, options)
    jobs=[get_job(data, x, options) for x in pairs]
    return jobs
    
################################################################################

def main(options):
    """
    Run
    """
    data=pyEigenstrat.load(options.eigenstrat)
    print("Building job list", file=sys.stderr)
    jobs=make_jobs(data, options)
    print("Detecting IBD", file=sys.stderr)
    if options.min_chunk>0:
        print("Restricting to chunks > "+str(options.min_chunk)+"Mb", file=sys.stderr)

    if options.ncore>1:
        print("Using "+str(options.ncore)+" cores", file=sys.stderr)
        pool=Pool(options.ncore)
        results=pool.map(estimate_sharing, jobs)
        pool.close()
    else:
        results=[estimate_sharing(job) for job in jobs]
    print("Interpreting results", file=sys.stderr)
    phibd_interpret.simple_autosomes(results)

################################################################################

def estimate_sharing(job):
    """
    the job is a dictionary containing one key "pair" which tells us which
    individuals are involved, and one key "chrX" for X in 1..24. Each of these 
    "chr" values is a dictionary with "pos" and "state" giving us the position of
    each marker, and whether these individuals share it, respectively. 
    
    We return, a dictionary containing the shared chunk length distributions 
    aggregated over the autosomes and X
    """
    chunks=[]
    het=[]
    hmms=[]
    for chrom in job["options"].auto:
        hmms.append(phibd_hmm.hmm2(job["pair"], job["chr"+chrom]["states"], job["chr"+chrom]["pos"]))
        
    multi_hmm=phibd_hmm.multi_hmm(hmms, tolerance=job["options"].tolerance, max_iters=job["options"].max_iters)
    multi_hmm.train()
    chunks=multi_hmm.get_chunks()
    min_length_b=job["options"].min_chunk*1e6
    lengths_filtered=np.array([[sum([z for z in x if z>min_length_b]) for x in y] for y in chunks])
    lengths=np.array([[sum(x) for x in y] for y in chunks])
    counts=np.array([[len(x) for x in y] for y in chunks])
    means=np.array([[np.mean(x) if len(x) else np.NaN for x in y ] for y in chunks])
    auto_total=np.sum(lengths)
    auto_total_filtered=np.sum(lengths_filtered)
    auto_state_total=np.sum(lengths, axis=0)
    auto_state_total_filtered=np.sum(lengths_filtered, axis=0)
    auto_state_proportions=np.sum(lengths, axis=0)/auto_total
    auto_state_proportions_filtered=auto_state_total_filtered/auto_total_filtered

    return {"pair":job["pair"],
            "p":multi_hmm.p,
            "auto_SNPs":sum([len(job["chr"+x]["states"]) for x in AUTOSOMES]),
            "auto_total":auto_total,
            "auto_state_total":auto_state_total,
            "auto_state_proportions":auto_state_proportions,
            "auto_chunks":chunks,
            "auto_lengths":lengths,
            "auto_counts":counts,
            "auto_means":means,
            "auto_total_filtered":auto_total_filtered,
            "auto_state_total_filtered":auto_state_total_filtered,
            "auto_state_proportions_filtered":auto_state_proportions_filtered,
            "auto_lengths_filtered":lengths_filtered,
            }
    
    
################################################################################

if __name__=="__main__":
    options=parse_options()
    main(options)

