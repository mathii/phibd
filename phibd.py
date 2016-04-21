# Estimating ibd and relatedness in pseudo-haploid samples

from __future__ import division, print_function
import sys, argparse, pyEigenstrat, itertools
import phibd_hmm, phibd_interpret
import numpy as np

################################################################################
#Default human-centric!

AUTOSOMES=["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11",
           "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21"]
CHRX=["23"]
CHRY=["24"]
CHROMS=AUTOSOMES+CHRX+CHRY

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

    return parser.parse_args()

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

def get_job(data, pair):
    """
    Get the match/mismatch and position matrix for chromosome each pair
    """
    p0=pair[0]
    i0=np.where(data.ind["IND"]==p0)[0][0]
    g0=data.geno()[:,i0]
    p1=pair[1]
    i1=np.where(data.ind["IND"]==p1)[0][0]
    g1=data.geno()[:,i1]
    
    job={"pair":pair}
    for chrom in  CHROMS:
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
    jobs=[get_job(data, x) for x in pairs]
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
    for chrom in AUTOSOMES:
        hmm=phibd_hmm.hmm2(job["pair"], job["chr"+chrom]["states"], job["chr"+chrom]["pos"])
        chunks.append(hmm.get_chunks())
        
    lengths=np.array([[sum(x) for x in y] for y in chunks])
    auto_total=np.sum(lengths)
    auto_state_total=np.sum(lengths, axis=0)
    auto_state_proportions=np.sum(lengths, axis=0)/auto_total

    return {"pair":job["pair"],
            "auto_SNPs":sum([len(job["chr"+x]["states"]) for x in AUTOSOMES]),
            "auto_total":auto_total,
            "auto_state_total":auto_state_total,
            "auto_state_proportions":auto_state_proportions,
            "auto_lengths":lengths}
    
    
################################################################################

if __name__=="__main__":
    options=parse_options()
    main(options)

