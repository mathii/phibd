#Functions for interpreting the distribution of chunks. 

from __future__ import division, print_function
import numpy as np

def simple_autosomes(results):
    """
    Just count the proportion that is not IBD0
    """
    IBD1s=[x["auto_state_proportions"][1] for x in results]
    IBD2s=[x["auto_state_proportions"][2] for x in results]
    proportions=[sum(x["auto_state_proportions"][1:]) for x in results]
    proportions_filtered=[sum(x["auto_state_proportions_filtered"][1:]) for x in results]
    pairs=[x["pair"][0]+"\t"+x["pair"][1] for x in results]
    means=[np.nanmean(x["auto_means"][:,1:])/1e6 for x in results]
    counts=[np.nansum(x["auto_counts"][:,1:]) for x in results]
    N_SNPs=[x["auto_SNPs"] for x in results]
    ps=[x["p"] for x in results]
    print("\t".join(["ID1", "ID2", "N_SNP", "p", "mean", "count", "IBD1", "IBD2", "IBD", "IBD_filtered"]))
          
    print("\n".join(["\t".join(["%s", "%d", "%1.4f", "%1.2f", "%d","%1.4f", "%1.4f","%1.4f", "%1.4f"])%x 
                     for x in zip(pairs, N_SNPs, ps, means, counts, IBD1s, IBD2s, proportions, proportions_filtered)]))    

################################################################################
#Map
  
interpret_map={"simple_auto":simple_autosomes}

