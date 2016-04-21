#Functions for interpreting the distribution of chunks. 

from __future__ import division, print_function

def simple_autosomes(results):
    """
    Just count the proportion that is not IBD0
    """
    proportions=[sum(x["auto_state_proportions"][1:]) for x in results]
    pairs=[x["pair"][0]+"\t"+x["pair"][1] for x in results]
    N_SNPs=[x["auto_SNPs"] for x in results]
    ps=[x["p"] for x in results]
   
    print("\n".join(["%s\t%d\t%1.4f\t%1.4f"%x for x in zip(pairs, N_SNPs, ps, proportions)]))    

################################################################################
#Map
  
interpret_map={"simple_auto":simple_autosomes}

