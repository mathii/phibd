#Read a .ind file and return a list of all pairs
#that are in the same population. Removes _[0-9]d_rel_xxx$
#tags which are our standard label for identified relatives. 

from __future__ import division, print_function
from collections import defaultdict
import fileinput, re, itertools

inds_by_pop=defaultdict(list)
strip_rel=re.compile("_[0-9]d_rel.*$")

for line in fileinput.input():
    bits=line.split()
    inds_by_pop[strip_rel.sub("", bits[2])].append(bits[0])
    
for group in inds_by_pop.values():
    if(len(group)>1):
        print("\n".join(["\t".join(x) for x in itertools.combinations(group, 2)]))