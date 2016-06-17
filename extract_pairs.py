#Read a .ind file and return a list of all pairs
#that are in the same population. Removes _[0-9]d_rel_xxx$
#tags which are our standard label for identified relatives. 

from __future__ import division, print_function
from collections import defaultdict
import fileinput, re, itertools

inds_by_pop=defaultdict(list)
strips=[re.compile("_[0-9]d_rel.*$"),
re.compile("_brother.*$"),
re.compile("_sister.*$"),
re.compile("_mother.*$"),
re.compile("_father.*$"),
re.compile("_son.*$"),
re.compile("_daughter.*$"),
re.compile("-[0-9]$")]

for line in fileinput.input():
    bits=line.split()
    pop=bits[2]
    for strip in strips:
        pop=strip.sub("", pop)
    inds_by_pop[pop].append(bits[0])
    
for group in inds_by_pop.values():
    if(len(group)>1):
        print("\n".join(["\t".join(x)+"\t"+group for x in itertools.combinations(group, 2)]))
