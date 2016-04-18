#!/Users/tstevens/venv/bin/python
from Bio.Restriction.Restriction_Dictionary import rest_dict
from Bio import SeqIO
from copy import copy
from itertools import combinations,permutations
import numpy as np
import os
import re
import sys

MIN_FRAG_LEN = 300
MIN_FRAG_RATIO = 0.3
MAX_FRAG_NUM = 4
MAX_ENZ = 2
FULL_ENZ_LIST = 'AatII,AscI,AsiSI,AvrII,BamHI,BbsI,BbvCI,BsaI,BsmBI,BtgZI,EcoRI,EcoRV,HindIII,KpnI,MluI,NcoI,NotI,PacI,PmeI,SacI,SacII,SalI,SbfI,SmaI,SpeI,SwaI,XbaI,XhoI,XmaI'.split(',')
# CutSmart, 37C Only, Use HF when available
ENZ_LIST = 'AatII,AscI,AsiSI,AvrII,BamHI,BbvCI,BsaI,BtgZI,EcoRI,EcoRV,HindIII,KpnI,NcoI,NotI,PacI,PmeI,SacI,SacII,SbfI,SmaI,SpeI,XbaI,XhoI,XmaI'.split(',')

class FragTable:
    '''Digestion fragment table for a plasmid sequence'''
    def __init__(self,fname):
        self.name = os.path.basename(fname)
        record = SeqIO.read(fname,'gb')
        self.seq = str(record.seq)
        self.length = len(self.seq)
        self.enz_table = dict(zip(range(1,MAX_FRAG_NUM+1),[[] for n in range(MAX_FRAG_NUM)]))
        self.build_table()
    def digest(self,enz_comb):
        '''Returns array of fragment lengths resulting from digestion'''
        indices = np.array([])
        for enz in enz_comb:
            site = rest_dict[enz]['site']
            indices = np.hstack([indices,[match.start() for match in re.finditer(site,self.seq)]])
        if len(indices)>0:
            return circular_diff(np.sort(indices),self.length)
        else:
            return indices
    def build_table(self):
        '''Populate MAX_FRAG_NUM entries of enzyme table'''
        for n_enz in range(1,MAX_ENZ+1):
            for enz_comb in combinations(ENZ_LIST, n_enz):
                frag_len_arr = self.digest(enz_comb)
                entry = len(frag_len_arr)
                if (0<entry<=MAX_FRAG_NUM):
                    if is_diff(frag_len_arr):
                        # Append as frozenset to allow sets of frozensets for permuted_intersection
                        self.enz_table[entry].append(frozenset(enz_comb))

def circular_diff(indices,length):
    '''Returns fragment lengths as difference between indices on circular string'''
    return np.diff(np.hstack([indices,[length+indices[0]]]))

def is_diff(frag_len_arr):
    '''Checks if all fragments are longer than MIN_FRAG_LEN and each band is +/- MIN_FRAG_RATIO'''
    if len(frag_len_arr)==1:
        return True
    elif np.any(frag_len_arr<MIN_FRAG_LEN):
        return False
    else:
        frag_len_arr.dtype = float
        frag_len_arr = np.sort(frag_len_arr)
        frag_len_diff = np.diff(np.hstack([[MIN_FRAG_LEN*MIN_FRAG_RATIO], frag_len_arr]))
        return np.all(frag_len_diff/frag_len_arr>=MIN_FRAG_RATIO)

def permuted_intersection(fragtable_list):
    '''Maps fragment numbers onto plasmids, returns enzyme intersection (if any).
    Ensures distinct fragment numbers for each plasmid using a common enzyme set.
    Perform all permutations to output fewest number of enzymes.'''
    candidates = []
    for perm in permutations(range(1,MAX_FRAG_NUM+1),len(fragtable_list)):
        # E.g. first permutation is (1,2,3)->(0,1,2)
        for f,p in enumerate(perm):
            if f==0:
                enz_set = set(fragtable_list[f].enz_table[p])
            else:
                enz_set.intersection_update(set(fragtable_list[f].enz_table[p]))
            if len(enz_set)==0:
                # No intersection, yield next permutation
                break
        if len(enz_set)>0:
            candidates.append(sorted(enz_set,key=len,reverse=True).pop())
    if len(candidates)>0:
        return sorted(candidates,key=len,reverse=True).pop()
    else:
        return None

def main(fname_list):
    '''USAGE: ./dxdigest_fragnum.py RMp1.ape RMp2.ape ...'''
    # Construct a FragTable for each input plasmid
    fragtable_list = []
    for fname in fname_list:
        fragtable_list.append(FragTable(fname))
    # Intersect tables and Digest plasmid with final set
    enz_comb = permuted_intersection(fragtable_list)
    if enz_comb is None:
        print('No {} enzymes distinguish '.format(MAX_ENZ)+' '.join([f.name for f in fragtable_list]))
    else:
        for fragtable in fragtable_list:
            print('\t'.join([fragtable.name, ' '.join(enz_comb),
                ' '.join(map(lambda(x):str(int(x)), np.sort(fragtable.digest(enz_comb))))]))

if __name__=='__main__':
    main(sys.argv[1:])
