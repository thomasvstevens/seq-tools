#!/Users/tvs/virtual-environments/venv2.7/bin/python

from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from random import random
import sys
import ipdb

RARE_FRAC = 0.1

class DiversifiedSeq:
    codon, codon_frac = {},{}
    def __init__(self,pep_str,codon_fname=None,len_pep_repeat=0,just_sample=False):
        self.pep = pep_str
        self.len_pep_repeat = len_pep_repeat
        if DiversifiedSeq.codon=={}:
            build_codon_back_table(codon_fname)
        # Initialize DP table and accumulator variables
        self.L = np.zeros([len(self.pep)*3]*2,dtype='int')
        self.lensub = 0
        self.substr = ''
        self.nt = ''
        # LCS
        if self.len_pep_repeat>0:
            self.back_trans_seed()
            self.fill_table(1,len(self.nt))
        if just_sample:
            self.weighted_sample()
            self.fill_table(1,len(self.nt))
        else:
            for p in range(self.len_pep_repeat,len(self.pep)):
               self.greedy_extend(p*3,self.pep[p])
    def back_trans_seed(self):
        # Make first repeat sequence by rotation through codon table
        repeat_pep = self.pep[:self.len_pep_repeat]
        repeat_nt = ''
        for aa in repeat_pep:
            p = self.codon[aa].pop()
            repeat_nt += p
            self.codon[aa].insert(0,p)
        self.nt = repeat_nt
    def fill_table(self,jl,ju):
        # Calculate DP entries between lower and upper bounds on i
        # For self-similarity, only need entries above diagonal (j>i)
        # Means max self-similarity is len-3
        # Diagonal would hold len of entire string
        # Work on a copy for greedy comparison
        for j in range(jl,ju):
            for i in range(j):
                if self.nt[i]==self.nt[j]:
                    if i==0:
                        self.L[i,j] = 1
                    else:
                        self.L[i,j] = self.L[i-1,j-1] + 1
                    if self.L[i,j]>self.lensub:
                        self.lensub = self.L[i,j]
                        # slice indices do not include endpoint_ext
                        self.substr = self.nt[j-self.lensub+1:j+1]
    def greedy_extend(self,n,aa):
        # Take the codon that does not extend lcsubstr
        # If all extend, take that which minimizes the lcsubstr
        clist = self.codon[aa]
        copies = []
        min_lensub = len(self.nt)
        min_idx = 0
        for idx, next_codon in enumerate(clist):
            sub_lcs = deepcopy(self)
            sub_lcs.nt = sub_lcs.nt + next_codon
            sub_lcs.fill_table(n,n+3)
            copies.append(sub_lcs)
            if sub_lcs.lensub < min_lensub:
                min_lensub = sub_lcs.lensub
                min_idx = idx
        # Populate with best
        self.nt = copies[min_idx].nt
        self.L = copies[min_idx].L
        self.lensub = copies[min_idx].lensub
        self.substr = copies[min_idx].substr
    def weighted_sample(self):
        for aa in self.pep:
            codon_index = next(index for index,value in enumerate(self.codon_frac[aa]) if random() < value)
            self.nt += self.codon[aa][codon_index]

def build_codon_back_table(codon_fname):
    codon_seq,codon_frac = {},{}
    df = pd.read_table(codon_fname,sep=' +',index_col=1,engine='python')
    for aa in df.index.unique():
        if df.loc[aa].ndim > 1:
            frac_sorted = df.loc[aa][df.loc[aa,'frac']>RARE_FRAC].sort('frac')
            # Exclude rare codons and sort list in Ascending fractional usage
            codon_seq[aa] = list(frac_sorted['codon'])
            cumulative = (frac_sorted['frac']/frac_sorted['frac'].sum()).cumsum()
            codon_frac[aa] = list(cumulative)
        else:
            codon_seq[aa] = [df.loc[aa,'codon']]
            codon_frac[aa] = [1.0]
    DiversifiedSeq.codon = codon_seq
    DiversifiedSeq.codon_frac = codon_frac

def main(pep_fa=None,
         codon_fname='codon_table.txt'):
    # Codon diversification for repetitive amino acid sequences
    # Modified longest common substring dynamic programming algorithm, O(n**2)
    # INPUT: fasta of peptide sequences, (optional) codon table flat file
    # OUTPUT: fasta of nucleotide sequences
    #########
    # If no input, run benchmark
    if pep_fa is None:
        repeat = 'GGG'
        #repeat = 'REPEAT'
        pep = repeat*20
        # Naive list rotation approach (use in predefined order)
        ds_rotate = DiversifiedSeq(pep,codon_fname,len_pep_repeat=len(pep))
        # Greedy algorithm (locally choose codon to minimize self-LCS)
        ds_greedy = DiversifiedSeq(pep,len_pep_repeat=len(repeat))
        # Draw samples to plot null distribution of self-LCS
        n_samples = 1000
        sampled_lensub = [DiversifiedSeq(pep,just_sample=True).lensub for i in range(n_samples)]
        plt.hist(sampled_lensub, bins=range(len(pep)), align='left', facecolor='c')
        y = np.linspace(*plt.ylim())
        plt.plot(np.tile(ds_rotate.lensub,len(y)),y,'k--',lw=2)
        plt.plot(np.tile(ds_greedy.lensub,len(y)),y,'r--',lw=2)
        plt.legend(['Naive (rotation)','Greedy (DP)'])
        plt.xlabel('Length of DNA repeat')
        plt.ylabel('Count')
        plt.xlim((0,40))
        plt.show()

if __name__=='__main__':
    main(*sys.argv[1:])
