#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Written for python 3, not tested under 2.
"""
"""
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = ""

N = "\n"
T = "\t"
# N = "<br/>
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils import MeltingTemp as mt
import math
from warnings import warn
import random

def deep_mutation_scan(region, section, target_temp=55, overlap_len=22, primer_range=(None,60), mutation='NNK'):
    """
    Designs primers for quikchange for deep mutation scanning.
    Based on the overlap principle of http://nar.oxfordjournals.org/content/43/2/e12.long
    In terms of calculating the melting temperature module is used.
    For now the calculations are based on everything after the NNK. The problem is the missing thermodynamics for the various mismatches.
    :param region: the sequence of a region including neighbouring parts and not only the seq of interest.
    :param section: the range of bases to make primers for as a list/tuple of two items or a slice
    :param target_temp: the temp threshold. Remember that Phusion HF buffer is secret, but they say the Tm is +3. For salt correction see below
    :param overlap_len: the length of the overlap of the two primers
    :param primer_range: the min and max len of the primer.
    :param mutation: str of the desired mutatated codon
    :return: a list of dictionaries with the following keys: base codon primer len_homology len_anneal len_primer homology_start homology_stop homology_Tm anneal_Tm

    Regarding salts. check out mt.salt_correction at http://biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-module.html#salt_correction
    if needed.
    """
    #
    # sanity
    if isinstance(region,Seq):
        pass
    elif isinstance(region,str):
        region=Seq(region, generic_dna)
    else:
        raise TypeError('Sequence is neither string or Seq.')
    region=region.ungap(' ') #spaces!
    if isinstance(section, slice):
        pass
    else:  #try if it is an iterable.
        section=slice(*section)
    gene=region[section]
    if len(gene)%3 !=0:
        warn('The length of the region is not a multiple of three: '+str(len(region)))
    #max size
    if not primer_range[0] or primer_range[0] < overlap_len:
        primer_range=(overlap_len,primer_range[1])
    # iterate across codons.
    geneball=[]
    for x in range(section.start,section.stop,3):
        start=x-int(overlap_len/2)
        stop=x+overlap_len-int(overlap_len/2)
        codon = {'codon': region[x:x+3],
                 'base':x,
                 'homology_start':start,
                 'homology_stop': stop,
                 'len_homology': overlap_len,
                 'homology_Tm': round(mt.Tm_NN(region[start:stop]),1)
                 }
        #iterate to find best fw primer.
        for dir in ('fw','rv'):
            for i in range(primer_range[0]-int(overlap_len/2)-3, primer_range[1]-int(overlap_len/2)-3):
                #the length of the annealing part of the primer is i+3 (the end of the mutated codon)
                #so the region prior to the mutation does not count: -int(overlap_len/2)-3

                # This cannot be done:
                #mut=region[start:x]+Seq('NNK')+region[x+3:x+i]
                #ori=region[start:x+i]
                #t= mt.Tm_NN(mut, c_seq=ori.complement())
                # ValueError: no thermodynamic data for neighbors 'GG/TA' available

                #this seems to pick the weakest
                # mut = region[start:x] + Seq('NNK') + region[x + 3:x + i]
                # t = mt.Tm_NN(mut)

                # so ignoring forepart
                if dir == 'fw':
                    mut=region[x+3:x+i]
                elif dir == 'rv':
                    mut = region[x - i:x].reverse_complement()
                else:
                    raise Exception

                t=mt.Tm_NN(mut)

                #check if the tms are good...
                if mut[-1].upper() in ['C','G'] and t>target_temp-1:
                    break
                elif t>target_temp+1:
                    break
            else:
                warn('Target temperature not met. {0}C > {1}C'.format(target_temp,t))

            '''
            # check if G:C nearby first.
            for j in (1,2,3):
                if mut[-j-1].upper() in ['C','G']:
                    mut2 = region[x + 3:x + i-j]
                    t2 = mt.Tm_NN(mut2)
                    if t2+ 1 > target_temp:
                        t=t2
                        mut=mut2
                        i=i-j
            '''
            if dir == 'fw':
                codon[dir+'_primer']=region[start:x].upper() + mutation.lower() + mut.upper()
            else: #dir == 'rv'
                codon[dir + '_primer'] = region[x+3:stop].reverse_complement().upper() + Seq(mutation).reverse_complement().lower() + mut.upper()
            codon[dir+'_len_primer']=len(codon[dir+'_primer'])
            codon[dir+'_anneal_Tm']=round(t,1)
            codon[dir+'_len_anneal']=i

        for i in range(primer_range[0] - int(overlap_len / 2) - 3, primer_range[1] - int(overlap_len / 2) - 3):

            t = mt.Tm_NN(mut)
            if t>target_temp:
                break

        geneball.append(codon)
    return geneball

def randomer(n):
    alphabet=['A','T','G','C']
    return ''.join([random.choice(alphabet) for x in range(n)])

def test():
    n = 30
    m=21
    query = randomer(n).lower() + randomer(m).upper() + randomer(n).lower()
    print('sequence',query)
    import csv
    w = csv.DictWriter(open('out.csv', 'w', newline=''),
                       fieldnames='base codon fw_primer rv_primer len_homology fw_len_anneal rv_len_anneal fw_len_primer rv_len_primer homology_start homology_stop homology_Tm fw_anneal_Tm rv_anneal_Tm'.split())
    w.writeheader()
    w.writerows(deep_mutation_scan(query, (n, n + m)))

if __name__ == "__main__":
    test()
