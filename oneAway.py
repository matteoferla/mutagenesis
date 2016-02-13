#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Written for python 3
__doc__ = '''
Wayne asked the following:
1. Read the first codon (say it’s "GTC" — Val).
2. Determine all 1-base changes that alter the amino acid (TTC - Phe; CTC — Leu; ATC — Ile; GCC — Ala; GAC — Asp; GGC — Gly).
3. Feed this amino acid set into a souped-up version of AA-Calculator to determine the optimal degenerate codon, or codons, to encode those amino acids (HTC + GVC).
4. Write this out. Move to the next codon.
5. At the end, output the set of mutagenic oligos and also some statistics (how many total variants? How many total oligos?)

Thinking about it some more, the best way to do that is to calculate the optimal mutagenic codon pools for each of the 64 codons, and stash them in a database. Then the program could retrieve pre-calculated codon pools.

So regardless of sequence is probably quicker just making a oneAway codon replacement table and then applying it to the sequence.
The codons can be split into:
* self
* accessible with one change
* inaccessible with one change
Due to the layout, synonymous codons behave differently.
Bottom-up or top-down?
The permutation can be made into a dictionary of AA and sets of codons and the "intersections" are the best.
One option is to do in the negative: unaccessible AA have to be removed from NNN.
How does one equate A + G = P? Doe biopython have it?
There are 15 bases options. 15^3 is only 3375 options.
So get the 3375 codons and make a set of AA made.
Given a codon's options, filter out the ones that give the inaccessible ones.
For the remaining find ones where all AA are made.
In case of a draw, combination of two or more.



Folder: Wayne Marsden
'''
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "23 / 01 / 16"

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import itertools as it
from collections import defaultdict, Counter
import pickle


N = "\n"
T = "\t"
#N = "<br/>

def initialise(avoidStop=False, avoidSelf=False):
    __doc__ = "Returns a dictionary with key = codon and value = dict with keys = degenarate codon and value its AA"
    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
    ft=standard_table.forward_table
    ft["TAA"]="*"
    ft["TGA"]="*"
    ft["TAG"]="*"
    degerated={"N":{"A","T","G","C"},
                "W":{"A","T"},
                "S":{"G","C"},
                "R":{"A","G"},
                "Y":{"T","C"},
                "M":{"A","C"},
                "K":{"G","T"},
                "B":{"T","G","C"},
                "D":{"A","T","G"},
                "H":{"A","T","C"},
                "V":{"A","G","C"},
                "A":{"A"},
                "T":{"T"},
                "G":{"G"},
                "C":{"C"}
              }
    #This dict contains sets of AA for each degenerate codon
    degdex={"".join(codon): {ft["".join(x)] for x in it.product(degerated[codon[0]],degerated[codon[1]],degerated[codon[2]])} for codon in it.product(degerated.keys(), repeat=3)}
    print("Is this dict of degenerate codons: set of AA correct? ",degdex)
    #{'RNR': {'V', 'I', 'R', 'K', 'M', 'T', 'E', 'G', 'A'}, 'YSV': {'P', 'W', 'C', 'R', 'S', '*'}, 'YRM':  etc.

    #For each codon split AA into acc and inacc.
    bases = "A T G C".split()
    aas = set("A R N D C Q E G H I L K M F P S T W Y V".split())  #removing the asterisk makes stop codons okay. but it seems to make little difference,
    if avoidStop:
        aas.add("*")
    accdex={"".join(codon) : {ft[c+codon[1]+codon[2]] for c in bases} | {ft[codon[0]+c+codon[2]] for c in bases} | {ft[codon[0]+codon[1]+c] for c in bases} for codon in it.product("ATGC",repeat=3)}
    for codon in accdex:
        if "*" in accdex[codon]:
            accdex[codon].remove("*")
        if avoidSelf:
            accdex[codon].remove(ft[codon])
    print("Is this dict of codons: set of accessible AA correct?",accdex)
    inaccdex={codon: aas-s for codon,s in accdex.items()}
    print("Is this dict of codons: set of inaccessible AA correct?",inaccdex)
    #remove self-AA from accdex: we really do not care if self-AA is in or not. But we don't want it on the inaccdex either.
    for codon in accdex:
        if ft[codon] != "*":
            accdex[codon].remove(ft[codon])
    #for each codon filter out inacc degcodons
    posdegdex={codon: {dc for dc in degdex if degdex[dc].isdisjoint(inaccdex[codon])} for codon in accdex}
    print("Is this dict of codons: set of poss deg codon correct?",posdegdex)

    #for each codon find if any contain all the accdex
    gooddex=dict()
    for codon in posdegdex:
        gooddex[codon]={dc for dc in posdegdex[codon] if accdex[codon].issubset(degdex[dc])}         #Actually nothing
        for n in range(2,10):
            if gooddex[codon]==set():
                gooddex[codon]={"+".join(dcball) for dcball in it.combinations(posdegdex[codon],n) if accdex[codon].issubset(set().union(*[degdex[dc] for dc in dcball]))}
    print("Is this dict of codons: set of covered deg codon correct?",gooddex)
    codex={codon: {dc: degdex[dc] for dc in list(gooddex[codon])[0].split("+")} for codon in gooddex}   #first goes in. Better picking can be done f need be.
    print("The counts of codons needed are: ", Counter([len(codex[codon].keys()) for codon in codex]))

    pickle.dump(codex,open("codon2degcodon.p","wb"))
    return codex  #codex is a dictionary keyed to each codon with a value of a dictionary that has the degenerate codons and the AA they represent.

def getGeneSeq(genome, symbol):
    __doc__ = "Return a gene from the genome given a symbol. Returns a seqFeature with a qualifier CDS with the DNA sequence"
    gene={gene.qualifiers["gene"][0]:gene for gene in genome.features if gene.type =="CDS"}[symbol]
    extr=gene.extract(genome)
    print(gene.qualifiers["gene"][0][0].capitalize()+gene.qualifiers["gene"][0][1:]+" sequence extracted correctly: ", extr.seq.translate(to_stop=True) == gene.qualifiers["translation"][0])
    gene.qualifiers["CDS"]=[str(extr.seq)]
    return gene

def generateMutants(seq):
    __doc__ = "seq is a custom mod of seq feature"
    codex=initialise(True)
    #codex=pickle.load(open("codon2degcodon.p","rb"))

    if type(seq) is not str:
        text=seq.qualifiers["gene"][0]
        seq=seq.qualifiers["CDS"][0]
    else:
        text=""
    geneball=[SeqRecord(Seq(seq,IUPAC.ambiguous_dna),"Wild_type","_",text+" (no mutations)")]
    for i in range(3, len(seq)-3, 3): #start codon is left alone...
        for dc,aas in codex[seq[i:i+3]].items():
            geneball.append(SeqRecord(Seq(seq[0:i]+dc+seq[i+3:],IUPAC.ambiguous_dna),"ID","_",text+" "+seq[i:i+3]+str(i)+" ("+CodonTable.unambiguous_dna_by_name["Standard"].forward_table[seq[i:i+3]]+str(int(i/3))+") mutated to "+dc+" ("+", ".join(aas)+")"))
    for gi,gene in enumerate(geneball):
        if gi >0:
            gene.id="mutant_"+str(gi)
    print("Number of sequences: ",len(geneball)-1)
    return geneball

def analysis():
    codex=initialise(avoidStop=True, avoidSelf=True)


if __name__ == "__main__":
    genome=SeqIO.read("NC_000913.gbk","gb") #Eco genome
    #triose phos isomerase is tpiA
    tpiA=getGeneSeq(genome, "tpiA")
    gb=generateMutants(tpiA)
    SeqIO.write(gb,"mutantlist.fa","fasta")

