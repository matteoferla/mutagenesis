#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Written for python 3

__author__ = 'Matteo'
__doc__ = ''''''

N = "\n"
T = "\t"
# N="<br/>"

from Bio.Seq import Seq
from Bio.Alphabet import NucleotideAlphabet
from Bio._py3k import basestring
import re, math

class Mutation():
    __doc__ = "arg: mutation string. seq if protein coding."

    def __init__(self, mutation, seq=None,forceDNA=False, coding=True):
        if forceDNA and mutation.find(">"):  # a hack...
            rex = re.match("(\D)(\d+)(\D)",mutation)
            if rex:
                mutation=str(rex.group(2))+str(rex.group(1))+">"+str(rex.group(3))
            else:
                MutationFormatError()
        rex=re.match("(\d+)(\w).(\w)",mutation)  #234A>T
        if rex:
            self.nfrom=rex.group(2)
            self.nto=rex.group(3)
            self.nnum=rex.group(1)
        elif re.match("(\D)(\d+)(\D)",mutation):
            raise Exception("mutations on DNA form protein is not yet implemented")
        else:
            raise MutationFormatError()
        if seq:
            pass
            ####HERE!
        if seq and coding:
            translation=seq.translate()
            r=math.floor()
            translation[]
            self.pfrom=None
            self.pto=None
            self.pnum=None
            self.codon=None
        else:
            self.pfrom=None
            self.pto=None
            self.pnum=None
            self.codon=None

    def __str__(self):
        return str(self.nnum)+self.nfrom+">"+self.nto

class MutationFormatError(Exception):
    message='''Error in the parsing a mutation notation.
    A mutation should be written as 123A>T for nucleotide or W45T for protein.
    If the method accepts multiple mutations, they should be separated with a space or as a list.
    '''
    def __init__(self, value=None):
        self.value=value
    def __str__(self):
        reply=self.message
        if self.value:
            reply += "Error raised due to "+str(self.value)
        return reply

class MutationDNASeq(Seq):
    __doc__ = '''A variant of the seq class, but with the method mutate which counts from one.
    method mutate(): mutates the sequence based on the  mutations, expressed as a string with spaces or list of strings.
        Different customs for nucleotide and protein are used to tell them apart:
        * 234A>T for DNA
        * A12F for protein, unless forceDNA is true.
        What if someone passes a Bio.Seq object at __init__?
        '''
    def __init__(self,data):  #For now only DNA.
        __doc__="This is just a copy of the Bio.Seq.Seq.__init__ method with the difference that it can only be nucleotide"
        if not isinstance(data, basestring):
            raise TypeError("The sequence data given to a Seq object should be a string (not another Seq object etc)")
        self.wt = data
        self._data = data
        self.alphabet = NucleotideAlphabet()  # Can only be nucleotide...
        self.mutations=[]
    def mutate(self, mutations,forceDNA=False):
        #reorganise and check mutations into a list of string
        if isinstance(mutations,str):
            mutations=mutations.split()
        if not isinstance(mutations, list):
            raise MutationFormatError()
        for mutation in mutations:
            if not isinstance(mutation, str):
                raise MutationFormatError()
        #parse mutations
        for mutation in mutations:
            if mutation.find(">") != -1: #DNA
                self.mutations.append(Mutation(mutation, self))
            else: #protein
                self.mutations.append(Mutation(mutation))

            rex=re.match("(\d+)(\w).(\w)",mutation)  #234A>T
            rexp = re.match("(\D)(\d+)(\D)",mutation)
            if rex:
                mutdex={"from": rex.group(2),"to": rex.group(3),"base":rex.group(1)}
            elif rexp:
                mutdex={"from": rex.group(1),"to": rex.group(3),"residue":rex.group(2)}
            else:
                raise MutationFormatError()

class mutationSpectrum():
    __doc__= '''Returns the mutational spectrum, an object with
    the mutation frequency,
    the base freq
    the SE of the mut freq
    the mutational rate

    '''
    def __init__(self, mutationList):
        raise Exception('CODE NOT WRITTEN')
    def __getitem__(self,direction): #as in spectrum["A>C"]? freq("A>C") better?
        raise Exception('CODE NOT WRITTEN')

def test():
    seq="ATTTTGGGGAATTTTGGGGA"
    print("Generate a mutationDNASeq instance: ", seq)
    print("slice 1:5", MutationDNASeq(seq))
    m="2A>T"
    print("Mutating "+str(m)+" ",MutationDNASeq(seq).mutate(m))
    print("Test complete")


if __name__ == "__main__":
    test()