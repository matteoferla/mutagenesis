#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Written for python 3

__author__ = 'Matteo'
__doc__ = ''''''

N = "\n"
T = "\t"
# N="<br/>"

from Bio.Seq import Seq
from Bio.Alphabet import NucleotideAlphabet
from Bio._py3k import basestring
import re, math

__doc__ = '''
Classes:
* mutation
* MutationFormatError
* MutationDNASeq
* mutationSpectrum


This is a partial rewrite of mutanalyst js code. As a result a lot of attribute names are in camelCase, following JS style as opposed to PEP8.
'''


class Mutation():
    __doc__ = '''
    Accepts as arguments:
    * a mutation string
    * (opt) Seq object
    * (opt) forceDNA flag (def. False), if DNA is used but with protein notation
    * (opt) coding flag (def. True) to use the ref sequence as coding.
    It has the following groups of attributes:
    * nucFrom, nucTo, nucNum: nucleotide from, to and number.
    * protFrom, protTo, protNum: protein from, to and number
    * codonFrom, codonTo: codon from and to
    * type: synonymous, non-synonymous and nonsense.
    It does not check whether the nucleotides are legitimate.

    Has also the method apply which returns a string where the mutation is applied to the Seq object (unchanged).
    '''

    def __init__(self, mutation, seq=None, forceDNA=False, coding=True):
        # TODO seq should be a weak reference.
        self.protFrom = None
        self.protTo = None
        self.protNum = None
        self.codonFrom = None
        self.codonTo = None
        self.nucFrom = None
        self.nucTo = None
        self.nucNum = None
        self.type = "ERROR"
        # deal with forceDNA flag
        if forceDNA:  # a hack...
            rex = re.match("(\D)(\d+)(\D)", mutation)
            if rex:
                mutation = str(rex.group(2)) + str(rex.group(1)) + ">" + str(rex.group(3))
            elif mutation.find(">") != -1:  # 234A>T
                print('forceDNA flag called even if DNA mutation given') # TODO warning
            else:
                MutationFormatError()
        rex = re.match("(\d+)(\w)\>(\w)", mutation)  # 234A>T
        rexprot = re.match("(\D)(\d+)(\D)", mutation) # A23T
        # NUCLEOTIDE
        if rex:
            self.nucFrom = rex.group(2)
            self.nucTo = rex.group(3)
            self.nucNum = int(rex.group(1))
            if seq:
                assert seq[self.nucNum - 1] == self.nucFrom, self.nucNum + " is " + seq[
                    self.nucNum - 1] + ", not " + self.nucFrom
            if seq and coding:
                translation = seq.translate()._data
                r = math.floor(self.nucNum / 3)
                self.protNum = r + 1
                self.codonFrom = seq[r * 3:r * 3 + 3]._data
                self.codonTo = seq[r * 3:self.nucNum - 1]._data + self.nucTo + seq[self.nucNum:r * 3 + 3]._data
                self.protFrom = translation[r]
                self.protTo = Seq(self.codonTo).translate()._data
                if self.protFrom == self.protTo:
                    self.type = "synonymous"
                elif self.protTo == "*":
                    self.type = "nonsense"
                else:
                    self.type = "non-synonymous"
        # PROTEIN
        elif rexprot:
            self.protFrom = rexprot.group(1)
            self.protTo = rexprot.group(3)
            self.protNum = int(rexprot.group(2))
            if self.protTo == self.protFrom:
                self.type = "synonymous"   #no questions asked.
            elif self.protTo == "*":
                self.type = "nonsense"
            else:
                self.type = "non-synonymous"
            if seq and coding:
                assert seq.translate()[self.protNum - 1] == self.protFrom, self.protNum + " is " + seq.translate()[self.protNum - 1] + ", not " + self.protFrom
                self.codonFrom = seq._data[(self.protNum -1) *3 : (self.protNum -1) *3 +3]
                raise  Exception('UNFINISHED')
                #find minimum mutations needed to get that mutation from given codon.
        else:
            raise MutationFormatError()
        print(self.__dict__)

    def apply(self, seq):
        return seq[:self.nucNum - 1]._data + self.nucTo + seq[self.nucNum:]._data

    def __str__(self):
        text = str(self.nucNum) + self.nucFrom + ">" + self.nucTo
        if self.protNum:
            text += " (" + self.type + ": " + self.protFrom + str(self.protNum) + self.protTo + ")"
        return text


class MutationFormatError(Exception):
    message = '''Error in the parsing a mutation notation.
    A mutation should be written as 123A>T for nucleotide or W45T for protein.
    If the method accepts multiple mutations, they should be separated with a space or as a list.
    '''

    def __init__(self, value=None):
        self.value = value

    def __str__(self):
        reply = self.message
        if self.value:
            reply += "Error raised due to " + str(self.value)
        return reply


class MutationDNASeq(Seq):
    __doc__ = '''A variant of the seq class, but with the method mutate and the attribute mutations.
    Also accepts a Seq object in addition to a string for the sequence.

method mutate(): mutates the sequence based on the  mutations, expressed as a string with spaces or list of strings.
Different customs for nucleotide and protein are used to tell them apart:
* 234A>T for DNA
* A12F for protein, unless forceDNA is true. This is not yet implemented and will require special coding.

The mutations list contains mutation objects.
'''

    def __init__(self, data):
        # TODO For now only DNA. In future translate and replace AA.
        __doc__ = '''This is just a copy of the Bio.Seq.Seq.__init__ method with the difference that
        * it can only be nucleotide
        * data can be string or Seq'''
        if isinstance(data, Seq):
            self._data = data._data
            # TODO assert if nucleotide
        elif not isinstance(data, basestring):
            raise TypeError("The sequence data given to a Seq object should be a string (not another Seq object etc)")
        else:
            self._data = data
        self.alphabet = NucleotideAlphabet()  # Can only be nucleotide...
        self.wt = data
        self.mutations = []

    def mutate(self, mutations, forceDNA=False):
        # reorganise and check mutations into a list of string
        if isinstance(mutations, str):
            mutations = mutations.split()
        if not isinstance(mutations, list):
            raise MutationFormatError()
        for mutation in mutations:
            if not isinstance(mutation, str):
                raise MutationFormatError()
        # parse mutations
        for mutation in mutations:
            mut = Mutation(mutation, self, forceDNA)
            self.mutations.append(mut)
            self._data = mut.apply(self)
        return self

class mutationTable:
    __doc__ = " ATGC^2 table. Basically a 2d enum"
    def __init__(self):
        self._data=[
                       [0, 0, 0, 0],
                       [0, 0, 0, 0],
                       [0, 0, 0, 0],
                       [0, 0, 0, 0]
                   ]
        #GO TO HERE...
        # figure out how the table gets in there in the first place

class mutationSpectrum:  # is this needed for Pedel?
    __doc__ = '''Returns the mutational spectrum, an object with
    the mutation frequency,
    the base freq
    the SE of the mut freq
    the mutational rate


    class method from sequences
    init from values

    A lot of this will be plagiarised from JS https://github.com/matteoferla/mutant_calculator/blob/master/mutationalBias.js
    The JS is...
    //**************************************************
    //the mutball object.
    //Commit and read were initially written as methods of mutball, but were moved out in order to quarantine interactions with document to the document.
    function mutagen() {
        var mutball = {}
        for (b in bases) {
            mutball["sum" + bases[b]] = "25";
        }
        for (b in ways) {
            mutball[ways[b]] = "0";
        }
        return mutball;
    }

    time for a rethink.


    '''
    types = ['TsOverTv', 'W2SOverS2W', 'W2N', 'S2N', 'W2S', 'S2W', 'ΣTs', 'Ts1', 'Ts2', 'ΣTv', 'TvW', 'TvN1', 'TvS',
             'TvN2'];
    bases = "A T G C".split()

    # ways=  #itertools... permutation?

    @classmethod
    def fromValues(cls,
                   source="loaded",
                   sequence="",
                   baseList="",
                   freqMean=0,
                   freqVar=0,
                   freqList=0,
                   mutTable=[
                       [0, 0, 0, 0],
                       [0, 0, 0, 0],
                       [0, 0, 0, 0],
                       [0, 0, 0, 0]
                   ],
                   mutFreq=[  # todo check what this is called in mutanalyst.js
                       [0, 0, 0, 0],
                       [0, 0, 0, 0],
                       [0, 0, 0, 0],
                       [0, 0, 0, 0]
                   ]
                   ):
        # this method dones not check if the values are correct.
        mut = cls.__new__(cls)
        # each individually as I should assert what they are
        mut.source = source
        mut.sequence = sequence
        mut.baseList = baseList
        mut.freqMean = freqMean
        mut.freqVar = freqVar
        mut.freqList = freqList

    @classmethod  # change to __init__
    def fromSeqs(cls, mutationList):  # class method
        # user gives a bunch of mutant sequences
        raise Exception('CODE NOT WRITTEN')
        for variant in mutationList:
            assert type(variant) is MutationDNASeq, str(variant) + " is not a instance of MutationDNASeq as expected."

    def __init__(self, ):
        raise Exception('CODE NOT WRITTEN')

    def __getitem__(self, direction):  # as in spectrum["A>C"]? freq("A>C") better?
        raise Exception('CODE NOT WRITTEN')


def test():
    seq = "ATGTTGGGGAATTTTGGGGAA"
    print("Generate a mutationDNASeq instance: ", seq)
    m = "M1L"
    print("Mutating " + str(m) + " ", MutationDNASeq(seq).mutate(m))
    #mutationSpectrum()
    print("Test complete")


if __name__ == "__main__":
    test()
