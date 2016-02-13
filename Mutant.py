__author__ = 'Matteo'
__doc__ = ''''''

N = "\n"
T = "\t"
# N="<br/>"

from Bio.Seq import Seq
from Bio.Alphabet import NucleotideAlphabet
from Bio._py3k import basestring
import re
class mutationDNASeq(Seq):
    __doc__ = '''A variant of the seq class, but with the function mutate which counts from one'''
    def __init__(self,data):  #For now only DNA.
        __doc__="This is just a copy of the Bio.Seq.Seq.__init__ method with the difference that it can only be nucleotide"
        if not isinstance(data, basestring):
            raise TypeError("The sequence data given to a Seq object should be a string (not another Seq object etc)")
        self._data = data
        self.alphabet = NucleotideAlphabet()  # Can only be nucleotide...
    def mutate(self, mutations,forceDNA=False):
        __doc__='''mutates the sequence based on the mutations list of strings (or string with spaces).
        Different customs for nucleotide and protein are used to tell them apart:
        * 234A>T for DNA
        * A12F for protein, unless forceDNA is true.
        '''
        if isinstance(mutations,str):
            mutations=mutations.split()
        if not isinstance(mutations, list):
            raise TypeError("The mutations given should be a list of strings or a string with spaces")
        for mutation in mutations:
            if not isinstance(mutation, str):
                raise TypeError("The mutation given should be a list of strings or a string with spaces")
            mutdex={}
            rex=re.match("(\d+)(\w).(\w)",mutation)  #234A>T
            rexp = re.match("(\D)(\d+)(\D)",mutation)
            if rex:
                mutdex={"from": rex.group(2),"to": rex.group(3),"base":rex.group(1)}
            elif rexp and forceDNA:
                mutdex={"from": rex.group(1),"to": rex.group(3),"base":rex.group(2)}
            elif rexp:
                mutdex={"from": rex.group(1),"to": rex.group(3),"residue":rex.group(2)}
            HERE!!!




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

class test():
    pass


if __name__ == "__main__":
    print(mutationDNASeq("ATTTTGGGGA")[1:5])