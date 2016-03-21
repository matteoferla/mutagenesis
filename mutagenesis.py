#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Written for python 3

from Bio.Seq import Seq
from Bio.Alphabet import NucleotideAlphabet
from Bio._py3k import basestring  # I am not overly sure what this does compared to str. I just copied.
from Bio.Data import CodonTable
import re
import math
import warnings
import itertools as it

__author__ = 'Matteo'
__version__ = "$Revision$"
# $Source$

N = "\n"
T = "\t"
# N="<br/>"



__doc__ = '''
Classes:
* mutation
* MutationFormatError
* MutationDNASeq
* mutationSpectrum
This is a partial rewrite of mutanalyst js code. As a result a lot of attribute names are in camelCase, following JS style as opposed to PEP8.
'''


class Mutation:
    __doc__ = '''
    Accepts as arguments:
    * a mutation string
    * (opt) Seq object
    * (opt) forceDNA flag (def. False), if DNA is used but with protein notation
    * (opt) coding flag (def. True) to use the ref sequence as coding.
    It has the following groups of attributes:
    * from_nuc, to_nuc, num_nuc: nucleotide from, to and number.
    * from_aa, to_aa, num_aa: protein from, to and number
    * from_codon, to_codon: codon from and to
    * type: synonymous, non-synonymous and nonsense, and frameshift
    * is_substitution: true if a substitution
    It does not check whether the nucleotides are legitimate.

    Has also the method apply which returns a string where the mutation is applied to the Seq object (unchanged).
    '''
    codon_codex = {
        'ATG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'YTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'AGG', 'M': 'ATG'},
        'ACC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'AGC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'TAA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TAA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'TTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'AGA', 'M': 'ATG'},
        'CAG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'CTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'CGG', 'M': 'ATG'},
        'CGC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'AGC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'TGC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'TCC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TGA', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'TCT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'TCT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'GAT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'AGT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'AAC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'AGC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'TAC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'TCC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'CCA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TRA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'CTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'CGA', 'M': 'ATG'},
        'TGA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TGA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'TTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'MGA', 'M': 'ATG'},
        'CAT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'AGT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'TTG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'TTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'AGG', 'M': 'ATG'},
        'GGA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TGA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'TTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'MGA', 'M': 'ATG'},
        'TTA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TRA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'TTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'AGA', 'M': 'ATG'},
        'ACG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'TTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'AGG', 'M': 'ATG'},
        'GTG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'YTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'AGG', 'M': 'ATG'},
        'AGC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'AGC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'AGR', 'M': 'ATG'},
        'GCG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'TTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'AGG', 'M': 'ATG'},
        'GTC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'AGC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'GAC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'AGC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'TGG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TGA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'TTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'MGG', 'M': 'ATG'},
        'AGG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'AGY', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'TTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'AGG', 'M': 'ATG'},
        'GGT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'AGT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'CCC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'TCC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'AGA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'AGY', 'Y': 'TAY',
                'V': 'GTA', '*': 'TGA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'TTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'AGA', 'M': 'ATG'},
        'GTT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'AGT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'TCG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'TTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'AGG', 'M': 'ATG'},
        'CAA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TAA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'CTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'CGA', 'M': 'ATG'},
        'CCG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'CTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'CGG', 'M': 'ATG'},
        'ACT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'AGT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'GGG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'TTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'MGG', 'M': 'ATG'},
        'ATT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'AGT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'GGC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'AGC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'GCC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'TCC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'GTA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TRA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'YTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'AGA', 'M': 'ATG'},
        'ACA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TRA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'TTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'AGA', 'M': 'ATG'},
        'CGG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'CTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'CGG', 'M': 'ATG'},
        'CCT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'TCT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'GAA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TAA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'TTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'AGA', 'M': 'ATG'},
        'ATA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TRA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'YTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'AGA', 'M': 'ATG'},
        'AAG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'TTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'AGG', 'M': 'ATG'},
        'TTC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'TCC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'TTR',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'TCA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TRA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'TTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'AGA', 'M': 'ATG'},
        'CTC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'AGC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'TAG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'TTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'AGG', 'M': 'ATG'},
        'CGA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TGA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'CTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'CGA', 'M': 'ATG'},
        'TTT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'TCT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'TTR',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'TGT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'TCT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TGA', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'TAT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'TCT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'TCC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'TCC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'AGT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'AGT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'AGR', 'M': 'ATG'},
        'CGT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'AGT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'AAA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TAA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'TTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'AGA', 'M': 'ATG'},
        'ATC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'AGC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'},
        'CTG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'CTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'CGG', 'M': 'ATG'},
        'CTA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TRA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'CTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'CGA', 'M': 'ATG'},
        'GCA': {'N': 'AAY', 'P': 'CCA', 'W': 'TGG', 'K': 'AAA', 'H': 'CAY', 'Q': 'CAA', 'S': 'TCA', 'Y': 'TAY',
                'V': 'GTA', '*': 'TRA', 'F': 'TTY', 'C': 'TGY', 'A': 'GCA', 'E': 'GAA', 'G': 'GGA', 'L': 'TTA',
                'D': 'GAY', 'I': 'ATA', 'T': 'ACA', 'R': 'AGA', 'M': 'ATG'},
        'GAG': {'N': 'AAY', 'P': 'CCG', 'W': 'TGG', 'K': 'AAG', 'H': 'CAY', 'Q': 'CAG', 'S': 'TCG', 'Y': 'TAY',
                'V': 'GTG', '*': 'TAG', 'F': 'TTY', 'C': 'TGY', 'A': 'GCG', 'E': 'GAG', 'G': 'GGG', 'L': 'TTG',
                'D': 'GAY', 'I': 'ATH', 'T': 'ACG', 'R': 'AGG', 'M': 'ATG'},
        'CTT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'AGT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'AAT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'AGT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'GCT': {'N': 'AAT', 'P': 'CCT', 'W': 'TGG', 'K': 'AAR', 'H': 'CAT', 'Q': 'CAR', 'S': 'TCT', 'Y': 'TAT',
                'V': 'GTT', '*': 'TAR', 'F': 'TTT', 'C': 'TGT', 'A': 'GCT', 'E': 'GAR', 'G': 'GGT', 'L': 'CTT',
                'D': 'GAT', 'I': 'ATT', 'T': 'ACT', 'R': 'CGT', 'M': 'ATG'},
        'CAC': {'N': 'AAC', 'P': 'CCC', 'W': 'TGG', 'K': 'AAR', 'H': 'CAC', 'Q': 'CAR', 'S': 'AGC', 'Y': 'TAC',
                'V': 'GTC', '*': 'TAR', 'F': 'TTC', 'C': 'TGC', 'A': 'GCC', 'E': 'GAR', 'G': 'GGC', 'L': 'CTC',
                'D': 'GAC', 'I': 'ATC', 'T': 'ACC', 'R': 'CGC', 'M': 'ATG'}}

    def __init__(self, mutation, seq=None, forceDNA=False, coding=True):
        # TODO frameshift.
        # regarding frameshifts and co. there are lots of notations (http://www.hgmd.cf.ac.uk/docs/mut_nom.html seems helpful).dels are marked with 76-78delACT or 76_78del 83^84insTG 76_77insT
        # I'll implement one first.
        # TODO check how unicode in code is handled when not on my machine... delta and omega would be cool.
        # TODO seq should be a weak reference.
        self.from_aa = None
        self.to_aa = None
        self.num_aa = None
        self.from_codon = None
        self.to_codon = None
        self.from_nuc = None
        self.to_nuc = None
        self.num_nuc = None
        self.is_substitution=False
        self.type = "ERROR"
        mutation = mutation.replace("_","-") #not implemented yet
        mutation = mutation.replace("del","\u0394") # \u0394 is uppercase delta
        rexprotsub = re.match("([A-Z])(\d+)([A-Z])", mutation)  # A23T
        rexnuclsub = re.match("(\d+)([A-Z])\>([A-Z])", mutation)  # 234A>T
        rexprotdel = re.match("([A-Z])(\d+)\u0394", mutation)  # A23del
        rexnucldel = re.match("(\d+)\u0394([A-Z]?)", mutation)  # 234delA
        rexprotmanydel = re.match("([A-Z])(\d+)\-([A-Z])(\d+)\u0394", mutation)  # A23-D24del
        rexnuclmanydel = re.match("(\d+)\-(\d+)\u0394([A-Z]+?)", mutation)  # 234-235delAT
        # deal with forceDNA flag
        if forceDNA:  # a hack...
            if rexprotsub:
                mutation = str(rexprotsub.group(2)) + str(rexprotsub.group(1)) + ">" + str(rexprotsub.group(3))
                rexnuclsub = re.match("(\d+)(\w)\>(\w)", mutation)
            elif rexprotdel:
                mutation = str(rexprotsub.group(2)) + "\u0394" + str(rexprotsub.group(1))
                rexnucldel = re.match("(\d+)\u0394(\w?)", mutation)  # 234delA
            elif mutation.find(">") != -1:  # 234A>T
                print('forceDNA flag called even if DNA mutation given')  # TODO warning
            else:
                MutationFormatError()
        # NUCLEOTIDE
        if rexnuclsub:
            self.is_substitution = True
            self.from_nuc = rexnuclsub.group(2)
            self.to_nuc = rexnuclsub.group(3)
            self.num_nuc = int(rexnuclsub.group(1))
            if seq:
                assert seq[self.num_nuc - 1] == self.from_nuc, str(self.num_nuc) + " is " + seq[
                    self.num_nuc - 1] + ", not " + self.from_nuc
            if seq and coding:
                translation = seq.translate()._data
                r = math.floor(self.num_nuc / 3)
                self.num_aa = r + 1
                self.from_codon = seq[r * 3:r * 3 + 3]._data
                self.to_codon = seq[r * 3:self.num_nuc - 1]._data + self.to_nuc + seq[self.num_nuc:r * 3 + 3]._data
                self.from_aa = translation[r]
                self.to_aa = Seq(self.to_codon).translate()._data
                if self.from_aa == self.to_aa:
                    self.type = "synonymous"
                elif self.to_aa == "*":
                    self.type = "nonsense"
                else:
                    self.type = "non-synonymous"
        elif rexnucldel: #rexnucldel = re.match("(\d+)\u0394(\w?)", mutation)  # 234delA
            self.from_nuc = rexnucldel.group(2)
            self.to_nuc = ''
            self.num_nuc = int(rexnucldel.group(1))
            if seq:
                if self.from_nuc:
                    assert seq[self.num_nuc - 1] == self.from_nuc, str(self.num_nuc) + " is " + seq[
                        self.num_nuc - 1] + ", not " + self.from_nuc
                else:
                    self.from_nuc = seq[self.num_nuc - 1]
            if seq and coding:
                translation = seq.translate()._data
                r = math.floor(self.num_nuc / 3)
                self.num_aa = r + 1
                self.from_codon = seq[r * 3:r * 3 + 3]._data
                self.to_codon = seq[r * 3:self.num_nuc - 1]._data + self.to_nuc + seq[self.num_nuc:r * 3 + 3]._data
                self.from_aa = translation[r]
                self.to_aa = Seq(self.to_codon).translate()._data  #TODO check if it is a frameshift
                self.type = "deletion"
        # PROTEIN
        elif rexprotsub:
            self.is_substitution = True
            self.from_aa = rexprotsub.group(1)
            self.to_aa = rexprotsub.group(3).replace("X","*")
            self.num_aa = int(rexprotsub.group(2))
            if self.to_aa == self.from_aa:
                self.type = "synonymous"  # no questions asked.
            elif self.to_aa == "*":
                self.type = "nonsense"
            else:
                self.type = "non-synonymous"
            if seq and coding:
                assert seq.translate()[self.num_aa - 1] == self.from_aa, str(self.num_aa) + " is " + seq.translate()[
                    self.num_aa - 1] + ", not " + self.from_aa
                self.from_codon = seq._data[(self.num_aa - 1) * 3: (self.num_aa - 1) * 3 + 3]
                self.to_codon = self.codon_codex[self.from_codon][self.to_aa]
                if self.from_aa == self.to_aa: #avoid raising errors...
                    self.from_nuc = self.from_codon[0]
                    self.to_nuc = self.from_nuc
                    self.num_nuc = self.num_aa * 3
                #crap. what if there are two or three mutations to make an aa change?
                diff=[i for i in range(3) if self.to_codon[i] != self.from_codon[i]]
                self.from_nuc = self.from_codon[diff[0]:diff[-1] + 1]
                self.to_nuc = self.to_codon[diff[0]:diff[-1] + 1]
                self.num_nuc = self.num_aa * 3 - 2 + diff[0]
        else:
            raise MutationFormatError(str(mutation))
        #TODO handle other cases

    def apply(self, seq):
        return seq[0:self.num_nuc - 1]._data + self.to_nuc + seq[self.num_nuc + len(self.from_nuc) - 1:]._data

    def __str__(self):
        text = str(self.num_nuc) + self.from_nuc + ">" + self.to_nuc
        if self.num_aa:
            text += " (" + self.type + ": " + self.from_aa + str(self.num_aa) + self.to_aa + ")"
        return text

    def shortform(self):
        return self.from_nuc+">"+self.to_nuc



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

It has the following arguments:
* The mutations list contains Mutation objects.
* alphabet, always NucleotideAlphabet()
* _data', a str accessible via str()
* wt, a str without the mutation
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
            mutations = mutations.replace(","," ").split()
        elif isinstance(mutations, Mutation):
            mutations=[mutations]
        elif not isinstance(mutations, list):
            raise MutationFormatError()
        for mutation in mutations:
            if not isinstance(mutation, str):
                raise MutationFormatError()
        # parse mutations
        for mutation in mutations:
            if isinstance(mutation, str):
                mut = Mutation(mutation, self, forceDNA)
            elif isinstance(mutations, Mutation):
                mut = mutation
            else:
                raise MutationFormatError()
            self.mutations.append(mut)
            self._data = mut.apply(self)
        return self


class MutationTable:
    __doc__ = '''ATGC^2 table. The values are accessed with a A>C notation.
    Due to the fact that an argument name must be a valid variable name "A>C" cannot be given as MutationTable(A>C=1), but has to be given as MutationTable({A>C: 1})
    To access a frequency, use instance["A>C"] notation'''
    _bases={"A":0,"T":1,"G":2,"C":3}

    def __init__(self,frequencies=None):
        #each inner list has all the changes of a single from_base
        self._data = [ #A T G C
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]
        ]
        if frequencies:
            for item in frequencies:
                if item.find(">") != -1:
                    self[item]=frequencies[item]

    def _parse_input(self,item):
        item = item.upper().replace("U","T")
        if item.find(">"):
            (frombase,tobase)=item.split(">")
        elif item.find("<"):
            (tobase,frombase)=item.split("<")
        else:
            raise ValueError('Only N>N or N<N forms are accepted.')
        return (self._bases[frombase],self._bases[tobase])

    def normalize(self, A=0.25, T=0.25, G=0.25, C=0.25): #returns a copy. I was not sure if to do it in place...
        #TODO make it so that a zero base freq does not cause a div by zero.
        freqs={"A": A, "T":T, "G":G, "C":C} #seems a bit circular...
        norm1 = {bfrom+">"+bto: self[bfrom+">"+bto]/freqs[bfrom] for bto in self._bases for bfrom in self._bases}
        norm2 = {d:norm1[d]/sum(norm1.values()) for d in norm1}
        return MutationTable(norm2)

    def __getitem__(self, item): #A>C or C<A... #TODO accept degenerate bases...
        (from_enum,to_enum)=MutationTable._parse_input(self,item)
        return self._data[from_enum][to_enum]

    def __setitem__(self, item, value):
        (from_enum,to_enum)=MutationTable._parse_input(self,item)
        self._data[from_enum][to_enum]  = value

    def to_dict(self):
        return {from_base+">"+to_base:self[from_base+">"+to_base] for from_base in self._bases for to_base in self._bases}

    def __str__(self):
        return str(self.to_dict())  #TODO Fix in future

class MutationSpectrum:  # is this needed for Pedel?
    __doc__ = '''Returns the mutational spectrum, an object with
    the mutation frequency,
    the base freq
    the SE of the mut freq
    the mutational rate


    class method from sequences
    init from values

    A lot of this will be plagiarised from JS https://github.com/matteoferla/mutant_calculator/blob/master/mutationalBias.js
    '''
    types = ['TsOverTv', 'W2SOverS2W', 'W2N', 'S2N', 'W2S', 'S2W', 'ΣTs', 'Ts1', 'Ts2', 'ΣTv', 'TvW', 'TvN1', 'TvS',
             'TvN2']
    _bases={"A":0,"T":1,"G":2,"C":3}

    # ways=  #itertools... permutation?

    def __init__(self, mutants):
        # user gives a bunch of mutant sequences
        self.source="inputted from mutants"
        # Mutational load
        self.freqMean=0
        self.freqVar=0
        self.freqList=0
        # mutational spectrum
        self.raw_table=MutationTable()
        self.seq=None
        self.base_count=None
        self.base_frequency={"A":0.25,"T":0.25,"G":0.25,"C":0.25}
        MutationSpectrum._process_mutations_from_mutants(self,mutants)
        MutationSpectrum._calculate_base_frequency(self)
        self.table=self.raw_table.normalize(**self.base_frequency)
        #mutation frequency

    def __getitem__(self, item):
        return self.table[item]

    def __setitem__(self, item, value):
        self.table[item] = value

    def __str__(self):
        return str(self.table)

    @classmethod
    def from_values(cls,
                   source="loaded",
                   sequence="",
                   mutations=[], #TODO fix this dangerous entry
                   freqMean=0,
                   freqVar=0,
                   freqList=0,
                   raw_table=MutationTable(),
                   table=MutationTable()  # todo check what this is called in mutanalyst.js
                   ):
        # todo this method dones not check if the values are correct.
        # each individually as I should assert what they are
        mut = cls.__new__(cls)
        mut.source = source
        mut.sequence = sequence
        mut.mutations = mutations
        mut.freqMean = freqMean
        mut.freqVar = freqVar
        mut.freqList = freqList
        return mut

    def add_mutations(self, mutations):
        for variant in mutations:
            if type(variant) is str:
                variant=Mutation(variant,forceDNA=True, seq=self.seq)
            assert type(variant) is Mutation, str(variant) + " is not a instance of Mutation as expected. Consider MutationSpectrum.from_seqs() or MutationSpectrum.from_values()."
            if variant.is_substitution:
                self.raw_table[variant.shortform()] +=1

    def _process_mutations_from_mutants(self, mutants):
        self.mutations=[]
        for variant in mutants:
            assert type(variant) is MutationDNASeq, str(variant) + " is not a instance of MutationDNASeq as expected."
            self.mutations.extend(variant.mutations)
            if not self.seq and variant.wt:
                self.seq = variant.wt
            elif variant.wt:
                assert self.seq == variant.wt," Mutants appear to not be variants of the same wt sequence"
        MutationSpectrum.add_mutations(self,self.mutations)

    def _calculate_base_frequency(self):
        if self.seq:
            self.base_count={b: self.seq.count(b) for b in self._bases}
            self.base_frequency={b: self.seq.count(b)/sum(self.base_count.values()) for b in self._bases}
            #sum(self.base_count.values()) and not len(self.seq) becuase the former loses the Ns Xs and weirdos that somehow were smuggled in.
        else:
            warnings.warn("Sequence not given, default (equal base frequency) used")


    @classmethod
    def from_mutation_list(cls, mutations, seq =None):  # class method
        self = MutationSpectrum([])
        self.seq=seq
        self.source = "inputted from mutations"
        self.add_mutations(mutations)
        return self


def mincodondist(codon,
                 aa):  # there must be a more elegant way, but this will do. Serine and stop are the problematic ones.
    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
    ft = standard_table.forward_table
    ft["TAA"] = "*"
    ft["TGA"] = "*"
    ft["TAG"] = "*"
    degerated = {"N": {"A", "T", "G", "C"},
                 "W": {"A", "T"},
                 "S": {"G", "C"},
                 "R": {"A", "G"},
                 "Y": {"T", "C"},
                 "M": {"A", "C"},
                 "K": {"G", "T"},
                 "B": {"T", "G", "C"},
                 "D": {"A", "T", "G"},
                 "H": {"A", "T", "C"},
                 "V": {"A", "G", "C"},
                 "A": {"A"},
                 "T": {"T"},
                 "G": {"G"},
                 "C": {"C"}
                 }
    bases = "A T G C".split()
    aas = set("A R N D C Q E G H I L K M F P S T W Y V *".split())
    hits = set()
    # one mutation away
    for x in reversed(range(3)):
        for b in bases:
            candidate = codon[0:x] + b + codon[x + 1:3]
            if b == codon[x] and ft[candidate] == aa:  # preference is given to mutation that does not change position.
                return candidate
            elif b != codon[x] and ft[candidate] == aa:
                hits.add(b)
        if hits:
            for d in degerated:
                if not hits.symmetric_difference(degerated[d]):
                    return codon[0:x] + d + codon[x + 1:3]
    # two mutations away, it will give a degenerate base only for second position, because the code works that way.
    for x1 in reversed(range(2)):
        for x2 in reversed(range(2)):
            if x1 >= x2:
                continue
            else:
                for b1 in bases:
                    for b2 in bases:
                        candidate = codon[0:x1] + b1 + codon[x1 + 1:x2] + b2 + codon[x2 + 1:3]
                        if ft[candidate] == aa:
                            hits.add(b2)
                    if hits:
                        for d in degerated:
                            if not hits.symmetric_difference(degerated[d]):
                                return codon[0:x1] + b1 + codon[x1 + 1:x2] + d + codon[x2 + 1:3]
    return {"A": "GCN", "R": "CGN", "N": "AAY", "D": "GAY", "C": "TGY", "Q": "CAR", "E": "GAR", "G": "GGN", "H": "CAY",
            "I": "ATH", "M": "ATG", "L": "CTN", "K": "AAR", "F": "TTY", "P": "CCN", "S": "TCN", "T": "ACN", "W": "TGG",
            "Y": "TAY", "V": "GTN", "*": "TAR"}[aa]


def generateCodonCodex():
    __doc__ = "To find the mutation from a codon to encode a different AA a pregenerate dictionary is needed. This is here for reference."
    aas = set("A R N D C Q E G H I L K M F P S T W Y V *".split())
    return {"".join(codon): {aa: mincodondist("".join(codon), aa) for aa in aas} for codon in
            it.product("ATGC", repeat=3)}

def test():
    seq = "ATGTTGGGGAATTTTGGGGAACCC"
    # print("Generate a mutationDNASeq instance: ", seq)
    m = "2T>G"
    print("Mutating " +seq +" " + m + " ", MutationDNASeq(seq).mutate(m))
    # mutationSpectrum()
    # print(mincodondist("ATG", "I"))  #ATH is correct answer
    #print(generateCodonCodex())  # ACG is correct answer
    # print("Test complete")
    print(MutationSpectrum([MutationDNASeq(seq).mutate(m)]))


if __name__ == "__main__":
    test()
