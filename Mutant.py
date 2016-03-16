#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Written for python 3

from Bio.Seq import Seq
from Bio.Alphabet import NucleotideAlphabet
from Bio._py3k import basestring  # I am not overly sure what this does compared to str. I just copied.
from Bio.Data import CodonTable
import re, math
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
    * type: synonymous, non-synonymous and nonsense.
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
        self.type = "ERROR"
        # deal with forceDNA flag
        if forceDNA:  # a hack...
            rex = re.match("(\D)(\d+)(\D)", mutation)
            if rex:
                mutation = str(rex.group(2)) + str(rex.group(1)) + ">" + str(rex.group(3))
            elif mutation.find(">") != -1:  # 234A>T
                print('forceDNA flag called even if DNA mutation given')  # TODO warning
            else:
                MutationFormatError()
        rex = re.match("(\d+)(\w)\>(\w)", mutation)  # 234A>T
        rexprot = re.match("(\D)(\d+)(\D)", mutation)  # A23T
        # NUCLEOTIDE
        if rex:
            self.from_nuc = rex.group(2)
            self.to_nuc = rex.group(3)
            self.num_nuc = int(rex.group(1))
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
        # PROTEIN
        elif rexprot:
            self.from_aa = rexprot.group(1)
            self.to_aa = rexprot.group(3)
            self.num_aa = int(rexprot.group(2))
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
            raise MutationFormatError()
        print(self.__dict__)

    def apply(self, seq):
        return seq[0:self.num_nuc - 1]._data + self.to_nuc + seq[self.num_nuc + len(self.from_nuc) - 1:]._data

    def __str__(self):
        text = str(self.num_nuc) + self.from_nuc + ">" + self.to_nuc
        if self.num_aa:
            text += " (" + self.type + ": " + self.from_aa + str(self.num_aa) + self.to_aa + ")"
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


class MutationTable:
    __doc__ = '''ATGC^2 table. The values are accessed with a A>C notation.
    Due to the fact that an argument name must be a valid variable name "A>C" cannot be given as MutationTable(A>C=1), but has to be given as MutationTable({A>C: 1})
    To access a frequency, use instance["A>C"] notation'''

    def __init__(self,frequencies):
        self._data = [ #A T G C
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]
        ]
        for item in frequencies:
            if item.find(">") != -1:
                self[item]=frequencies[item]


    def __getitem__(self, item): #A>C  TODO allow backwards assignment alla R
        bases={"A":0,"T":1,"G":2,"C":3}
        (frombase,tobase)=item.upper().replace("U","T").split(">")
        return self._data[bases[frombase]][bases[tobase]]

    def __setitem__(self, item, value):
        bases={"A":0,"T":1,"G":2,"C":3}
        (frombase,tobase)=item.upper().replace("U","T").split(">")
        self._data[bases[frombase]][bases[tobase]] = value

class MutationSpectrum:  # is this needed for Pedel?
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
             'TvN2']
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
    seq = "ATGTTGGGGAATTTTGGGGAA"
    # print("Generate a mutationDNASeq instance: ", seq)
    m = "G3*"
    #print("Mutating " +seq +" " + m + " ", MutationDNASeq(seq).mutate(m))
    # mutationSpectrum()
    # print(mincodondist("ATG", "I"))  #ATH is correct answer
    #print(generateCodonCodex())  # ACG is correct answer
    # print("Test complete")
    print(MutationTable({"A>G": 2})["A>G"])


if __name__ == "__main__":
    test()
