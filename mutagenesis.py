#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Written for python 3

import itertools as it
import numpy as np
# import scipy as sp
from scipy import misc, optimize, special
import math
import re
import warnings

from Bio.Alphabet import NucleotideAlphabet
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio._py3k import basestring  # I am not overly sure what this does compared to str. I just copied.

__author__ = 'Matteo'
__version__ = "N/A"

# This sloppiness will be refactored out once finished.
N = "\n"
T = "\t"
# N="<br/>"
"""Classes:
* mutation
* MutationFormatError
* MutationDNASeq
* mutationSpectrum
This is a partial rewrite of mutanalyst js code. As a result a lot of attribute names are in camelCase, following JS style as opposed to PEP8.
"""


class Mutation:
    """Accepts as arguments:
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
    """
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
        self.is_substitution = False
        self.type = "ERROR"
        mutation = mutation.replace("_", "-")  # not implemented yet
        mutation = mutation.replace("del", "\u0394")  # \u0394 is uppercase delta
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
                warnings.warn('forceDNA flag called even if DNA mutation given')
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
                r = math.floor((self.num_nuc - 1) / 3)
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
        elif rexnucldel:  # rexnucldel = re.match("(\d+)\u0394(\w?)", mutation)  # 234delA
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
                r = math.floor((self.num_nuc - 1) / 3)
                self.num_aa = r + 1
                self.from_codon = seq[r * 3:r * 3 + 3]._data
                self.to_codon = seq[r * 3:self.num_nuc - 1]._data + self.to_nuc + seq[self.num_nuc:r * 3 + 3]._data
                self.from_aa = translation[r]
                self.to_aa = Seq(self.to_codon).translate()._data  # TODO check if it is a frameshift
                self.type = "deletion"
        # PROTEIN
        elif rexprotsub:
            self.is_substitution = True
            self.from_aa = rexprotsub.group(1)
            self.to_aa = rexprotsub.group(3).replace("X", "*")
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
                if self.from_aa == self.to_aa:  # avoid raising errors...
                    self.from_nuc = self.from_codon[0]
                    self.to_nuc = self.from_nuc
                    self.num_nuc = self.num_aa * 3
                # crap. what if there are two or three mutations to make an aa change?
                diff = [i for i in range(3) if self.to_codon[i] != self.from_codon[i]]
                self.from_nuc = self.from_codon[diff[0]:diff[-1] + 1]
                self.to_nuc = self.to_codon[diff[0]:diff[-1] + 1]
                self.num_nuc = self.num_aa * 3 - 2 + diff[0]
        else:
            raise MutationFormatError(str(mutation))
            # TODO handle other cases

    def apply(self, seq):
        return seq[0:self.num_nuc - 1]._data + self.to_nuc + seq[self.num_nuc + len(self.from_nuc) - 1:]._data

    def __str__(self):
        text = str(self.num_nuc) + self.from_nuc + ">" + self.to_nuc
        if self.num_aa:
            text += " (" + self.type + ": " + self.from_aa + str(self.num_aa) + self.to_aa + ")"
        return text

    def shortform(self):
        return self.from_nuc + ">" + self.to_nuc


class MutationFormatError(Exception):
    """Error in the parsing a mutation notation.
    A mutation should be written as 123A>T for nucleotide or W45T for protein.
    If the method accepts multiple mutations, they should be separated with a space or as a list.
    """

    def __init__(self, value=None):
        self.value = value

    def __str__(self):
        reply = self.message
        if self.value:
            reply += "Error raised due to " + str(self.value)
        return reply


class MutationDNASeq(Seq):
    """A variant of the seq class, but with the method mutate and the attribute mutations.
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
    """

    def __init__(self, data):  # TODO For now only DNA. In future translate and replace AA.
        """This is just a copy of the Bio.Seq.Seq.__init__ method with the difference that
        * it can only be nucleotide
        * data can be string or Seq"""
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
        """
        Adds mutations to a seq. See Mutation class for format of mutation.
        :param mutations: comma or space separated string, list of Mutation instances or single Mutation istance.
        :param forceDNA: boolean seen with Mutation class.
        :return: Self
        """
        # reorganise and check mutations into a list of string
        if isinstance(mutations, str):
            mutations = mutations.replace(",", " ").split()
        elif isinstance(mutations, Mutation):
            mutations = [mutations]
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

    def variant(self, mutations, forceDNA=False):
        """
        Returns a copy of the sequence (wt only, no preëxisting mutations) mutated with the given mutations.
        :param mutations: same as mutate method.
        :param forceDNA:
        :return: new instance
        """
        return MutationDNASeq(self.wt).mutate(mutations, forceDNA)

    def variants(self, mutations, forceDNA=False):
        if isinstance(mutations,str):
            mutations=list(filter(len,mutations.replace(";","\n").replace("\r","\n").split("\n")))
        pass
    #HERE


    def to_regular_seq(self):
        """
        Returns a Seq() object, where the sequence is that with the mutation added (`_data`), not the wild type (`wt`).
        :return: Seq instance.
        """
        return Seq(self._data, alphabet=NucleotideAlphabet())


class MutationTable:
    """ATGC^2 table. The values are accessed with a A>C notation.
    Due to the fact that an argument name must be a valid variable name "A>C" cannot be given as MutationTable(A>C=1), but has to be given as MutationTable({A>C: 1})
    To access a frequency, use instance["A>C"] notation"""
    _bases = {"A": 0, "T": 1, "G": 2, "C": 3}

    def __init__(self, frequencies=None):
        # each inner list has all the changes of a single from_base
        self._data = [  # A T G C
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0]
        ]
        if frequencies:
            for item in frequencies:
                if item.find(">") != -1:
                    self[item] = frequencies[item]

    def _parse_input(self, item):
        item = item.upper().replace("U", "T")
        if item.find(">") == 1:
            (frombase, tobase) = item.split(">")
        elif item.find("<") == 1:
            (tobase, frombase) = item.split("<")
        else:
            raise ValueError('Only N>N or N<N forms are accepted.')
        return (self._bases[frombase], self._bases[tobase])

    def normalize(self, A=0.25, T=0.25, G=0.25, C=0.25):  # returns a copy. I was not sure if to do it in place...
        # TODO make it so that a zero base freq does not cause a div by zero.
        freqs = {"A": A, "T": T, "G": G, "C": C}  # seems a bit circular...
        norm1 = {bfrom + ">" + bto: self[bfrom + ">" + bto] / freqs[bfrom] for bto in self._bases for bfrom in
                 self._bases}
        norm2 = {d: norm1[d] / sum(norm1.values()) for d in norm1}
        return MutationTable(norm2)

    def __getitem__(self, item):  # A>C or C<A... #TODO accept degenerate bases...
        (from_enum, to_enum) = self._parse_input(item)
        return self._data[from_enum][to_enum]

    def __setitem__(self, item, value):
        (from_enum, to_enum) = self._parse_input(item)
        self._data[from_enum][to_enum] = value

    def to_dict(self):
        return {from_base + ">" + to_base: self[from_base + ">" + to_base] for from_base in self._bases for to_base in
                self._bases}

    def __str__(self):
        return str(self.to_dict())  # TODO Fix in future. Ordereddictionary

    def __iter__(self):
        for o in self._bases:
            for i in self._bases:
                yield o + ">" + i

    @staticmethod
    def ibase():
        for o in "A T G C".split():
            for i in "A T G C".split():
                yield (o, i)


class MutationLoad:
    """Calculates the mutational load of a pool of mutants.
    It assumes it is from a simple amplification reaction, however a sum of two ind. Poisson distributions is a Poisson, so if two samples were pooled the lamdba would be the same.


                    freqMean=0,
                    freqVar=0,
                    freqList=0,
    """

    def __init__(self, mutants, pcr_efficiency=None, pcr_cycles=32):
        self.mutation_tally = []
        self.seq = None
        for variant in mutants:
            assert type(variant) is MutationDNASeq, str(variant) + " is not a instance of MutationDNASeq as expected."
            self.mutation_tally.append(len(variant.mutations))
            if not self.seq and variant.wt:
                self.seq = variant.wt
            elif variant.wt:
                assert self.seq == variant.wt, " Mutants appear to not be variants of the same wt sequence"
        tally = np.array(self.mutation_tally)  # make np earlier still?
        self.mutation_frequency = np.bincount(tally)  # overkill?
        self.mean = NumSEM(tally.mean(), tally.std() / np.sqrt(len(tally)))
        self.arange = np.arange(0, tally.max() + 1)
        parameters, cov_matrix = optimize.curve_fit(poisson, self.arange, self.mutation_frequency)
        self.lamb = NumSEM(parameters[0], np.sqrt(np.diag(cov_matrix)))
        if pcr_efficiency:
            self.pcr_efficiency = pcr_efficiency
            self.pcr_cycles = pcr_cycles
            pcr = pcr_distribution_factory(self.pcr_efficiency, self.pcr_cycles)
            parameters, cov_matrix = optimize.curve_fit(pcr, self.arange, self.mutation_frequency)
            self.pcr = NumSEM(parameters[0], np.sqrt(np.diag(cov_matrix)))


def pcr_distribution_factory(efficiency, cycles):
    def pcr_distribution(k, mu):
        # fnx copied from Drummond 2005
        # except that the author is a casinista when it comes to naming variables.
        # there <m_{nt}> means both Poisson-lambda (here mu) and m = Poisson-k (here k)
        # k is used as the summation index, here i.
        # p = (1 + h)^{-n}\sum_{i=0}^{n}\binom{n}{i}h^k\frac{(ix)^m e^{-ik}}{k!}
        # x = (mu (1 + h))/(n h)
        x = mu * (1 + efficiency) / (cycles * efficiency)
        return (1 + efficiency) ** (-cycles) * np.sum(
            [special.binom(cycles, i) * efficiency ** i * ((i * x) ** k) * np.exp(-i * x) / misc.factorial(k) for i in
             np.arange(0, cycles + 1)])

    return pcr_distribution


def poisson(k, lamb):
    return (lamb ** k / misc.factorial(k)) * np.exp(-lamb)


class MutationSpectrum:  # is this needed for Pedel?
    """Returns the mutational spectrum of a pool of mutants, an object that can be accessed with the "N>N" or "N<N" index notation,
    __E.g.__ `spectrum_instance["A>T"]`, like MutationTable, but has also derived keys such as `spectrum_instance["S>N"]`.
    I has the following attributes:

    * `raw_table`, A MutationTable (int) with the raw values
    * `seq`, the wt sequences of the mutants (has to be common)
    * `base_count`, the base count of seq
    * `base_freq`, the base frequency
    * `table`, the normalised table (MutationTable of floats)
    * `avg_table`, the normalised table averaged based on strand equivalent mutations
    * `var_table`, the variance of the above
    * `se_table`, the standard error
    * `_data`, where the dictionary (not MutationTable) of NumSEM values is stored.


    class method from sequences
    init from values

    A lot of this will be plagiarised from JS https://github.com/matteoferla/mutant_calculator/blob/master/mutationalBias.js
    """
    _bases = {"A": 0, "T": 1, "G": 2, "C": 3}

    # ways=  #itertools... permutation?

    def __init__(self, mutants):
        # user gives a bunch of mutant sequences
        self.source = "inputted from mutants"
        # mutational spectrum
        self.raw_table = MutationTable()
        self.seq = None
        self.base_count = None
        self.base_frequency = {"A": 0.25, "T": 0.25, "G": 0.25, "C": 0.25}
        self._process_mutations_from_mutants(mutants)
        self._calculate_base_frequency()
        self.refresh()

    def refresh(self):
        self.table = self.raw_table.normalize(**self.base_frequency)
        self._calculate_se_table()
        self._calculate_advanced()

    @staticmethod
    def _sanify_input(item):
        item = item.upper().replace("U", "T")
        subdex = {"TvW": "W>W", "Ts1": "W>S:Ts", "Ts2": "S>W:Ts", "W>S:Tv": "TvN1", "TvS": "S>S",
                  "S>W:Tv": "TvN2"}  # substitution dictionary
        if item in subdex:
            item = subdex[item]
        item.replace("sum", "Σ")
        item.replace("∑", "Σ")  # greek vs. maths. Unicode is hard to embrace, ae?
        item.replace("Transversions", "Tv")
        item.replace("Transions", "Ts")
        if item.find("<") == 1:
            if item.find(":") == -1:
                (tobase, frombase) = item.split("<")
                item = frombase + ">" + tobase
            else:
                (direction, qual) = item.split(":")
                (tobase, frombase) = direction.split("<")
                item = frombase + ">" + tobase + ":" + qual
        return item

    def __getitem__(self, item):
        item = MutationSpectrum._sanify_input(item)
        return self._data[item]

    def __setitem__(self, item, value):
        item = MutationSpectrum._sanify_input(item)
        self._data[item] = value

    def __str__(self):
        return str("\n".join([str(x) + ":\t" + str(y) for x, y in self._data.items()]))

    @classmethod
    def from_values(cls,
                    source="loaded",
                    sequence="",
                    mutations=None,
                    raw_table=MutationTable(),
                    table=MutationTable()  # todo check what this is called in mutanalyst.js
                    ):
        # todo this method dones not check if the values are correct.
        # each individually as I should assert what they are
        mut = cls.__new__(cls)
        mut.source = source
        mut.sequence = sequence
        if mutations:
            mut.mutations = mutations
        else:
            mut.mutations = []
        raise Exception('Incomplete')
        return mut

    def add_mutations(self, mutations):
        for variant in mutations:
            if type(variant) is str:
                variant = Mutation(variant, forceDNA=True, seq=self.seq)
            assert type(variant) is Mutation, str(
                variant) + " is not a instance of Mutation as expected. Consider MutationSpectrum.from_seqs() or MutationSpectrum.from_values()."
            if variant.is_substitution:
                self.raw_table[variant.shortform()] += 1

    def _process_mutations_from_mutants(self, mutants):
        self.mutations = []
        for variant in mutants:
            assert type(variant) is MutationDNASeq, str(variant) + " is not a instance of MutationDNASeq as expected."
            self.mutations.extend(variant.mutations)
            if not self.seq and variant.wt:
                self.seq = variant.wt
            elif variant.wt:
                assert self.seq == variant.wt, " Mutants appear to not be variants of the same wt sequence"
        self.add_mutations(self.mutations)

    def _calculate_base_frequency(self):
        if self.seq:
            self.base_count = {b: self.seq.count(b) for b in self._bases}
            self.base_frequency = {b: self.seq.count(b) / sum(self.base_count.values()) for b in self._bases}
            # sum(self.base_count.values()) and not len(self.seq) becuase the former loses the Ns Xs and weirdos that somehow were smuggled in.
        else:
            warnings.warn("Sequence not given, default (equal base frequency) used")

    def _calculate_se_table(self):
        """
        Calculates the average and se. Fills the `_data` attribute (a dictionary) with the basic spectrum as NumSEM instances.
        Despite `_data` not being a MutationTable, its members can be accessed by indexing of the instance with the "N>N" or "N<N" notation.

        >>> spectro = MutationSpectrum([mutant1, mutant2])
        >>> print(spectro["A>T"])
        0.4375±0.21875

        The derived fields (_e.g._ "W>N") are calculated in _calculate_advanced.
        :return:
        """
        comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
        self.avg_table = MutationTable(
            {b1 + ">" + b2: (self.table[b1 + ">" + b2] + self.table[comp[b1] + ">" + comp[b2]]) / 2 for (b1, b2) in
             MutationTable.ibase()})
        # JS mutball.compTable[h+x][g+y]=muts[h+x][g+y]/2 +muts[h+!x][g+!y]/2;//it actually switches the diagonal. But it's zero.
        self.var_table = MutationTable({b1 + ">" + b2: sum(
            [(self.avg_table[x] - self.table[x]) ** 2 for x in [b1 + ">" + b2, comp[b1] + ">" + comp[b2]]]) / 2 for
                                        (b1, b2) in MutationTable.ibase()})
        self.se_table = MutationTable({d: math.sqrt(self.var_table[d] / 2) for d in self.var_table})  # n, not n-1?
        a = self.avg_table.to_dict()
        s = self.se_table.to_dict()
        self._data = {k: NumSEM(a[k], s[k]) for k in a}

    def _calculate_advanced(self):  # UNFINISHED
        """
        Calculates the derived fields:
        * W>S
        * W>W
        * W>N
        * S>S
        * S>W
        * S>N
        * W>S/S>W
        * W>S:Ts
        * S>W:Ts
        * W>S:Tv
        * S>W:Tv
        * ΣTs
        * ΣTv
        * ΣTs/ΣTv
        :return:
        """
        self["W>S"] = self["A>G"] + self["T>G"] + self["A>C"] + self["T>C"]
        self["W>W"] = self["A>T"] + self["T>A"]  # TvW in JS
        self["S>W"] = self["G>A"] + self["G>T"] + self["C>A"] + self["C>T"]
        self["S>S"] = self["G>C"] + self["C>G"]
        self["W>N"] = self["W>S"] + self["S>S"]
        self["S>N"] = self["S>W"] + self["W>W"]
        self["W>S/S>W"] = self["W>S"] / self["S>W"]
        self["W>S:Ts"] = self["A>G"] + self["T>C"]  # Ts1 in JS
        self["S>W:Ts"] = self["G>A"] + self["C>T"]  # Ts2 in JS
        self["ΣTs"] = self["W>S:Ts"] + self["S>W:Ts"]
        self["W>S:Tv"] = self["A>C"] + self["T>G"]
        self["S>W:Tv"] = self["G>T"] + self["C>A"]
        self["ΣTv"] = self["W>S:Tv"] + self["S>W:Tv"] + self["S>S"] + self["W>W"]
        self["ΣTs/ΣTv"] = self["ΣTs"] / self["ΣTv"]

    @classmethod
    def from_mutation_list(cls, mutations, seq=None):  # class method
        self = MutationSpectrum([])
        self.seq = seq
        self.source = "inputted from mutations"
        self.add_mutations(mutations)
        return self


class NumSEM:
    """
    A class to handle numbers with SEM.
    I am not sure why there is nothing that does this and whether this is the best way of doing.
    For now errors are propagated parametrically based on the maths I discuss in [mutanalyst](http://www.mutanalyst.com/)
    """
    option_space_between_numbers = True

    def __init__(self, num, sem, df=2):
        """
        :param num: the number (mean)
        :param sem: the standard error
        :param df: the number of samples used to determine the SE
        :return: an object with the three inputs as _num, _sem and _df attributes.
        """
        self._num = float(num)
        self._sem = float(sem)
        self._df = int(df)

    def __str__(self):
        """
        The lowest digit to show is the first digit of the error. So math.floor(math.log10(se))
        :return:
        """
        if self.option_space_between_numbers:
            s = " "
        else:
            s = ""
        if self._sem == 0:

            return str(self._num) + s + "±" + s + "Ind."
        else:
            sig = -math.floor(math.log10(self._sem))
            if sig < 0:
                sig = 0
            txt = "{:." + str(sig) + "f}" + s + "±" + s + "{:." + str(sig) + "f}"
            return txt.format(self._num, self._sem)

    def __add__(self, other):
        """
        Returns the addition of either two NumSEM objects or a NumSEM and a float/int
        :param other: NumSEM or int or float
        :return: a new NumSEM instance where the variance is based on the Binaymé rule if both NumSEM.
        """
        if type(other) is NumSEM:
            v = self._var() + other._var()
            df = NumSEM._df(self._df, other._df)
            return NumSEM(self._num + other._num, math.sqrt(v / df), df)
        else:  # assume int or float
            return NumSEM(self._num + other, self._sem, self._df)

    def __sub__(self, other):
        """
        Returns the subtraction of either two NumSEM objects or a NumSEM and a float/int
        :param other: NumSEM or int or float
        :return: a new NumSEM instance where the variance is based on the Binaymé rule (var summed) if both NumSEM.
        """
        if type(other) is NumSEM:
            v = self._var() + other._var()
            df = NumSEM._df(self._df, other._df)
            return NumSEM(self._num - other._num, math.sqrt(v / df), df)
        else:  # assume int or float
            return NumSEM(self._num - other, self._sem, self._df)

    def __mul__(self, other):
        """
        Returns the mutiplication of either two NumSEM objects or a NumSEM and a float/int
        :param other: NumSEM or int or float
        :return: a new NumSEM instance where, if both NumSEM, the variance is the var(x)/mean(x)^2 + var(y)/mean(y)^2.
        The latter formula stems from the first order Taylor approximation of the König–Huygens theorem applied to a function:
        > Var(f(x))\approx [f'(x)]^2 · Var(x)
        And converting Var(xy) to Var(e^ln(xy)) and solving.
        > Var(e^ln(xy)) = (e^ln(xy))^2 · Var(ln(xy)) = x^2 · y^2 · (Var(ln(x))+Var(ln(y))) = x^2 · y^2 · (Var(x)/x^2+Var(y)/y^2) _etc._
        """
        if type(other) is NumSEM:
            v = self._var() * other._num ** 2 + other._var() * self._num ** 2
            df = NumSEM._df(self._df, other._df)
            return NumSEM(self._num * other._num, math.sqrt(v / df), df)

    def __truediv__(self, other):
        """
        Same principle as multiplication
        :param other: NumSEM or int or float
        :return: a new NumSEM instance where, if both NumSEM, the variance is the var(x)/mean(y)^2 + var(y)*mean(x)^2/mean(y)^4.
        """
        if type(other) is NumSEM:
            v = self._var() / other._num ** 2 + other._var() * self._num ** 2 / other._num ** 4
            df = NumSEM._df(self._df, other._df)
            return NumSEM(self._num / other._num, math.sqrt(v / df), df)

    def _var(self):
        return self._sem ** 2 * self._df

    @staticmethod
    def _df(a, b):
        return min(a,
                   b)  # I am unsure if min is best, hence the method. I assume that the worst case scenario is the smallest.
        # Also it is not degrees of freedom but sample size...


def mincodondist(codon, aa):
    # there must be a more elegant way, but this will do. Serine and stop are the problematic ones.
    """
    This method is called solely by generateCodonCodex, which is here just to show how the codoncodex dictionary was made.
    :param codon: a string of three characters corresponding to the starting codon
    :param aa: a string of a letter corresponsponding to the target AA
    :return: the degenerate codon needed to reach that AA from that codon.
    """
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
    """To find the mutation from a codon to encode a different AA a pregenerate dictionary is needed.
    This is here for reference."""
    aas = set("A R N D C Q E G H I L K M F P S T W Y V *".split())
    return {"".join(codon): {aa: mincodondist("".join(codon), aa) for aa in aas} for codon in
            it.product("ATGC", repeat=3)}


##### Wayne's stuff

def glue_completeness(vsize, lsize):
    """
    Find the library completeness given:
    :param vsize: Number of possible variants
    :param lsize: Size of library
    :return: completeness as [0,1] number
    """
    return 1 - math.exp(-lsize / vsize)


def glue_library(vsize, completeness):
    """
    Find the required library size given:
    :param vsize: Number of possible variants
    :param completeness: completeness as [0,1] number
    :return: library size required
    """
    if completeness > 1:  # percentage...
        completeness /= 100
    return -vsize * math.log(1 - completeness)


def glue_probability(vsize, probability):
    """
    Find the req library size given:
    :param vsize: Number of possible variants
    :param probability: probability the library is 100% complete
    :return: library size required
    """
    return -vsize * math.log(1 - math.exp(math.log(probability) / vsize))


def glue(vsize, lsize=None, completeness=None, probability=None):  # TODO UNFINISHED
    """
    Wrapper for the glue functions. Depending on what is given it will figure it out.
    :param vsize: Number of possible variants
    :param lsize: library size required
    :param completeness: completeness as [0,1] number
    :param probability: probability the library is 100% complete
    :return: dictionary
    """
    if not vsize:
        raise Exception("Why would you not know that?")
    if not lsize:
        if not completeness and not probability:
            pass
    raise Exception('NOT FINISHED. I need to think when certain features would be useful')


def pedel(lsize, seq_len, mps, dist_fx=poisson):
    """
    For poisson distribution use poisson (default).
    For pcr distribution, use first pcr_distribution_factory(efficiency, cycles) to obtain a function specific to those parameters.
    >>> pedel(1e6,len(wt),4,pcr_distribution_factory(0.4, 32))
    :param lsize: library size
    :param seq_len: sequence length (bases)
    :param mps: mutation per seq.
    :param dist_fx: function of the distribution
    :return: number of distinct sequences in library
    """

    def Lx(x):
        return poisson(x, mps) * lsize

    def Cx(x):
        """
        Number of unique variants in libray with x mutations
        :param x: n muts
        :return: C_x
        """
        if x == 0:
            return 1  # wt
        # eq 6.
        return Vx(x) * (1 - math.exp(-Lx(x) / Vx(x)))

    def Vx(x):
        """
        Number of possible variants with x muts
        :param x: number of mutations
        :return: V_x
        """
        if x == 0:
            return 1  # wt
        # eq. 5
        return special.binom(seq_len, x) * 3 ** x

    xu = (mps + math.log(0.1 / lsize)) / math.log(mps / (2 * seq_len))
    xl = (mps + math.log(3 / lsize)) / math.log(mps / (3 * seq_len))
    s1 = math.floor(xl)
    s2 = math.ceil(xu)
    # note that range(min,max) does not include max.
    C = 1 + \
        sum([Vx(x) for x in range(1, s1 + 1)]) + \
        sum([Cx(x) for x in range(s1, s2)]) + \
        lsize * (1 - sum([dist_fx(x, mps) for x in range(0, s2)]))
    return C


def driver(lsize, seq_len, cross, positions, observable=True):
    """
    DRIVeR.
    :param lsize: library size
    :param seq_len: seq length
    :param cross: mean number of crossovers
    :param positions: positions of variable bases
    :param observable: boolean for observable (True) or all (False) crossovers
    :return: number of possible seqs, expected number of distinct seqs, mean num of actual cross, mean num of obs cross.
    """
    pass


class LibraryStatistics():
    def __init__(self, lsize, mload, mspectrum):
        pass

    '''

    '''


def test():
    wt = MutationDNASeq(
        "ATGAACACAGACGACATTCTGTTTTCTTACGGAGAAGAAGACATTCCTTTGAAGGCGCTGTCGTTTCCCATCTTCGAAACGACGAATTTCTACTTCGACAGTTTCGACGAGATGTCGAAAGCCCTCAGAAACGGAGACTACGAATTCGTTTACAAAAGAGGAAGTAATCCCACAACGAGACTGGTGGAGAAGAAACTCGCAGCGCTTGAAGAGTGTGAAGATGCCCGCCTCGTTGCCTCTGGAATGAGCGCCATTTCGCTTTCCATCCTTCATTTCCTCAGCTCGGGAGACCACGTCGTGTGTGTGGACGAGGCTTACTCCTGGGCGAAAAAGTTCTTCAACTACCTTTCAAAGAAGTTCGATATAGAAGTCAGCTACGTTCCTCCCGACGCGGAAAGAATAGTCGAAGCCATCACGAAGAAGACGAAGCTCATCTACCTCGAAAGTCCCACGAGTATGAGAATGAAAGTGATCGATATAAGAAAGGTCACAGAAGCGGCAGGAGAACTCAAGATAAAAACCGTCATAGACAACACCTGGGCGTCGCCGATCTTTCAAAAACCAAAGCTTCTGGGAGTGGATGTGGTGGTCCACTCTGCGACGAAGTACATCTCAGGACACGGAGACGTGATGGCAGGAGTGATCGCAGGAGACGTCGAAGATATGAAGAACATCTTCGTGGATGAATACAAAAACATCGGACCGGTTCTCTCGCCCATAGAAGCCTGGCTCATCTTGAGAGGTCTTAGAACGCTGGAACTCCGTATGAAAAAGCACTACGAAAACGCTCTTGTGGTGTCTGACTTCCTCATGGATCACCCGAAGGTCCTCGAGGTGAACTACCCGATGAATCCAAGATCACCGCAGTACGAACTCGCTTCCTCTCAGATGAGCGGTGGCTCAGGACTGATGAGCTTCAGGCTGAAAACGGACAGCGCAGAGAAAGTCAAAGAGTTCGTCGAAAGTCTGAGGGTTTTCAGGATGGCTGTGAGCTGGGGAAGTCACGAGAACCTTGTTGTTCCAAGGGTGGCTTATGGAGACTGCCCGAAAAAAGACGTGAACCTGATAAGAATCCATGTGGGTCTCGGAGATCCAGAAAAGCTCGTGGAAGATCTGGATCAGGCACTCAAAAAGATTTAA")
    mutations='''G286A T306C A687T T880C
WT
WT
C372T A923G
G832A
A650C
A720C
C449A A556G A814G A853G C826T C1059T
G793A
A1048G
G726T
G982C
C53T G678T T798A
A224T A447G A586T C964T
A185T
T860A
A979T A1005T
C453T
WT
'''
    wt.variants(mutations,forceDNA=True)
    return None
    mutball = [MutationDNASeq(seq).mutate(m),
               MutationDNASeq(seq).mutate("3G>T"),
               MutationDNASeq(seq).mutate(["1A>G", "2T>C"])]
    spectro = MutationSpectrum(mutball)
    print(spectro)
    # load = MutationLoad(mutball, 0.6, 32)
    print("Test complete")


def test_wayne():
    def comp(n, a, b):
        print("{0}= found:{1} exp:{2}".format(n, a, b))

    comp("Glue completeness", glue_completeness(1e6, 3e6), 0.9502)
    comp("Glue library", glue_library(1e6, 0.95), 2.996e6)
    comp("Glue prob", glue_probability(1e6, 0.95), 1.679e7)
    comp("Pedel", pedel(1e7, 1000, 4, poisson), 8.872e6)


if __name__ == "__main__":
    test()
    # TODO fix PCR distribution
    # TODO pretty print MutationalSpectum
    # TODO add the following:
    # glue-it
    # codoncalculator
    # aa-calculator
    # Pedel
    # pedel-AA
    # driver
