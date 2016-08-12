# Mutants!
Biopython is lacking when it comes to mutagenesis. Hopefully this module will help!
This module will be the basis of a revamped Pedel server, you can see the proto-pre-alpha
at [pedel2-matteoferla.rhcloud.com](pedel2-matteoferla.rhcloud.com).

> This is pre-release. You can use some features of it, but some might not work or might change in future.
> Also the code contains a lot of features that are not documented here. This is a just a MD version of the actual docstring, so best check there for an upto date details.

# Usage
This Python module allow an easy way to add mutations (either nucleic acid or protein) to a nucleic acid sequence.
It also performs several features present in the [mutanalyst site](http://www.mutanalyst.com)
namely the calculation of the mutational spectrum and mutational load of a set of mutated sequences using mathemagical
steps (_i.e._ just pretend the maths is magic) to get the most accurate values.

The module uses the "123A>T" style of annotating mutations at the nucleic acid level.

```
>>> seq =
>>> MutationDNASeq("ATGTTGGGGAATTTTGGGGAACCC").mutate("3G>T")
#spectro = MutationSpectrum([MutationDNASeq(seq).mutate("1A>T"), MutationDNASeq(seq).mutate("3G>T"),MutationDNASeq(seq).mutate("4T>A")])
```
For now the help is a munged pydoc man file... For a more upto date version download the code andimport
```
>>> import mutagenesis
>>> help(mutagenesis)
```

# Description
Classes:

* mutation
* MutationFormatError
* MutationDNASeq
* mutationSpectrum

This is a partial rewrite of mutanalyst js code. As a result a lot of attribute names are in camelCase, following JS style as opposed to PEP8.

# CLASSES
## Mutation
class Mutation(builtins.object)
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

Methods defined here:

__init__(self, mutation, seq=None, forceDNA=False, coding=True)
    Initialize self.  See help(type(self)) for accurate signature.

__str__(self)
    Return str(self).

apply(self, seq)

shortform(self)

----------------------------------------------------------------------
Data descriptors defined here:

__dict__
dictionary for instance variables (if defined)

__weakref__
list of weak references to the object (if defined)

----------------------------------------------------------------------
Data and other attributes defined here:

codon_codex = {'AAA': {'*': 'TAA', 'A': 'GCA', 'C': 'TGY', 'D': 'GAY',...
    
# MutationDNASeq

class MutationDNASeq(Bio.Seq.Seq)
A variant of the seq class, but with the method mutate and the attribute mutations.
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

Method resolution order:
    MutationDNASeq
    Bio.Seq.Seq
    builtins.object

Methods defined here:

__init__(self, data)
    Create a Seq object.
    
Arguments:
*  seq - Sequence, required (string)
*  alphabet - Optional argument, an Alphabet object from Bio.Alphabet

You will typically use Bio.SeqIO to read in sequences from files as
SeqRecord objects, whose sequence will be exposed as a Seq object via
the seq property.

However, will often want to create your own Seq objects directly:

> from Bio.Seq import Seq
> from Bio.Alphabet import IUPAC
> my_seq = Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF",
...              IUPAC.protein)
> my_seq
Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF', IUPACProtein())
> print(my_seq)
MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
> my_seq.alphabet
IUPACProtein()

mutate(self, mutations, forceDNA=False)

----------------------------------------------------------------------
Methods inherited from Bio.Seq.Seq:

__add__(self, other)
Add another sequence or string to this sequence.

If adding a string to a Seq, the alphabet is preserved:

> from Bio.Seq import Seq
> from Bio.Alphabet import generic_protein
> Seq("MELKI", generic_protein) + "LV"
Seq('MELKILV', ProteinAlphabet())

When adding two Seq (like) objects, the alphabets are important.
Consider this example:

> from Bio.Seq import Seq
> from Bio.Alphabet.IUPAC import unambiguous_dna, ambiguous_dna
> unamb_dna_seq = Seq("ACGT", unambiguous_dna)
> ambig_dna_seq = Seq("ACRGT", ambiguous_dna)
> unamb_dna_seq
Seq('ACGT', IUPACUnambiguousDNA())
> ambig_dna_seq
Seq('ACRGT', IUPACAmbiguousDNA())

If we add the ambiguous and unambiguous IUPAC DNA alphabets, we get
the more general ambiguous IUPAC DNA alphabet:

> unamb_dna_seq + ambig_dna_seq
Seq('ACGTACRGT', IUPACAmbiguousDNA())

However, if the default generic alphabet is included, the result is
a generic alphabet:

> Seq("") + ambig_dna_seq
Seq('ACRGT', Alphabet())

You can't add RNA and DNA sequences:

> from Bio.Alphabet import generic_dna, generic_rna
> Seq("ACGT", generic_dna) + Seq("ACGU", generic_rna)
Traceback (most recent call last):
   ...
TypeError: Incompatible alphabets DNAAlphabet() and RNAAlphabet()

You can't add nucleotide and protein sequences:

> from Bio.Alphabet import generic_dna, generic_protein
> Seq("ACGT", generic_dna) + Seq("MELKI", generic_protein)
Traceback (most recent call last):
   ...
TypeError: Incompatible alphabets DNAAlphabet() and ProteinAlphabet()

__contains__(self, char)
Implements the 'in' keyword, like a python string.

e.g.

> from Bio.Seq import Seq
> from Bio.Alphabet import generic_dna, generic_rna, generic_protein
> my_dna = Seq("ATATGAAATTTGAAAA", generic_dna)
> "AAA" in my_dna
True
> Seq("AAA") in my_dna
True
> Seq("AAA", generic_dna) in my_dna
True

Like other Seq methods, this will raise a type error if another Seq
(or Seq like) object with an incompatible alphabet is used:

> Seq("AAA", generic_rna) in my_dna
Traceback (most recent call last):
   ...
TypeError: Incompatible alphabets DNAAlphabet() and RNAAlphabet()
> Seq("AAA", generic_protein) in my_dna
Traceback (most recent call last):
   ...
TypeError: Incompatible alphabets DNAAlphabet() and ProteinAlphabet()

__eq__(self, other)
Compare the sequence to another sequence or a string (README).

Historically comparing Seq objects has done Python object comparison.
After considerable discussion (keeping in mind constraints of the
Python language, hashes and dictionary support), Biopython now uses
simple string comparison (with a warning about the change).

Note that incompatible alphabets (e.g. DNA to RNA) will trigger a
warning.

During this transition period, please just do explicit comparisons:

> from Bio.Seq import Seq
> from Bio.Alphabet import generic_dna
> seq1 = Seq("ACGT")
> seq2 = Seq("ACGT")
> id(seq1) == id(seq2)
False
> str(seq1) == str(seq2)
True

The new behaviour is to use string-like equality:

> from Bio.Seq import Seq
> from Bio.Alphabet import generic_dna
> seq1 == seq2
True
> seq1 == "ACGT"
True
> seq1 == Seq("ACGT", generic_dna)
True

__getitem__(self, index)
Returns a subsequence of single letter, use my_seq[index].

__hash__(self)
Hash for comparison.

See the __cmp__ documentation - this has changed from past
versions of Biopython!

__le__(self, other)
Less than or equal, see __eq__ documentation.

__len__(self)
Returns the length of the sequence, use len(my_seq).

__lt__(self, other)
Less than, see __eq__ documentation.

__ne__(self, other)
Not equal, see __eq__ documentation.

__radd__(self, other)
Adding a sequence on the left.

If adding a string to a Seq, the alphabet is preserved:

> from Bio.Seq import Seq
> from Bio.Alphabet import generic_protein
> "LV" + Seq("MELKI", generic_protein)
Seq('LVMELKI', ProteinAlphabet())

Adding two Seq (like) objects is handled via the __add__ method.

__repr__(self)
Returns a (truncated) representation of the sequence for debugging.

__str__(self)
Returns the full sequence as a python string, use str(my_seq).

Note that Biopython 1.44 and earlier would give a truncated
version of repr(my_seq) for str(my_seq).  If you are writing code
which need to be backwards compatible with old Biopython, you
should continue to use my_seq.tostring() rather than str(my_seq).

back_transcribe(self)
Returns the DNA sequence from an RNA sequence. New Seq object.

> from Bio.Seq import Seq
> from Bio.Alphabet import IUPAC
> messenger_rna = Seq("AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG",
...                     IUPAC.unambiguous_rna)
> messenger_rna
Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG', IUPACUnambiguousRNA())
> messenger_rna.back_transcribe()
Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG', IUPACUnambiguousDNA())

Trying to back-transcribe a protein or DNA sequence raises an
exception:

> my_protein = Seq("MAIVMGR", IUPAC.protein)
> my_protein.back_transcribe()
Traceback (most recent call last):
   ...
ValueError: Proteins cannot be back transcribed!

complement(self)
Returns the complement sequence. New Seq object.

> from Bio.Seq import Seq
> from Bio.Alphabet import IUPAC
> my_dna = Seq("CCCCCGATAG", IUPAC.unambiguous_dna)
> my_dna
Seq('CCCCCGATAG', IUPACUnambiguousDNA())
> my_dna.complement()
Seq('GGGGGCTATC', IUPACUnambiguousDNA())

You can of course used mixed case sequences,

> from Bio.Seq import Seq
> from Bio.Alphabet import generic_dna
> my_dna = Seq("CCCCCgatA-GD", generic_dna)
> my_dna
Seq('CCCCCgatA-GD', DNAAlphabet())
> my_dna.complement()
Seq('GGGGGctaT-CH', DNAAlphabet())

Note in the above example, ambiguous character D denotes
G, A or T so its complement is H (for C, T or A).

Trying to complement a protein sequence raises an exception.

> my_protein = Seq("MAIVMGR", IUPAC.protein)
> my_protein.complement()
Traceback (most recent call last):
   ...
ValueError: Proteins do not have complements!

count(self, sub, start=0, end=9223372036854775807)
Non-overlapping count method, like that of a python string.

This behaves like the python string method of the same name,
which does a non-overlapping count!

Returns an integer, the number of occurrences of substring
argument sub in the (sub)sequence given by [start:end].
Optional arguments start and end are interpreted as in slice
notation.

Arguments:
*  sub - a string or another Seq object to look for
*  start - optional integer, slice start
*  end - optional integer, slice end

e.g.

> from Bio.Seq import Seq
> my_seq = Seq("AAAATGA")
> print(my_seq.count("A"))
5
> print(my_seq.count("ATG"))
1
> print(my_seq.count(Seq("AT")))
1
> print(my_seq.count("AT", 2, -1))
1

HOWEVER, please note because python strings and Seq objects (and
MutableSeq objects) do a non-overlapping search, this may not give
the answer you expect:

> "AAAA".count("AA")
2
> print(Seq("AAAA").count("AA"))
2

An overlapping search would give the answer as three!

endswith(self, suffix, start=0, end=9223372036854775807)
Does the Seq end with the given suffix?  Returns True/False.

This behaves like the python string method of the same name.

Return True if the sequence ends with the specified suffix
(a string or another Seq object), False otherwise.
With optional start, test sequence beginning at that position.
With optional end, stop comparing sequence at that position.
suffix can also be a tuple of strings to try.  e.g.

> from Bio.Seq import Seq
> my_rna = Seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
> my_rna.endswith("UUG")
True
> my_rna.endswith("AUG")
False
> my_rna.endswith("AUG", 0, 18)
True
> my_rna.endswith(("UCC", "UCA", "UUG"))
True

find(self, sub, start=0, end=9223372036854775807)
Find method, like that of a python string.

This behaves like the python string method of the same name.

Returns an integer, the index of the first occurrence of substring
argument sub in the (sub)sequence given by [start:end].

Arguments:
*  sub - a string or another Seq object to look for
*  start - optional integer, slice start
*  end - optional integer, slice end

Returns -1 if the subsequence is NOT found.

e.g. Locating the first typical start codon, AUG, in an RNA sequence:

> from Bio.Seq import Seq
> my_rna = Seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
> my_rna.find("AUG")
3

lower(self)
Returns a lower case copy of the sequence.

This will adjust the alphabet if required. Note that the IUPAC alphabets
are upper case only, and thus a generic alphabet must be substituted.

> from Bio.Alphabet import Gapped, generic_dna
> from Bio.Alphabet import IUPAC
> from Bio.Seq import Seq
> my_seq = Seq("CGGTACGCTTATGTCACGTAG*AAAAAA", Gapped(IUPAC.unambiguous_dna, "*"))
> my_seq
Seq('CGGTACGCTTATGTCACGTAG*AAAAAA', Gapped(IUPACUnambiguousDNA(), '*'))
> my_seq.lower()
Seq('cggtacgcttatgtcacgtag*aaaaaa', Gapped(DNAAlphabet(), '*'))

See also the upper method.

lstrip(self, chars=None)
Returns a new Seq object with leading (left) end stripped.

This behaves like the python string method of the same name.

Optional argument chars defines which characters to remove.  If
omitted or None (default) then as for the python string method,
this defaults to removing any white space.

e.g. print(my_seq.lstrip("-"))

See also the strip and rstrip methods.

reverse_complement(self)
Returns the reverse complement sequence. New Seq object.

> from Bio.Seq import Seq
> from Bio.Alphabet import IUPAC
> my_dna = Seq("CCCCCGATAGNR", IUPAC.ambiguous_dna)
> my_dna
Seq('CCCCCGATAGNR', IUPACAmbiguousDNA())
> my_dna.reverse_complement()
Seq('YNCTATCGGGGG', IUPACAmbiguousDNA())

Note in the above example, since R = G or A, its complement
is Y (which denotes C or T).

You can of course used mixed case sequences,

> from Bio.Seq import Seq
> from Bio.Alphabet import generic_dna
> my_dna = Seq("CCCCCgatA-G", generic_dna)
> my_dna
Seq('CCCCCgatA-G', DNAAlphabet())
> my_dna.reverse_complement()
Seq('C-TatcGGGGG', DNAAlphabet())

Trying to complement a protein sequence raises an exception:

> my_protein = Seq("MAIVMGR", IUPAC.protein)
> my_protein.reverse_complement()
Traceback (most recent call last):
   ...
ValueError: Proteins do not have complements!

rfind(self, sub, start=0, end=9223372036854775807)
Find from right method, like that of a python string.

This behaves like the python string method of the same name.

Returns an integer, the index of the last (right most) occurrence of
substring argument sub in the (sub)sequence given by [start:end].

Arguments:
*  sub - a string or another Seq object to look for
*  start - optional integer, slice start
*  end - optional integer, slice end

Returns -1 if the subsequence is NOT found.

e.g. Locating the last typical start codon, AUG, in an RNA sequence:

> from Bio.Seq import Seq
> my_rna = Seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
> my_rna.rfind("AUG")
15

rsplit(self, sep=None, maxsplit=-1)
Right split method, like that of a python string.

This behaves like the python string method of the same name.

Return a list of the 'words' in the string (as Seq objects),
using sep as the delimiter string.  If maxsplit is given, at
most maxsplit splits are done COUNTING FROM THE RIGHT.
If maxsplit is omitted, all splits are made.

Following the python string method, sep will by default be any
white space (tabs, spaces, newlines) but this is unlikely to
apply to biological sequences.

e.g. print(my_seq.rsplit("*",1))

See also the split method.

rstrip(self, chars=None)
Returns a new Seq object with trailing (right) end stripped.

This behaves like the python string method of the same name.

Optional argument chars defines which characters to remove.  If
omitted or None (default) then as for the python string method,
this defaults to removing any white space.

e.g. Removing a nucleotide sequence's polyadenylation (poly-A tail):

> from Bio.Alphabet import IUPAC
> from Bio.Seq import Seq
> my_seq = Seq("CGGTACGCTTATGTCACGTAGAAAAAA", IUPAC.unambiguous_dna)
> my_seq
Seq('CGGTACGCTTATGTCACGTAGAAAAAA', IUPACUnambiguousDNA())
> my_seq.rstrip("A")
Seq('CGGTACGCTTATGTCACGTAG', IUPACUnambiguousDNA())

See also the strip and lstrip methods.

split(self, sep=None, maxsplit=-1)
Split method, like that of a python string.

This behaves like the python string method of the same name.

Return a list of the 'words' in the string (as Seq objects),
using sep as the delimiter string.  If maxsplit is given, at
most maxsplit splits are done.  If maxsplit is omitted, all
splits are made.

Following the python string method, sep will by default be any
white space (tabs, spaces, newlines) but this is unlikely to
apply to biological sequences.

e.g.

> from Bio.Seq import Seq
> my_rna = Seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
> my_aa = my_rna.translate()
> my_aa
Seq('VMAIVMGR*KGAR*L', HasStopCodon(ExtendedIUPACProtein(), '*'))
> my_aa.split("*")
[Seq('VMAIVMGR', HasStopCodon(ExtendedIUPACProtein(), '*')), Seq('KGAR', HasStopCodon(ExtendedIUPACProtein(), '*')), Seq('L', HasStopCodon(ExtendedIUPACProtein(), '*'))]
> my_aa.split("*", 1)
[Seq('VMAIVMGR', HasStopCodon(ExtendedIUPACProtein(), '*')), Seq('KGAR*L', HasStopCodon(ExtendedIUPACProtein(), '*'))]

See also the rsplit method:

> my_aa.rsplit("*", 1)
[Seq('VMAIVMGR*KGAR', HasStopCodon(ExtendedIUPACProtein(), '*')), Seq('L', HasStopCodon(ExtendedIUPACProtein(), '*'))]

startswith(self, prefix, start=0, end=9223372036854775807)
Does the Seq start with the given prefix?  Returns True/False.

This behaves like the python string method of the same name.

Return True if the sequence starts with the specified prefix
(a string or another Seq object), False otherwise.
With optional start, test sequence beginning at that position.
With optional end, stop comparing sequence at that position.
prefix can also be a tuple of strings to try.  e.g.

> from Bio.Seq import Seq
> my_rna = Seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG")
> my_rna.startswith("GUC")
True
> my_rna.startswith("AUG")
False
> my_rna.startswith("AUG", 3)
True
> my_rna.startswith(("UCC", "UCA", "UCG"), 1)
True

strip(self, chars=None)
Returns a new Seq object with leading and trailing ends stripped.

This behaves like the python string method of the same name.

Optional argument chars defines which characters to remove.  If
omitted or None (default) then as for the python string method,
this defaults to removing any white space.

e.g. print(my_seq.strip("-"))

See also the lstrip and rstrip methods.

tomutable(self)
Returns the full sequence as a MutableSeq object.

> from Bio.Seq import Seq
> from Bio.Alphabet import IUPAC
> my_seq = Seq("MKQHKAMIVALIVICITAVVAAL",
...              IUPAC.protein)
> my_seq
Seq('MKQHKAMIVALIVICITAVVAAL', IUPACProtein())
> my_seq.tomutable()
MutableSeq('MKQHKAMIVALIVICITAVVAAL', IUPACProtein())

Note that the alphabet is preserved.

tostring(self)
Returns the full sequence as a python string (DEPRECATED).

You are now encouraged to use str(my_seq) instead of
my_seq.tostring().

transcribe(self)
Returns the RNA sequence from a DNA sequence. New Seq object.

> from Bio.Seq import Seq
> from Bio.Alphabet import IUPAC
> coding_dna = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG",
...                  IUPAC.unambiguous_dna)
> coding_dna
Seq('ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG', IUPACUnambiguousDNA())
> coding_dna.transcribe()
Seq('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG', IUPACUnambiguousRNA())

Trying to transcribe a protein or RNA sequence raises an exception:

> my_protein = Seq("MAIVMGR", IUPAC.protein)
> my_protein.transcribe()
Traceback (most recent call last):
   ...
ValueError: Proteins cannot be transcribed!

translate(self, table='Standard', stop_symbol='*', to_stop=False, cds=False)
Turns a nucleotide sequence into a protein sequence. New Seq object.

This method will translate DNA or RNA sequences, and those with a
nucleotide or generic alphabet.  Trying to translate a protein
sequence raises an exception.

Arguments:
* table - Which codon table to use?  This can be either a name
  (string), an NCBI identifier (integer), or a CodonTable
  object (useful for non-standard genetic codes).  This
  defaults to the "Standard" table.
* stop_symbol - Single character string, what to use for terminators.
  This defaults to the asterisk, "*".
* to_stop - Boolean, defaults to False meaning do a full translation
  continuing on past any stop codons (translated as the
  specified stop_symbol).  If True, translation is
  terminated at the first in frame stop codon (and the
  stop_symbol is not appended to the returned protein
  sequence).
* cds - Boolean, indicates this is a complete CDS.  If True,
  this checks the sequence starts with a valid alternative start
  codon (which will be translated as methionine, M), that the
  sequence length is a multiple of three, and that there is a
  single in frame stop codon at the end (this will be excluded
  from the protein sequence, regardless of the to_stop option).
  If these tests fail, an exception is raised.

e.g. Using the standard table:

> coding_dna = Seq("GTGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
> coding_dna.translate()
Seq('VAIVMGR*KGAR*', HasStopCodon(ExtendedIUPACProtein(), '*'))
> coding_dna.translate(stop_symbol="@")
Seq('VAIVMGR@KGAR@', HasStopCodon(ExtendedIUPACProtein(), '@'))
> coding_dna.translate(to_stop=True)
Seq('VAIVMGR', ExtendedIUPACProtein())

Now using NCBI table 2, where TGA is not a stop codon:

> coding_dna.translate(table=2)
Seq('VAIVMGRWKGAR*', HasStopCodon(ExtendedIUPACProtein(), '*'))
> coding_dna.translate(table=2, to_stop=True)
Seq('VAIVMGRWKGAR', ExtendedIUPACProtein())

In fact, GTG is an alternative start codon under NCBI table 2, meaning
this sequence could be a complete CDS:

> coding_dna.translate(table=2, cds=True)
Seq('MAIVMGRWKGAR', ExtendedIUPACProtein())

It isn't a valid CDS under NCBI table 1, due to both the start codon and
also the in frame stop codons:

> coding_dna.translate(table=1, cds=True)
Traceback (most recent call last):
    ...
TranslationError: First codon 'GTG' is not a start codon

If the sequence has no in-frame stop codon, then the to_stop argument
has no effect:

> coding_dna2 = Seq("TTGGCCATTGTAATGGGCCGC")
> coding_dna2.translate()
Seq('LAIVMGR', ExtendedIUPACProtein())
> coding_dna2.translate(to_stop=True)
Seq('LAIVMGR', ExtendedIUPACProtein())

NOTE - Ambiguous codons like "TAN" or "NNN" could be an amino acid
or a stop codon.  These are translated as "X".  Any invalid codon
(e.g. "TA?" or "T-A") will throw a TranslationError.

NOTE - Does NOT support gapped sequences.

NOTE - This does NOT behave like the python string's translate
method.  For that use str(my_seq).translate(...) instead.

ungap(self, gap=None)
Return a copy of the sequence without the gap character(s).

The gap character can be specified in two ways - either as an explicit
argument, or via the sequence's alphabet. For example:

> from Bio.Seq import Seq
> from Bio.Alphabet import generic_dna
> my_dna = Seq("-ATA--TGAAAT-TTGAAAA", generic_dna)
> my_dna
Seq('-ATA--TGAAAT-TTGAAAA', DNAAlphabet())
> my_dna.ungap("-")
Seq('ATATGAAATTTGAAAA', DNAAlphabet())

If the gap character is not given as an argument, it will be taken from
the sequence's alphabet (if defined). Notice that the returned sequence's
alphabet is adjusted since it no longer requires a gapped alphabet:

> from Bio.Seq import Seq
> from Bio.Alphabet import IUPAC, Gapped, HasStopCodon
> my_pro = Seq("MVVLE=AD*", HasStopCodon(Gapped(IUPAC.protein, "=")))
> my_pro
Seq('MVVLE=AD*', HasStopCodon(Gapped(IUPACProtein(), '='), '*'))
> my_pro.ungap()
Seq('MVVLEAD*', HasStopCodon(IUPACProtein(), '*'))

Or, with a simpler gapped DNA example:

> from Bio.Seq import Seq
> from Bio.Alphabet import IUPAC, Gapped
> my_seq = Seq("CGGGTAG=AAAAAA", Gapped(IUPAC.unambiguous_dna, "="))
> my_seq
Seq('CGGGTAG=AAAAAA', Gapped(IUPACUnambiguousDNA(), '='))
> my_seq.ungap()
Seq('CGGGTAGAAAAAA', IUPACUnambiguousDNA())

As long as it is consistent with the alphabet, although it is redundant,
you can still supply the gap character as an argument to this method:

> my_seq
Seq('CGGGTAG=AAAAAA', Gapped(IUPACUnambiguousDNA(), '='))
> my_seq.ungap("=")
Seq('CGGGTAGAAAAAA', IUPACUnambiguousDNA())

However, if the gap character given as the argument disagrees with that
declared in the alphabet, an exception is raised:

> my_seq
Seq('CGGGTAG=AAAAAA', Gapped(IUPACUnambiguousDNA(), '='))
> my_seq.ungap("-")
Traceback (most recent call last):
   ...
ValueError: Gap '-' does not match '=' from alphabet

Finally, if a gap character is not supplied, and the alphabet does not
define one, an exception is raised:

> from Bio.Seq import Seq
> from Bio.Alphabet import generic_dna
> my_dna = Seq("ATA--TGAAAT-TTGAAAA", generic_dna)
> my_dna
Seq('ATA--TGAAAT-TTGAAAA', DNAAlphabet())
> my_dna.ungap()
Traceback (most recent call last):
   ...
ValueError: Gap character not given and not defined in alphabet

upper(self)
Returns an upper case copy of the sequence.

> from Bio.Alphabet import HasStopCodon, generic_protein
> from Bio.Seq import Seq
> my_seq = Seq("VHLTPeeK*", HasStopCodon(generic_protein))
> my_seq
Seq('VHLTPeeK*', HasStopCodon(ProteinAlphabet(), '*'))
> my_seq.lower()
Seq('vhltpeek*', HasStopCodon(ProteinAlphabet(), '*'))
> my_seq.upper()
Seq('VHLTPEEK*', HasStopCodon(ProteinAlphabet(), '*'))

This will adjust the alphabet if required. See also the lower method.

----------------------------------------------------------------------
Data descriptors inherited from Bio.Seq.Seq:

__dict__
dictionary for instance variables (if defined)

__weakref__
list of weak references to the object (if defined)

## MutationFormatError    
class MutationFormatError(builtins.Exception)
Common base class for all non-exit exceptions.

Method resolution order:
MutationFormatError
builtins.Exception
builtins.BaseException
builtins.object

Methods defined here:

__init__(self, value=None)
Initialize self.  See help(type(self)) for accurate signature.

__str__(self)
Return str(self).

----------------------------------------------------------------------
Data descriptors defined here:

__weakref__
list of weak references to the object (if defined)

----------------------------------------------------------------------
Data and other attributes defined here:

message = 'Error in the parsing a mutation notation.\n    A ...ould be...

----------------------------------------------------------------------
Methods inherited from builtins.Exception:

__new__(*args, **kwargs) from builtins.type
Create and return a new object.  See help(type) for accurate signature.

----------------------------------------------------------------------
Methods inherited from builtins.BaseException:

__delattr__(self, name, /)
Implement delattr(self, name).

__getattribute__(self, name, /)
Return getattr(self, name).

__reduce__(...)
helper for pickle

__repr__(self, /)
Return repr(self).

__setattr__(self, name, value, /)
Implement setattr(self, name, value).

__setstate__(...)

with_traceback(...)
Exception.with_traceback(tb) --
set self.__traceback__ to tb and return self.

----------------------------------------------------------------------
Data descriptors inherited from builtins.BaseException:

__cause__
exception cause

__context__
exception context

__dict__

__suppress_context__

__traceback__

args
## MutationSpectrum

class MutationSpectrum(builtins.object)
Returns the mutational spectrum, an object with
the mutation frequency,
the base freq
the SE of the mut freq
the mutational rate


class method from sequences
init from values

A lot of this will be plagiarised from JS https://github.com/matteoferla/mutant_calculator/blob/master/mutationalBias.js

Methods defined here:

__getitem__(self, item)

__init__(self, mutants)
Initialize self.  See help(type(self)) for accurate signature.

__setitem__(self, item, value)

__str__(self)
Return str(self).

add_mutations(self, mutations)

----------------------------------------------------------------------
Class methods defined here:

from_mutation_list(mutations, seq=None) from builtins.type

from_values(source='loaded', sequence='', mutations=None, freqMean=0, freqVar=0, freqList=0, raw_table=<mutagenesis.MutationTable object at 0x102c2b400>, table=<mutagenesis.MutationTable object at 0x102c2dbe0>) from builtins.type

----------------------------------------------------------------------
Data descriptors defined here:

__dict__
dictionary for instance variables (if defined)

__weakref__
list of weak references to the object (if defined)

----------------------------------------------------------------------
Data and other attributes defined here:

types = ['TsOverTv', 'W2SOverS2W', 'W2N', 'S2N', 'W2S', 'S2W', 'ΣTs', ...
 
## MutationTable
class MutationTable(builtins.object)
ATGC^2 table. The values are accessed with a A>C notation.
Due to the fact that an argument name must be a valid variable name "A>C" cannot be given as MutationTable(A>C=1), but has to be given as MutationTable({A>C: 1})
To access a frequency, use instance["A>C"] notation

Methods defined here:

__getitem__(self, item)

__init__(self, frequencies=None)
Initialize self.  See help(type(self)) for accurate signature.

__iter__(self)

__setitem__(self, item, value)

__str__(self)
Return str(self).

normalize(self, A=0.25, T=0.25, G=0.25, C=0.25)

to_dict(self)

----------------------------------------------------------------------
Static methods defined here:

ibase()

----------------------------------------------------------------------
Data descriptors defined here:

__dict__
dictionary for instance variables (if defined)

__weakref__
list of weak references to the object (if defined)

# FUNCTIONS

generateCodonCodex()

mincodondist(codon, aa)

test()

# DATA

N = '\n'
T = '\t'

# VERSION

$Revision$

# AUTHOR

Matteo

# FILE

/Users/matfer/PycharmProjects/mutagenesis/mutagenesis.py


