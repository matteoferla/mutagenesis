# mutagenesis
Revamp of the Pedel server. Work in progress.

## Mutants!
Biopython is terrible when it comes to mutagenesis. Hopefully this module will help!

Firstly, it has the class `MutationDNASeq`. This behaves like a normal Seq object except that it has a method called mutate and the attribute mutations.
The instantiation argument is a string like for a seq object or a seq object itself.
The method `mutate()` mutates the sequence based on the  mutations, expressed as a string with spaces or list of strings.
Different customs for nucleotide and protein are used to tell them apart:
* 234A>T for DNA
* A12F for protein, unless forceDNA is true. This is not yet implemented and will require special coding.

The mutations list contains mutation objects.