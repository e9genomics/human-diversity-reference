# Human diversity reference (DivRef) bundle

This resource bundle contains a custom human diversity reference (DivRef) composed of short base sequences around common variation in the
human population, with a tool for remapping coordinates from this custom reference back to GRCh38.

DivRef is constructed by computing empirical phased haplotypes within 25 BPs over 0.5% allele frequency from the 
[Human Genome Diversity Panel (HGDP)](https://www.internationalgenome.org/data) using the [phased Hail 
dataset provided by the gnomAD team at the Broad Institute](https://gnomad.broadinstitute.org/downloads), 
merged with single variants over 0.5% AF from the gnomAD v4.1.0 summary release. The base context including
and surrounding these variants and haplotypes is generated from the GRCh38 reference sequence, and each variant
and haplotype is exported to a FASTA file.

These FASTA files are intended to be used with tools for guide design and off-target nomination that already accept
reference sequences as FASTA files.

The haplotype-only version of the resource contains 1,795,661 sequences overlapping 2+ variants.
The gnomAD + HGDP merge contains 36,343,619 sequences including single variants and haplotypes.

For documentation on using the bundled distributed on Zenodo, including a changelog, see [the bundle README](./bundle/README.md).

---

From the team at [E9 Genomics](https://e9genomics.com/).
