# Human diversity reference (DivRef) bundle

*Version 1.1, 2025-02-12*

This resource bundle contains a custom human diversity reference (DivRef) composed of short base sequences around common variation in the
human population, with a tool for remapping coordinates from this custom reference back to GRCh38.

DivRef is constructed from two data resources, the 
[Human Genome Diversity Panel (HGDP)](https://www.internationalgenome.org/data) and [gnomAD 4.1](https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1/).
Using the [phased HGDP Hail dataset provided by the gnomAD team at the Broad Institute](https://gnomad.broadinstitute.org/downloads), we compute groups of common (> 0.5%) phased haplotypes in close 
proximity (<25 base pairs) and merge with the common variants in the gnomAD 4.1 release. **Note that given the sample size of the HGDP dataset,
there is limited power to detect haplotypes under 1%; we have included haplotypes discovered between 0.5% and 1%, but would expect to find more
haplotypes in that frequency range with a larger dataset.**
We then look up the 25-base-pair sequence context around these variants and haplotypes using the GRCh38 reference
sequence, and export all such sequences as FASTA reference files.

These FASTA files are intended to be used with tools for guide design and off-target nomination that already accept
reference sequences in FASTA format.

The haplotype-only version of the resource contains 2026924 sequences overlapping 2+ variants.
The gnomAD + HGDP merge contains 36574883 sequences including single variants and haplotypes.

For full details, see the GitHub repository: [human-diversity-reference](https://github.com/e9genomics/human-diversity-reference)

## Changelog

- v1.1, 2025-02-11:
  - Updated haplotype filtering algorithm to improve discovery of haplotypes <= 1% AF, which added ~200k haplotypes.
  - Fixed bug in `remap_divref.py` script

## Contents

This bundle includes:

- A README file (you're reading it)
- License files for scripts (CC-BY-4.0-LICENSE) and data (LICENSE)
- Two sets of FASTA files containing sequence context around common human variation:
  - `DivRef-v1.1.haplotypes.fasta`, which contains all sequences around 2 or more phased common variants in HGDP, provided for convenience for applications that do not require gnomAD variants.
  - 24 FASTA files `DivRef-v1.1.haplotypes_gnomad_merge.*.fasta`, one for each chromosome, which is a superset of `DivRef-v1.1.haplotypes.fasta`, including sequences surrounding single variants from gnomAD 4.1. 
- A DuckDB index used for fast coordinate liftover. A single index is provided for both resources; a smaller index file for haplotype-only applications is available on request.
- Two compressed TSV files with all included sequences and HGDP / gnomAD frequency information in a more readable format. One is HGDP-only, the other is merged with single variants from gnomAD.
- A remapping script (`remap_divref.py`)

## License & usage restrictions

The software (scripts, etc.) in this project is MIT-licensed with no restrictions on either commercial or non-commercial usage with no attribution requirement.
The data resources are licensed under Creative Commons CC BY 4.0, which permits both commercial and non-commercial usage, and requires attribution to E9 Genomics and this resource in redistributions and derivative works.

For the full license terms, see the license files in the bundle directory or [GitHub repository](https://github.com/e9genomics/human-diversity-reference).

## Dependencies

The reindexing script depends on the *amazing* Python environment management tool `uv` -- [install here](https://docs.astral.sh/uv/).

The script declares its dependencies inline (reference documentation [here](https://docs.astral.sh/uv/guides/scripts/#declaring-script-dependencies)).

## FASTA indices and dictionaries

[SAMTools](https://www.htslib.org/download/) is used to generate FASTA indices and dictionaries. Indices and dictionaries
are not included in the bundle in order to save space; the large number of "contigs" in DivRev mean that
the index and dictionary files are about as large as the FASTA itself.

To create indices for a FASTA file with SAMTools, run:
```bash
$SAMTOOLS faidx $FASTA
```

To create a FASTA dictionary, run:
```bash
$SAMTOOLS dict $FASTA -o ${FASTA}.dict -u NA
```

Note the inclusion of the `"-u NA"` argument, which prevents the repetition of the full file URI for each line of the
dictionary and helps control dictionary file size.

## Using the Remapping Tool

The bundle includes a tool for remapping coordinates from DivRef space back to GRCh38. Currently only TSV files resembling
CALITAS outputs are supported, but others can be added on request.

### CALITAS Output Remapping

To remap coordinates from CALITAS output:

```bash
uv run remap_divref.py calitas input.tsv output.tsv
```

Optional parameters:

- `-i INDEX_PATH`: Path to the DuckDB index file (optional: not necessary when the index database is in the same directory as the script)
- `-s SEPARATOR`: Input/output file separator (default: tab)

The input file must contain these columns:

- `chromosome`: DivRev sequence ID
- `coordinate_start`: Start position in DivRev coordinates
- `coordinate_end`: End position in DivRev coordinates

The output will contain the original columns plus:

- `divref_sequence_id`, `divref_start`, `divref_end`: Original DivRev coordinates
- `chromosome`, `coordinate_start`, `coordinate_end`: Remapped reference genome coordinates (in same positions)
- `all_variants`: All variants included in sequence, in the format `chr:start:ref:alt`, comma-delimited
- `variants_involved`: Variants overlapping the alignment interval, in the format `chr:start:ref:alt`, comma-delimited
- `n_variants_involved`: Number of variants overlapping the alignment interval
- `max_pop`: Population with the highest frequency for this variant/haplotype
- `variant_source`: Origin of the sequence ("HGDP_haplotype" or "gnomAD_variant")
- `popmax_empirical_AF`: The empirical allele frequency of the variant/haplotype in the population with highest frequency. This is the frequency of the entire haplotype, not the subset of variants overlapping the alignment interval.
- `popmax_empirical_AC`: The allele count of the variant/haplotype in the population with highest frequency. This is the count of the entire haplotype, not the subset of variants overlapping the alignment interval.
- `fraction_phased`: For haplotypes, the ratio of the observed haplotype frequency to the frequency of the rarest component variant (always 1.0 for single gnomAD variants).
- `estimated_gnomad_AF`: Estimated gnomAD frequency, computed by multiplying the gnomAD AF of the rarest component variant by `fraction_phased`
- `population_frequencies_json`: JSON string listing the population frequencies for each variant, in the format 
   `{"afr":[0.68271,0.30024],"amr":[0.53207,0.47267],"eas":[0.83507,0.83443],"nfe":[0.46066,0.42511],"sas":[0.57101,0.536]}`.

⚠️ Note: The remapping tool computes the variants overlapped by an alignment region, which might be fewer variants than are included in the
full haplotype. The `fraction_phased` and `popmax_empirical_AF` are computed for the entire haplotype, not just the variants overlapping the
alignment interval, so the `popmax_empirical_AF` might be an *underestimate* of the frequency of the variants overlapping the alignment interval.
This may be improved in a future version of DivRef.

### Haplotype filtering algorithm details

The haplotypes included in the resource are generated and filtered according to the following criteria:

- Individuals in HGDP are annotated with continental ancestry using [gnomAD labels](https://gnomad.broadinstitute.org/data).
- Genotypes at variants with less than 0.5% AF in each individual's annotated ancestry are removed.
- Sets of phased alleles are computed over 100-base-pair windows of the genome.
- The phase ratio is computed by dividing the empirical AF of the haplotype by the empirical AF of the rarest component variant.
- The gnomAD haplotype frequency is estimated by multiplying the gnomAD AF of the rarest component variant by the phase ratio.
- Haplotypes with estimated gnomAD frequency < 0.5% are removed.

### Questions / Feedback

Please contact Tim Poterba (tpoterba@e9genomics.com) with questions or feedback!

---

From the team at [E9 Genomics](https://e9genomics.com/).
