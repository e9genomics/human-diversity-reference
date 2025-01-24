# Human diversity reference (DivRef) bundle

This resource bundle contains a custom human diversity reference (DivRef) composed of short base sequences around common variation in the
human population, with a tool for remapping coordinates from this custom reference back to GRCh38.

DivRef is constructed from two data resources, the 
[Human Genome Diversity Panel (HGDP)](https://www.internationalgenome.org/data) and [gnomAD 4.1](https://gnomad.broadinstitute.org/news/2024-04-gnomad-v4-1/).
Using the [phased HGDP Hail dataset provided by the gnomAD team at the Broad Institute](https://gnomad.broadinstitute.org/downloads), we compute groups of common (> 0.5%) phased haplotypes in close 
proximity (<25 base pairs) and merge with the common variants in the gnomAD 4.1 release. 
We then look up the 25-base-pair sequence context around these variants and haplotypes using the GRCh38 reference
sequence, and export all such sequences as FASTA reference files.

These FASTA files are intended to be used with tools for guide design and off-target nomination that already accept
reference sequences in FASTA format.

The haplotype-only version of the resource contains 1795661 sequences overlapping 2+ variants.
The gnomAD + HGDP merge contains 36343619 sequences including single variants and haplotypes.

For full details, see the GitHub repository: [human-diversity-reference](https://github.com/e9genomics/human-diversity-reference)

## Contents

This bundle includes:

- A README file (you're reading it)
- License files for scripts and data
- Two directories, `DivRef` and `DivRef_just_haplotypes`, which each contain:
  - FASTA files(s) containing sequence context around common human variation
    - The `DivRef` folder contains 24 FASTA files, one for each chromosome: `DivRef.haplotypes_gnomad_merge.*.fasta`
    - The `DivRef_just_haplotypes` folder contains one such FASTA file, `DivRef.haplotypes.fasta`, only with sequences 
      around >=2 phased variants in HGDP. This is a subset of the sequences found in the `DivRef` folder, and is
      provided for convenience and efficiency for applications that only require looking at the haplotypes and not
      gnomAD variants.
  - A DuckDB index used for fast coordinate liftover
  - a compressed TSV file with all included sequences and HGDP / gnomAD frequency information in a more readable format
  - A remapping script (`remap_divref.py`)

## License & usage restrictions

The software (scripts, etc.) in this project is MIT-licensed with no restrictions on either commercial or non-commercial usage.
The data resources are licensed under Creative Commons CC BY 4.0, which permits either commercial or non-commercial usage, but requires attribution to E9 Genomics and this resource in redistributions.

For the full license terms, see the license files in the bundle directory or [GitHub repository]().

## Dependencies

The reindexing script depends on the *amazing* Python environment management tool `uv`  -- [install here](https://docs.astral.sh/uv/).

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

- `divref_chromosome`, `divref_start`, `divref_end`: Original DivRev coordinates
- `chromosome`, `coordinate_start`, `coordinate_end`: Remapped reference genome coordinates (in same positions)
- `variants_involved`: Variants included in the called region, in the format `chr:start:ref:alt` delimited by commas
- `n_variants_involved`: Number of variants in the region
- `min_gnomad_popmax`: Minimum gnomAD popmax frequency of included variants, derived by first computing the maximum
  allele frequency in any population per variant, then taking the minimum across variants.
- `min_gnomad_popmax_pop`: Population with the minimum gnomAD popmax frequency


### Questions / Feedback

Please contact Tim Poterba (tpoterba@e9genomics.com) with questions or feedback!

---

From the team at [E9 Genomics](https://e9genomics.com/).
