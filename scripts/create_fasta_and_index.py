import os

import duckdb
import hail as hl
import polars
import typer

app = typer.Typer()


def to_hashable_items(d):
    return tuple(sorted(d.items()))


def get_haplo_sequence(context_size, variants):
    sorted_variants = hl.sorted(variants, key=lambda x: x.locus.position)
    min_variant = sorted_variants[0]
    max_variant = sorted_variants[-1]
    min_pos = min_variant.locus.position
    max_pos = max_variant.locus.position
    max_variant_size = hl.len(max_variant.alleles[0])
    full_context = hl.get_sequence(min_variant.locus.contig, min_pos, before=context_size,
                                   after=(max_pos - min_pos + max_variant_size + context_size - 1),
                                   reference_genome='GRCh38')

    # translate a locus position into an index in the string
    # (min_pos - index_translation) should equal context size
    index_translation = min_pos - context_size

    def get_chunk_until_next_variant(i):
        v = sorted_variants[i]
        variant_size = hl.len(v.alleles[0])
        reference_buffer_size = hl.if_else(i == hl.len(sorted_variants) - 1,
                                           context_size,
                                           sorted_variants[i + 1].locus.position - (v.locus.position + variant_size))
        start = v.locus.position - index_translation + variant_size
        return v.alleles[1] + full_context[start:start + reference_buffer_size]

    return (full_context[:context_size] +
            hl.delimit(hl.range(hl.len(sorted_variants)).map(lambda i: get_chunk_until_next_variant(i)), ''))


def variant_distance(v1, v2):
    # 1:1:A:T and 1:3:A:T have distance 1 (1 base in between)
    # 1:1:AA:T and 1:3:A:T have distance 0

    return v2.locus.position - v1.locus.position - hl.len(v1.alleles[0])


def split_haplotypes(ht, window_size):
    # generate indices where the distance between variants is greater than the window size

    breakpoints = hl.range(1, hl.len(ht.variants)).filter(
        lambda i: (i == 0) | (variant_distance(ht.variants[i - 1], ht.variants[i]) >= window_size))

    def get_range(i):
        start_index = hl.if_else(i == 0, 0, breakpoints[i - 1])
        end_index = hl.if_else(i == hl.len(breakpoints), hl.len(ht.variants), breakpoints[i])
        return hl.range(start_index, end_index)

    split_hap_indices = hl.range(0, hl.len(breakpoints) + 1).map(get_range).filter(lambda r: hl.len(r) > 1)
    ht = ht.annotate(haplotype_indices=split_hap_indices)
    ht = ht.explode('haplotype_indices')
    ht = ht.annotate(
        haplotype=ht.haplotype_indices.map(lambda i: ht.haplotype[i]),
        variants=ht.haplotype_indices.map(lambda i: ht.variants[i]),
        gnomad_freqs=ht.haplotype_indices.map(lambda i: ht.gnomad_freqs[i]),
    )
    return ht.drop('haplotype_indices')


@app.command()
def main(
        haplotypes_table_path: str = typer.Option(default=..., help="haplotypes table path"),
        gnomad_va_file: str = typer.Option(default=..., help="gnomAD computed variant frequencies"),
        reference_fasta: str = typer.Option(default=..., help="fasta path"),
        window_size: int = typer.Option(default=..., help="Base window size"),
        output_base: str = typer.Option(default=..., help="Output base path"),
        merge: bool = typer.Option(default=False, help="Merge gnomad variants"),
        frequency_cutoff: float = typer.Option(default=0.005, help="Frequency cutoff for gnomAD truncation"),
        split_contigs: bool = typer.Option(default=False, help="Split contigs"),
):
    """
    Process VCF files with gnomAD annotations and output filtered results.
    """
    # Initialize Hail
    hl.init()
    #
    ht = hl.read_table(haplotypes_table_path).key_by()
    va = hl.read_table(gnomad_va_file)
    pops_legend = va.pops.collect()[0]

    va.describe()
    ht.describe()

    hl.get_reference('GRCh38').add_sequence(reference_fasta)

    ht.describe()

    print(f'haplotype table contains {ht.count()} unique haplotypes above frequency threshold')
    ht = split_haplotypes(ht, window_size)
    ht = ht.key_by('haplotype').distinct().key_by().drop('haplotype')
    ht = ht.annotate(source='HGDP_haplotype')
    print(
        f'after splitting at window size {window_size}, haplotype table contains {ht.count()} unique haplotypes above frequency threshold')

    # now add in gnomAD variants if merging
    if merge:
        va = va.rename({'pop_freqs': 'gnomad_freqs'})
        va = va.key_by()
        va = va.select(max_pop=hl.argmax(va.gnomad_freqs.map(lambda x: hl.max(x.AF))),
                       max_pop_freq=hl.max(va.gnomad_freqs.map(lambda x: hl.max(x.AF))),
                       max_pop_observed=hl.max(va.gnomad_freqs.map(lambda x: hl.max(x.AF))),
                       all_pop_freqs=hl.range(hl.len(va.gnomad_freqs)).map(lambda i: hl.struct(
                           pop=i,
                           observed=va.gnomad_freqs[i].AC,
                           frequency=va.gnomad_freqs[i].AF,
                       )),
                       source='gnomAD_variant',
                       variants=[hl.struct(locus=va.locus, alleles=va.alleles)],
                       gnomad_freqs=[va.gnomad_freqs])
        va = va.filter(va.max_pop_freq >= frequency_cutoff)
        ht = ht.union(va, unify=True)

    ht.describe()

    ht = ht.add_index()
    ht = ht.annotate(sequence=get_haplo_sequence(window_size, ht.variants))
    ht = ht.annotate(variant_strs=ht.variants.map(lambda x: hl.variant_str(x)))

    ht = ht.annotate(sequence_length=hl.len(ht.sequence),
                     sequence_id=ht.idx,
                     n_variants=hl.len(ht.variants),
                     ).drop('idx')

    file_suffix = '.haplotypes' if not merge else '.haplotypes_gnomad_merge'

    ht = ht.checkpoint('/tmp/' + f'{file_suffix}.ht', overwrite=True)

    ht.select('sequence',
              'sequence_length',
              'sequence_id',
              'n_variants',
              'max_pop_freq',
              'max_pop_observed',
              hgdp_max_pop=hl.literal(pops_legend)[ht.max_pop],
              variants=hl.delimit(ht.variant_strs, ','),
              **{f'gnomAD_AF_{pop}': hl.delimit(
                  ht.gnomad_freqs.map(lambda x: hl.format("%.5f", x[i].AF)),
                  ','
              ) for i, pop in enumerate(pops_legend)}
              ).export(output_base + f'{file_suffix}.tsv.bgz')

    df = polars.read_csv(output_base + f'{file_suffix}.tsv.bgz', separator='\t',
                         schema_overrides={'sequence_id': polars.String})
    if split_contigs:
        df = df.with_columns(contig=df['variants'].str.split(':').list.get(0))

        for chr in df['contig'].unique().to_list():

            typer.echo('creating FASTA for chromosome ' + chr)
            df2 = df.filter(df['contig'] == chr)

            with open(output_base + f'{file_suffix}.{chr}.fasta', 'w') as f:
                for sequence, sequence_id in df2.select('sequence', 'sequence_id').iter_rows():
                    f.write(f'>{sequence_id}\n{sequence}\n')
    else:
        typer.echo('creating FASTA')
        with open(output_base + f'{file_suffix}.fasta', 'w') as f:
            for sequence, sequence_id in df.select('sequence', 'sequence_id').iter_rows():
                f.write(f'>{sequence_id}\n{sequence}\n')

    duckdb_file = output_base + f'{file_suffix}.index.duckdb'
    if os.path.exists(duckdb_file):
        os.remove(duckdb_file)
    con = duckdb.connect(output_base + f'{file_suffix}.index.duckdb')
    con.execute("CREATE TABLE haplotypes AS SELECT * FROM df")
    con.execute("CREATE INDEX idx_sequence_id ON haplotypes(sequence_id)")

    # create a single window_size value
    con.execute(f"CREATE TABLE window_size AS SELECT {window_size} AS window_size")

    # create a single pops legend value with the list of population identifiers in order
    con.execute(f"CREATE TABLE pops_legend AS SELECT {pops_legend} AS pops_legend")

    con.close()


if __name__ == "__main__":
    app()
