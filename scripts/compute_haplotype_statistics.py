import hail as hl
import typer
from pydantic import BaseModel

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
        window_sizes: str = typer.Option(default=..., help="Base window sizes, delimited by comma"),
        frequency_cutoffs: str = typer.Option(default=..., help="Frequency thresholds, delimited by comma"),
        output_base: str = typer.Option(default=..., help="Output file path"),
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

    class HGDPResult(BaseModel):
        frequency_cutoff: float
        window_size: int
        hgdp_haplotype_count: int

    class GnomADResult(BaseModel):
        frequency_cutoff: float
        gnomad_variant_count: int

    hgdp_results = []
    gnomad_results = []

    for frequency_str in frequency_cutoffs.split(','):
        frequency = float(frequency_str)
        for window_size_str in window_sizes.split(','):
            window_size = int(window_size_str)
            ht2 = ht.filter(ht.max_pop_freq >= frequency)
            ht2 = split_haplotypes(ht2, window_size)
            n_unique_haplotypes = ht2.key_by('haplotype').distinct().key_by().count()
            print(f'for frequency {frequency} and window size {window_size}, found {n_unique_haplotypes} haplotypes')
            hgdp_results.append(HGDPResult(
                frequency_cutoff=frequency,
                window_size=window_size,
                hgdp_haplotype_count=n_unique_haplotypes))

        gnomad_count = va.filter(hl.max(va.pop_freqs.map(lambda x: hl.max(x.AF))) >= frequency).count()
        print(f'for frequency {frequency}, found {gnomad_count} gnomAD variants')
        gnomad_results.append(GnomADResult(
            frequency_cutoff=frequency,
            gnomad_variant_count=gnomad_count))

    with open(f'{output_base}.hgdp.tsv', 'w') as f:
        f.write('frequency\twindow_size\thgdp_haplotype_count\n')
        for result in hgdp_results:
            f.write(f'{result.frequency_cutoff}\t{result.window_size}\t{result.hgdp_haplotype_count}\n')

    with open(f'{output_base}.gnomad.tsv', 'w') as f:
        f.write('frequency\tgnomad_variant_count\n')
        for result in gnomad_results:
            f.write(f'{result.frequency_cutoff}\t{result.gnomad_variant_count}\n')


if __name__ == "__main__":
    app()
