#!/usr/bin/env python

import hail as hl
import typer

app = typer.Typer()


def to_hashable_items(d):
    return tuple(sorted(d.items()))


@app.command()
def main(
        vcfs_path: str = typer.Option(default=..., help="vcfs_path"),
        gnomad_va_file: str = typer.Option(default=..., help="gnomAD computed variant frequencies"),
        gnomad_sa_file: str = typer.Option(default=..., help="gnomAD HGDP sample metadata"),
        output_ht: str = typer.Option(default=..., help="Output path"),
):
    """
    Process VCF files with gnomAD annotations and output filtered results.
    """
    # Initialize Hail
    hl.init()
    #
    gnomad_sa = hl.read_table(gnomad_sa_file)
    gnomad_va = hl.read_table(gnomad_va_file)
    mt = hl.import_vcf(vcfs_path, reference_genome='GRCh38', min_partitions=64)
    mt = mt.select_rows().select_cols()
    mt = mt.annotate_rows(freq=gnomad_va[mt.row_key].pop_freqs)
    mt = mt.filter_rows(hl.is_defined(mt.freq))
    mt = mt.annotate_rows(max_pop_freq=hl.max(mt.freq.map(lambda x: hl.max(x.AF))))

    freq_thresholds = [0, 0.0001, 0.001, 0.005, 0.01, 0.5, 0.1]
    from pprint import pprint
    pprint(mt.aggregate_rows(hl.struct(
        **{f'n_sites_above_{x}': hl.agg.count_where(mt.max_pop_freq > x) for x in freq_thresholds},
    )))

    mt = mt.annotate_cols(pop=gnomad_sa[mt.col_key].pop)
    mt = mt.annotate_cols(counts=hl.struct(
        **{f'n_sites_above_{x}':
               hl.agg.count_where(mt.GT.is_non_ref() & (mt.max_pop_freq > x)) for x in freq_thresholds},
    ))

    mt.cols().write(output_ht, overwrite=True)

if __name__ == "__main__":
    app()
