#!/usr/bin/env python
import hail as hl
import typer

app = typer.Typer()

GNOMAD_SITES_TABLE_PATH = 'gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_variant_annotations.ht'
GNOMAD_HGDP_SAMPLE_DATA_PATH = 'gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.ht'
HGDP_REPRESENTED_POPS = ['afr', 'amr', 'eas', 'sas', 'nfe']


def to_hashable_items(d):
    return tuple(sorted(d.items()))


@app.command()
def main(
        output_file_va: str = typer.Argument(..., help="Output table path for variant frequencies"),
        output_file_sa: str = typer.Argument(..., help="Output table path for sample metadata"),
        freq_threshold: float = typer.Option(0.001, help="Frequency threshold for filtering variants")
):
    """
    Process VCF files with gnomAD annotations and output filtered results.
    """
    # Initialize Hail
    hl.init()

    va = hl.read_table(GNOMAD_SITES_TABLE_PATH)

    freq_meta = va.globals.gnomad_freq_meta[:15].collect()[0]

    map_to_index = {to_hashable_items(x): i for i, x in enumerate(freq_meta)}

    pop_indices = []
    for pop in HGDP_REPRESENTED_POPS:
        idx = map_to_index.get(to_hashable_items({'group': 'adj', 'pop': pop}))
        if idx is None:
            raise ValueError(f"Population {pop} not found in gnomAD frequency metadata")
        pop_indices.append(idx)

    # some empty filters are {}, some are NA
    va = va.filter(hl.coalesce(hl.len(va.filters) == 0, True))

    # pick out only the pop freqs
    va = va.select_globals(pops=HGDP_REPRESENTED_POPS)
    va = va.select(pop_freqs=hl.literal(pop_indices).map(lambda i: va.gnomad_freq[i]))

    # drop sites under `freq_threshold` in all pops
    va = va.filter(hl.any(lambda x: x.AF > freq_threshold, va.pop_freqs))

    va.naive_coalesce(64).write(output_file_va, overwrite=True)

    sa = hl.read_table(GNOMAD_HGDP_SAMPLE_DATA_PATH).select_globals()
    sa = sa.select(pop=sa.gnomad_population_inference.pop)
    sa.naive_coalesce(4).write(output_file_sa, overwrite=True)


if __name__ == "__main__":
    app()
