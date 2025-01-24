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
        window_size: int = typer.Option(default=..., help="Base window size"),
        freq_threshold: float = typer.Option(default=..., help="Frequency threshold for keeping haplotypes"),
        output_base: str = typer.Option(default=..., help="Output base path"),
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

    pop_legend = gnomad_va.globals.pops.collect()[0]
    pop_ints = {pop: i for i, pop in enumerate(pop_legend)}
    mt = mt.annotate_cols(pop_int=hl.literal(pop_ints).get(gnomad_sa[mt.col_key].pop))
    mt = mt.filter_cols(hl.is_defined(mt.pop_int))

    mt = mt.add_row_index().add_col_index()

    mt = mt.annotate_rows(
        pops_and_ids_left=hl.agg.filter(
            mt.GT[0] != 0, hl.agg.collect(hl.struct(pop=mt.pop_int, sample=mt.col_idx))),
        pops_and_ids_right=hl.agg.filter(
            mt.GT[1] != 0, hl.agg.collect(hl.struct(pop=mt.pop_int, sample=mt.col_idx))))
    ht = mt.rows().select('freq', 'pops_and_ids_left', 'pops_and_ids_right', 'row_idx')

    ht = ht.checkpoint('/tmp/xht.ht', overwrite=True)

    def get_haplotypes(ht, windower_f, idx):
        new_locus = windower_f(ht.locus)
        ht = ht.annotate(new_locus=new_locus)

        def agg_haplos(arr):
            flat = hl.agg.explode(lambda elt: hl.agg.collect(elt.annotate(row_idx=ht.row_idx)), arr)
            pop_grouped = hl.group_by(lambda x: x.pop, flat)
            pop_to_haps_array = pop_grouped.map_values(
                lambda arr_per_pop: hl.array(
                    hl.array(hl.group_by(
                        lambda inner_elt: inner_elt.sample,
                        arr_per_pop))
                    # remove single variants
                    .filter(lambda sample_and_records: hl.len(sample_and_records[1]) > 1)
                    # after this map, the inner array per pop is a list of observed haplotypes: array<array<int64>>
                    .map(
                        lambda sample_and_records:
                        hl.sorted(sample_and_records[1].map(lambda e: e.row_idx)))
                    # group by the identity now
                    .group_by(lambda x: x)
                    .map_values(lambda arr: hl.len(arr)))
            )

            return pop_to_haps_array

        ht_grouped = ht.group_by('new_locus').aggregate(
            row_map=hl.dict(hl.agg.collect((ht.row_idx, ht.row.select('locus', 'alleles', 'freq')))),
            left_haplos=agg_haplos(ht.pops_and_ids_left),
            right_haplos=agg_haplos(ht.pops_and_ids_right)
        )
        # left_haplos and right_haplos are array<(pop_id, array<(sample_id, array<row_idx>)>)>

        count_by_pop = hl.literal(mt.aggregate_cols(hl.agg.counter(mt.pop_int)))

        # operates on inner array
        def collapse_haplos_across_samples(pop, arr):
            # assumes all AN == 2 * N_samples
            freqs = arr.map(
                lambda t: hl.struct(haplotype=t[0], pop=pop, observed=t[1], frequency=(t[1] / count_by_pop[pop] / 2)))
            freqs = freqs.filter(lambda x: x.frequency > freq_threshold)
            return hl.array(freqs)

        # compute pop => (haplotype, freq) filtered by freq_threshold
        ht_grouped = ht_grouped.transmute(
            all_haplos=hl.literal(list(pop_ints.values())).flatmap(
                lambda pop: collapse_haplos_across_samples(pop, hl.array([ht_grouped.left_haplos,
                                                                          ht_grouped.right_haplos]).flatmap(
                    lambda d: d.get(pop)))))

        def get_haplotype_summary(a):
            # a is an array of struct(haplotype, pop, frequency)
            a_sorted = hl.sorted(a, key=lambda x: x.frequency, reverse=True)
            return dict(max_pop=a_sorted[0].pop,
                        max_pop_freq=a_sorted[0].frequency,
                        max_pop_observed=a_sorted[0].observed,
                        all_pop_freqs=a_sorted.map(lambda x: x.drop('haplotype')))

        ht_grouped = ht_grouped.transmute(
            all_haplos=hl.array(hl.group_by(lambda x: x.haplotype, ht_grouped.all_haplos))
            .map(lambda t: hl.struct(haplotype=t[0], **get_haplotype_summary(t[1])))
        )

        hte = ht_grouped.explode('all_haplos')
        hte = hte.key_by().drop('new_locus')

        def get_variant(row_idx):
            r = hte.row_map[row_idx]
            return r.select('locus', 'alleles')

        def get_gnomad_freq(row_idx):
            r = hte.row_map[row_idx]
            return r.freq

        hte = hte.select(**hte.all_haplos,
                         variants=hte.all_haplos.haplotype.map(lambda row_idx: get_variant(row_idx)),
                         gnomad_freqs=hte.all_haplos.haplotype.map(lambda row_idx: get_gnomad_freq(row_idx))
                         )

        typer.echo(f"Writing {output_base}.{idx}.ht...")
        return hte.checkpoint(f"{output_base}.{idx}.ht", overwrite=True)

    # compute staggered windows
    window1 = get_haplotypes(ht, lambda locus: locus - (locus.position % window_size), 1)
    window2 = get_haplotypes(ht, lambda locus: locus - ((locus.position + window_size // 2) % window_size), 2)

    htu = window1.union(window2)
    htu.describe()
    typer.echo(f"Writing final {output_base}.ht...")

    htu.key_by('haplotype').distinct().naive_coalesce(64).write(f"{output_base}.ht", overwrite=True)


if __name__ == "__main__":
    app()
