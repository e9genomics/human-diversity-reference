#!/usr/bin/env python

import hail as hl
import typer

app = typer.Typer(pretty_exceptions_enable=False)


def to_hashable_items(d):
    return tuple(sorted(d.items()))


@app.command()
def main(
    vcfs_path: str = typer.Option(default=..., help="vcfs_path"),
    gnomad_va_file: str = typer.Option(
        default=..., help="gnomAD computed variant frequencies"
    ),
    gnomad_sa_file: str = typer.Option(default=..., help="gnomAD HGDP sample metadata"),
    window_size: int = typer.Option(default=..., help="Base window size"),
    freq_threshold: float = typer.Option(
        default=..., help="Frequency threshold for keeping haplotypes"
    ),
    output_base: str = typer.Option(default=..., help="Output base path"),
    temp_dir: str = typer.Option(default="/tmp", help="Temporary directory"),
):
    """
    Process VCF files with gnomAD annotations and output filtered results.
    """
    # Initialize Hail
    hl.init()
    #
    gnomad_sa = hl.read_table(gnomad_sa_file)
    gnomad_va = hl.read_table(gnomad_va_file)
    gnomad_va = gnomad_va.filter(
        hl.max(gnomad_va.pop_freqs.map(lambda x: x.AF)) >= freq_threshold
    )
    mt = hl.import_vcf(vcfs_path, reference_genome="GRCh38", min_partitions=64)
    mt = mt.select_rows().select_cols()
    mt = mt.annotate_rows(freq=gnomad_va[mt.row_key].pop_freqs)
    mt = mt.filter_rows(hl.is_defined(mt.freq))

    pop_legend = gnomad_va.globals.pops.collect()[0]
    pop_ints = {pop: i for i, pop in enumerate(pop_legend)}
    mt = mt.annotate_cols(pop_int=hl.literal(pop_ints).get(gnomad_sa[mt.col_key].pop))
    mt = mt.filter_cols(hl.is_defined(mt.pop_int))

    mt = mt.add_row_index().add_col_index()

    mt.rows().describe()

    # don't consider sites where the variant is not above a frequency threshold in gnomAD
    mt = mt.filter_entries(mt.freq[mt.pop_int].AF >= freq_threshold)

    mt = mt.annotate_rows(
        pops_and_ids_left=hl.agg.filter(
            mt.GT[0] != 0, hl.agg.collect(hl.struct(pop=mt.pop_int, sample=mt.col_idx))
        ),
        pops_and_ids_right=hl.agg.filter(
            mt.GT[1] != 0, hl.agg.collect(hl.struct(pop=mt.pop_int, sample=mt.col_idx))
        ),
        frequencies_by_pop=hl.agg.group_by(mt.pop_int, hl.agg.call_stats(mt.GT, 2)),
    )
    ht = mt.rows().select(
        "freq",
        "pops_and_ids_left",
        "pops_and_ids_right",
        "row_idx",
        "frequencies_by_pop",
    )

    # ht = ht.checkpoint(f"{temp_dir}/xht.ht", overwrite=True)

    def get_haplotypes(ht, windower_f, idx):
        new_locus = windower_f(ht.locus)
        ht = ht.annotate(new_locus=new_locus)

        def agg_haplos(arr):
            flat = hl.agg.explode(
                lambda elt: hl.agg.collect(elt.annotate(row_idx=ht.row_idx)), arr
            )
            pop_grouped = hl.group_by(lambda x: x.pop, flat)
            pop_to_haps_array = pop_grouped.map_values(
                lambda arr_per_pop: hl.array(
                    hl.array(
                        hl.group_by(lambda inner_elt: inner_elt.sample, arr_per_pop)
                    )
                    # remove single variants
                    .filter(
                        lambda sample_and_records: hl.len(sample_and_records[1]) > 1
                    )
                    # after this map, the inner array per pop is a list of observed haplotypes: array<array<int64>>
                    .map(
                        lambda sample_and_records: hl.sorted(
                            sample_and_records[1].map(lambda e: e.row_idx)
                        )
                    )
                    # group by the identity now
                    .group_by(lambda x: x)
                    .map_values(lambda arr: hl.len(arr))
                )
            )

            return pop_to_haps_array


        ht_grouped = ht.group_by("new_locus").aggregate(
            row_map=hl.dict(
                hl.agg.collect((ht.row_idx, ht.row.select("locus", "alleles", "freq", "frequencies_by_pop")))
            ),
            left_haplos=agg_haplos(ht.pops_and_ids_left),
            right_haplos=agg_haplos(ht.pops_and_ids_right),
        )
        # left_haplos and right_haplos are array<(pop_id, array<(sample_id, array<row_idx>)>)>

        # operates on inner array
        def collapse_haplos_across_samples(pop, arr1, arr2):
            # assumes all AN == 2 * N_samples
            flat = hl.array([arr1, arr2]).flatmap(lambda x: x.get(pop))
            
            def map_haplo_group(t):
                # t is (haplotype, array<(haplotype, n)>)
                haplotype = t[0]
                n_observed = hl.sum(t[1].map(lambda x: x[1]))
                component_variant_frequencies = haplotype.map(
                    lambda x: ht_grouped.row_map[x].frequencies_by_pop[pop]
                )
                min_AN = hl.min(component_variant_frequencies.map(lambda x: x.AN))
                return hl.struct(
                    haplotype=haplotype,
                    pop=pop,
                    empirical_AC=n_observed,
                    min_variant_frequency=hl.min(
                        component_variant_frequencies.map(lambda x: x.AF[1])
                    ),
                    empirical_AF=n_observed / min_AN,
                )

            return hl.array(hl.group_by(lambda x: x[0], flat)).map(map_haplo_group)

        # compute pop => (haplotype, freq) filtered by freq_threshold
        ht_grouped = ht_grouped.annotate(
            all_haplos=hl.literal(list(pop_ints.values())).flatmap(
                lambda pop: collapse_haplos_across_samples(pop, ht_grouped.left_haplos, ht_grouped.right_haplos)))

        def get_haplotype_summary(a):
            # a is an array of struct(haplotype, pop, frequency)
            a_sorted = hl.sorted(a, key=lambda x: x.empirical_AF, reverse=True)
            return dict(
                max_pop=a_sorted[0].pop,
                max_empirical_AF=a_sorted[0].empirical_AF,
                max_empirical_AC=a_sorted[0].empirical_AC,
                min_variant_frequency=a_sorted[0].min_variant_frequency,
                all_pop_freqs=a_sorted.map(lambda x: x.drop("haplotype")),
            )

        ht_grouped = ht_grouped.transmute(
            all_haplos=hl.array(
                hl.group_by(lambda x: x.haplotype, ht_grouped.all_haplos)
            ).map(lambda t: hl.struct(haplotype=t[0], **get_haplotype_summary(t[1])))
        )

        hte = ht_grouped.explode("all_haplos")
        hte = hte.key_by().drop("new_locus")
        hte.describe()

        def get_variant(row_idx):
            r = hte.row_map[row_idx]
            return r.select("locus", "alleles")

        def get_gnomad_freq(row_idx):
            r = hte.row_map[row_idx]
            return r.freq

        hte = hte.select(
            **hte.all_haplos,
            variants=hte.all_haplos.haplotype.map(lambda row_idx: get_variant(row_idx)),
            gnomad_freqs=hte.all_haplos.haplotype.map(
                lambda row_idx: get_gnomad_freq(row_idx)
            ),
        )

        # Select the haplotype with the highest empirical_AF across all populations
        hte = hte.group_by("haplotype").aggregate(**hl.sorted(hl.agg.collect(hte.row.drop('haplotype')), key=lambda row: -row.max_empirical_AF)[0])

        typer.echo(f"Writing {output_base}.{idx}.ht...")
        return hte.checkpoint(f"{output_base}.{idx}.ht", overwrite=True)

    # compute staggered windows
    window1 = get_haplotypes(
        ht, lambda locus: locus - (locus.position % window_size), 1
    )
    window2 = get_haplotypes(
        ht, lambda locus: locus - ((locus.position + window_size // 2) % window_size), 2
    )

    htu = window1.union(window2)
    htu.describe()
    typer.echo(f"Writing final {output_base}.ht...")

    htu.key_by("haplotype").naive_coalesce(64).write(
        f"{output_base}.ht", overwrite=True
    )


if __name__ == "__main__":
    app()
