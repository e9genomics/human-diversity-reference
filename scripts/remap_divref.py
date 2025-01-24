# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "duckdb",
#     "numpy",
#     "pandas",
#     "pydantic",
#     "tqdm",
#     "typer",
# ]
# ///
from pathlib import Path
from typing import Optional

import duckdb
import numpy as np
import pandas as pd
import typer
from pydantic import BaseModel
from tqdm import tqdm

app = typer.Typer()

POPULATION_LEGEND = {'afr': 0, 'amr': 1, 'eas': 2, 'sas': 3, 'nfe': 4}


class Variant(BaseModel):
    chromosome: str
    position: int
    reference: str
    alternate: str

    def render(self):
        return f"{self.chromosome}:{self.position}:{self.reference}:{self.alternate}"


def intervals_overlap(start1, end1, start2, end2):
    return start1 < end2 and start2 < end1


class ReferenceMapping(BaseModel):
    chromosome: str
    start: int
    end: int
    min_gnomad_popmax: Optional[float]
    min_gnomad_popmax_pop: Optional[int]
    variants_involved: list[Variant]
    first_variant_index: Optional[int]
    last_variant_index: Optional[int]

    def variants_involved_str(self):
        return ','.join([v.render() for v in self.variants_involved])


class Haplotype(BaseModel):
    sequence_id: str
    sequence: str
    sequence_length: int
    n_variants: int
    variants: str
    gnomad_AF: str

    _variants: Optional[list[Variant]] = None

    def parsed_variants(self) -> list[Variant]:
        if self._variants is not None:
            return self._variants
        vs = []
        for v_str in self.variants.split(','):
            chrom, pos, ref, alt = v_str.strip().split(':')
            vs.append(Variant(chromosome=chrom, position=int(pos), reference=ref, alternate=alt))
        self._variants = vs
        return vs

    def parse_gnomad_AF(self) -> list[list[float]]:
        return [list(map(lambda x: 0 if x == 'null' else float(x), x.split('|'))) for x in self.gnomad_AF.split(',')]

    def contig(self):
        variants = self.parsed_variants()
        return variants[0].chromosome

    def reference_mapping(self, start: int, end: int, context_size: int) -> ReferenceMapping:
        vs = self.parsed_variants()

        # translate a locus position into an index in the string
        # as examples:
        #   2-6 for 1:500:AAA:T with context window 0 should be 502-503
        #   2-6 for 1:500:AAA:T with context window 2 should be 500-501

        # translate variants into [start, end) intervals in 0-indexed haplotype sequence space
        variant_intervals = []

        # update index_translation based on variant size as we go
        index_translation = vs[0].position - context_size
        for i, v in enumerate(vs):
            v_start = v.position - index_translation
            v_end = v_start + len(v.alternate)
            index_translation += len(v.reference) - len(v.alternate)
            variant_intervals.append((v_start, v_end))

        first_variant_index = None
        last_variant_index = None
        for i, (v_start, v_end) in enumerate(variant_intervals):
            if intervals_overlap(start, end, v_start, v_end):
                if first_variant_index == None:
                    first_variant_index = i
                last_variant_index = i

        def translate_coordinate_to_ref(coord: int, sign: int) -> int:
            # Either the coordinate is contained within a variant interval, or it isn't
            # if it is, (1) return the start of the variant if sign<0 or end of variant if sign>0
            # if it isn't, (2) return add the distance from the previous variant interval end to that variant's position
            # unless (3) the coordinate is before than the first variant, in which case we translate from the first variant's position

            first_variant_start = variant_intervals[0][0]
            if coord < first_variant_start:
                # path (3)
                return vs[0].position - (first_variant_start - coord)

            last_smaller_variant = 0
            for i, (v_start, v_end) in enumerate(variant_intervals):

                # if contained in an interval, path (1)
                if v_start <= coord < v_end:

                    # if the coordinate is contained in a variant interval, return the start or end based on the sign
                    if sign < 0:
                        return vs[i].position
                    else:
                        return vs[i].position + len(vs[i].reference)

                if v_start > coord:
                    break

                last_smaller_variant = i

            # if we're here, we know that the coordinate is not contained in any variant interval
            v = vs[last_smaller_variant]
            v_end = v.position + len(v.reference)
            return v_end + (coord - variant_intervals[last_smaller_variant][1])

        reference_coord_start = translate_coordinate_to_ref(start, -1)
        reference_coord_end = translate_coordinate_to_ref(end, 1)

        gnomad_afs = self.parse_gnomad_AF()

        min_gnomad_popmax = None
        min_gnomad_popmax_pop = None
        if first_variant_index is not None:
            for i in range(first_variant_index, last_variant_index + 1):
                local_min_popmax_pop = np.argmax(gnomad_afs[i])
                local_min_gnomad_popmax = gnomad_afs[i][local_min_popmax_pop]
                if min_gnomad_popmax is None or local_min_gnomad_popmax < min_gnomad_popmax:
                    min_gnomad_popmax = local_min_gnomad_popmax
                    min_gnomad_popmax_pop = local_min_popmax_pop

        rm = ReferenceMapping(
            chromosome=self.contig(),
            start=reference_coord_start,
            end=reference_coord_end,
            min_gnomad_popmax=min_gnomad_popmax,
            min_gnomad_popmax_pop=min_gnomad_popmax_pop,
            variants_involved=vs[first_variant_index:last_variant_index + 1] if first_variant_index is not None else [],
            first_variant_index=first_variant_index,
            last_variant_index=last_variant_index,
        )

        return rm


def get_index_path(path: Optional[Path]):
    # if None, look in the same directory as this file
    if path is None:
        import os
        # walk files in the same directory as this script
        for root, dirs, files in os.walk(Path(__file__).parent):
            for file in files:
                if file.endswith(".duckdb"):
                    path = Path(root) / file
                    break
    conn = duckdb.connect(path)
    return conn


@app.command(name='calitas')
def calitas(
        input_path: Path = typer.Argument(..., help="Path to the CALITAS output file"),
        output_path: Path = typer.Argument(..., help="Path to the remapped output file"),
        index_path: Optional[Path] = typer.Option(None, '-i', help="Path to the FASTA index file"),
        sep: str = typer.Option('\t', '-s', help="Separator in the file"),
        batch_size: int = typer.Option(25000, '-b', help="Number of rows to process in each batch"),
):
    conn = get_index_path(index_path)

    df = pd.read_csv(input_path, sep=sep)
    chrom_field = 'chromosome'
    start_field = 'coordinate_start'
    end_field = 'coordinate_end'

    # if any of these fields are absent, error
    if not all(x in df.columns for x in (chrom_field, start_field, end_field)):
        raise ValueError(f"Required fields not found in the input file: {chrom_field}, {start_field}, {end_field}")

    if df[chrom_field].dtype != str:
        df[chrom_field] = df[chrom_field].astype(str)

    window_size: int = conn.execute(f"SELECT * FROM window_size").fetchone()[0]
    print(f"Remapping DivRev coordinates, using window size {window_size}")
    pops_legend = conn.execute(f"SELECT * FROM pops_legend").fetchone()[0]

    # Lists to store results
    contigs = []
    starts = []
    ends = []
    variants_involved = []
    n_variants_involved = []
    min_gnomad_popmax = []
    min_gnomad_popmax_pop = []

    # Process in batches
    for batch_start in tqdm(range(0, len(df), batch_size)):
        batch_end = min(batch_start + batch_size, len(df))

        batch_df = df.iloc[batch_start:batch_end]
        batch_hap_ids = batch_df[chrom_field].tolist()

        # Query database for batch
        results = conn.execute("""
            SELECT * FROM haplotypes 
            WHERE haplotypes.sequence_id IN (SELECT unnest($1::STRING[]))
            """, [batch_hap_ids]).fetchall()

        columns = [desc[0] for desc in conn.description]
        id_to_hap = {}
        for row in results:
            hap = Haplotype(**dict(zip(columns, row)))
            id_to_hap[hap.sequence_id] = hap

        # Process each row in the batch
        for _, df_row in batch_df.iterrows():
            start = df_row[start_field]
            end = df_row[end_field]
            hap_id = df_row[chrom_field]

            strand = df_row['strand']
            padded_target = df_row['padded_target']
            target = df_row['unpadded_target_sequence']

            padded_len_adj = len(padded_target) - len(target)

            # account for PAM in protospacer sequence
            if strand == '+':
                end += padded_len_adj
            else:
                start -= padded_len_adj

            hap = id_to_hap[hap_id]
            rm = hap.reference_mapping(start, end, window_size)

            # Append results to lists
            contigs.append(rm.chromosome)
            starts.append(rm.start)
            ends.append(rm.end)
            variants_involved.append(rm.variants_involved_str())
            n_variants_involved.append(len(rm.variants_involved))
            min_gnomad_popmax.append(rm.min_gnomad_popmax)
            min_gnomad_popmax_pop.append(
                pops_legend[rm.min_gnomad_popmax_pop] if rm.min_gnomad_popmax_pop is not None else None
            )

    # Update DataFrame with results, maintaining original structure
    df['divref_chromosome'] = df[chrom_field]
    df['divref_start'] = df[start_field]
    df['divref_end'] = df[end_field]
    df[chrom_field] = contigs
    df[start_field] = starts
    df[end_field] = ends
    df['genome_build'] = 'GRCh38'
    df['variants_involved'] = variants_involved
    df['n_variants_involved'] = n_variants_involved
    df['min_gnomad_popmax'] = min_gnomad_popmax
    df['min_gnomad_popmax_pop'] = min_gnomad_popmax_pop

    df.to_csv(output_path, sep=sep, index=False)

if __name__ == "__main__":
    app()
