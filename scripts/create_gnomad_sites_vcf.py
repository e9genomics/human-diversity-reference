import hail as hl
import typer

app = typer.Typer()


@app.command()
def main(
        path: str = typer.Argument(..., help="Path to the gnomAD sites table"),
        vcf_path: str = typer.Argument(..., help="Path to the VCF file output"),
        min_popmax: float = typer.Option(..., help="Minimum gnomAD popmax frequency"),
):
    hl.init()

    # Load the gnomAD sites table
    ht = hl.read_table(path)

    filt = ht.filter(hl.max(ht.pop_freqs[0].map(lambda x: x.AF)) >= min_popmax)
    pops = ht.pops.collect()[0]

    filt = filt.annotate(info=hl.struct(
        **{
            f"{pop}_{fd}": filt.pop_freqs[i][fd] for i, pop in enumerate(pops) for fd in filt.pop_freqs[i]
        }
    ))
    hl.export_vcf(filt, vcf_path)


if __name__ == "__main__":
    app()
