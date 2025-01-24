import tqdm
import typer

app = typer.Typer()


@app.command(help="Rewrite FASTA file without alt contigs (keeps only chr1-22, X, Y, MT)")
def main(
        fasta_path: str = typer.Argument(default=..., help="Path to the FASTA file"),
        output_path: str = typer.Argument(default=..., help="Output path"),
):
    contigs_to_keep = {f'chr{i}' for i in range(1, 23)} | {'chrX', 'chrY', 'chrMT'}

    with open(fasta_path) as f, open(output_path, 'w') as out:
        header_line = None
        for line in tqdm.tqdm(f):
            if line[0] == '>':
                contig = line.split()[0][1:]
                keep = contig in contigs_to_keep
                if keep:
                    out.write(line)
            else:
                if keep:
                    out.write(line)


if __name__ == "__main__":
    app()
