import typer 
from pathlib import Path
app = typer.Typer()

def read_fastq(filepath):
    with open(filepath,"r") as file:
        lines = [line.strip() for line in file.readlines()]

    batch = list(zip(*[iter(lines)] * 4))

    seq_info = [(tups[1], tups[3]) for tups in batch]

    return seq_info

def read_count(seq_info):
    read_count = len(seq_info)

    return read_count

def gc_count(seq_info): 
    sequence = [seq[0] for seq in seq_info]
    gc_count = [round((seq.count("G")+seq.count("C"))/len(seq) * 100, 2) for seq in sequence]

    return sum(gc_count)/len(gc_count)

def mean_quality(seq_info):
    phred = [tups[1] for tups in seq_info]
    quality = []

    for asci in phred:
        for qual in asci:
            quality.append(ord(qual) -33 )
    return sum(quality)/ len(quality)


def parse_vcf(filepath):
    with open(filepath, "r") as file:
        vcf_lines = [line.strip() for line in file.readlines()]
        variants = [line.split() for line in vcf_lines if not line.startswith("#")]
    
    return variants


def filter_by_qual(variants, min_qual):
    filtered_qual = [lists for lists in variants if float(lists[5]) > min_qual]

    return filtered_qual

@app.command()
def fastq_stats(filepath: Path = typer.Argument(..., exists=True, help="Path to FASTQ file")):
    seq_info = read_fastq(filepath)
    count = read_count(seq_info)
    gc = gc_count(seq_info)
    quality = mean_quality(seq_info)

    print(f"Total number of reads: {count}")
    print(f"Mean GC percentage : {gc}")
    print(f"Mean quality : {quality}")

@app.command()
def vcf_filter(filepath: str , min_qual : float = 30.0):
    variants = parse_vcf(filepath)
    filtered = filter_by_qual(variants, min_qual)

    for variant in filtered:
        print("\t".join(variant))


if __name__ == "__main__":
    app()

