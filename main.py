import pandas as pd

import argparse
import libgal2


# Create the argument parser
parser = argparse.ArgumentParser(description="Example script with long arguments")

parser.add_argument(
    "--input", type=str, default="KO_mcr_rose_union.tsv", help="Input file"
)

parser.add_argument(
    "--dir",
    type=str,
    default="/home/antony/development/data/modules/genome",
    help="Path to databases",
)

parser.add_argument(
    "--out",
    type=str,
    default="KO_mcr_rose_union_annotated.tsv",
    help="Path to the output file",
)

# Add long arguments
parser.add_argument("--genome", type=str, default="mm10", help="Genome e.g. grch38")
parser.add_argument(
    "--closest", type=int, default=5, help="How many closest genes to add"
)
parser.add_argument("--annotation", action="store_true", help="Enable annotation")
parser.add_argument(
    "--blacklist", action="store_true", help="Enable blacklist annotation"
)

parser.add_argument("--tads", action="store_true", help="Enable TAD annotation")

parser.add_argument(
    "--centromere",
    action="store_true",
    help="Enable centromere annotation",
)

parser.add_argument(
    "--telomere", action="store_true", help="Enable telomere annotation"
)

parser.add_argument(
    "--promoter-lim",
    type=int,
    nargs=2,
    default=libgal2.PROMOTER_LIM,
    help="Promoter limits",
)

parser.add_argument(
    "--max-dist", type=int, default=2000000, help="Maximum distance to consider"
)

parser.add_argument("--verbose", action="store_true", help="Verbose mode")

# Parse the command-line arguments
args, unknown = parser.parse_known_args()

input_file = args.input
output_file = args.out
dir = args.dir
db_file = f"{dir}/gtf_{args.genome}.db"
blacklist_db_file = f"{dir}/blacklist_{args.genome}.db"
centromere_db_file = f"{dir}/centromere_{args.genome}.db"
telomere_db_file = f"{dir}/telomere_{args.genome}.db"

# tad_db_file = f"/home/antony/development/data/modules/genome/tads_{args.genome}.db"
# tad_db_file = f"{dir}/tads_{args.genome}.db"
tad_db_file = f"{dir}/tads_grch38.db"


closest_n = args.closest
verbose = args.verbose
annotate_mode = True  # args.annotation

blacklist_mode = args.blacklist


centromere_mode = True  # args.centromere
telomere_mode = True  # args.telomere
tad_mode = True  # args.tads


promoter_lim = args.promoter_lim
max_distance = args.max_dist

prom_header = (
    f"Relative To Gene (prom=-{promoter_lim[0]/1000}/+{promoter_lim[1]/1000} kb)"
)


print(input_file, output_file, db_file, closest_n, verbose, annotate_mode)

df_query = pd.read_csv(
    input_file,
    sep="\t",
    header=0,
    index_col=0,
    keep_default_na=False,
)

if annotate_mode:
    annotation = libgal2.DataframeAnnotation(
        closest_n,
        max_distance,
        promoter_lim,
    )
    annotation.open(db_file, df_query)
    df_query = annotation.annotate_genes()

    if closest_n > 0:
        df_query = annotation.annotate_closest_genes(closest_n)
        df_query.to_csv(output_file, sep="\t", header=True, index=True)

    annotation.close()

if blacklist_mode:
    blacklist = libgal2.BlacklistAnnotation()
    blacklist.open(blacklist_db_file)
    blacklist.annotate_df(df_query)
    blacklist.close()

if centromere_mode:
    centromere = libgal2.CentromereAnnotation()
    centromere.open(centromere_db_file)
    centromere.annotate_df(df_query)
    centromere.close()

if telomere_mode:
    telomere = libgal2.TelomereAnnotation()
    telomere.open(telomere_db_file)
    telomere.annotate_df(df_query)
    telomere.close()

if tad_mode:
    tad_annotation = libgal2.TADAnnotation()
    tad_annotation.open(tad_db_file)
    tad_ann = tad_annotation.annotate_df(df_query)
    tad_annotation.close()


df_query.to_csv(output_file, sep="\t", header=True, index=True)

print(f"Finished writing to {output_file}.")
