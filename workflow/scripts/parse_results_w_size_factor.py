# This script summarizes the species level
# counts for each of the short read datasets.
# The script outputs two tables: one with the
# raw counts and one with the a library size
# adjustment. We've used counts per million
# (CPM) to adjust for differences in sequencing
# depth (samples as rows, species as columns)
# The script is  meant to run on the output of the
# genome mapping portion of the CAO pipeline

# Conda env: envs/sum_gen_map.yml


import os
import re
import pandas as pd
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("--counts", required=True, nargs="+", help="List of counts files")
parser.add_argument("--logs", required=True, nargs="+", help="List of fastp logs")
parser.add_argument("--raw_counts", required=True, type=str, help="Raw counts output")
parser.add_argument(
    "--size_adjusted", required=True, type=str, help="Size adjusted counts"
)
parser.add_argument("--lib_sizes", required=True, type=str, help="Library sizes")
args = parser.parse_args()

# Paths to count results and fastp logs (for testing)
# res_dir = "result/"
# fastp_logs = "logs/fastp/"

# On Uppmax
# res_dir = '/crex/proj/uppstore2017149/BioDiversityAnalysisCAO/LTS-BiodiversityAnalysisCAO/results/genome_mappings/target_species/counts/'
# fastp_logs = '/crex/proj/uppstore2017149/BioDiversityAnalysisCAO/LTS-BiodiversityAnalysisCAO/results/logs/fastp/'

# Empty list for sample names
dfidx = []

# Empty list for dictionaries (each has the
# result per sample)
data = []
# dirs = os.listdir(res_dir)
for file in args.counts:
    if file.endswith("species_counts.txt"):
        sample = re.sub("_species_counts.txt", "", os.path.basename(file))
        dfidx.append(sample)
        # Go through each results file, adding
        # the species and counts to a dictionary
        # for the sample
        fh = open(file, "r")
        dic = {}
        for line in fh:
            line = line.strip().split(",")
            sp = str(line[0])
            count = int(line[1])
            dic[sp] = count
        # Add dictionary to list of dictionaries
        data.append(dic)
        fh.close()


# Make a dataframe using the list of sample
# names as the index and the list of
# dictionaries
df = pd.DataFrame(data, index=dfidx)

# Convert 'NA' to zero, then convert all
# to integers
df = df.fillna(0)
df = df.astype(int)

# Write to file
df.to_csv(args.raw_counts, sep=",")


############################################
# Adjust for library size differences. Use
# counts per million (CPM) and the read counts
# after fastp filtering.

# Go through each sample file. Add the name of the
# sample as the dictionary key and the library read
# count as value to lib_sizes dictionary.

lib_sizes = {}
# dir2s = os.listdir(fastp_logs)
for file in args.logs:
    if file.endswith(".fastp.log"):
        sample = re.sub(".fastp.log", "", os.path.basename(file))
        fh = open(file, "r")
        for line in fh:
            if line.startswith("reads passed filter:"):
                l = line.strip().split(":")
                c = int(l[1].strip())
                lib_sizes[sample] = c
        fh.close()
# Make lib_sizes dictionary into df
sizes_df = pd.DataFrame.from_dict(lib_sizes, orient="index", columns=["lib_size"])

# list of species columns
sp_cols = list(df.columns)

# Leftjoin then get counts per million
df_sz = df.merge(sizes_df, how="left", left_index=True, right_index=True)

# Get counts per million, then write to file
cpm_df = (df_sz[sp_cols].div(df_sz["lib_size"], axis="index")) * 1000000
cpm_df.to_csv(args.size_adjusted, sep=",")

# Write lib sizes to file
sizes_df.to_csv(args.lib_sizes, sep=",")
