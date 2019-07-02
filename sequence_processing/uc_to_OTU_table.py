#!/usr/bin/env python

#Matti Ruuskanen Jan 2018
#Makes an OTU table from a clustered *.uc file.
#Possibility to remove OTUs from the OTU table with less than a specified number of reads.
#If given a matching fasta file of cluster seeds, renames the sequences and removes
#OTUs with less than the user specified minimum number of reads also from the fasta file.

#check Python version
import sys
if sys.version_info[0] < 3:
    raise "Must be using Python 3"

#import needed packages
import numpy as np
import pandas as pd
import argparse
import os
import re

#parse arguments given to the script
parser = argparse.ArgumentParser()
parser.add_argument("-c", "--clusters", dest = "clusters", help = "cluster file name")
parser.add_argument("-o", "--out", dest = "out", default = os.getcwd(), help = "output folder [default = current working directory]")
parser.add_argument("-f", "--fasta", dest = "fasta", help = "fasta file name containing cluster seeds")
parser.add_argument("-m", "--min", dest = "min", default = 0, help = "minimum number of reads in OTU to be included [default = disabled]", type = int)
args = parser.parse_args()

#stop execution if no input file is given
if args.clusters is None:
    print("At least one argument must be supplied (input file).")
    parser.print_help()
    sys.exit(1)

#change execution location to the specified (or default) output directory
os.chdir(args.out)

#read in the input .uc file and remove useless columns
df = pd.read_table(args.clusters, header = None)
df = df.drop(df.columns[[5,6]], axis = 1)
df.ix[:,1] = "OTU_" + df.ix[:,1].astype(str)


#remove clusters with less than specified number of counts
cSizes = df.loc[df[0] == "C"]
cSizes = cSizes.iloc[:, [1,2,6]]
cSizes.columns = ["OTU", "count", "fastaheader"]
if args.min > 0:
    toRemove = cSizes[cSizes["count"] <= args.min]
    df = df.loc[df[0] != "C"]
    df = df.iloc[:, [1, 6]]
    df.columns = ["OTU", "Sample"]
    df = df[~df["OTU"].isin(toRemove["OTU"])]
    print("Removed " + str(len(toRemove)) + " OTUs with a read count of " + str(args.min) + " or below.\n")
else:
    df = df.loc[df[0] != "C"]
    df = df.iloc[:, [1, 6]]
    df.columns = ["OTU", "Sample"]

df = df.sort_values(by=["OTU"])
df.loc[:, "Sample"] = [re.sub("(.*)_.*", "\\1", x) for x in df["Sample"]]
df = df.groupby(["Sample", "OTU"]).size().reset_index()
nReads = df[0].sum()
df = df.pivot(index="Sample", columns="OTU", values=0)
df[np.isnan(df)] = 0
nOTUs = len(df.columns)
nSamples = len(df.index)

print("Saved a table of " + str(nSamples) + " samples with " + str(nReads) + " counts of " + str(nOTUs) + " OTUs.\n")

df.to_csv(path_or_buf="otu_table.txt", sep="\t")

if args.fasta != None:
    from Bio import SeqIO
    record_dict = SeqIO.index(args.fasta, "fasta")
    headers = pd.DataFrame({"orig" : [key for key in record_dict.keys()], "fastaheader" : [key.split(';')[0] for key in record_dict.keys()]})
    toRemove = cSizes[cSizes["count"] <= args.min]
    toKeep = headers[~headers["fastaheader"].isin(toRemove["fastaheader"])]
    toKeep = pd.merge(toKeep, cSizes, on="fastaheader") 
    fileName = str(re.sub("(.+)\\.f.*", "\\1", args.fasta)) + "_renamed.fasta"
    with open(fileName, "w") as out_handle:
        i = 0
        for id in toKeep["orig"]:
            record = record_dict[id]
            record.id = toKeep.ix[i,2]
            record.description = toKeep.ix[i,2]
            SeqIO.write(record, out_handle, "fasta")
            i += 1
