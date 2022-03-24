import sys

input_file = snakemake.input[0]
output_file = snakemake.output[0]
chromosome = snakemake.wildcards["chr"]

fhandle = open(input_file, "r")
fout = open(output_file, "w")
for line in fhandle:
    line = line.replace("\n","")
    chromosome = line.split("\t")[123]
    if (chromosome == chromosome) | (chromosome == "CHROM"):
            fout.write(line + "\n")

fhandle.close()
fout.close()
