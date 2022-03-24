input_file = snakemake.input["file"]
output_file = snakemake.output[0]
outliers = snakemake.input["outliers"]

fhandle = open(input_file)
fout = open(output_file, "w")

bad_sample_file = open(outliers, "r")
bad_samples = bad_sample_file.readlines()
bad_samples_stripped = [s.strip() for s in bad_samples]  

print(bad_samples_stripped)

for line in fhandle:
    line = line.replace("\n","")
    sample = line.split("\t")[138]
    print(sample)
    if sample in bad_samples_stripped:
            continue
    else:
            fout.write(line + "\n")

fhandle.close()
fout.close()
