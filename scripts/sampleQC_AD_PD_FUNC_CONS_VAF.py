#!/usr/bin/python
import sys
import re
import csv 

fhandle = open(snakemake.input[0])
fout = open(snakemake.output[0], 'w')

bad_sample_file = open(snakemake.params["TCGA_lowqual_samples"], "r")
bad_samples = bad_sample_file.readlines()
bad_samples_stripped = [s.strip() for s in bad_samples]  

with open(snakemake.input[0], 'r') as infile:
    reader = csv.DictReader(infile, delimiter = "\t")
    header = reader.fieldnames

sample_index = header.index("SAMPLE")
ID_index = header.index("Otherinfo3")
AD_index = header.index("AD")
DP_index = header.index("DP")


i = 1
for line in fhandle:
    line = line.replace("\n","")
    if i == 1:
      fout.write(line + "\t" +  "VAF" + "\n")
    else:
      sample = line.split("\t")[sample_index][0:15]
      ID = line.split("\t")[ID_index]
      AD = line.split("\t")[AD_index].split(",")
      DP = line.split("\t")[DP_index]
      sample_type_index = sample[13]
      print(sample_type_index)
 
      if ((sample in bad_samples_stripped) or (sample_type_index != "1")):
        continue    
      else:   
        if (len(AD) == 2) & (DP != "."):
          VAF = float(AD[1])/(float(AD[0]) + float(AD[1]))
          if (VAF >= 0.1) & (float(DP) >= 10):
              VAFchr = str(VAF)
              print(AD)
              fout.write(line + "\t" +  VAFchr + "\n")
    i +=1
    print(i)
    
fhandle.close()
fout.close()
