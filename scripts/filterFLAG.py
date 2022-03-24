fhandle = open(snakemake.input[0])
fout = open(snakemake.output[0], 'w')
i = 1
for line in fhandle:
    line = line.replace("\n","")
    line_split = line.split("\t")
    if i == 1:
        TCGA_pass_index = line_split.index("FILTER_ID")
        gnomAD_pass_index = line_split.index("FILTER_gnomad")
        gnomAD_wgs_pass_index = line_split.index("FLAG_gnomad_WGS")
        GQ_index = line_split.index("GQ")
        fout.write(line + "\n")

        print(TCGA_pass_index)
        print(gnomAD_pass_index)
        print(gnomAD_wgs_pass_index)
        print(GQ_index)


    else:
        TCGA_pass = line_split[TCGA_pass_index]
        gnomAD_pass = line_split[gnomAD_pass_index]
        gnomAD_wgs_pass = line_split[gnomAD_wgs_pass_index]
        GQ = int(line_split[GQ_index])

        if ((gnomAD_pass == "PASS") and (TCGA_pass == "PASS") and (GQ >= 20)):
            fout.write(line + "\n")
        if ((TCGA_pass == "PASS") and (gnomAD_pass == "") and (gnomAD_wgs_pass in ["PASS", ""]) and (GQ >= 20)):
            fout.write(line + "\n")
    
    i += 1
fhandle.close()
fout.close()
