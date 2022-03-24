rule get_GT_raw:
    input: f"{OUTDIR}all.chr.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.{{chr}}.tsv"
    output: temp(f"{OUTDIR}all.chr.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.GT.header.{{chr}}.tsv")
    resources:
        mem_mb=16000,
        walltime=90
    shell: "cut -f 1,138,139 {input} > {output}"


rule merge_GT_files_raw:
    input: expand("{OUTDIR}all.chr.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.GT.header.{chr}.tsv", chr = CHRS, OUTDIR = OUTDIR)
    output: temp(f"{OUTDIR}all.chr.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.GT.header.tsv")
    resources:
        mem_mb=64000,
        walltime=90
    shell: "cat {input} | sort | uniq | sort -r > {output}"

rule get_Nvariants:
    input: f"{OUTDIR}all.chr.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.GT.header.tsv"
    output: f"{OUTDIR}N_variants_df.tsv"
    resources:
        mem_mb=32000,
        walltime=90
    script: "get_Nvariants.R"
