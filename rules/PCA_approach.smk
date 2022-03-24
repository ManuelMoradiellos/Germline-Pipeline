rule filter_common:
    input: f"{OUTDIR}all.chr.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.SNV.TCGA_AF.header.{{chr}}.tsv"
    output: f"{OUTDIR}all.chr.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.SNV.TCGA_AF.common_{{threshold}}.header.{{chr}}.tsv"
    params:
        threshold = f"{{threshold}}"
    resources:
        mem_mb=128000,
        walltime=180
    script: "filter_common.R"


rule get_GT:
    input: f"{OUTDIR}all.chr.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.SNV.TCGA_AF.common.header.{{chr}}.tsv"
    output: temp(f"{OUTDIR}all.chr.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.SNV.TCGA_AF.common.GT.header.{{chr}}.tsv")
    resources:
        mem_mb=16000,
        walltime=90
    shell: "cut -f 1,139,140 {input} > {output}"


rule merge_GT_files:
    input: expand("{OUTDIR}all.chr.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.SNV.TCGA_AF.common.GT.header.{chr}.tsv", chr = CHRS, OUTDIR = OUTDIR)
    output: f"{OUTDIR}all.chr.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.SNV.TCGA_AF.common.GT.header.tsv"
    resources:
        mem_mb=64000,
        walltime=90
    shell: "cat {input} | sort | uniq | sort -r > {output}"


rule get_IQR_thresholds:
    input: f"{OUTDIR}all.chr.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.SNV.TCGA_AF.common.GT.header.tsv"
    output: f"{OUTDIR}IQR_samples.txt"
    resources:
        mem_mb=64000,
        walltime=90
    script: "get_IQR_samples.R"


rule filter_IQR:
    input:
        file = f"{OUTDIR}all.chr.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.SNV.TCGA_AF.common.header.{{chr}}.tsv",
        outliers = f"{OUTDIR}IQR_samples.txt"
    output: f"{OUTDIR}all.chr.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.SNV.TCGA_AF.common.35IQR.header.{{chr}}.tsv"
    resources:
        mem_mb=64000,
        walltime=90
    script: "filter_outliers.py"


rule get_PCA:
    input: f"{OUTDIR}all.chr.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.SNV.TCGA_AF.common_{{threshold}}.header.{{chr}}.tsv"
    output:
        PCA3levs = temp(f"{OUTDIR}PCA.input.3levels.{{chr}}.common_{{threshold}}.tsv"),
        PCA2levs = temp(f"{OUTDIR}PCA.input.2levels.{{chr}}.common_{{threshold}}.tsv")
    params: 
        samples_ID = config["samples_ID"]
    resources:
        mem_mb=64000,
        walltime=90
    script: "getPCA_input.R"


rule paste_PCA2levs:
    input: expand("{OUTDIR}PCA.input.2levels.{chr}.common_{{threshold}}.tsv", chr = CHRS, OUTDIR = OUTDIR)
    output: f"{OUTDIR}PCA.input.2levels.common_{{threshold}}.tsv"
    resources:
        mem_mb=16000,
        walltime=90
    shell: "cat {OUTDIR}PCA_header.txt {input} > {output}"


rule paste_PCA3levs:
    input: expand("{OUTDIR}PCA.input.3levels.{chr}.common_{{threshold}}.tsv", chr = CHRS, OUTDIR = OUTDIR)
    output: f"{OUTDIR}PCA.input.3levels.common_{{threshold}}.tsv"
    resources:
        mem_mb=16000,
        walltime=90
    shell: "cat {OUTDIR}PCA_header.txt {input} > {output}"

