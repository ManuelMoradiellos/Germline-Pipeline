## GERMLINE PIPELINE SNAKEMAKE FILE

# Here you will find the main rules for the TCGA germline pipeline (PCAs not included).
# Scripts are placed in scrips/
# The other smk files are in smk/
# Neccesary files can be found in files/

# Warning messages obtained from tdido's bollito code.
class ansitxt:
    RED = '\033[31m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

def warning(msg):
    print(f"\n{ansitxt.BOLD}{ansitxt.RED}{msg}{ansitxt.ENDC}\n",file=sys.stderr)

try:
    configfile: "config.yaml"
except WorkflowError:
    warning("ERROR: config.yaml does not exist or is incorrectly formatted. Please see the README file for details. Quitting now.")
    sys.exit(1)


# Define parameters.
OUTDIR = config["outdir"]
CHRS = config["CHRS"]
common_threshold = config["common_threshold"]
rare_threshold = config["rare_threshold"]
print(rare_threshold)

TCGA_lowqual_samples = config["TCGA_lowqual_samples"]
doPCA = config["doPCA"]
pathogenic = config["pathogenic"]
pext_dir = config["PEXT"]["PEXT_dir"]

rule all:
    input:
        expand("{OUTDIR}rare_{rare_threshold}/{patho_type}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.VEP.{patho_type}.rare_{rare_threshold}.header.vcf",  OUTDIR = OUTDIR, rare_threshold = rare_threshold, patho_type = pathogenic),
        expand("{OUTDIR}rare_{rare_threshold}/{patho_type}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.PEXT.{patho_type}.2nd.rare_{rare_threshold}.header.tsv",  OUTDIR = OUTDIR, rare_threshold = rare_threshold, patho_type = pathogenic)


#rule divide_by_chr:
#    input: f"{OUTDIR}germline.header.tsv"
#    output: f"{OUTDIR}germline.header.{{chr}}.tsv"
#    resources:
#        mem_mb=8000,
#        walltime=90
#    script: "script/divide_chr.py"

rule rename:
    input: f"/storage/scratch01/users/lgjimeno/tcga_pipeline/final_mapping/out/final_mapping/{{chr}}.tsv"
    output: temp(f"{OUTDIR}germline.header.{{chr}}.tsv")
    resources:
        mem_mb=4000,
        walltime=90
    shell: "cp {input} {output}"

rule first_filter:
    input: f"{OUTDIR}germline.header.{{chr}}.tsv"
    output: temp(f"{OUTDIR}germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.header.{{chr}}.tsv")
    params: 
        TCGA_lowqual_samples = config["TCGA_lowqual_samples"]
    resources:
        mem_mb=4000,
        walltime=90
    script: "script/sampleQC_AD_PD_FUNC_CONS_VAF.py"



# FLAG pass // PASS-PASS / LowQual-PASS / PASS-. /  TCGA-gnomAD
rule FLAGpass:
    input: f"{OUTDIR}germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.header.{{chr}}.tsv"
    output: f"{OUTDIR}germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.{{chr}}.tsv"
    resources:
        mem_mb=28000,
        walltime=480
    script: "script/filterFLAG.py"
    
"""
rule paste_header:
    input: f"{OUTDIR}germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.{{chr}}.tsv"
    output: f"{OUTDIR}germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.{{chr}}.tsv"
    resources:
        mem_mb=4000,
        walltime=90
    shell: f"cat {OUTDIR}/header.tsv {{input}} > {{output}}"
"""

rule filter_SNV:
    input: f"{OUTDIR}germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.{{chr}}.tsv"
    output: f"{OUTDIR}germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.SNV.header.{{chr}}.tsv"
    resources:
        mem_mb=4000,
        walltime=90
    script: "script/filter_SNVs.R"


rule get_AFs:
    input: f"{OUTDIR}germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.header.{{chr}}.tsv"
    output: f"{OUTDIR}germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.header.{{chr}}.tsv"
    params:
        pop_info = config["pop_info"],
        pop_counts = config["pop_counts"]
    resources:
        mem_mb=240000,
        walltime=180
    script: "script/getTCGA_AFs.R"


rule filter_by_function:
    input: f"{OUTDIR}germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.header.{{chr}}.tsv"
    output: temp(f"{OUTDIR}germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.header.{{chr}}.tsv")
    params:
        func_refgene_vector = config["functional_filter"]
    resources:
        mem_mb=160000,
        walltime=180
    script: "script/filter_by_function.R"

rule filter_rare:
    input: f"{OUTDIR}germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.header.{{chr}}.tsv"
    output: f"{OUTDIR}rare_{{rare_threshold}}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.rare_{{rare_threshold}}.header.{{chr}}.tsv"
    params:
        rare_threshold = f"{{rare_threshold}}"
    resources:
        mem_mb=160000,
        walltime=180
    script: "script/filter_rare.R"


rule filter_pathogenic:
    input: f"{OUTDIR}rare_{{rare_threshold}}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.rare_{{rare_threshold}}.header.{{chr}}.tsv"
    output: f"{OUTDIR}rare_{{rare_threshold}}/{{patho_type}}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.{{patho_type}}.rare_{{rare_threshold}}.header.{{chr}}.tsv"
    params:
        patho_value = f"{{patho_type}}",
        delmis_method = config["delmis_method"]
    resources:
        mem_mb=24000,
        walltime=180
    script: "script/filter_pathogenic.R"

rule filter_pathogenic_2nd:
    input: f"{OUTDIR}rare_{{rare_threshold}}/{{patho_type}}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.{{patho_type}}.rare_{{rare_threshold}}.header.{{chr}}.tsv"
    output: f"{OUTDIR}rare_{{rare_threshold}}/{{patho_type}}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.{{patho_type}}.2nd.rare_{{rare_threshold}}.header.{{chr}}.tsv"
    params:
        patho_value = f"{{patho_type}}",
        VEP_file = config["VEP_pLOF"]
    resources:
        mem_mb=24000,
        walltime=180
    conda: "envs/R.yaml"
    script: "script/filter_pathogenic2.R"


rule merge_chrs_1st:
    input: expand("{OUTDIR}rare_{{rare_threshold}}/{{patho_type}}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.{{patho_type}}.rare_{{rare_threshold}}.header.{chr}.tsv", chr = CHRS, OUTDIR = OUTDIR)
    output: f"{OUTDIR}rare_{{rare_threshold}}/{{patho_type}}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.{{patho_type}}.rare_{{rare_threshold}}.header.tsv"
    resources:
        mem_mb=16000,
        walltime=180
    shell: "head -n 1 {input[1]} > {output}; tail -n +2 -q {input} >> {output}"

rule merge_chrs_2nd:
    input: expand("{OUTDIR}rare_{{rare_threshold}}/{{patho_type}}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.{{patho_type}}.2nd.rare_{{rare_threshold}}.header.{chr}.tsv", chr = CHRS, OUTDIR = OUTDIR)
    output: f"{OUTDIR}rare_{{rare_threshold}}/{{patho_type}}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.{{patho_type}}.2nd.rare_{{rare_threshold}}.header.tsv"
    resources:
        mem_mb=16000,
        walltime=180
    shell: "head -n 1 {input[1]} > {output}; tail -n +2 -q {input} >> {output}"


rule filter_PEXT_step1:
    input: f"{OUTDIR}rare_{{rare_threshold}}/{{patho_type}}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.{{patho_type}}.2nd.rare_{{rare_threshold}}.header.{{chr}}.tsv"
    output: temp(f"{OUTDIR}rare_{{rare_threshold}}/{{patho_type}}/germline.PEXT.target_variants.{{patho_type}}.2nd.rare_{{rare_threshold}}.header.{{chr}}.txt")
    resources:
        mem_mb=4000,
        walltime=20
    shell:
        "cut -f 1 {input} | sort | uniq > {output}; echo 'chrom_pos_ref_alt_.' >> {output}"

rule filter_PEXT_step2:
    input: f"{pext_dir}PEXT.{{chr}}.header.txt"
    output: f"{pext_dir}PEXT.ID.{{chr}}.header.txt"
    resources:
        mem_mb=4000,
        walltime=20
    shell:
        "awk -f script/get_PEXT_ID.awk {input} > {output}"


rule filter_PEXT_step3:
    input: 
        target_variants = f"{OUTDIR}rare_{{rare_threshold}}/{{patho_type}}/germline.PEXT.target_variants.{{patho_type}}.2nd.rare_{{rare_threshold}}.header.{{chr}}.txt",
        PEXT_file = f"{pext_dir}PEXT.ID.{{chr}}.header.txt"
    output: 
        PEXT_file_output = temp(f"{pext_dir}PEXT.filtered.{{patho_type}}.2nd.rare_{{rare_threshold}}.{{chr}}.header.txt")
    resources:
        mem_mb=4000,
        walltime=20
    shell:
        "grep -f {input.target_variants} {input.PEXT_file} > {output.PEXT_file_output}"


rule filter_PEXT_step4:
    input:
        patho_input = f"{OUTDIR}rare_{{rare_threshold}}/{{patho_type}}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.{{patho_type}}.2nd.rare_{{rare_threshold}}.header.{{chr}}.tsv",
        PEXT_file_input = f"{pext_dir}PEXT.filtered.{{patho_type}}.2nd.rare_{{rare_threshold}}.{{chr}}.header.txt"
    output: f"{OUTDIR}rare_{{rare_threshold}}/{{patho_type}}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.PEXT.{{patho_type}}.2nd.rare_{{rare_threshold}}.header.{{chr}}.tsv"
    params:
        PEXT_expression_threshold = config["PEXT"]["PEXT_expression_threshold"],
        PEXT_tissue_equivalences = config["PEXT"]["PEXT_equivalences"],
        TCGA_metadata = config["TCGA_metadata"]
    resources:
        mem_mb=4000,
        walltime=20
    script:
        "script/PEXT_filter.R"


rule merge_chrs_PEXT:
    input: expand("{OUTDIR}rare_{{rare_threshold}}/{{patho_type}}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.PEXT.{{patho_type}}.2nd.rare_{{rare_threshold}}.header.{chr}.tsv", chr = CHRS, OUTDIR = OUTDIR)
    output: f"{OUTDIR}rare_{{rare_threshold}}/{{patho_type}}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.PEXT.{{patho_type}}.2nd.rare_{{rare_threshold}}.header.tsv"
    resources:
        mem_mb=4000,
        walltime=180
    shell: "head -n 1 {input[1]} > {output}; tail -n +2 -q {input} >> {output}"




rule to_vep:
    input: f"{OUTDIR}rare_{{rare_threshold}}/{{patho_type}}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.{{patho_type}}.rare_{{rare_threshold}}.header.tsv"
    output: f"{OUTDIR}rare_{{rare_threshold}}/{{patho_type}}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.nohead.{{patho_type}}.rare_{{rare_threshold}}.header.vcf"
    params:
        rare_threshold = f"{{rare_threshold}}"
    resources:
        mem_mb=8000,
        walltime=180
    script: "script/to_VCF.R"

rule to_vep_2:
    input: f"{OUTDIR}rare_{{rare_threshold}}/{{patho_type}}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.nohead.{{patho_type}}.rare_{{rare_threshold}}.header.vcf"
    output: f"{OUTDIR}rare_{{rare_threshold}}/{{patho_type}}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.{{patho_type}}.rare_{{rare_threshold}}.header.vcf"
    params:
        rare_threshold = f"{{rare_threshold}}"
    resources:
        mem_mb=8000,
        walltime=180
    conda: 
        "envs/bcftools.yaml"
    shell: 
        "cat files/header_vcf.txt {input} > {output}"

rule to_vep_3:
    input: f"{OUTDIR}rare_{{rare_threshold}}/{{patho_type}}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.{{patho_type}}.rare_{{rare_threshold}}.header.vcf"
    output: f"{OUTDIR}rare_{{rare_threshold}}/{{patho_type}}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.sorted.{{patho_type}}.rare_{{rare_threshold}}.header.vcf"
    params:
        rare_threshold = f"{{rare_threshold}}"
    resources:
        mem_mb=8000,
        walltime=180
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools sort {input} > {output}"

include: "smk/vep.smk"
#include: "smk/N_variants.smk"
#include: "smk/PCA_appoach.smk"
