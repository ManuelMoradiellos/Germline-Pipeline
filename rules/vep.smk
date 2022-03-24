rule get_vep_cache:
    output:
        directory("vep/cache")  # if not working try resources/vep/cache
    params:
        species="homo_sapiens",
        build="GRCh37",
        release="104"
    threads: 1
    resources:
        mem=16000,
        walltime=720
    log:
        f"log/vep/cache.log"
    cache: True  # save space and time with between workflow caching (see docs)
    wrapper:
        "v0.86.0/bio/vep/cache"

rule download_vep_plugins:
    output:
        directory("vep/plugins")
    params:
        release=100
    wrapper:
        "v0.86.0/bio/vep/plugins"

rule annotate_variants:
    input:
        calls=f"{OUTDIR}rare_{{rare_threshold}}/{{patho_type}}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.sorted.{{patho_type}}.rare_{{rare_threshold}}.header.vcf",  # .vcf, .vcf.gz or .bcf
        cache="vep/cache",  # can be omitted if fasta and gff are specified
        plugins="resources/vep/plugins",
        # optionally add reference genome fasta
        # fasta="genome.fasta",
        # fai="genome.fasta.fai", # fasta index
        # gff="annotation.gff",
        # csi="annotation.gff.csi", # tabix index
        # add mandatory aux-files required by some plugins if not present in the VEP plugin directory specified above.
        # aux files must be defined as following: "<plugin> = /path/to/file" where plugin must be in lowercase
        # revel = path/to/revel_scores.tsv.gz
    output:
        calls= f"{OUTDIR}rare_{{rare_threshold}}/{{patho_type}}/germline.TCGAqc_filt.GQ.DP.VAF.AD.BCD.PASS.TCGA_AF.FUNC.VEP.{{patho_type}}.rare_{{rare_threshold}}.header.vcf",  # .vcf, .vcf.gz or .bcf
        stats= f"{OUTDIR}rare_{{rare_threshold}}/{{patho_type}}/stats.{{patho_type}}.rare_{{rare_threshold}}.header.html"
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=["LoFtool", "SpliceAI", "LOFTEE"],
        extra="--everything",  # optional: extra arguments
    log:
        f"logs/vep/annotate.{{patho_type}}.rare_{{rare_threshold}}.log"
    threads: 8
    resources:
        mem_mb=16000,
        walltime=1400
    conda:
        "../envs/vep.yaml"
    script:
        "../script/vep_ann.py"

