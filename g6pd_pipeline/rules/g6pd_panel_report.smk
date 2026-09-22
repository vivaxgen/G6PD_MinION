from ngs_pipeline.rules import pkg

include: pkg("ngs_pipeline::msf/panel_varcall_lr.smk")

# rule full_report:
#     input:
#         f"{outdir}/merged_genetic_report.tsv"

# rule merge_g6pd_report:
#     input:
#         expand(f'{outdir}/samples/{{sample}}/genetic_report.tsv', sample=read_files.samples())
#     output:
#         f'{outdir}/merged_genetic_report.tsv',
#     run:
#         import pandas as pd
#         dfs = []
#         for infile in input:
#             df = pd.read_table(infile, sep="\t")
#             dfs.append(df)
#         merged_df = pd.concat(dfs, ignore_index=True, axis=0)
#         merged_df.to_csv(output[0], sep="\t", index=False)

rule check_multiple_missense:
    input:
        bam = f"{outdir}/samples/{{sample}}/maps/mapped-final.bam",
        idx = f"{outdir}/samples/{{sample}}/maps/mapped-final.bam.bai",
        phase_info = get_abspath(config.get('phase_info')),
    output:
        tsv = f"{outdir}/samples/{{sample}}/vcfs/multiple_missense_report.tsv"
    log:
        f"{outdir}/samples/{{sample}}/logs/multiple_missense_report.log"
    params:
        min_qual = config.get("min_read_qual", 10),
        min_mapq = config.get("min_mapq", 30)
    shell:
        """
        ngs-pl construct-pseudo-haplotypes --min_qual {params.min_qual} \
        --min_mapq {params.min_mapq} --variants_list {input.phase_info} \
        --single --sample {wildcards.sample} --log {log} -o {output.tsv} {input.bam}
        """

rule gen_g6pd_report:
    threads: 1
    input:
        vcf = f"{outdir}/samples/{{sample}}/vcfs/variants.vcf.gz",
        vcf_idx = f"{outdir}/samples/{{sample}}/vcfs/variants.vcf.gz.tbi",
        missenses = f"{outdir}/samples/{{sample}}/vcfs/multiple_missense_report.tsv",
        variant_info = variant_info
    output:
        tsv = f"{outdir}/samples/{{sample}}/genetic_report.tsv"
    params:
        min_var_qual = config.get('min_variant_qual', 10),
        min_depth = config.get('report_calling_mindepth', 20),
    shell:
        """
        ngs-pl generate-g6pd-panel-report --infofile {input.variant_info} \
        --outfile {output.tsv} \
        --mindepth {params.min_depth} --min_var_qual {params.min_var_qual} \
        --multiple_missenses {input.missenses} {input.vcf}
        """

ruleorder: gen_g6pd_report > gen_report