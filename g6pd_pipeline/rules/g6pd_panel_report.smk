from ngs_pipeline.rules import pkg

include: pkg("ngs_pipeline::msf/panel_varcall_lr.smk")
include: "set_variant_gt.smk"

## Required files for configs:
# Hard fail if not found
amplicon_bed_file = get_abspath(config.get("amplicon_bed"))
phase_info_file = get_abspath(config.get('phase_info'))

rule full_details_report:
    input:
        f"{outdir}/merged_genetic_report.pre.tsv",
        f"{outdir}/amplicon_coverage.tsv"

rule merge_coverage_report:
    input:
        tsv = expand(f"{outdir}/samples/{{sample}}/amplicon_coverage.tsv",
                     sample=read_files.samples()),
    output:
        tsv = f"{outdir}/amplicon_coverage.tsv",
    run:
        import pandas as pd
        dfs = [pd.read_table(f) for f in input.tsv]
        merged_df = pd.concat(dfs, ignore_index=True)
        pivoted = merged_df.pivot(index="sample", columns="Amplicon_name", values=["coverage", "meandepth"])
        pivoted = pivoted.swaplevel(axis=1).sort_index(axis=1, level=0)
        pivoted.columns = [f"{amp}_{stat}" for amp, stat in pivoted.columns]
        pivoted = pivoted.reset_index()
        pivoted.to_csv(output.tsv, sep="\t", index=False)
        
rule amplicon_coverage_report:
    input:
        bam = f"{outdir}/samples/{{sample}}/maps/mapped-final.bam",
        bai = f"{outdir}/samples/{{sample}}/maps/mapped-final.bam.bai",
        amplicon_target = amplicon_bed_file
    output:
        depth_coverage = f"{outdir}/samples/{{sample}}/amplicon_coverage.tsv",
    run:
        import pandas as pd
        import numpy as np
        from io import StringIO
        markers = pd.read_table(input.amplicon_target, header=None)
        if markers.shape[1] < 4:
            markers.loc[:, 3] = markers.apply(lambda x: f"{x[0]}:{x[1]}-{x[2]}", axis=1)
        markers.columns = ["Chr", "Start", "End", "Amplicon_name"]
        all_results = []

        rows = []
        for chrom, g in markers.groupby("Chr", sort=False):
            starts = g["Start"].to_numpy()
            ends = g["End"].to_numpy()
            names = g["Amplicon_name"].to_numpy()
            idx = g.index.to_numpy()
            # every start/end is a potential breakpoint; between consecutive
            # breakpoints, coverage (how many markers span that stretch) is constant
            breakpoints = np.unique(np.concatenate([starts, ends]))
            for seg_start, seg_end in zip(breakpoints[:-1], breakpoints[1:]):
                covering = np.where((starts <= seg_start) & (ends >= seg_end))[0]
                if covering.size == 1:          # keep only stretches with exactly 1 marker
                    i = covering[0]
                    rows.append((chrom, seg_start, seg_end, names[i], idx[i]))
        frag = pd.DataFrame(
            rows, columns=["Chr", "Start", "End", "Amplicon_name", "_marker_idx"]
        )
        # re-merge fragments that are still contiguous (same marker, no gap),
        # so a marker only gets split when a real overlap carved it in two
        frag = frag.sort_values(["_marker_idx", "Start"]).reset_index(drop=True)
        prev_end = frag.groupby("_marker_idx")["End"].shift()
        frag["_group"] = (frag["Start"] != prev_end).groupby(frag["_marker_idx"]).cumsum()
        markers_to_test = (
            frag.groupby(["_marker_idx", "_group"], as_index=False)
                .agg(Chr=("Chr", "first"),
                     Start=("Start", "min"),
                     End=("End", "max"),
                     Amplicon_name=("Amplicon_name", "first"))
                .drop(columns=["_marker_idx", "_group"])
                .sort_values(["Chr", "Start"])
                .reset_index(drop=True)
        )

        markers_to_test["region"] = markers_to_test["Chr"] + ":" + markers_to_test["Start"].astype(str) + "-" + markers_to_test["End"].astype(str)
        for marker in markers_to_test["region"]:
            temp = pd.read_table(StringIO(shell(f"samtools coverage -H -r {marker} {input.bam}", read= True)), header=None, names = ["rname", "startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq"])
            temp["sample"] = wildcards.sample
            temp["region"] = marker
            all_results.append(temp)
        all_results = pd.concat(all_results)
        full_result = markers_to_test.merge(all_results, left_on="region", right_on="region", how="outer").drop("region", axis=1)
        full_result.to_csv(output.depth_coverage, sep="\t", index=False)
        

if config.get("neg_control", None):
    for neg_sample in config["neg_control"]:
        if neg_sample not in read_files.samples():
            raise ValueError(f"Negative control sample '{neg_sample}' not found in input samples.")
            exit(1)

    checkpoint determine_highest_negative_depth:
        input:
            expand(f"{outdir}/samples/{{sample}}/genetic_report.pre.full_details.tsv", sample=config["neg_control"])
        output:
            f"{outdir}/highest_negative_depth.txt"
        params:
            current_mindepth = config.get('report_calling_mindepth', 10),
            required_min_multiplication = 1.25,
        run:
            import pandas as pd
            import math

            highest_depth = max(max(
                pd.read_table(report, usecols=["DP"])["DP"].max()
                for report in input
            )*params.required_min_multiplication, params.current_mindepth)
            with open(output[0], "w") as depth_file:
                depth_file.write(f"{math.ceil(highest_depth)}\n")
else:
    checkpoint null_highest_negative_depth:
        output:
            f"{outdir}/highest_negative_depth.txt"
        params:
            current_mindepth = config.get('report_calling_mindepth', 10)
        run:
            with open(output[0], "w") as depth_file:
                depth_file.write(f"{params.current_mindepth}\n")


rule check_multiple_missense:
    input:
        bam = f"{outdir}/samples/{{sample}}/maps/mapped-final.bam",
        idx = f"{outdir}/samples/{{sample}}/maps/mapped-final.bam.bai",
        phase_info = phase_info_file,
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

use rule merge_vcfs as merge_vcfs_gt_set with:
    input:
        vcfs = expand(f"{outdir}/samples/{{sample}}/vcfs/variants.setgt.vcf.gz",
                      sample=read_files.samples()),
        idx = expand(f"{outdir}/samples/{{sample}}/vcfs/variants.setgt.vcf.gz.csi",
                     sample=read_files.samples()),
    output:
        vcf = f"{outdir}/merged_variants.setgt.vcf.gz",

use rule merge_vcfs as merge_norm_vcfs_gt_set with:
    input:
        vcfs = expand(f"{outdir}/samples/{{sample}}/vcfs/variants.setgt.norm.vcf.gz",
                      sample=read_files.samples()),
        idx = expand(f"{outdir}/samples/{{sample}}/vcfs/variants.setgt.norm.vcf.gz.csi",
                     sample=read_files.samples()),
    output:
        vcf = f"{outdir}/merged_variants.setgt.norm.vcf.gz",

rule link_final_vcf:
    localrule: True
    input:
        vcf = f"{outdir}/merged_variants.setgt.norm.vcf.gz",
    output:
        vcf = f"{outdir}/merged.vcf.gz",
    shell:
        """
        ln -srf {input.vcf} {output.vcf}
        """

rule link_final_report:
    localrule: True
    input:
        tsv = f"{outdir}/samples/{{sample}}/genetic_report.norm.tsv"
    output:
        tsv = f"{outdir}/samples/{{sample}}/genetic_report.tsv"
    shell:
        """
        ln -srf {input.tsv} {output.tsv}
        """

rule gen_g6pd_report:
    threads: 1
    input:
        vcf = f"{outdir}/samples/{{sample}}/vcfs/variants.setgt.vcf.gz",
        vcf_idx = f"{outdir}/samples/{{sample}}/vcfs/variants.setgt.vcf.gz.tbi",
        missenses = f"{outdir}/samples/{{sample}}/vcfs/multiple_missense_report.tsv",
        variant_info = variant_info,
        negative_depth = [],
    output:
        multiext(f"{outdir}/samples/{{sample}}/genetic_report.pre", tsv=".tsv", tsv_full=".full_details.tsv")
    params:
        min_var_qual = config.get('min_variant_qual', 10),
        min_depth = config.get('report_calling_mindepth', 10),
    shell:
        """
        ngs-pl generate-g6pd-panel-report --infofile {input.variant_info} \
        --outfile {output.tsv} \
        --mindepth {params.min_depth} --min_var_qual {params.min_var_qual} \
        --multiple_missenses {input.missenses} {input.vcf}
        """

use rule gen_g6pd_report as gen_norm_g6pd_report with:
    input:
        vcf = f"{outdir}/samples/{{sample}}/vcfs/variants.setgt.norm.vcf.gz",
        vcf_idx = f"{outdir}/samples/{{sample}}/vcfs/variants.setgt.norm.vcf.gz.tbi",
        missenses = f"{outdir}/samples/{{sample}}/vcfs/multiple_missense_report.tsv",
        variant_info = variant_info,
        negative_depth = normalized_depth_file,
    output:
        multiext(f"{outdir}/samples/{{sample}}/genetic_report.norm", tsv=".tsv", tsv_full=".full_details.tsv")
    params:
        min_var_qual = config.get('min_variant_qual', 10),
        min_depth = normalized_minimum_depth,

use rule merge_report as merge_norm_report with:
    input:
        tsv = expand(f"{outdir}/samples/{{sample}}/genetic_report.tsv",
                     sample=read_files.samples()),
    output:
        tsv = f"{outdir}/merged_genetic_report.tsv",

use rule merge_report as merge_non_norm_report with:
    input:
        tsv = expand(f"{outdir}/samples/{{sample}}/genetic_report.pre.tsv",
                     sample=read_files.samples()),
    output:
        tsv = f"{outdir}/merged_genetic_report.pre.tsv",

ruleorder: merge_norm_report > merge_report
ruleorder: link_final_report > gen_report
ruleorder: link_final_vcf > merge_vcfs