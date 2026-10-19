# clair3 variant caller for G6PD single-sample (discovery) flow
# use as varcaller_wf: g6pd_pipeline::ssf_varcall_clair3_lr.smk

from ngs_pipeline.rules import pkg

# clair3.smk evaluates target_variants_vcf while defining its rules (at include time),
# so it must exist before the include; empty means no --vcf_fn, ie. discovery calling
target_variants_vcf = ""


include: pkg("ngs_pipeline::varcaller/clair3.smk")

rule clair3_model_file:
    # ssf flow has no rule generating reads/model-{idx}.txt (msf uses check_clair3_models)
    localrule: True
    input:
        "reads/raw-{idx}_R0.fastq.gz"
    output:
        "reads/model-{idx}.txt"
    run:
        from ngs_pipeline.cmds.check_clair3_models import read_first_line, process_header
        model = process_header(read_first_line(input[0]))
        if not model:
            raise ValueError(f"cannot determine Clair3 model from FASTQ header of {input[0]}")
        with open(output[0], "w") as f:
            f.write(model + "\n")
