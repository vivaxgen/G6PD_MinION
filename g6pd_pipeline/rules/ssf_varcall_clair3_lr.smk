from ngs_pipeline.rules import pkg

target_variants_vcf = ""
include: pkg("ngs_pipeline::varcaller/clair3.smk")

rule clair3_model_file:
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