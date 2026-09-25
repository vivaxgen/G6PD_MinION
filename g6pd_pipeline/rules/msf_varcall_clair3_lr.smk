import ngs_pipeline.rules
include: pkg("ngs_pipeline::varcaller/clair3.smk")

#NGS_PIPELINE_BASE = config['NGS_PIPELINE_BASE']

ruleorder: link_variants > clair3_symlink > index_tbi
ruleorder: clair3_g6pd > clair3 > index_tbi

use rule clair3 as clair3_g6pd with:
    params:
        contig = "--contig X",

rule link_variants:
    input:
        "<sp>vcfs/clair3/merge_output.vcf.gz"
    output: 
        final = "<sp>vcfs/variants.vcf.gz",
        tbi = "<sp>vcfs/variants.vcf.gz.tbi",
    shell:
        '''
        ln -srf {input} {output.final}
        tabix -f {output.final}
        '''


# EOF
