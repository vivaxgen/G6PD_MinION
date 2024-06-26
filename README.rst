G6PD Amplicon Sequencing Pipeline
=================================

This is ONT-based amplicon sequencing pipeline for detecting and diagnosing of human G6PD variants.
The pipeline leverages on `vivaxGEN ngs-pipeline <https://github.com/vivaxgen/ngs-pipeline>`_ to provide most of its functionality. 


Quick Installation
------------------

To install the pipeline and all of its dependencies (including reference sequence), run the following command from your terminal/shell::

    "${SHELL}" <(curl -L https://raw.githubusercontent.com/vivaxgen/G6PD_MinION/main/install.sh)

Make sure that the installation completes sucessfully.
Take a note of the activation script once the installation process finished.


Quick Start
-----------

* Source the activation script that was mentioned after the installation process finished, such as::

    source SOME_DIRECTORY/bin/activate.sh

* Go to the parent directory of where the sample files reside.
  The sample files must be in compressed FASTQ format (fastq.gz).
  Note that the sample file names (without the extension) would be taken as sample code, eg. my-sample-1.fastq.gz would yield to sample code of my-sample-1. Symbolic links are also fine.

* To run the panel variant calling, which will check only known variants as
  listed in 
  `g6pd-variant-info.tsv <https://github.com/vivaxgen/G6PD_MinION/blob/main/refs/g6pd-variant-info.tsv>`_,
  execute the following command (assuming the compressed FASTQ files are under
  my_data directory and the output be written to my_output directory)::

    ngs-pl run-g6pd-panel-variant-caller -o my_output my_data/*.fastq.gz

  The list of variant of each sample is written in merged_report.tsv file
  (as well as merged_report.xlsx) in the output directory.
  The following symbols are used to indicate the mutation of each variant:

  ===== =====================================================================
   \-   The interested variant is not found, i.e., if the interested variant is a mutation, the allele in this position is still the same as reference, **or** if the interested variant is the reference, then an alternate allele is called.
   \+   The interested variant is found, i.e., if the interested variant is a mutation, the allele in this position is the alternate base, **or** if the interested variant is the reference, allele in this position is same as the reference.
   -/+  Heterozygote mutation, both the interested allele and non-interested allele (which could be reference or another alternate allele) is at this position
   ?    No known state of mutation, either because low quality position or no reads were mapped
  ===== =====================================================================

* To run the discovery variant calling, execute the following command::

    ngs-pl run-g6pd-discovery-variant-caller -o my_output my_data/*.fastq.gz

  Note that the default is using mpileup variant calling.
  To use different variant calller, such as freebayes or clair3, use the --caller argument::

    ngs-pl run-g6pd-discovery-variant-caller --caller freebayes -o my_output my_data/*.fastq.gz

* To increase the number of parallel processes of the analysis, use -j option, eg in a properly setup HPC system::

    ngs-pl run-g6pd-discovery-variant-caller -j 512 -o my_output my_data/*.fastq.gz


Introduction to the pipeline
----------------------------

[to be written]


Detailed Usage
--------------

[to be written]


