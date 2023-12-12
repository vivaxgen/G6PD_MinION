G6PD Amplicon Sequencing Pipeline
=================================

This is ONT-based amplicon sequencing pipeline for detecting and diagnosing of human G6PD variants.


Quick Installation
------------------

To install the pipeline and all of its dependencies (including reference sequence), run the following command from your terminal/shell::

    "${SHELL}" <(curl -L https://raw.githubusercontent.com/vivaxgen/install/main/G6PD-pl.sh)

Take a note of the activation script once the installation process finished.


Quick Start
-----------

* Source the activation script that was mentioned after the installation process finished, such as::

    source SOME_DIRECTORY/bin/activate.sh

* Go to the parent directory of where the sample files reside.
  The sample files must be in compressed FASTQ format (fastq.gz).
  Note that the sample file names (without the extension) would be taken as sample code, eg. my-sample-1.fastq.gz would yield to sample code of my-sample-1. Symbolic links are also fine.

* To run the panel variant calling, execute the following command (assuming the compressed FASTQ files are under my_data directory)::

    ngs-pl run-g6pd-panel-variant-caller -o my_output my_data/*.fastq.gz

  The list of mutation of each sample is written in merged_report.tsv file in the output directory.

* To run the discovery variant calling, execute the following command::

    ngs-pl run-g6pd-discovery-variant-caller -o my_output my_data/*.fastq.gz

  Note that the default is using mpileup variant calling.
  To use different variant calller, such as clair3, use the --snakefile option to direct the correct snakefile::

    ngs-pl run-g6pd-discovery-variant-caller --snakefile discovery_clair3_varcall_lr.smk -o my_output my_data/*.fastq.gz

* To increase the number of parallel processes of the analysis, use -j option, eg in a properly setup HPC system::

    ngs-pl run-g6pd-discovery-variant-caller -j 512 -o my_output my_data/*.fastq.gz


Introduction to the pipeline
----------------------------

[to be written]


Detailed Usage
--------------

[to be written]


