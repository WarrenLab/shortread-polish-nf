# shortread-polish-nf
Nextflow pipeline for polishing an assembly with short reads and freebayes

## Introduction
When using error-prone long reads (e.g., PacBio CLR) to assemble a genome, it is
often necessary to use short reads (e.g., Illumina) to error-correct ("polish")
the assembly. This repository contains a [nextflow](http://nextflow.io) pipeline
implementing the [Vertebrate Genomes Project][vgp] best practices for using
freebayes to polish an assembly (see [their github][vgp-github] for more
information).

## Requirements
### Data
* Short reads
* An assembly created from erroneous long reads

### Nextflow
* [nextflow](https://www.nextflow.io/) — can be installed with the command
  `curl -s https://get.nextflow.io | bash`

### Other software
* [minimap2](https://github.com/lh3/minimap2)
* [freebayes](https://github.com/freebayes/freebayes)
* [samtools](http://www.htslib.org/)
* [bcftools](http://www.htslib.org/)

This pipeline is set up to use mamba to create an environment with these
programs in it, but you can always install them yourself or use modules. See
configuration section below.

## Configuration
Nextflow configuration is handled by the file `nextflow.config` in the directory
where you're running the nextflow command. The configuration file in this
repository is for running the pipeline on the lewis cluster at Mizzou using
SLURM and mamba, but you can adjust it to use any batch or cloud system you
want. Check out the [nextflow docs](https://www.nextflow.io/docs/edge/index.html)
for more information.

### Lewis-specific stuff
Nextflow needs a filesystem where locking is allowed for keeping track of which
jobs are running, but not to actually store the data or temporary files you're
creating. On Lewis, HTC allows locking but is slow and HPC does not allow
locking but is fast. To take advantage of the best of both worlds, run this
pipeline from within a project directory in HTC, but set the environment
variable `$NXF_WORK` to point to an empty directory on HPC.

## Running
To download the pipeline and run it on your assembly, just run the command:
```
nextflow run WarrenLab/shortread-polish-nf \
    --reference unpolished_assembly.fa \
    --sra SRX1234567
```
This will download short reads from SRA, align all the reads to your reference,
and then use the alignments to correct the assembly.

## TODO
* Write help message

[vgp]: https://vertebrategenomesproject.org/
[vgp-github]: https://github.com/VGP/vgp-assembly/tree/master/pipeline/freebayes-polish
