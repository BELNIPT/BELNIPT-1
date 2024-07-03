#!/usr/bin/env nextflow

// Script parameters
params.genome_dir = "genome.fa"

log.info """\
Using genome directory : $params.genome_dir
"""

process mapping {

  conda '/conda/envs/mapping'
  publishDir "$params.genome_dir/bowtie2/v2.4.1", mode: 'copy'

  input:
    path genome from params.genome_dir
    
  output:
    file 'genome*' into genome_index

  script:
  """
  cp $genome/genome.fa .
  bowtie2-build genome.fa genome
  """
}
