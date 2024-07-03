#!/usr/bin/env nextflow

// Script parameters
params.genome_path = "/tmp/ref"
params.samplename = "sample"
params.path = ""
params.read1 = "R1.fastq.gz"
params.read2 = "R2.fastq.gz"

log.info """\
Using genome : $params.genome_path
Sample : $params.samplename
Path : $params.path
Read1 : $params.read1
Read2 : $params.read2
"""

process mapping {

  conda '/conda/envs/mapping'
  publishDir "$params.path", mode: 'copy', pattern: '*.{err,fastq.metrics,bam}'

  input:
    path genome from params.genome_path
    path read1 from params.read1
    path read2 from params.read2
    val name from params.samplename
    
  output:
    path "${name}.bam" into bam_file
    path "${name}.err" into mapping_stats
    path "${name}.fastq.metrics" into fastq_metrics

  script:
  """
  bowtie2 --threads 16 \
		--fast-local \
		--rg-id ${name} --rg ID:${name} --rg SM:${name} --rg PL:ILLUMINA --rg CN:AZDelta --rg LB:${name}  \
		-x $genome/bowtie2/v2.4.1/genome \
		-1 $read1 -2 $read2 2> ${name}.err | \
	elprep filter /dev/stdin ${name}.bam \
		--sorting-order coordinate \
		--nr-of-threads 16 \
		--mark-duplicates \
		--mark-optical-duplicates ${name}".fastq.metrics"
  """
}

process samtools {

  conda '/conda/envs/samtools'

  input:
    path bam from bam_file
    val name from params.samplename

  output:
    path "${name}.bam.bai" into bai_file

  script:
  """
  samtools index $bam ${name}.bam.bai
  """
}

process wisecondor {

  conda '/conda/envs/wisecondor'
  publishDir "$params.path", mode: 'copy' 

  input:
    path bam from bam_file
    path bai from bai_file
    val name from params.samplename

  output:
    path "${name}.npz" into npz_file

  script:
  """
  #generate npz
  WisecondorX convert $bam ${name}.npz
  """

}

process seqff {

  conda '/conda/envs/seqff'
  publishDir "$params.path", mode: 'copy'

  input:
    path bam from bam_file
    path bai from bai_file
    val name from params.samplename

  output:
    path "${name}_seqff.csv" into seqff_csv

  script:
  """
  echo \$PWD
  workdir=\$PWD
  samtools view -S -F 0x4 -f 0x1 -F 0x400 $bam > seqff_withoutHeader.sam
  cd /seqff
  Rscript seqff.r /seqff \$workdir/seqff_withoutHeader.sam \$workdir \$workdir/${name}_seqff.csv sam
  """

}
