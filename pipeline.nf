#!/usr/bin/env nextflow

// Script parameters
params.genome_path = "/tmp/ref"
params.sample_path = "/tmp/in"
params.sample_code = "sample"
params.pipeline_name = "nipt"
params.pipeline_version = "v1.0.0"

log.info """\
Directories
-----------
Genome          : $params.genome_path
Sample path     : $params.sample_path
Sample Information
------------------
  * Sample      : $params.sample_code
Pipeline
--------
  * name        : $params.pipeline_name
  * version     : $params.pipeline_version
"""

output_path = "$params.sample_path" + "/" + "$params.pipeline_name" + "/" + "$params.pipeline_version"


process combine_fastq {

  conda '/conda/envs/mapping'

  input:
    path sample_path from params.sample_path
    val name from params.sample_code
    
  output:
    path "${name}.R1.fastq.gz" into fastq1
    path "${name}.R2.fastq.gz" into fastq2

  script:
  """
  zcat $sample_path/fastq/${name}*_R1_001.fastq.gz >> ${name}.R1.fastq
  zcat $sample_path/fastq/${name}*_R2_001.fastq.gz >> ${name}.R2.fastq
  gzip ${name}.R1.fastq
  gzip ${name}.R2.fastq
  """
}

process mapping {

  conda '/conda/envs/mapping'
  publishDir "$output_path", mode: 'copy', pattern: '*.{err,fastq.metrics}'

  input:
    path genome from params.genome_path
    path fastq1 from fastq1
    path fastq2 from fastq2
    val name from params.sample_code
    
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
		-1 $fastq1 -2 $fastq2 2> ${name}.err | \
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
    val name from params.sample_code

  output:
    path "${name}.bam.bai" into bai_file

  script:
  """
  samtools index $bam ${name}.bam.bai
  """
}

process yu_length {

  conda '/conda/envs/samtools'
  publishDir "$output_path", mode: 'copy'

  input:
    path bamfile from bam_file
    val name from params.sample_code

  output:
    path "${name}.yu_ff.tsv" into yuFF_file

  script:
  """
  samtools view -S -F 0x4 -f 0x1 -F 0x400 -q 30 $bamfile | awk 'BEGIN{long=0; short=0}{if(\$9>0){if(\$9>=100 && \$9<=150){short=short+1}if(\$9>=163 && \$9<=169){long=long+1}}}END{print "FF\t"(15.625*(short/long)-9.053)}' > $name".yu_ff.tsv" 
  """
}

process seqff_samtools {

  conda '/conda/envs/samtools'

  input:
    path bam from bam_file
    path bai from bai_file
    val name from params.sample_code

  output:
    path "seqff_withoutHeader.sam" into seqff_sam

  script:
  """
  samtools view -S -F 0x4 -f 0x1 -F 0x400 $bam > seqff_withoutHeader.sam
  """
}

process seqff {

  conda '/conda/envs/seqff'
  publishDir "$output_path", mode: 'copy'

  input:
    path bam from bam_file
    path bai from bai_file
    val name from params.sample_code
    path seqff_sam from seqff_sam

  output:
    path "${name}_seqff.csv" into seqff_csv

  script:
  """
  echo \$PWD
  workdir=\$PWD
  # samtools view -S -F 0x4 -f 0x1 -F 0x400 $bam > seqff_withoutHeader.sam
  cd /seqff
  # Rscript seqff.r /seqff \$workdir/seqff_withoutHeader.sam \$workdir \$workdir/${name}_seqff.csv sam
  Rscript seqff.r /seqff \$workdir/$seqff_sam \$workdir \$workdir/${name}_seqff.csv sam
  """
}

process wisecondor_npz {

  conda '/conda/envs/wisecondor'
  publishDir "$output_path", mode: 'copy'

  input:
    path genome from params.genome_path
    path bam from bam_file
    path bai from bai_file
    val name from params.sample_code
    val pipeline_name from params.pipeline_name
    val pipeline_version from params.pipeline_version

  output:
    path "${name}.npz" into npz_file

  script:
  """
  #generate npz
  WisecondorX convert $bam ${name}.npz
  """
}

process wisecondor_gender {

  conda '/conda/envs/wisecondor'
  publishDir "$output_path", mode: 'copy'

  input:
    path genome from params.genome_path
    val name from params.sample_code
    path npz_file from npz_file
    val pipeline_name from params.pipeline_name
    val pipeline_version from params.pipeline_version

  output:
    path "${name}_gender.txt" into wisecondor_gender

  script:
  """
  #Gender
  WisecondorX gender ${npz_file} ${genome}/$pipeline_name/$pipeline_version/reference_output.npz > ${name}"_gender.txt"
  """
}


process wisecondor_predict {

  conda '/conda/envs/wisecondor'
  publishDir "$output_path", mode: 'copy'

  input:
    path genome from params.genome_path
    val name from params.sample_code
    path npz_file from npz_file
    val pipeline_name from params.pipeline_name
    val pipeline_version from params.pipeline_version

  output:
    path "${name}_aberrations.bed" into aberrations_bed
    path "${name}_bins.bed" into bins_bed
    path "${name}_chr_statistics.txt" into chr_statistics
    path "${name}_segments.bed" into segements_bed
    path "${name}.plots" into plots
    path "${name}.plots/genome_wide.png" into genome_wide_plot

  script:
  """
  #Predict
  WisecondorX predict ${name}.npz ${genome}/$pipeline_name/$pipeline_version/reference_output.npz ${name} --bed --plot
  sleep 5
  """
}