# BELNIPT IQC pipeline

This pipeline was developed using Nextflow WSL1, using Nextflow 21.04.3 and miniconda3-py39_4.10.3 on a RHEL 8.
All needed conda environments and tools are available in the conda_yml directory.

## Cite:

Martens, G. et al. Integration of Trophoblastic Mosaic Ratio to Enhance Post-Test Counseling in Genome-Wide NIPT Screening. (paper submitted)

## Create needed reference files
### Create reference for bowtie

Make sure the hg19 reference genome is available in /ref as genome.fa

```bash
nextflow prepare_bowtie.nf --genome_dir /ref/
```

### Prepare reference samples

This assumes that your reference samples are placed in /data. Each sample has 1 R1.fastq.gz and 1 R2.fastq.gz file (```bash zcat $sample"_R1_001.fastq.gz" | gzip >> $sample".R1.fastq.gz"```).
Data will be stored in the same directory as the fastq files.

```bash
nextflow prepare_samples_for_refset.nf --genome_path /ref/ --samplename $sample --path /data/ --read1 /data/$sample".R1.fastq.gz" --read2 /data/$sample".R2.fastq.gz"
```

### Prepare the refset

The needed reference will be created, started from the data generated in the previous step, and the result will be stored in the reference directory.

```bash
nextflow prepare_refest.nf --input /data/ --output /ref/
```

## Running a sample

Samples are ran one-by-one. Data needs to be placed in the /sample directory: /sample/fastq for the fastq files. Results will be placed in /sample/nipt/v1.0.0/.
The seqff directory found in this project must be stored in /seqff.

```bash
nextflow pipeline.nf --genome_path /ref/ --sample_path /sample --sample_code $sample
```

## Final used metrics:

The used metrics in the paper can be found:

- reads mapping once: calculated from $sample.err file: *aligned concordantly exactly 1 time* + *aligned discordantly 1 time* + 2 * *aligned exactly 1 time*
- seqFF: found in $sample_seqff.csv
- Yu FF: found in $sample.yu_ff.tsv

## Used Tools:

- [Bowtie2](https://github.com/BenLangmead/bowtie2)  
- [elprep](https://github.com/ExaScience/elprep)  
- [samtools](https://github.com/samtools/samtools)  
- [WiseCondorX](https://github.com/CenterForMedicalGeneticsGhent/WisecondorX)  
- [SeqFF](https://obgyn.onlinelibrary.wiley.com/doi/10.1002/pd.4615)  