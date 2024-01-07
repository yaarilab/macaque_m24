The pre-processing pipeline for the R24 macaque sequences.

### The library protocol is as such:

The 5' of read 1 is comprised of (i) a given number of random bases for a sample (4-8) (ii) followed by an Illumina index  (iii) followed by the RACE primer sequence.

The UMI in the header is comprised of the 4-8 random bases. The TRIM comprises of (i),(ii), and (iii) which is retained in the header after trimming the sequence.



### Input files:

1. Two fastq file of paird sequencing
2. sample_name
3. run_FastQC (yes/no)



### Output files:

1. {sampleName}.fasta
2. log tab file for each steps
3. report for each steps
4. pipeline_statistic.csv - table of pass and fail reads for each step.


### Pipeline container:

* Docker: immcantation/suite:4.4.0



### Sequence processing steps:


**1. FastQC**

> fastqc

**2. Paired-end assembly**

> AssemblePairs

**3. Quality control**

> FilterSeq




