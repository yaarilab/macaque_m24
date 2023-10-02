The processing pipeline for the R24 macaque sequences.

The library protocol is as such:

The 5' of read 1 is comprised of (i) a given number of random bases for a sample (4-8) (ii) followed by an Illumina index  (iii) followed by the RACE primer sequence.

The UMI in the header is comprised of the 4-8 random bases that we are using with CollapseSeq for deduplication. The TRIM comprises of (i),(ii), and (iii) which is retained in the header after trimming the sequence.

