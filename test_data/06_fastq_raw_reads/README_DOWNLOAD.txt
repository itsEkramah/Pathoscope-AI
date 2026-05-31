
================================================================================
                    BIOLOGICAL SEQUENCING DATA DOWNLOAD GUIDE
================================================================================
To test the pipeline with real biological sequencing datasets, you can download 
a raw FASTQ file from the NCBI Sequence Read Archive (SRA) or European Nucleotide Archive (ENA).

Recommended Dataset: Escherichia Phage Lambda (SRR11246031)
- Platform: Illumina MiSeq paired-end sequencing
- Purpose: High-coverage, small-genome test case ideal for local development pipelines.

INSTRUCTIONS:
1. Open your terminal or a browser.
2. Download R1 (Forward) and R2 (Reverse) fastq files:
   Forward: https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR112/031/SRR11246031/SRR11246031.fastq.gz
   (Optional: You can download and run paired-end QC!)

3. Save the downloaded gzipped FASTQ file directly into the target folder:
   Path: D:\FUCTIONAL_GENOMICS_PROJECT\test_data\06_fastq_raw_reads

4. The pipeline coordinator and preprocessor will automatically detect the gzipped
   fastq format (.fastq.gz) and execute quality control and adapter trimming using
   either the native 'fastp' binary or the high-fidelity pure-Python fallback parser.
================================================================================
