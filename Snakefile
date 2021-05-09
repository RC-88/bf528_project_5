
rule all:
	input: 'P0_1_tophat/accepted_hits.bam'

rule download_samples:
	output: 'SRR1727914.1'
	shell: 'wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/SRR1727914/SRR1727914.1'

rule renamed_sample:
	input: 'SRR1727914.1'
	output: 'P0_1.sra'
	shell: 'mv {input} {output}'

rule exact_fastq:
	input: 'P0_1.sra'
	output:
		fastq1='P0_1_1.fastq',
		fastq2='P0_1_2.fastq'
	shell: 'fastq-dump -I --split-files {input}'

rule run_fastqc:
	input: 
	        FASTQ_1='P0_1_1.fastq',
        	FASTQ_2='P0_1_2.fastq',
		OUTDIR='fastq_output'
	output: 'P0_1_1.fastqc.html', 'P0_1_2.fastqc.html'
	shell: 'fastqc -o {input.OUTDIR} --noextract {input.FASTQ_1} {input.FASTQ_2}'

rule run_samtools:
	input:
        	fastq1='P0_1_1.fastq',
        	fastq2='P0_1_2.fastq',
        	reference='/project/bf528/project_2/reference',
		annot='/project/bf528/project_2/reference/annot/mm9.gtf',
		outdir='P0_1_tophat'
	output: 'P0_1_tophat/accepted_hits.bam'
	shell: 'tophat -r 200 -G {input.annot} --segment-length=20 --segment-mismatches=1 --no-novel-juncs -o {input.outdir} -p 16 {input.reference} {input.fastq1} {input.fastq2}'
		
