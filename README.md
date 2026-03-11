https://www.pacb.com/wp-content/uploads/Application-note-Bioinformatics-tools-for-full-length-isoform-sequencing.pdf
https://www.youtube.com/watch?v=-4C_G5IOhyk
https://pacb.my.salesforce.com/sfc/p/70000000IVif/a/PF0000008HN3/ISVWd.KJLJtqmowyMrlJ_ZwduWwx4zfwmx69eNYsFdQ
https://www.pacb.com/wp-content/uploads/SMRT-Link-Kinnex-full-length-RNA-troubleshooting-guide.pdf


# Human reference 
https://downloads.pacbcloud.com/public/dataset/MAS-Seq/REF-pigeon_ref_sets/Human_hg38_Gencode_v39/

# isoform_CWOW
Long and short read isoform differences among human Lewy body dementia brains

Output form TGen:
The results are in unmapped BAM file format (u.bam) and the workflow used to process the results from the sequencer is located here https://github.com/tgen/tgs_readgeneration?tab=readme-ov-file#pacbio-revio. Steps 1-3 were followed from PacBio's workflow https://isoseq.how/clustering/cli-workflow.html.


Yes, we do have a pipeline that you can use to identify novel isoforms. Since you are working with full-length Iso-seq data you would start with the https://isoseq.how/clustering/cli-workflow.html (which includes a step for merging SMRT cells) and then follow with the https://isoseq.how/classification/workflow.html. The tool for isoform classification is pigeon, which is based off SQANTI3: https://isoseq.how/classification/categories.

Do you know if the TGen sequencing core did any type of preliminary processing of the data, specifically demultiplexing at the Kinnex barcoded adapters level and deconcatenating the arrays? If not, you might actually need to run a couple of steps before starting the Iso-seq pipeline. If you are not sure, please send me the list of files that you received from the core, along with the header of one of your BAM files (samtools view -H <file_name>.bam)

# Samtools header 
-lima for demultiplexing at the barcoded adapters level
-skera for deconcatenating the Kinnex arrays into S-reads
-lima for primer removal and demultiplexing
-isoseq refine for removal of the poly-A tails

#-------
conda activate isoseq 


ls /tgen_labs/jfryer/cores/tgen/r84132_20250721_215156/1_D01/CWOW_NA07150/CWOW_NA07150_1_BR_LBD_C1_PKFLR_A38419_SMRT3337_ACAGTC_L004.hifi.u.bam /tgen_labs/jfryer/cores/tgen/r84132_20250726_004614/1_C01/CWOW_NA07150/CWOW_NA07150_1_BR_LBD_C1_PKFLR_A38419_SMRT3340_ACAGTC_L003.hifi.u.bam > CWOW_NA07150.fofn


#-------- TGen
1. lima for demultiplexing at the barcoded adapters level
2. skera for deconcatenating the Kinnex arrays into S-reads
3. lima for primer removal and demultiplexing
4. isoseq refine for removal of the poly-A tails

#-------- Me
merged the SMRT cells
isoseq cluster2
aligned via pbmm2.


Prepare reference files for pigeon


isoseq collapse --do-not-collapse-extra-5exons NA15-031_CONTROL_mapped.bam NA15-031_CONTROL_collapsed.gff
pigeon sort NA15-031_CONTROL_collapsed.gff
pigeon classify NA15-031_CONTROL_collapsed.gff.sorted ../GRCh38/gencode.v38.annotation.sorted.gtf ../GRCh38/GRCh38.primary_assembly.genome.fa




#-------- SQANTI3
wget https://github.com/ConesaLab/SQANTI3/releases/download/v5.5/SQANTI3_v5.5.zip
mkdir sqanti3
unzip SQANTI3_v5.5.zip -d sqanti3

The SQANTI3 tool is designed to enable quality control and filtering of long read-defined transcriptomes, which are often rich in artifacts and false-positive isoforms. Therefore, a good curation of the transcriptome is indispensable to proceed with FIT analysis and produce valid, biologically sound conclusions/hypothesis.

the files head NA18-285_LBD_collapsed_sorted_classification.filtered_lite_classification.txt contain the counts 
After running the pigeon filter pipeline step, the file that contains the counts for each isoform is named with the extension _filtered.classification.txt.

This file contains isoform annotation information, which includes columns for full-length (FL) read counts. Specifically, the columns FL.<Sample> and fl_assoc provide the full-length read counts for each isoform.


#------- FLAIR https://flair.readthedocs.io/en/latest/modules.html
conda activate flair

Run flair correct and collapse individually on each sample, then combine your transcriptomes using flair combine. This method will be faster and easier, but you may miss some low-expression transcripts.
Your other option is to merge your samples before running FLAIR. If using PacBio reads, be careful doing this, as the reads may not have unique IDs. You may need to label each read with its sample ID to keep the read IDs unique. You can either merge the FASTA/FASTQ files before running FLAIR (simplest, recommended).

Modules must be run in order (align, correct, collapse).

** If you want to compare multiple samples, there are two primary ways of doing this:
Run FLAIR align, correct, and collapse (or FLAIR transcriptome) on each sample separately (better for large sets of samples)