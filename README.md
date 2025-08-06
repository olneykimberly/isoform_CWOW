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
