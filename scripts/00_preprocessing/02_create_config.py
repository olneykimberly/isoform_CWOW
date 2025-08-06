#!/usr/bin/python3

# create a new output file
outfile = open('config.json', 'w')

# get all sample names
allSamples = list()
numSamples = 0

with open('simple_kbase_name.tsv', 'r') as infile:
    for line in infile:
        numSamples += 1

        line = line.replace("\t", "_")
        split = line.split()
        sampleAttributes = split[0].split('_')  # CWOW_NA15333_NA15-333_SMRT3337_1_D01_CWOW_NA15333_1_BR_CONTROL_C1_PKFLR_A38414_SMRT3337_ACAGTC_L004
                                                #  E1_BR.FCHVC2VDRXY_L1_R1_ITAAGTGGT-CTTAAGCC.fastq.gz pigID_tissue.sequencer_lane_read_X-X.fastq.gz
        # create a shorter sample name
        stemName = sampleAttributes[2] + '_' + sampleAttributes[10] # CWOW_NA15333_CONTROL
        allSamples.append(stemName)

# create header and write to outfile
header = '''{{
    "Commment_Input_Output_Directories": "This section specifies the input and output directories for scripts",
    "merged_SMRTcells" : "/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/merged_SMRTcells/",
    "cluster2" : "/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/cluster2/",
    "mapped" : "/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/mapped/",

    "Comment_Reference" : "This section specifies the location of the human, Genocode reference genome",
    "ref_fa" : "/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/GRCh38/GRCh38.primary_assembly.genome.fa",

    "Comment_Sample_Info": "The following section lists the samples that are to be analyzed",
    "sample_names": {0},
'''
outfile.write(header.format(allSamples))

# config formatting
counter = 0
with open('simple_kbase_name.tsv', 'r') as infile:
    for line in infile:
        counter += 1

        # make naming consistent, we will rename using only underscores (no hyphens)
        line = line.replace("\t", "_")
        split = line.split()
        sampleAttributes = split[0].split('_')  # project_uniqueNum_1_tissue_group_XX_XX_sequencer_adapter_lane_read_001

        base = '/tgen_labs/jfryer/cores/tgen/r84132_20250721_215156/1_D01/' + sampleAttributes[6] + '_' + sampleAttributes[7] + '/' + sampleAttributes[6] + '_' + sampleAttributes[7] +  '_' + sampleAttributes[8] + '_' + sampleAttributes[9] + '_' + sampleAttributes[10] +'_' + sampleAttributes[11] +'_' + sampleAttributes[12] +'_' + sampleAttributes[13] +'_' + sampleAttributes[14] +'_' + sampleAttributes[15] +'_' + sampleAttributes[16]
        sampleName1 = base 
        sampleName2 = sampleName1.replace("SMRT3337_ACAGTC_L004", "SMRT3340_ACAGTC_L003")
        sampleName2 = sampleName2.replace("/tgen_labs/jfryer/cores/tgen/r84132_20250721_215156/1_D01/", "/tgen_labs/jfryer/cores/tgen/r84132_20250726_004614/1_C01/")
        sampleInfo = split[0]

        # create a shorter sample name
        stemName = sampleAttributes[2] + '_' + sampleAttributes[10] # NAID_GROUP
        shortName1 = stemName + '_SMRT3337_L004'
        shortName2 = stemName + '_SMRT3340_L003'

        # break down fastq file info
        # @A00127:312:HVNLJDSXY:2:1101:2211:1000
        # @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>
        sampleInfo = sampleInfo.split('_')
        instrument = sampleInfo[0]
        runNumber = sampleInfo[1]
        flowcellID = sampleInfo[2]

        lane = sampleInfo[3]
        ID = stemName  # ID tag identifies which read group each read belongs to, so each read group's ID must be unique
        SM = stemName  # Sample
        PU = flowcellID + "." + lane  # Platform Unit
        LB = stemName

        out = '''
    "{0}":{{
        "bam1": "{1}",
        "bam2": "{2}",
        "bam_SMRT1": "{3}",
        "bam_SMRT2": "{4}",
        "merge" : "{5}",
        "ID": "{5}",
        "SM": "{6}",
        "PU": "{7}",
        "LB": "{8}",
        "PL": "PacBio"
        '''
        outfile.write(out.format(stemName, sampleName1, sampleName2, shortName1, shortName2, stemName, stemName, PU, LB))
        if (counter == numSamples):
            outfile.write("}\n}")
        else:
            outfile.write("},\n")
outfile.close()
