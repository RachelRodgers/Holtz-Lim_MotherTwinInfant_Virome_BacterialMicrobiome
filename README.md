# Holtz Lab: Analysis of Mother-Twin Infant Virome & Bacterial Microbiome

----
## Bacterial Microbiome Analysis

Analysis  consists of five scripts:

### QC & Pre-processing
1. TwinMom16S_Quality_ErrorRate_Checking.Rmd
2. TwinMom16S_dada2.R

### Data Cleaning
3. TwinMom16S_Pre-processing.Rmd

### Analysis
4. TwinMom16S_Analysis
5. TwinMom16S_Shared_RelAbd_Plots.Rmd

----
## Virome

Raw sequencing data was processed through VirusSeeker Virome v0.063 (https://wupathlabs.wustl.edu/virusseeker/). Phage reads and corresponding blastX files were used to generate .RMA files for each sample in MEGAN6 Community Edition v6.10.5 using the script **run_blast2rma.sh**.  RMA files were used to generate the Compare file using absolute counts and ignoring unclassified reads.  A BIOM1 file was then created using the following steps:

1. Opened compare file
2. Selected virus node
3. Uncollapsed subtree below virus node
4. Selected subtree
5. File > Export > BIOM1 Format... > (export taxa at official ranks only?) > No
6. Saved .biom file

The resulting biom file was utilized in all viral analyses.
