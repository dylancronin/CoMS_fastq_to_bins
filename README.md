# CoMS_fastq_to_bins
- Creating bins from reads  
- Original talk: https://tinyurl.com/ReadsToMags  
- Most recent set of [slides](https://github.com/dylancronin/CoMS_fastq_to_bins/blob/main/CoMS-20250304.pptx) for talk given on 03/04/25

# Table of Contents

1. [The Data](#the-data)  
     - [CAMI Datasets](#cami-datasets)  

2. [Read Quality Control](#read-quality-control)  
   - [bbtools for Read QC and General Formatting](#bbtools-for-read-qc-and-general-formatting)  
     - [Installation](#installation)  
     - [Running bbduk](#running-bbduk)  
     - [Output Interpretation](#output-interpretation)  
   - [Read Quality Visualization](#read-quality-visualization)  
     - [FastQC Installation](#fastqc-installation)  
     - [Running FastQC](#running-fastqc)  
     - [MultiQC for Summary Reports](#multiqc-for-summary-reports)  

3. [Assembly](#assembly)  
   - [MEGAHIT Assembly](#megahit-assembly)  
     - [Installation](#installation)  
     - [Running MEGAHIT](#running-megahit)  
     - [Output Files](#output-files)  

4. [Binning](#binning)   
   - [Read Mapping](#read-mapping)  

5. [Initial Binning Methods](#initial-binning-methods)  
   - [MetaBAT2](#metabat2)  
   - [MaxBin2](#maxbin2)  
   - [CONCOCT](#concoct)  

6. [Ensemble Binning (Combining Multiple Binners)](#ensemble-binning-combining-multiple-binners)  

7. [Bin Refinement](#bin-refinement)  
   - [CheckM for Bin Quality Control](#checkm-for-bin-quality-control)  
   - [RefineM for Further Cleaning](#refinem-for-further-cleaning)  

8. [Simpler Workflow - Aviary](#simpler-workflow---aviary)  
   - [Installation](#installation)  
   - [Setting Up Read Files](#setting-up-read-files)  
   - [Running Aviary](#running-aviary)  


## The Data
[CAMI Datasets](https://cami-challenge.org/datasets/)
Specifically targeting the "High Complexity" fastq files for all 5 samples.

## Read Quality Control

Reads that come back from sequencing runs are not perfect, so they typically require some correction so that errors do not propagate to additional data products from the sequencing reads.

### bbtools for read QC and general formatting

JGI provides a set of scripts that can handle many different tasks you would want to perform on read files called bbtools https://bitbucket.org/berkeleylab/jgi-bbtools/src/master/docs/guides/. Read QC typically is performed through two tools [bbduk](https://bitbucket.org/berkeleylab/jgi-bbtools/src/master/docs/guides/BBDukGuide.txt) and [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). Both tools work and can process your reads fairly quickly, but bbduk can work fairly faster, so that is what we will use today. 

Installation through conda (other options available through the above links). Version 38.51 is used here, but newer versions are now available as well.
```
conda create --name=bbmap-38.51
conda activate bbmap-38.51
conda install -c bioconda bbmap=38.51
```

Run bbduk on the read sets for each sample separately. Likely for many samples, you would want to write a for loop in bash or your favorite language to execute these commands.
```
bbduk.sh in=RH_S001__insert_270.fq out1=RH_S001__insert_270.clean.1.fq.gz out2=RH_S001__insert_270.clean.2.fq.gz minlength=51 qtrim=rl maq=10 maxns=0 trimq=10 ref=/bbmap_38.51/resources/adapters.fa ow=t
bbduk.sh in=RH_S002__insert_270.fq out1=RH_S002__insert_270.clean.1.fq.gz out2=RH_S002__insert_270.clean.2.fq.gz minlength=51 qtrim=rl maq=10 maxns=0 trimq=10 ref=/bbmap_38.51/resources/adapters.fa ow=t
bbduk.sh in=RH_S003__insert_270.fq out1=RH_S003__insert_270.clean.1.fq.gz out2=RH_S003__insert_270.clean.2.fq.gz minlength=51 qtrim=rl maq=10 maxns=0 trimq=10 ref=/bbmap_38.51/resources/adapters.fa ow=t
bbduk.sh in=RH_S004__insert_270.fq out1=RH_S004__insert_270.clean.1.fq.gz out2=RH_S004__insert_270.clean.2.fq.gz minlength=51 qtrim=rl maq=10 maxns=0 trimq=10 ref=/bbmap_38.51/resources/adapters.fa ow=t
bbduk.sh in=RH_S005__insert_270.fq out1=RH_S005__insert_270.clean.1.fq.gz out2=RH_S005__insert_270.clean.2.fq.gz minlength=51 qtrim=rl maq=10 maxns=0 trimq=10 ref=/bbmap_38.51/resources/adapters.fa ow=t
```
Notice that this command can take interleaved reads (forward and reverse are both contained in one file) and output into paired format (forward and reverse separated). Reads will typically be in paired format, but the same structure of the command still applies. This one command packages in both read quality control and trimming, and the associated stats for read removal and such are listed below.
```
Forcing interleaved input because paired output was specified for a single input file.
0.036 seconds.
Initial:
Memory: max=51693m, total=51693m, free=50614m, used=1079m

Added 2970 kmers; time:         0.053 seconds.
Memory: max=51693m, total=51693m, free=48996m, used=2697m

Input is being processed as paired
Started output streams: 0.358 seconds.
Processing time:                206.542 seconds.

Input:                          99811870 reads          14971780500 bases.
Contaminants:                   0 reads (0.00%)         0 bases (0.00%)
QTrimmed:                       45800763 reads (45.89%)         306744010 bases (2.05%)
Low quality discards:           85516 reads (0.09%)     12412925 bases (0.08%)
Total Removed:                  86038 reads (0.09%)     319156935 bases (2.13%)
Result:                         99725832 reads (99.91%)         14652623565 bases (97.87%)

Time:                           206.955 seconds.
Reads Processed:      99811k    482.29k reads/sec
Bases Processed:      14971m    72.34m bases/sec
```
As you can see, most reads are retained in the sample with bases removed from the trimming. Generally, removing obvious errors and very low quality sequences is good practice, but using overly strong QC can change your resulting assemblies and MAG generation.

### Read quality visualization

Installation
```
`conda install -c bioconda fastqc`
```

Running fastqc on just one sample's forward cleaned reads.
```
fastqc RH_S001__insert_270.clean.1.fq.gz
```
Command line output:
```
Started analysis of RH_S001__insert_270.clean.1.fq.gz
Approx 5% complete for RH_S001__insert_270.clean.1.fq.gz
Approx 10% complete for RH_S001__insert_270.clean.1.fq.gz
Approx 15% complete for RH_S001__insert_270.clean.1.fq.gz
Approx 20% complete for RH_S001__insert_270.clean.1.fq.gz
Approx 25% complete for RH_S001__insert_270.clean.1.fq.gz
Approx 30% complete for RH_S001__insert_270.clean.1.fq.gz
Approx 35% complete for RH_S001__insert_270.clean.1.fq.gz
Approx 40% complete for RH_S001__insert_270.clean.1.fq.gz
Approx 45% complete for RH_S001__insert_270.clean.1.fq.gz
Approx 50% complete for RH_S001__insert_270.clean.1.fq.gz
Approx 55% complete for RH_S001__insert_270.clean.1.fq.gz
Approx 60% complete for RH_S001__insert_270.clean.1.fq.gz
Approx 65% complete for RH_S001__insert_270.clean.1.fq.gz
Approx 70% complete for RH_S001__insert_270.clean.1.fq.gz
Approx 75% complete for RH_S001__insert_270.clean.1.fq.gz
Approx 80% complete for RH_S001__insert_270.clean.1.fq.gz
Approx 85% complete for RH_S001__insert_270.clean.1.fq.gz
Approx 90% complete for RH_S001__insert_270.clean.1.fq.gz
Approx 95% complete for RH_S001__insert_270.clean.1.fq.gz
Analysis complete for RH_S001__insert_270.clean.1.fq.gz
```
Output:
HTML file with quality information. Generally good for a sanity check on your data in case there are any obvious issues. 

Note: Also check out [MultiQC](https://multiqc.info/) for generating summary reports of many samples at once.

## Assembly

Two of the more popular tools for assembly are MetaSPAdes and MEGAHIT. Either tool should suit your needs. Generally speaking, MetaSPAdes has shown great performance for assembly, but it is also VERY compute intensive, meaning one will likely need a large amount of memory, cores/threads, and CPU hours. MEGAHIT, on the other hand, performs well but is less compute intensive. 

For this workflow, I have chosen megahit due to the compute needs associated with MetaSPAdes.

Installation
```
conda create --name=megahit-1.2.9
conda activate megahit-1.2.9
conda install -c bioconda megahit
```
Assembly (one example sample)
```
megahit -1 RH_S001__insert_270.clean.1.fq.gz -2 RH_S001__insert_270.clean.2.fq.gz -o RH_S001__insert_270_megahit_out -t 20
```
Output
```
checkpoints.txt  
done  
final.contigs.fa  
intermediate_contigs  
log  
options.json
```
Where final.contigs.fa is ultimately what we care most about, which is the final assembly fasta file output from megahit.

## Binning
Binning more generally relies on two things:
1) Abundances of contigs across samples
2) Tetranucleotide frequencies/genomic features

Requirements for binning based on the above:
Assembly (Contigs/Scaffolds) and BAM file(s)

### Read mapping

For each assembly that you have, you want to map as many related samples as possible to the assembly to maximize the amount of data to better resolve the clustering. For this, we need to do read-mapping to determine the coverages of contigs across samples. Read-mapping generates bam files which are compressed sam files that contain information about the reads mapping to the various contigs. CoverM is a  wrapper to do the read-mapping and some filtering for us. It leverages Minimap2-sr for the actual read-mapping.

Installation
```
conda create --name=coverm-0.6.1
conda activate coverm-0.6.1
conda install -c bioconda coverm=0.6.1
```
Read map every sample against the 1st  sample assembly
```
coverm contig -r RH_S001__insert_270_megahit_out/final.contigs.fa -t 20 \
--bam-file-cache-directory bam_files \
-c RH_S001__insert_270.clean.1.fq.gz RH_S001__insert_270.clean.2.fq.gz \
RH_S002__insert_270.clean.1.fq.gz RH_S002__insert_270.clean.2.fq.gz \
RH_S003__insert_270.clean.1.fq.gz RH_S003__insert_270.clean.2.fq.gz \
RH_S004__insert_270.clean.1.fq.gz RH_S004__insert_270.clean.2.fq.gz \
RH_S005__insert_270.clean.1.fq.gz RH_S005__insert_270.clean.2.fq.gz
```
This will generate five bam files for each set of paired reads that were mapped:
```
RH_S001__insert_270_megahit_out.final.contigs.fna.RH_S001__insert_270.clean.1.fq.gz.bam
RH_S001__insert_270_megahit_out.final.contigs.fna.RH_S002__insert_270.clean.1.fq.gz.bam
RH_S001__insert_270_megahit_out.final.contigs.fna.RH_S003__insert_270.clean.1.fq.gz.bam
RH_S001__insert_270_megahit_out.final.contigs.fna.RH_S004__insert_270.clean.1.fq.gz.bam
RH_S001__insert_270_megahit_out.final.contigs.fna.RH_S005__insert_270.clean.1.fq.gz.bam
```

After generation of the bam files, create sorted and indexed bam files using samtools for input into other programs that require this formatting 
Sort the bam files
```
cd bam_files
for bam in *.bam ; do 
	samtools sort -@ 20 -o ${bam::-4}.sorted.bam $bam ; 
	done
```
Index the bam files
```
for bam in *.sorted.bam ; do 
	samtools index $bam;
	done
cd ..
```



### Initial Binning
**Metabat2**
https://bitbucket.org/berkeleylab/metabat/src/master/
Installation
```
conda create --name=metabat2
conda activate metabat2
conda install -c bioconda/label/cf201901 metabat2
```
Generate metabat coverage
```
coverm contig -b bam_files/*.sorted.bam -m metabat -t 20 > metabat_coverage.tsv
```
Run metabat2
```
metabat2 -i ../RH_S001__insert_270_megahit_out/final.contigs.fa -a metabat_coverage.tsv -o metabat2 -t 20
```
**Maxbin2**
https://sourceforge.net/projects/maxbin2/
Installation
```
conda create --name=maxbin-2.2.7
conda activate maxbin-2.2.7
conda install -c bioconda maxbin2
```
Leverage coverm to calculate the trimmed mean coverage estimates (can use other coverage estimates other than trimmed_mean)
```
cd bam_files
for bam in *.sorted.bam ; do 
	coverm contig -t 20 -m trimmed_mean -b $bam > ../${bam::-4}.trimmed_mean.tsv ; 
	done
```
Create a file called *abund_list.txt* pointing to each tsv generated above that would look like this if the tsv files are located in the same directory:
```
RH_S001__insert_270_megahit_out.final.contigs.fna.RH_S001__insert_270.clean.1.fq.gz.trimmed_mean.tsv
RH_S001__insert_270_megahit_out.final.contigs.fna.RH_S004__insert_270.clean.1.fq.gz.trimmed_mean.tsv
RH_S001__insert_270_megahit_out.final.contigs.fna.RH_S002__insert_270.clean.1.fq.gz.trimmed_mean.tsv
RH_S001__insert_270_megahit_out.final.contigs.fna.RH_S005__insert_270.clean.1.fq.gz.trimmed_mean.tsv
RH_S001__insert_270_megahit_out.final.contigs.fna.RH_S003__insert_270.clean.1.fq.gz.trimmed_mean.tsv
```
Run Maxbin2 algorithm with a minimum contig length of 2500bp
```
mkdir maxbin2
run_MaxBin.pl  -contig ../RH_S001__insert_270_megahit_out/final.contigs.fa -abund_list abund_list.txt -out maxbin2/maxbin2 -min_contig_length 2500 -thread 20
```
**CONCOCT**
https://concoct.readthedocs.io/en/latest/
CONCOCT requires that you first cut up your contigs into chunks before coverage estimation
```
cut_up_fasta.py ../RH_S001__insert_270_megahit_out/final.contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
```
Generate CONCOCT coverage table based on the sorted bam files generated earlier
```
concoct_coverage_table.py contigs_10K.bed bam_files/*.sorted.bam > concoct_coverage_table.tsv
```
Run CONCOCT algorithm and merge the cut contigs from step 1 
```
concoct -l 2500 --threads 20 --composition_file contigs_10K.fa --coverage_file concoct_coverage_table.tsv -b concoct/

merge_cutup_clustering.py concoct/clustering_gt2500.csv > concoct/clustering_merged.csv
```
### Ensemble binning
Ensemble binning takes the output bins from the standalone "initial binners" and uses the information across bin sets to build a more high quality superset of bins. DAS Tool and MetaWRAP are the two most common at the moment. For this tutorial we will be using the MetaWRAP bin_refinement module.

**MetaWRAP**
https://github.com/bxlab/metaWRAP

For the ensemble binning, we are using the metawrap bin_refinement module. This module can accept the input from three different binners, so we are using the bins generated from the "Initial Binners" above -- metabat2, maxbin2, concoct -- as input into the ensemble binning module from metawrap.

Installation
```
conda create --name=metawrap-1.3.2
conda activate metawrap-1.3.2
conda install -y -c ursky metawrap-mg=1.3.2
```
Run metawrap bin_refinement using the output from the three "Initial binners"
```
metawrap bin_refinement -A metabat2/ -B maxbin2/fasta_bins/ \
-C concoct/fasta_bins/ -t 20 -o metawrap_bin_refinement/
```

Result:
By default, metawrap filters for bins that are >70% complete and <10% contaminated based on [CheckM](https://github.com/Ecogenomics/CheckM). CheckM uses lineage specific single-copy genes to estimate completeness and contamination of a particular genome. Roughly speaking, if a genome bin has more than one single-copy gene, that adds to the contamination because it was only expected once. Then, the completeness is determined by how many single-copy genes the genome bin has versus how many there "should" be for that lineage. 


## Bin Refinement

**Automated Bin QC -- RefineM**

RefineM will remove any potential contaminants and outliers based on tetranucleotide frequencies and abundances. It also has several other functions such as removing contigs that do not match taxonomically, but that will not be utilized here.

Install [RefineM](https://github.com/dparks1134/RefineM)
```
conda create --name=refinem-0.1.2
conda activate refinem-0.1.2
conda install refinem
```
Generate the scaffold/contig stats for the contigs in the bins generated (e.g. generates coverage profiles and tetranucleotide frequencies)
```
refinem scaffold_stats -x fa -c 20 \
../RH_S001__insert_270_megahit_out/final.contigs.fa \
metawrap_bin_refinement/metawrap_70_10_bins/ \
scaffold_stats bam_files/*.sorted.bam
```
Determine contigs that are outliers in the dataset
```
refinem outliers scaffold_stats/scaffold_stats.tsv scaffold_stats_outliers
```
Filter out outliers from bins
```
refinem filter_bins -x fa \
metawrap_bin_refinement/metawrap_70_10_bins/ \
scaffold_stats_outliers/outliers.tsv \
metawrap_bin_refinement/metawrap_70_10_bins_refinem_filtered/
```

# Simpler Workflow - Aviary (https://github.com/rhysnewell/aviary)
```
conda create -n aviary -c bioconda aviary
conda activate aviary
```
See https://github.com/rhysnewell/aviary for installing databases as well.

We need to first make sure our readsets aren't interleaved. They are by default when downloading the data. reformat.sh does the trick for this.
```
conda activate bbmap-38.51
reformat.sh in=RH_S001__insert_270.fq out1=RH_S001__insert_270.1.fq.gz out2=RH_S001__insert_270.2.fq.gz
reformat.sh in=RH_S002__insert_270.fq out1=RH_S002__insert_270.1.fq.gz out2=RH_S002__insert_270.2.fq.gz
reformat.sh in=RH_S003__insert_270.fq out1=RH_S003__insert_270.1.fq.gz out2=RH_S003__insert_270.2.fq.gz
reformat.sh in=RH_S004__insert_270.fq out1=RH_S004__insert_270.1.fq.gz out2=RH_S004__insert_270.2.fq.gz
reformat.sh in=RH_S005__insert_270.fq out1=RH_S005__insert_270.1.fq.gz out2=RH_S005__insert_270.2.fq.gz
```

Then just point aviary to wherever these read files are stored.
```
mkdir reads/
mv *.1.fq.gz reads/
mv *.2.fq.gz reads/
```
```
conda activate aviary
aviary recover -1 reads/*.1.fq.gz -2 reads/*.2.fq.gz --output aviary_output/ --max_threads 12 --n_cores 24
```
And that's it! Assuming you have your databases set up correctly, you will have CheckM, GTDB, among other results.


