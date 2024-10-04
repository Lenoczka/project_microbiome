#!/bin/bash

# Microbiome in colorectal cancer https://www.ebi.ac.uk/ena/browser/view/PRJEB57580
#The 16S rRNA amplicon library was amplified with 341F/805R primers (V3-V4) (CCTACGGGNGGCWGCAG (17), GGACTACHVGGGTATCTAAT(20)). 2x300 Illumina MiSeq. 
#805-341+1=465 library size, 600-465=135 overlap

mkdir /home/stud50/project_microbiome
mkdir /home/stud50/project_microbiome/raw_fastq
#conda create -n microbiome
#conda activate microbiome
#conda install bioconda::fastq-dl
#fastq-dl -a PRJEB57580 #download .fastq.gz #raw data from ENA bioproject

gunzip *.fastq.gz
mkdir fastqc

conda activate blastim_py3
fastqc -o ./fastqc ./*.fastq #ls | wc -l control for number of files
multiqc .
#multiqc_report.html #view multiqc summary of raw .fastq quality report 
#trimming of quality and primers (they were not trimmed befofe deposition in ENA)
#341F
CCTACGGGNGGCWGCAG (17)
805R
GGACTACHVGGGTATCTAAT (20)

conda activate blastim_py3
mkdir /home/stud50/project_microbiome/raw_fastq/trimmed/

#341F primer trimming, quality trimming -q 25 did not use 
parallel --dry-run -j 8 'cutadapt -u 17 -o "/home/stud50/project_microbiome/raw_fastq/trimmed/{/}" "{}"' ::: $(ls /home/stud50/project_microbiome/raw_fastq/*_1.fastq)

#805R primer trimming, quality trimming -q 25 did not use 
parallel -j 8 'cutadapt -u 20 -o "/home/stud50/project_microbiome/raw_fastq/trimmed/{/}" "{}"' ::: $(ls /home/stud50/project_microbiome/raw_fastq/*_2.fastq)

conda activate microbiome
#conda install bioconda::pear
#conda install parallel
mkdir /home/stud50/project_microbiome/raw_fastq/trimmed/pear

#merge paired-end reads
nano pear.sh 
#                                                 
#!/bin/bash
sample_ID=$(ls /home/stud50/project_microbiome/raw_fastq/trimmed/ | cut -d_ -f1 | sort | uniq)
for file in $sample_ID;
do
    pear -j 4 --min-overlap 50 -q 25 -f /home/stud50/project_microbiome/raw_fastq/trimmed/${file}_1.fastq \
-r /home/stud50/project_microbiome/raw_fastq/trimmed/${file}_2.fastq \
-o /home/stud50/project_microbiome/raw_fastq/trimmed/pear/$file;
done

#for single sample
pear -j 4 --min-overlap 50 -q 25 -f /home/stud50/project_microbiome/raw_fastq/trimmed/ERR10493148_1.fastq \
-r /home/stud50/project_microbiome/raw_fastq/trimmed/ERR10493148_2.fastq \
-o /home/stud50/project_microbiome/raw_fastq/trimmed/pear/ERR10493148


#create tmux session for continue process after technical dissconection from server
tmux new -s pear
bash pear.sh  #take times
ctrl +B, D #exit from tmux session
#tmux ls
#tmux a -t pear #enter to session again
#ps x
#tmux kill-session -t o
 

#"Turn on" environment with QIIME2, so it becomes reacheable
conda activate qiime2-amplicon-2024.5


# create manifest txt file with pathways to the input fastq files - microbiome samples folder, fastq format - one sapmple per one fastq file (merged with pear)
mkdir qiime2

#location cd /home/stud50/project_microbiome/raw_fastq/trimmed/pear/
ls *assembled.fastq | cut -d. -f1,2 | sort | uniq > /home/stud50/project_microbiome/qiime2/manifest1 #add string for sample-id future column 

ls /home/stud50/project_microbiome/raw_fastq/trimmed/pear/*assembled.fastq | sort > /home/stud50/project_microbiome/qiime2/manifest2 #add string for absolute-filepath future column 

#join manifest1 and manifest2 files
paste -d"	" manifest1 manifest2 > manifest

#add strings for future direction column
awk 'BEGIN { OFS = "	" } { print $0, "forward" }' manifest >  manifest3

#add columns sample-id,absolute-filepath,direction to manifest file
awk 'BEGIN { print "sample-id	absolute-filepath	direction" } { print }' manifest3 > manifest_

#Concatenate fastq into qza archive
#15 mins
qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path manifest_ --output-path ./seq.qza --input-format SingleEndFastqManifestPhred33V2 

#FILTER AND CLUSTERIZE FASTQ READS WITH DADA2 algorythms
mkdir denoise_out
#4 hours
qiime dada2 denoise-single --i-demultiplexed-seqs ./seq.qza --p-trunc-len 300 --p-trunc-q 10 --o-table "denoise_out/feature_table.qza" --o-representative-sequences "denoise_out/rep_seqs.qza" --o-denoising-stats "denoise_out/stats.qza" 

#Taxonomy object V3-V4 16S
##Install Silva reference data base

qiime rescript get-silva-data \
--p-version "138.1" \
--p-target "SSURef_NR99" \
--o-silva-sequences silva-138.1-ssu-rna-seqs.qza \
--o-silva-taxonomy silva-138.1-ssu-nr99-tax.qza

qiime rescript reverse-transcribe --i-rna-sequences silva-138.1-ssu-rna-seqs.qza --o-dna-sequences silva-138.1-ssu-dna-seqs.qza

qiime rescript cull-seqs --i-sequences silva-138.1-ssu-dna-seqs.qza --o-clean-sequences silva-138.1-ssu-dna-seqs-clean.qza

qiime rescript filter-seqs-length-by-taxon --i-sequences silva-138.1-ssu-dna-seqs-clean.qza  --i-taxonomy silva-138.1-ssu-nr99-tax.qza --o-filtered-seqs silva-138.1-ssu-nr99-seqs-filt.qza --o-discarded-seqs ./discarded --p-labels Archaea Bacteria Eukaryota  --p-min-lens 900 1200 1400

qiime rescript dereplicate  --p-threads 4 --i-sequences silva-138.1-ssu-nr99-seqs-filt.qza --i-taxa silva-138.1-ssu-nr99-tax.qza --p-mode "uniq" --o-dereplicated-sequences silva138.1-ssu-nr99-tax-derep-uniq.qza --o-dereplicated-taxa ./dereplicated

##Train full 16S classifier
qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads silva138.1-ssu-nr99-tax-derep-uniq.qza --o-classifier silva138.1-ssu-nr99-full16S-classifier.qza --i-reference-taxonomy dereplicated.qza

#Extract classifier only for V3-V4 CCTACGGGNGGCWGCAG, GGACTACHVGGGTATCTAAT
qiime feature-classifier extract-reads \
  --i-sequences silva138.1-ssu-nr99-tax-derep-uniq.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GGACTACHVGGGTATCTAAT \
  --p-min-length 250 \
  --p-max-length 500 \
  --o-reads v3-v4.qza

qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads v3-v4.qza --o-classifier V3-V4-classifier.qza --i-reference-taxonomy dereplicated.qza


ln -s /home/stud50/project_microbiome/qiime2/V3-V4-classifier.qza classifier
ln -s /home/stud50/project_microbiome/qiime2/denoise_out/feature_table.qza feature_table
ln -s /home/stud50/project_microbiome/qiime2/denoise_out/rep_seqs.qza rep_seqs_qza

# EXPORT REPRESENTATIVE SEQUENCES AS .fasta
qiime tools export --input-path rep_seqs_qza --output-path "denoise_out/"

ln -s /home/stud50/project_microbiome/qiime2/denoise_out/dna-sequences.fasta rep_seqs


# EXPORT DENOISING STATS
qiime tools export --input-path "denoise_out/stats.qza" --output-path "denoise_out/"

ln -s /home/stud50/project_microbiome/qiime2/denoise_out/stats.tsv dada_stats


# CLASSIFY REPRESENTATIVE SEQUENCES WITH --p-confidence {0.94}
mkdir taxonomy
qiime feature-classifier classify-sklearn --i-classifier classifier --i-reads rep_seqs_qza --o-classification "taxonomy/taxonomy.qza" --p-confidence 0.94


# EXPORT OUTPUT TABLE
qiime tools export --input-path "taxonomy/taxonomy.qza" --output-path "taxonomy"


#RAREFY and calculate ALPHA-diversity

mkdir rarefied
qiime feature-table rarefy --i-table "denoise_out/feature_table.qza" --p-sampling-depth 10000 --o-rarefied-table "rarefied/otus_rar_10K.qza"
qiime diversity alpha --i-table "rarefied/otus_rar_10K.qza" --p-metric "chao1" --o-alpha-diversity "rarefied/alpha_chao.qza"
qiime tools export --input-path "rarefied/alpha_chao.qza" --output-path "rarefied/alpha_chao.tsv" --output-format "AlphaDiversityFormat"


#make readable taxa tables
mkdir otus

qiime tools export --input-path "rarefied/otus_rar_10K.qza" --output-path "otus"

qiime taxa collapse --i-table "rarefied/otus_rar_10K.qza" --i-taxonomy "taxonomy/taxonomy.qza" --p-level 7 --o-collapsed-table "otus/collapse_7.qza"

qiime tools export --input-path "otus/collapse_7.qza" --output-path "otus/summarized_taxa"


#Convert bacteria-sample abundance tables to human readable formats

biom convert -i "otus/summarized_taxa/feature-table.biom" -o "otus/summarized_taxa/otu_table_L7.txt" --to-tsv

biom convert -i "otus/feature-table.biom" -o "otus/summarized_taxa/otu_table.txt" --to-tsv #without taxonomy table, only ASV names 


