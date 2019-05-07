#!/bin/bash

#SBATCH --job-name="SX238_MIRA"                       # name of the job submitted
#SBATCH -p short                                        # name of the queue you are submitting to
#SBATCH -n 20                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -N 1
#SBATCH --mem=60G                                      # memory for the job
#SBATCH -t 48:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH --mail-user=julestrachsel@gmail.com             # will receive an email when job starts, ends or fails
#SBATCH --mail-type=BEGIN,END,FAIL                      # will receive an email when job starts, ends or fails
#SBATCH -o "stdout.%j.%N"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N"                               # optional but it prints our standard error



# ENTER COMMANDS HERE:

module load java                                        # for bbtools (I'm not loading the scinet version of bbtools because I have installed these in my home directory)
module load mira/4.9.6                                  # MIRA
module load fastx                                       # for renaming reads so MIRA doesnt complain (MIRA hates long read names)
module load spades                                      # spades for comparison
module load bbmap                                       # for those without a personal installation

##########################

# Requirements:
# you need all these things in the same directory  

# 1) Your PE reads named like this : samplename_R1.fq.gz and samplename_R2.fq.gz
# 2) A template MIRA manifest file (TEMPmanifest.conf)
# 3) This template SLURM script
# 4) a samples.txt file of your sample names ( everything before _R1.fq.gz )

# you use this template to make SLURM scripts for each genome you assemble using the following bash loop:

# while read line
# do
# cat MIRAassemblyPipeline.sh | sed "s/SX238/$line/g" > "$line"_MIRA.SLURM
# done < samples.txt


##########################

# DISCLAIMER #
# I make no claims that these are best practices for assembling with MIRA, suggestions welcome.

##########################


# make the mira manifest from a template #

cat TEMPmanifest.conf | sed 's/ISOLATE/SX238/g' > SX238_manifest.conf


# This script is an altered version of the assemblyPipeline.sh script that ships with bbtools

# --- Setup ---

# Interleave reads
reformat.sh in1=SX238_R1.fq.gz in2=SX238_R2.fq.gz out=SX238_reads.fq.gz

#Link the interleaved input file as "SX238_temp.fq.gz"
rm SX238_temp.fq.gz; ln -s SX238_reads.fq.gz SX238_temp.fq.gz

# --- Preprocessing ---

#Remove optical duplicates
clumpify.sh in=SX238_temp.fq.gz out=SX238_clumped.fq.gz dedupe optical
rm SX238_temp.fq.gz; ln -s SX238_clumped.fq.gz SX238_temp.fq.gz

#Remove low-quality regions
filterbytile.sh in=SX238_temp.fq.gz out=SX238_filtered_by_tile.fq.gz
rm SX238_temp.fq.gz; ln -s SX238_filtered_by_tile.fq.gz SX238_temp.fq.gz

#Trim adapters.  Optionally, reads with Ns can be discarded by adding "maxns=0" and reads with really low average quality can be discarded with "maq=8".
bbduk.sh in=SX238_temp.fq.gz out=SX238_trimmed.fq.gz ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=70 ref=/home/julian.trachsel/bbmap/resources/adapters.fa forcetrimright=299 ordered
rm SX238_temp.fq.gz; ln -s SX238_trimmed.fq.gz SX238_temp.fq.gz

#Remove synthetic artifacts and spike-ins by kmer-matching.
bbduk.sh in=SX238_temp.fq.gz out=SX238_filtered.fq.gz k=31 ref=/home/julian.trachsel/bbmap/resources/sequencing_artifacts.fa.gz,/home/julian.trachsel/bbmap/resources/phix174_ill.ref.fa.gz ordered cardinality
rm SX238_temp.fq.gz; ln -s SX238_filtered.fq.gz SX238_temp.fq.gz

#Error-correct phase 1
bbmerge.sh in=SX238_temp.fq.gz out=SX238_ecco.fq.gz ecco mix vstrict ordered ihist=SX238_ihist_merge1.txt
rm SX238_temp.fq.gz; ln -s SX238_ecco.fq.gz SX238_temp.fq.gz

#Error-correct phase 2
clumpify.sh in=SX238_temp.fq.gz out=SX238_eccc.fq.gz ecc passes=4 reorder
rm SX238_temp.fq.gz; ln -s SX238_eccc.fq.gz SX238_temp.fq.gz

#Error-correct phase 3
#Low-depth reads can be discarded here with the "tossjunk", "tossdepth", or "tossuncorrectable" flags.
#For very large datasets, "prefilter=1" or "prefilter=2" can be added to conserve memory.
tadpole.sh in=SX238_temp.fq.gz out=SX238_ecct.fq.gz ecc k=62 ordered tossjunk=t
rm SX238_temp.fq.gz; ln -s SX238_ecct.fq.gz SX238_temp.fq.gz

#Merge
#This phase handles overlapping reads,
#and also nonoverlapping reads, if there is sufficient coverage and sufficiently short inter-read gaps
#For very large datasets, "prefilter=1" or "prefilter=2" can be added to conserve memory.
bbmerge-auto.sh in=SX238_temp.fq.gz out=SX238_merged.fq.gz outu=SX238_unmerged.fq.gz strict rem ordered ihist=SX238_ihist_merge.txt adapter1=CTGTCTCTTATACACATCT adapter2=CTGTCTCTTATACACATCT k=62 extend2=50 ecct
rm SX238_temp.fq.gz; ln -s SX238_merged.fq.gz SX238_temp.fq.gz


#Quality-trim the unmerged reads.  TURNED THIS OFF BECAUSE NOT USNIG NON-MERGED READS
#bbduk.sh in=SX238_unmerged.fq.gz out=SX238_qtrimmed.fq.gz qtrim=r trimq=10 minlen=150 ordered

This is attempting to normalize the merged reads, targeting 80x coverage throwing out those below 5x
bbnorm.sh in=SX238_merged.fq.gz out=SX238_merge_normalized.fq.gz target=80 min=5
rm SX238_temp.fq.gz; ln -s SX238_merge_normalized.fq.gz SX238_temp.fq.gz

pigz -dc SX238_merge_normalized.fq.gz > SX238.fastq

fastx_renamer -n COUNT -i SX238.fastq > SX238_rn.fastq

# --- Assembly ---
##Assemble with MIRA
#

# this runs MIRA with the manifest file that was created above
mira SX238_manifest.conf > SX238_mira.log

# this runs spades for comparison
spades.py -s SX238_temp.fq.gz -o SX238 --careful --assembler-only -t 20

