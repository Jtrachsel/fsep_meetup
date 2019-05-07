#!/bin/bash

#SBATCH --job-name="REPLACE_MIRA"                       # name of the job submitted
#SBATCH -p short                                        # name of the queue you are submitting to
#SBATCH -n 40                                           # number of cores/tasks in this job
#SBATCH -N 1
#SBATCH --mem=125G                                      # memory for the job
#SBATCH -t 48:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH --mail-user=julestrachsel@gmail.com             # will receive an email when job starts, ends or fails
#SBATCH --mail-type=BEGIN,END,FAIL                      # will receive an email when job starts, ends or fails
#SBATCH -o "stdout.%j.%N"                               # standard out %j=job number %N=node name
#SBATCH -e "stderr.%j.%N"                               # optional but it prints our standard error



# ENTER COMMANDS HERE:

module load java                                        # for bbtools
module load mira/4.9.6                                  # MIRA
module load fastx                                       # for renaming reads so MIRA doesnt complain 
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
# cat MIRAassemblyPipeline.sh | sed "s/REPLACE/$line/g" > "$line"_MIRA.SLURM
# done < samples.txt


##########################

# DISCLAIMER #
# I make no claims that these are best practices for assembling with MIRA, suggestions welcome.

##########################


# make the mira manifest from a template #

cat TEMPmanifest.conf | sed 's/ISOLATE/REPLACE/g' > REPLACE_manifest.conf


# This script is an altered version of the assemblyPipeline.sh script that ships with bbtools

# --- Setup ---

# Interleave reads
reformat.sh in1=REPLACE_R1.fq.gz in2=REPLACE_R2.fq.gz out=REPLACE_reads.fq.gz

#Link the interleaved input file as "REPLACE_temp.fq.gz"
rm REPLACE_temp.fq.gz; ln -s REPLACE_reads.fq.gz REPLACE_temp.fq.gz

# --- Preprocessing ---

#Remove optical duplicates
clumpify.sh in=REPLACE_temp.fq.gz out=REPLACE_clumped.fq.gz dedupe optical
rm REPLACE_temp.fq.gz; ln -s REPLACE_clumped.fq.gz REPLACE_temp.fq.gz

#Remove low-quality regions
filterbytile.sh in=REPLACE_temp.fq.gz out=REPLACE_filtered_by_tile.fq.gz
rm REPLACE_temp.fq.gz; ln -s REPLACE_filtered_by_tile.fq.gz REPLACE_temp.fq.gz

#Trim adapters.  Optionally, reads with Ns can be discarded by adding "maxns=0"
# and reads with really low average quality can be discarded with "maq=8"
bbduk.sh in=REPLACE_temp.fq.gz out=REPLACE_trimmed.fq.gz ktrim=r k=23 mink=11 hdist=1 tbo tpe minlen=70 ref=/software/7/apps/bbtools/37.02/resources/adapters.fa forcetrimright=299 ordered
rm REPLACE_temp.fq.gz; ln -s REPLACE_trimmed.fq.gz REPLACE_temp.fq.gz

#Remove synthetic artifacts and spike-ins by kmer-matching.
bbduk.sh in=REPLACE_temp.fq.gz out=REPLACE_filtered.fq.gz k=31 ref=/software/7/apps/bbtools/37.02/resources/sequencing_artifacts.fa.gz,/software/7/apps/bbtools/37.02/resources/phix174_ill.ref.fa.gz ordered cardinality
rm REPLACE_temp.fq.gz; ln -s REPLACE_filtered.fq.gz REPLACE_temp.fq.gz

#Error-correct phase 1
bbmerge.sh in=REPLACE_temp.fq.gz out=REPLACE_ecco.fq.gz ecco mix vstrict ordered ihist=REPLACE_ihist_merge1.txt
rm REPLACE_temp.fq.gz; ln -s REPLACE_ecco.fq.gz REPLACE_temp.fq.gz

#Error-correct phase 2
clumpify.sh in=REPLACE_temp.fq.gz out=REPLACE_eccc.fq.gz ecc passes=4 reorder
rm REPLACE_temp.fq.gz; ln -s REPLACE_eccc.fq.gz REPLACE_temp.fq.gz

#Error-correct phase 3
#Low-depth reads can be discarded here with the "tossjunk", "tossdepth", or "tossuncorrectable" flags.
#For very large datasets, "prefilter=1" or "prefilter=2" can be added to conserve memory.
tadpole.sh in=REPLACE_temp.fq.gz out=REPLACE_ecct.fq.gz ecc k=62 ordered tossjunk=t
rm REPLACE_temp.fq.gz; ln -s REPLACE_ecct.fq.gz REPLACE_temp.fq.gz

#Merge
#This phase handles overlapping reads,
#and also nonoverlapping reads, if there is sufficient coverage and sufficiently short inter-read gaps
#For very large datasets, "prefilter=1" or "prefilter=2" can be added to conserve memory.
bbmerge-auto.sh in=REPLACE_temp.fq.gz out=REPLACE_merged.fq.gz outu=REPLACE_unmerged.fq.gz strict rem ordered ihist=REPLACE_ihist_merge.txt adapter1=CTGTCTCTTATACACATCT adapter2=CTGTCTCTTATACACATCT k=62 extend2=50 ecct
rm REPLACE_temp.fq.gz; ln -s REPLACE_merged.fq.gz REPLACE_temp.fq.gz


#Quality-trim the unmerged reads.  TURNED THIS OFF BECAUSE NOT USNIG NON-MERGED READS
#bbduk.sh in=REPLACE_unmerged.fq.gz out=REPLACE_qtrimmed.fq.gz qtrim=r trimq=10 minlen=150 ordered

# This is attempting to normalize the merged reads, targeting 80x coverage throwing out those below 5x
bbnorm.sh in=REPLACE_merged.fq.gz out=REPLACE_merge_normalized.fq.gz target=80 min=5
rm REPLACE_temp.fq.gz; ln -s REPLACE_merge_normalized.fq.gz REPLACE_temp.fq.gz

pigz -dc REPLACE_merge_normalized.fq.gz > REPLACE.fastq

fastx_renamer -n COUNT -i REPLACE.fastq > REPLACE_rn.fastq

# --- Assembly ---
##Assemble with MIRA
#

# this runs MIRA with the manifest file that was created above
mira REPLACE_manifest.conf > REPLACE_mira.log

# this runs spades for comparison
spades.py -s REPLACE_temp.fq.gz -o REPLACE_spades --careful --assembler-only -t 20

