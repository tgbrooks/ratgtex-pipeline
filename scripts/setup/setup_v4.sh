set -e

#########################
## Prepare annotations ##
#########################

cd ref

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/323/735/GCF_036323735.1_GRCr8/GCF_036323735.1_GRCr8_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/323/735/GCF_036323735.1_GRCr8/GCF_036323735.1_GRCr8_genomic.gtf.gz

gunzip GCF_036323735.1_GRCr8_genomic.fna.gz
gunzip GCF_036323735.1_GRCr8_genomic.gtf.gz

python3 ../scripts/setup/prepare_refseq.py \
    --map sequence_report.tsv \
    --gtf-in GCF_036323735.1_GRCr8_genomic.gtf \
    --gtf-out GCF_036323735.1_GRCr8_genomic.chr.gtf \
    --fasta-in GCF_036323735.1_GRCr8_genomic.fna \
    --fasta-out GCF_036323735.1_GRCr8_genomic.chr.fa

cd ..

###################
## Set up for QC ##
###################

# NOTE: could not run this section

while read tissue; do
    mkdir -p v4/${tissue}
    wget -O v4/${tissue}/rat_ids.txt "https://ratgtex.org/data/v4/rat_ids/rat_ids.${tissue}.v4.txt"
    # Fastq-map now done by the setup_sra.py script
done <tissues.dup.txt

mkdir -p geno/original
cat v4/*/rat_ids.txt | sort | uniq >geno/original/rfids_for_geno_11_2.txt

## Transfer round11_2_vcf_ids.txt from TSCC. These IDs have extra info around RFIDs and must be parsed.

## Make lists of RFIDs and corresponding VCF IDs to extract genotypes
python scripts/setup/genotype_ids.py \
    --rfids geno/original/rfids_for_geno_11_2.txt \
    --vcf-ids geno/original/round11_2_vcf_ids.txt \
    --out-tsv geno/original/id_map_to_extract_11_2.tsv

## Transfer output to TSCC and run the commands documented in comments in genotypes_v4.sh

#######################
## Set up for Pantry ##
#######################

# NOTE: more phenos are used for other steps, so I stopped trying and gave up here
# we may just not run that part of the pipeline
while read tissue; do
    mkdir -p v4/$tissue/phenos
    mkdir -p v4/$tissue/phenos/output
    mkdir -p v4/$tissue/pheast/intermediate
    mkdir -p v4/$tissue/pheast
    mkdir -p v4/$tissue/pheast/intermediate/covar
    # NOTE: I am not confident about the sources of these files!
    # norm or unnorm?
    wget -O v4/$tissue/phenos/output/alt_polyA.bed.gz "https://ratgtex.org/data/v4/phenos/phenos.$tissue.alt_polyA.unnorm.v4_rn8.bed.gz"
    wget -O v4/$tissue/pheast/intermediate/covar/alt_polyA.covar.tsv "https://ratgtex.org/data/v4/covar/covar.Liver.alt_polyA.v4.tsv"
    tabix -p bed v4/$tissue/phenos/output/alt_polyA.bed.gz
    cp scripts/config_phenos.yml v4/$tissue/phenos/config.yml
    # (Then update read length in config.yml)
    # NOTE: the above appears to need to be done manually?
    echo "USER MUST MANUALLY update read lengths in v4/$tissue/phenos/config.yml"
    mkdir -p v4/$tissue/phenos/intermediate/bam
done <tissues.dup.txt

## Link FASTQ and BAM files

# NOTE: we now download from SRA instead
# ln -s ../../fastq/Adipose_Liver/Fat v4/Adipose/fastq
# ln -s ../../fastq/BLA_NAcc2_PL2 v4/BLA/fastq
# ln -s ../../fastq/Brain v4/Brain/fastq
# ln -s ../../fastq/Eye v4/Eye/fastq
# ln -s ../../fastq/IC_NAcc4_pVTA2 v4/IC/fastq
# ln -s ../../fastq/IL_LHb_NAcc_OFC_PL v4/IL/fastq
# ln -s ../../fastq/IL_LHb_NAcc_OFC_PL v4/LHb/fastq
# ln -s ../../fastq/Adipose_Liver/Liver v4/Liver/fastq
# ln -s ../../fastq/IL_LHb_NAcc_OFC_PL v4/NAcc1/fastq
# ln -s ../../fastq/BLA_NAcc2_PL2 v4/NAcc2/fastq
# ln -s ../../fastq/NAcc3_PL3_pVTA v4/NAcc3/fastq
# ln -s ../../fastq/IC_NAcc4_pVTA2 v4/NAcc4/fastq
# ln -s ../../fastq/IL_LHb_NAcc_OFC_PL v4/OFC/fastq
# ln -s ../../fastq/IL_LHb_NAcc_OFC_PL v4/PL1/fastq
# ln -s ../../fastq/BLA_NAcc2_PL2 v4/PL2/fastq
# ln -s ../../fastq/NAcc3_PL3_pVTA v4/PL3/fastq
# ln -s ../../fastq/NAcc3_PL3_pVTA v4/pVTA1/fastq
# ln -s ../../fastq/IC_NAcc4_pVTA2 v4/pVTA2/fastq
# ln -s ../../fastq/RMTg v4/RMTg/fastq

# while read tissue; do
#     while read sample; do
#         ln -s ../../../star_out/$sample.Aligned.sortedByCoord.out.bam v4/$tissue/phenos/intermediate/bam/$sample.bam
#         ln -s ../../../star_out/$sample.Aligned.sortedByCoord.out.bam.bai v4/$tissue/phenos/intermediate/bam/$sample.bam.bai
#     done <v4/$tissue/rat_ids.txt
# done <tissues.dup.txt

# ## After running phenotyping
#
# bash scripts/setup_merged_tissues.sh
#
# while read tissue; do
#     ## Make combined phenotype table
#     cd v4/$tissue/phenos/
#     bash scripts/combine_modalities.sh alt_polyA alt_TSS expression isoforms splicing stability
#     cd ../../../
#
#     rsync -av ~/tools/Pantry/pheast/ v4/$tissue/pheast --exclude input --exclude intermediate --exclude output --exclude .snakemake
#     cp scripts/config_pheast.yml v4/$tissue/pheast/config.yml
# done <tissues.txt
#
# ## After Brain plink files are generated
#
# mkdir -p ref/LDREF
# for chrom in {1..20}; do
#     plink --bfile v4/Brain/geno --chr $chrom --make-bed --out ref/LDREF/Brain_v4.chr$chrom
# done
# cat ref/LDREF/Brain_v4.chr{1..20}.bim | cut -f2 >ref/LD_SNPs.txt

#####################
# Setup containers
####################

# ncbi/sra-tools does not provide bash so we can't use it with snakemake
# and it also breaks the fakeroot command, so we have to disable that while building it
mkdir -p images/
apptainer build --ignore-subuid --ignore-fakeroot-command images/sratools.sif containers/sratools.def

# General bioinformatics container (STAR, samtools, bcftools, htslib, plink2,
# plink, gatk4, beagle, python + libs) used by the alignment, QC, genotype
# processing, and phasing rules. Replaces the conda environment.
apptainer build images/bioinfo.sif containers/bioinfo.def

# GPU container for the tensorQTL rules. Only needed when the QTL steps are
# re-enabled in the Snakefile; build on a machine set up for the target CUDA.
# apptainer build --ignore-subuid --ignore-fakeroot-command images/tensorqtl.sif containers/tensorqtl.def
