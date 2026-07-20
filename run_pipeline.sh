module load miniconda
eval "$(/appl/miniconda/bin/conda shell.bash hook)"
conda activate ratgtex

bsub -e logs/snakemake.err \
    -o logs/snakemake.out \
    snakemake \
    -c 100 \
    -j 100 \
    --executor lsf \
    --default-resources lsf_project=ratgtex lsf_queue=normal mem_mb=4000 \
    --resources sra_downloads=3 \
    --use-singularity --singularity-args "-B /project/itmatlab/" \
    "$@"
