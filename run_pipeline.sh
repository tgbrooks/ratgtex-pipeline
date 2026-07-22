bsub -e logs/snakemake.err \
    -o logs/snakemake.out \
    uv run snakemake \
    -c 100 \
    -j 100 \
    --executor lsf \
    --default-resources lsf_project=ratgtex lsf_queue=rhel9 mem_mb=4000 \
    --resources sra_downloads=3 \
    --use-singularity --singularity-args "-B /project/itmatlab/" \
    "$@"
