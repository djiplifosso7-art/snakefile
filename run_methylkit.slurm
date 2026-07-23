#!/bin/bash
#SBATCH --job-name=methylkit_all
#SBATCH --output=logs/methylkit_all_%j.out
#SBATCH --error=logs/methylkit_all_%j.err
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G

set -euo pipefail

PROJECT="/env/cns/proj/projet_DWS/Pipeline_methylation"
cd "$PROJECT"

mkdir -p logs \
  results/methylkit/input \
  results/methylkit/diff \
  results/methylkit/tables \
  results/methylkit/plots \
  results/methylkit/logs

module load micromamba 2>/dev/null || true
micromamba activate smk9_modkit

module load bioconductor 2>/dev/null || true

echo "===== ENVIRONNEMENT ====="
echo "Conda env: ${CONDA_DEFAULT_ENV:-no_conda_env}"

echo "Snakemake:"
which snakemake
snakemake --version

echo "Rscript:"
which Rscript
Rscript --version

echo "===== TEST PACKAGES R ====="
Rscript -e 'library(readr); library(dplyr); library(tidyr); library(ggplot2); library(methylKit); cat("Packages R + methylKit OK\n")'

echo "===== DRY-RUN FIGURE FINALE ====="

snakemake -s snakefile -np \
  results/methylkit/tables/methylkit_hyper_hypo_counts_all.tsv \
  results/methylkit/plots/methylkit_hyper_hypo_all_comparisons.png \
  --rerun-triggers mtime \
  --latency-wait 300

echo "===== LANCEMENT FIGURE FINALE ====="

snakemake -s snakefile --cores 1 -p \
  results/methylkit/tables/methylkit_hyper_hypo_counts_all.tsv \
  results/methylkit/plots/methylkit_hyper_hypo_all_comparisons.png \
  --rerun-incomplete \
  --rerun-triggers mtime \
  --latency-wait 300

echo "===== FIN ====="

ls -lh results/methylkit/tables/methylkit_hyper_hypo_counts_all.tsv
ls -lh results/methylkit/plots/methylkit_hyper_hypo_all_comparisons.png

echo "===== APERÇU TABLE FINALE ====="
head results/methylkit/tables/methylkit_hyper_hypo_counts_all.tsv

snakemake -s snakefile \
  -j 1 \
  -p \
  results/global_methylation/global_methylation_from_beds.tsv \
  results/global_methylation/global_methylation_summary.tsv \
  results/global_methylation/global_methylation_violin.png \
  --rerun-triggers mtime \
  --latency-wait 300
