#!/usr/bin/env bash
set -euo pipefail

# One-time cluster setup: create directories and download shared data.
# Expects BIOINF_DATA_DIR, BIOINF_RUNS_DIR, BIOINF_GNOMAD_DIR in ~/.secrets.env

source "$HOME/.secrets.env"

echo "Creating directories..."
mkdir -p "$BIOINF_DATA_DIR" "$BIOINF_RUNS_DIR" "$BIOINF_GNOMAD_DIR"

echo "Downloading ClinVar VCF..."
cd "$BIOINF_DATA_DIR"
wget -nc https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget -nc https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

echo "Setup complete."
echo "  Data dir:   $BIOINF_DATA_DIR"
echo "  Runs dir:   $BIOINF_RUNS_DIR"
echo "  gnomAD dir: $BIOINF_GNOMAD_DIR"
