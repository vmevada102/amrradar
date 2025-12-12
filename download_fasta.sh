#!/bin/bash

# Check if accession number is provided
if [ -z "$1" ]; then
  echo "Usage: $0 <accession_number>"
  exit 1
fi

ACCESSION=$1
OUTPUT_DIR="results"
OUTPUT_FILE="$OUTPUT_DIR/sequence.fasta"

# Create results directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Download FASTA file from NCBI using efetch
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACCESSION}&rettype=fasta&retmode=text" -o "$OUTPUT_FILE"

# Confirm download
if [ -s "$OUTPUT_FILE" ]; then
  echo "FASTA file for accession $ACCESSION saved to $OUTPUT_FILE"
else
  echo "Failed to download FASTA file for accession $ACCESSION"
  exit 1
fi