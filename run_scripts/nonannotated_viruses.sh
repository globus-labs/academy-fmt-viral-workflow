#!/bin/bash

# Input files

ANNOTATION_FILE="/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/06_annotate/annotated_AvRCv1.Merged_PredictedHosts.csv"
HITS_FILE="/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/05D_mergeblast/AVrC_allrepresentatives.fasta/clusterRes_rep_seq.fasta.txt"

# Output files

ANNOTATION_CONTIGS="/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/nonannotated_viral/annotation_contigs.txt"
HITS_CONTIGS="/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/nonannotated_viral/hits_contigs.txt"
ANNOTATION_CONTIGS_SORTED="/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/nonannotated_viral/annotation_contigs_sorted.txt"
HITS_CONTIGS_SORTED="/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/nonannotated_viral/hits_contigs_sorted.txt"
UNIQUE_HITS_CONTIGS="/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/nonannotated_viral/unique_hits_contigs.txt"

mkdir -p /xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/nonannotated_viral

# Step 1: Extract contig IDs from the annotation file

echo "Extracting contig IDs from annotation file..."
cut -d',' -f1 "$ANNOTATION_FILE" | tail -n +2 > "$ANNOTATION_CONTIGS"

# Step 2: Extract contig IDs from the hits file (second column) and remove duplicates

echo "Extracting contig IDs from hits file and removing duplicates..."
cut -f2 "$HITS_FILE" | sort | uniq > "$HITS_CONTIGS"

# Step 3: Sort the contig ID lists

echo "Sorting contig ID lists..."
sort "$ANNOTATION_CONTIGS" -o "$ANNOTATION_CONTIGS_SORTED"
sort "$HITS_CONTIGS" -o "$HITS_CONTIGS_SORTED"

# Step 4: Find contig IDs that are unique to the hits file (contig IDs only in the hits file but not in the annotation file)

echo "Finding unique contig IDs in the hits file..."
comm -23 "$HITS_CONTIGS_SORTED" "$ANNOTATION_CONTIGS_SORTED" > "$UNIQUE_HITS_CONTIGS"

# Step 5: Summary

echo "Summary of results:"
echo "Total contig IDs in annotation file: $(wc -l < "$ANNOTATION_CONTIGS")"
echo "Total contig IDs in hits file: $(wc -l < "$HITS_CONTIGS")"
echo "Total unique contig IDs in hits file: $(wc -l < "$UNIQUE_HITS_CONTIGS")"
echo "Unique contig IDs saved in $UNIQUE_HITS_CONTIGS"


