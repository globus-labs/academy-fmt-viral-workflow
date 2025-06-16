#!/bin/bash


# Input files

FASTA_FILE="/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/query/clusterRes_rep_seq.fasta"
HITS_FILE="/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/05D_mergeblast/AVrC_allrepresentatives.fasta/clusterRes_rep_seq.fasta.txt"

# Output files

ALL_NODES="/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/novel_viral/all_nodes.txt"
HITS_NODES="/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/novel_viral/hits_nodes.txt"
ALL_NODES_SORTED="/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/novel_viral/all_nodes_sorted.txt"
HITS_NODES_SORTED="/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/novel_viral/hits_nodes_sorted.txt"
NOVEL_NODES="/xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/novel_viral/novel_nodes.txt"

mkdir -p /xdisk/bhurwitz/virus_hunting/kolodisner/viral_detection_pipeline/results/novel_viral

# Step 1: Extract node IDs from the FASTA file

echo "Extracting node IDs from FASTA file..."
grep '^>' "$FASTA_FILE" | cut -d' ' -f1 | sed 's/^>//' > "$ALL_NODES"

# Step 2: Extract node IDs from the hits file

echo "Extracting node IDs from hits file..."
cut -f1 "$HITS_FILE" | sort | uniq > "$HITS_NODES"

# Step 3: Sort the node ID lists

echo "Sorting node ID lists..."
sort "$ALL_NODES" -o "$ALL_NODES_SORTED"
sort "$HITS_NODES" -o "$HITS_NODES_SORTED"

# Step 4: Find novel nodes (nodes only in the first list)

echo "Finding novel nodes..."
comm -23 "$ALL_NODES_SORTED" "$HITS_NODES_SORTED" > "$NOVEL_NODES"

# Step 5: Summary

echo "Summary of results:"
echo "Total nodes in FASTA file: $(wc -l < "$ALL_NODES")"
echo "Total nodes in hits file: $(wc -l < "$HITS_NODES")"
echo "Total novel nodes: $(wc -l < "$NOVEL_NODES")"

echo "Novel nodes saved in $NOVEL_NODES"


