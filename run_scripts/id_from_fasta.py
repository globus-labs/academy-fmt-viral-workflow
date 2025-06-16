#!/usr/bin/env python3
"""
Author : Me  <bhurwitz@junonia.hpc.arizona.edu>
Date   : 2024-07-03
Purpose: get ids from fasta
"""

import argparse
from Bio import SeqIO

# --------------------------------------------------
def get_args():
    """Get command-line arguments"""

    parser = argparse.ArgumentParser(
        description='get ids from fasta',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-c',
                        '--contigs',
                        help='A contigs file',
                        metavar='str',
                        default=None)

    parser.add_argument('-d',
                        '--dvf',
                        help='An output file from deepvirfinder',
                        metavar='FILE',
                        type=argparse.FileType('rt'),
                        default=None)
    
    parser.add_argument('-o',
                        '--outfile',
                        help='Output file',
                        metavar='FILE',
                        type=argparse.FileType('wt'),
                        default='dvf.fasta')    

    return parser.parse_args()


# --------------------------------------------------
def main():
    """Make a jazz noise here"""

    args = get_args()

    ids = read_ids_from_file(args.dvf)
    print_fasta_sequences(args.contigs, ids, args.outfile)


def read_ids_from_file(dvf):
    # Read the IDs from the given file and return them as a list
    ids = []
    for line in dvf:
        name, length, score, pvalue = line.split("\t")
        ids.append(name)
    ids = ids[1:]
    return ids

def print_fasta_sequences(contigs, ids, outfile): 
    # Use SeqIO to parse the FASTA file
    for record in SeqIO.parse(contigs, "fasta"):
        # Check if the record's ID is in the list of IDs to print
        if record.id in ids:
            # Print the ID and sequence in FASTA format
            SeqIO.write(record , outfile , "fasta")
   

# --------------------------------------------------
if __name__ == '__main__':
    main()
