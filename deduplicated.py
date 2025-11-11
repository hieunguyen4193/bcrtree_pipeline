#####------------------------------------------------------------------------#####
# This function was forked from the function deduplicated.py in gctree 
# https://github.com/matsengrp/gctree
#####------------------------------------------------------------------------#####
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import MultipleSeqAlignment
from collections import defaultdict, Counter
import os
import sys

#####------------------------------------------------------------------------#####
##### draft
#####------------------------------------------------------------------------#####
# aln_file = "./example/m11_IGHV9-3*01_IGHJ2*01_45.aln.fasta"
# root = "GL" # sequence ID of the root sequence (germline V-J gene sequences)
# frame=None
# id_abundances=False
# outputdir = "./example/output"
# output_name = None

#####------------------------------------------------------------------------#####
##### MAIN
#####------------------------------------------------------------------------#####
def get_parser(argv):
    parser = argparse.ArgumentParser(
        description=(
            " Deduplicate sequences in a fasta file, write to stdout in phylip format,"
            " and create a few other files (see arguments). An additional sequence"
            " representing the outgroup/root must be included"
            " (even if one or more observed sequences are identical to it)."
            " forked from gctree https://github.com/matsengrp/gctree and modified by work@hieunguyen.de"
            " In the CRC1382-BCRTree analysis project, sequence IDs in the input FASTA file must be integers"
            " representing the abundance of the sequence. We will keep this information as input to the downstream"
            " tree generation process."
        )
    )
    parser.add_argument(
        "--input",
        type=str,
        help="Fasta file with less than or equal to 10 characters unique header ID. "
        "Because dnapars will name internal nodes by integers a node name must include"
        "at least one non number character (unless using the id_abundances option).",
    )
    parser.add_argument(
        "--id_abundances",
        action="store_true",
        help="flag to interpret integer ids in input as abundances",
    )
    parser.add_argument(
        "--root",
        type=str,
        default="root",
        help=(
            "ID of the root sequence in fasta file. This ID will be"
            " used as the unique id for the root sequence, and any"
            " observed sequences matching the root sequence."
        ),
    )
    parser.add_argument(
        "--frame", type=int, default=None, help="codon frame", choices=(0, 1, 2, 3)
    )
    parser.add_argument(
        "--output", type=str, default=".", help="path to folder storing output"
    )
    parser.add_argument(
        "--output_name", type=str, default="test", help="identical name (sample ID) for all output files"
    )
    return parser.parse_args(argv)


def parse_sequence_file(aln_file, root, frame = None, id_abundances = False):
    ##### detect input file format
    if aln_file.endswith("fasta") or aln_file.endswith("fa"):
        aln_format = "fasta"
    elif aln_file.endswith("phylip") or aln_file.endswith("phy"):
        aln_format = "phylip"
    else:
        raise ValueError("unrecognized alignment file type: " + aln_file)
    aln = AlignIO.read(aln_file, aln_format)
    
    sequence_length = aln.get_alignment_length()
    if frame != 0:
        start = frame - 1
        end = start + 3 * ((sequence_length - start) // 3)
    else:
        start = 0
        end = sequence_length
    
    seqs_unique_counts = defaultdict(list) # by storing sequences as keys in the defaultdict, we keep
    # unique sequence only. We cound their abundances and assign to the value. 
    root_seq = None
    for seq in aln:
        seqstr = str(seq.seq)[start:end].upper()
        if seq.id == root:
            root_seq = seqstr
            if seqstr not in seqs_unique_counts:
                seqs_unique_counts[seqstr] = (
                        []
                )  # no observed root unless we see it elsewhere
        elif seq.id.isdigit() and id_abundances: # count the sequence abundance
            # if id is just an integer, assume it represents count of that sequence
            # sum the abundances of several sequences list.extend.
            seqs_unique_counts[seqstr].extend([seq.id for _ in range(int(seq.id))])
        else:
            seqs_unique_counts[seqstr].append(seq.id)
    
    # create the root seq and loop to count the number of appearances for each seq. 
    # that's why we use defaultdict here
    new_aln = MultipleSeqAlignment([SeqRecord(Seq(root_seq), id=root)])
    counts = {
        root: len(seqs_unique_counts[root_seq])
        }  # Add the count for the root sequence
    id_map = {root: [x for x in seqs_unique_counts[root_seq] if x != root]}
    id_map_seq = {root: [x for x in seqs_unique_counts[root_seq] if x != root]}
    
    del seqs_unique_counts[root_seq]  # Now delete the root so it does not appear twice
    for i, seq in enumerate(seqs_unique_counts, 1):
        new_id = "seq" + str(i)
        new_aln.append(SeqRecord(Seq(seq), id=new_id))
        counts[new_id] = len(seqs_unique_counts[seq])
        id_map[new_id] = seqs_unique_counts[seq]
        id_map_seq[new_id] = seq
    return new_aln, counts, id_map, id_map_seq

def main(args):
    ##### get input args via argparse
    aln_file = args.input
    root = args.root
    frame = args.frame
    id_abundances = args.id_abundances
    outputdir = args.output
    output_name = args.output_name

    os.system("mkdir -p {}".format(outputdir))
    ##### parse the input squence files, change sequence ID and count their abundances
    new_aln, counts, id_map, id_map_seq = parse_sequence_file(aln_file = aln_file,
                                                              root = root, 
                                                              frame = frame, 
                                                              id_abundances = id_abundances)
    
    ##### write the new sequence alignment files to outputdir/output_name phylip format file. 
    if output_name is None:
        input_name = str(aln_file).split("/")[-1]
        output_name = "".join(input_name.split(".")[0])

    print("output phylip format sequence file ...")
    outputfile = open(os.path.join(outputdir, "{}.phylip".format(output_name)), 'w')
    print(format(new_aln, "phylip"), file = outputfile)
    outputfile.close()

    print("output abundance file ...")
    ##### write abundance files
    with open(os.path.join(outputdir, "{}.abundance.csv".format(output_name)), "w") as f:
        for seqID, count in counts.items():
            print("{},{}".format(seqID, count), file=f)

    ##### write id_map file
    print("output id map file ...")
    with open(os.path.join(outputdir, "{}.id_map.csv".format(output_name)), "w") as f:
        for seq_id, cell_ids in id_map.items():
            print("{},{}".format(seq_id, ":".join(cell_ids)), file=f)

    with open(os.path.join(outputdir, "{}.id_map_seq.csv".format(output_name)), "w") as f:
        for seq_id, cell_ids in id_map_seq.items():
            print("{},{}".format(seq_id, cell_ids), file=f)
            
if __name__ == '__main__':
    main(get_parser(sys.argv[1:]))



