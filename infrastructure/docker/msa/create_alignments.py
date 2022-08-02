# Original Copyright 2021 DeepMind Technologies Limited
# Modification Copyright 2022 # Copyright 2021 AlQuraishi Laboratory
# Modifications Copyright 2022 Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

import os
import logging
import argparse
import re
from msa.alignment_runner import AlignmentRunner

logging.basicConfig()
logger = logging.getLogger(__file__)
logger.setLevel(level=logging.INFO)

def main(args):
 # Create the output directory
    os.makedirs(args.output_dir, exist_ok=True)

    output_dir_base = args.output_dir
    if not os.path.exists(output_dir_base):
        os.makedirs(output_dir_base)
    alignment_dir = os.path.join(output_dir_base, "alignments")

    for fasta_file in list_files_with_extensions(args.fasta_dir, (".fasta", ".fa")):
        # Gather input sequences
        with open(os.path.join(args.fasta_dir, fasta_file), "r") as fp:
            data = fp.read()
    
        tags, seqs = parse_fasta(data)

        precompute_alignments(tags, seqs, alignment_dir, args)

def precompute_alignments(tags, seqs, alignment_dir, args):
    for tag, seq in zip(tags, seqs):
        tmp_fasta_path = os.path.join(args.output_dir, f"tmp_{os.getpid()}.fasta")
        with open(tmp_fasta_path, "w") as fp:
            fp.write(f">{tag}\n{seq}")

        local_alignment_dir = os.path.join(alignment_dir, tag)
        logger.info(f"Generating alignments for {tag}...")
        if not os.path.exists(local_alignment_dir):
            os.makedirs(local_alignment_dir)

        alignment_runner = AlignmentRunner(
            jackhmmer_binary_path=args.jackhmmer_binary_path,
            hhblits_binary_path=args.hhblits_binary_path,
            hhsearch_binary_path=args.hhsearch_binary_path,
            uniref90_database_path=args.uniref90_database_path,
            mgnify_database_path=args.mgnify_database_path,
            bfd_database_path=args.bfd_database_path,
            uniclust30_database_path=args.uniclust30_database_path,
            pdb70_database_path=args.pdb70_database_path,
            use_small_bfd=args.use_small_bfd,
            no_cpus=args.cpus,
        )
        alignment_runner.run(
            tmp_fasta_path, local_alignment_dir
        )

        # Remove temporary FASTA file
        os.remove(tmp_fasta_path)

def parse_fasta(data):
    data = re.sub('>$', '', data, flags=re.M)
    lines = [
        l.replace('\n', '')
        for prot in data.split('>') for l in prot.strip().split('\n', 1)
    ][1:]
    tags, seqs = lines[::2], lines[1::2]

    tags = [t.split()[0] for t in tags]

    return tags, seqs

def list_files_with_extensions(dir, extensions):
    return [f for f in os.listdir(dir) if f.endswith(extensions)]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "fasta_dir", type=str,
        help="Path to directory containing FASTA files, one sequence per file"
    )
    parser.add_argument(
        "--output_dir", type=str, default=os.getcwd(),
        help="""Name of the directory in which to output the prediction""",
    )
    parser.add_argument(
        "--cpus", type=int, default=4,
        help="""Number of CPUs with which to run alignment tools"""
    )
    parser.add_argument(
        '--uniref90_database_path', type=str, default=None,
    )
    parser.add_argument(
        '--mgnify_database_path', type=str, default=None,
    )
    parser.add_argument(
        '--pdb70_database_path', type=str, default=None,
    )
    parser.add_argument(
        '--uniclust30_database_path', type=str, default=None,
    )
    parser.add_argument(
        '--bfd_database_path', type=str, default=None,
    )
    parser.add_argument(
        '--jackhmmer_binary_path', type=str, default='/usr/bin/jackhmmer'
    )
    parser.add_argument(
        '--hhblits_binary_path', type=str, default='/usr/bin/hhblits'
    )
    parser.add_argument(
        '--hhsearch_binary_path', type=str, default='/usr/bin/hhsearch'
    )
    parser.add_argument(
        '--kalign_binary_path', type=str, default='/usr/bin/kalign'
    )    
    parser.add_argument(
        '--use_small_bfd', type=str, default=None,
    )        
    args = parser.parse_args()

    main(args)