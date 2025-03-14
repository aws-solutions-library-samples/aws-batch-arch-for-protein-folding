# Original Copyright 2022 Facebook, Inc. and its affiliates.
# Modifications Copyright 2022 Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: MIT-0

import argparse
import json
import logging
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import os
import pyfastx
import torch
from transformers import AutoTokenizer, EsmForProteinFolding
from tqdm import tqdm

logging.basicConfig(
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%m/%d/%Y %H:%M:%S",
    level=logging.INFO,
)


def plot_pae(pae, output) -> None:
    fig = plt.figure()
    ax = fig.add_subplot()
    hcls_cmap = LinearSegmentedColormap.from_list(
        "hclscmap", ["#FFFFFF", "#007FAA", "#005276"]
    )
    _ = plt.imshow(pae, vmin=0.0, vmax=pae.max(), cmap=hcls_cmap)
    ax.set_title("Predicted Aligned Error")
    ax.set_xlabel("Scored residue")
    ax.set_ylabel("Aligned residue")
    fig.savefig(output)
    plt.close(fig)
    return None


def predict_structures(
    seqs: list,
    pretrained_model_name_or_path: str = "facebook/esmfold_v1",
    output_dir: str = "output",
    num_recycles: int = 4
):

    device = "cuda:0" if torch.cuda.is_available() else "cpu"

    tokenizer = AutoTokenizer.from_pretrained(pretrained_model_name_or_path)
    model = EsmForProteinFolding.from_pretrained(pretrained_model_name_or_path).to(
        device
    )

    logging.info(f"Predicting structures for {len(seqs)} sequences")
    for n, seq in tqdm(
        enumerate(seqs),
        desc=f"Generating structures",
    ):
        logging.info(f"Sequence {n+1} of {len(seqs)}")
        metrics = {"name": seq.name, "sequence": seq.seq, "sequence_length": len(seq.seq)}
        inputs = tokenizer(seq.seq, return_tensors="pt", add_special_tokens=False).to(
            device
        )

        inputs['num_recycles'] = num_recycles
        with torch.inference_mode():
            outputs = model(**inputs)

        output = {key: value.cpu() for key, value in outputs.items()}
        pdb_string = model.output_to_pdb(output)[0]
        output_dir = os.path.join(args.output_dir, str(n))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        output_file = os.path.join(output_dir, seq.name + ".pdb")
        header_str = f"REMARK   1\nREMARK   1 {seq.name}\n"
        with open(output_file, "w") as f:
            f.write(header_str)
            f.write(pdb_string)
        metrics.update(
            {
                "mean_plddt": round(torch.mean(output["plddt"]).item(), 3),
                "ptm": round(output["ptm"].item(), 3),
                "max_predicted_aligned_error": round(
                    output["max_predicted_aligned_error"].item(), 3
                ),
            }
        )
        torch.save(output, os.path.join(output_dir, seq.name + ".pt"))
        pae = output["predicted_aligned_error"]
        plot_pae(pae[0], os.path.join(output_dir, seq.name + ".png"))
        with open(os.path.join(output_dir, seq.name + ".json"), "w") as f:
            json.dump(metrics, f)
            f.write("\n")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "fasta",
        help="Path to input fasta file with sequences to process",
        type=str,
    )
    parser.add_argument(
        "--model-dir",
        help="ESMFold model to use",
        default="facebook/esmfold_v1",
        type=str,
    )
    parser.add_argument(
        "--pdb",
        help="(Optional) Path to output dir",
        default="output",
        type=str,
    )
    parser.add_argument(
        "--num-recycles",
        type=int,
        default=None,
        help="Number of recycles to run. Defaults to number used in training (4).",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=None,
        help="Chunks axial attention computation to reduce memory usage from O(L^2) to O(L). "
        "Equivalent to running a for loop over chunks of of each dimension. Lower values will "
        "result in lower memory usage at the cost of speed. Recommended values: 128, 64, 32. "
        "Default: None.",
    )
    parser.add_argument("--cpu-only", help="CPU only", action="store_true")
    parser.add_argument("--cpu-offload", help="Enable CPU offloading", action="store_true")

    args = parser.parse_args()
    seqs = [seq for seq in pyfastx.Fasta(args.fasta)]

    predict_structures(
        seqs,
        args.model_dir,
        args.pdb,
        args.num_recycles
    )
