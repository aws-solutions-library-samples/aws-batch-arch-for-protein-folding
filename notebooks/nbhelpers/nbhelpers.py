# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""
Helper functions for the AWS-Alphafold notebook.
"""
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
from string import ascii_uppercase, ascii_lowercase
import py3Dmol
from nbhelpers import parsers, utils
from IPython import display
from ipywidgets import GridspecLayout
from ipywidgets import Output
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBList import PDBList
from Bio.PDB.Polypeptide import PPBuilder
import requests
import subprocess
import boto3
import json
from datetime import datetime


def download_pdb_file(pdb_code, output_dir, file_format="pdb"):
    pdb_code = str.upper(pdb_code)
    pdbl = PDBList()
    os.makedirs(output_dir, exist_ok=True)
    ent_filename = pdbl.retrieve_pdb_file(
        pdb_code=pdb_code, file_format=file_format, pdir=output_dir, overwrite=True
    )
    pdb_filename = os.path.join(output_dir, pdb_code + ".pdb")
    os.rename(ent_filename, pdb_filename)
    return pdb_filename


def download_fasta_file(pdb_code, output_dir):
    pdb_code = str.upper(pdb_code)
    os.makedirs(output_dir, exist_ok=True)
    r = requests.get(f"https://www.rcsb.org/fasta/entry/{pdb_code}")
    fasta_filename = os.path.join(output_dir, pdb_code + ".fasta")
    with open(fasta_filename, "wb") as f:
        f.write(r.content)
    return fasta_filename


def get_pdb_data(pdb_list):

    pdbl = PDBList()
    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    writer = PDBIO()
    ppb = PPBuilder()
    os.makedirs("data/pdb", exist_ok=True)
    os.makedirs("data/fasta", exist_ok=True)

    for id in pdb_list:
        structure_id, chain_id = id.split("_")
        print(f"structure ID is {structure_id}")
        print(f"chain ID is {chain_id}")

        filename = pdbl.retrieve_pdb_file(
            pdb_code=structure_id, file_format="pdb", pdir="data"
        )
        structure = parser.get_structure(structure_id, filename)
        chain = structure[0][chain_id]
        writer.set_structure(chain)
        pdb_file = "data/pdb/" + structure_id + "_" + chain_id + ".pdb"
        writer.save(pdb_file)
        os.remove(filename)
        seq_record = SeqRecord(
            seq=Seq(ppb.build_peptides(chain)[0].get_sequence()),
            id=id,
            description="",
            annotations={"molecule_type": "protein"},
        )
        fasta_file = "data/fasta/" + structure_id + "_" + chain_id + ".fasta"
        with open(fasta_file, "w") as f_out:
            SeqIO.write(seq_record, f_out, "fasta")


def list_files_by_extension(extension="", dir="."):
    """List files in dir with a given extension"""
    absdir = os.path.abspath(dir)
    return [
        os.path.join(absdir, f) for f in os.listdir(absdir) if f.endswith(extension)
    ]


# Color bands for visualizing plddt
PLDDT_BANDS = [
    (0, 50, "#FF7D45"),
    (50, 70, "#FFDB13"),
    (70, 90, "#65CBF3"),
    (90, 100, "#0053D6"),
]

def reduce_stockholm_file(sto_file):
    """Read in a .sto file and parse format it into a numpy array of the
    same length as the first (target) sequence
    """
    msa = AlignIO.read(sto_file, "stockholm")
    msa_arr = np.array([list(rec) for rec in msa])
    return msa_arr[:, msa_arr[0, :] != "-"]


def plot_msa_array(msa_arr, id=None):

    total_msa_size = len(msa_arr)

    if total_msa_size > 1:
        aa_map = {res: i for i, res in enumerate("ABCDEFGHIJKLMNOPQRSTUVWXYZ-")}
        msa_arr = np.array([[aa_map[aa] for aa in seq] for seq in msa_arr])
        plt.figure(figsize=(12, 3))
        plt.title(
            f"Per-Residue Count of Non-Gap Amino Acids in the MSA for Sequence {id}"
        )
        plt.plot(np.sum(msa_arr != aa_map["-"], axis=0), color="black")
        plt.ylabel("Non-Gap Count")
        plt.yticks(range(0, total_msa_size + 1, max(1, int(total_msa_size / 3))))

        return plt

    else:
        print("Unable to display MSA of length 1")
        return None


def plot_msa_folder(msa_folder, id=None):
    combined_msa = None
    with os.scandir(msa_folder) as it:
        for obj in it:
            obj_path = os.path.splitext(obj.path)
            if "pdb_hits" not in obj_path[0] and obj_path[1] == ".sto":
                msa_arr = reduce_stockholm_file(obj.path)
                if combined_msa is None:
                    combined_msa = msa_arr
                else:
                    combined_msa = np.concatenate((combined_msa, msa_arr), axis=0)
    if combined_msa is not None:
        print(f"Total number of aligned sequences is {len(combined_msa)}")
        plot_msa_array(combined_msa, id).show()
        return None
    else:
        return None


def plot_msa_output_folder(path, id=None):
    """Plot MSAs in a folder that may have multiple chain folders"""
    plots = []
    monomer = True
    with os.scandir(path) as it:
        for obj in it:
            if obj.is_dir():
                monomer = False
                plot_msa_folder(obj.path, id + " " + obj.name)
        if monomer:
            plot_msa_folder(path, id)
    return None


def plot_plddt_legend():
    """Plots the legend for pLDDT."""

    thresh = [
        "Very low (pLDDT < 50)",
        "Low (70 > pLDDT > 50)",
        "Confident (90 > pLDDT > 70)",
        "Very high (pLDDT > 90)",
    ]

    colors = [x[2] for x in PLDDT_BANDS]

    plt.figure(figsize=(2, 2))
    for c in colors:
        plt.bar(0, 0, color=c)
    plt.legend(thresh, frameon=False, loc="center", fontsize=20)
    plt.xticks([])
    plt.yticks([])
    ax = plt.gca()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    plt.title("Model Confidence", fontsize=20, pad=20)
    return plt


def get_best_alphafold2_model(ranking_debug_file):
    with open(ranking_debug_file, "r") as f:
        rankings = json.load(f)
    return rankings["order"][0]


def pdb_plot(pdb_path, show_sidechains=True):
    """Create a plot of PDB structures"""

    pkl_files = list_files_by_extension(extension=".pkl", dir=pdb_path)
    if len(pkl_files) > 1:  # AlphaFold
        best_model = get_best_alphafold2_model(
            os.path.join(pdb_path, "ranking_debug.json")
        )
        result = np.load(
            os.path.join(pdb_path, f"result_{best_model}.pkl"), allow_pickle=True
        )
        best_pdb_file = os.path.abspath(
            os.path.join(pdb_path, f"relaxed_{best_model}.pdb")
        )

    else:  # OpenFold
        result = np.load(pkl_files[0], allow_pickle=True)
        best_pdb_file = list_files_by_extension(extension="_relaxed.pdb", dir=pdb_path)[
            0
        ]

    pae_outputs = (
        result["predicted_aligned_error"],
        result["max_predicted_aligned_error"],
    )

    # Construct multiclass b-factors to indicate confidence bands
    # 0=very low, 1=low, 2=confident, 3=very high

    with open(best_pdb_file) as f:
        best_pdb = f.read()

    banded_b_factors = []
    for plddt in result["plddt"]:
        for idx, (min_val, max_val, _) in enumerate(PLDDT_BANDS):
            if plddt >= min_val and plddt <= max_val:
                banded_b_factors.append(idx)
                break
    if "final_atom_mask" in result:
        banded_b_factors = (
            np.array(banded_b_factors)[:, None] * result["final_atom_mask"]
        )
    else:
        banded_b_factors = (
            np.array(banded_b_factors)[:, None]
            * result["structure_module"]["final_atom_mask"]
        )

    to_visualize_pdb = utils.overwrite_b_factors(best_pdb, banded_b_factors)

    # Color the structure by per-residue pLDDT
    color_map = {i: bands[2] for i, bands in enumerate(PLDDT_BANDS)}
    view = py3Dmol.view(width=800, height=600)
    view.addModelsAsFrames(to_visualize_pdb)
    style = {"cartoon": {"colorscheme": {"prop": "b", "map": color_map}}}

    if show_sidechains:
        style["stick"] = {}
    view.setStyle({"model": -1}, style)
    view.zoomTo()

    grid = GridspecLayout(1, 2)
    out = Output()
    with out:
        view.show()
    grid[0, 0] = out

    out = Output()
    with out:
        plot_plddt_legend().show()
    grid[0, 1] = out

    display.display(grid)

    # Display pLDDT and predicted aligned error (if output by the model).
    if pae_outputs:
        num_plots = 2
    else:
        num_plots = 1

    plt.figure(figsize=[8 * num_plots, 6])
    plt.subplot(1, num_plots, 1)
    plt.plot(result["plddt"])
    plt.title("Predicted LDDT")
    plt.xlabel("Residue")
    plt.ylabel("pLDDT")

    if num_plots == 2:
        plt.subplot(1, 2, 2)
        pae, max_pae = pae_outputs
        plt.imshow(pae, vmin=0.0, vmax=max_pae, cmap="Greens_r")
        plt.colorbar(fraction=0.046, pad=0.04)
        plt.title("Predicted Aligned Error")
        plt.xlabel("Scored residue")
        plt.ylabel("Aligned residue")

    return plt


def run_tmscore(pdb1, pdb2):

    cmd = [
        "TMscore",
        "-seq",
        pdb1,
        pdb2,
    ]
    output = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    parse_float = lambda x: float(x.split("=")[1].split()[0])
    o = {}
    for line in output.stdout.split("\n"):
        line = line.rstrip()
        if line.startswith("RMSD"):
            o["rms"] = parse_float(line)
        if line.startswith("TM-score"):
            o["tms"] = parse_float(line)
        if line.startswith("GDT-TS-score"):
            o["gdt"] = parse_float(line)
    return o


def get_batch_logs(logStreamName, boto_session=boto3.session.Session()):

    """
    Retrieve and format logs for batch job.
    """

    logs_client = boto_session.client("logs")

    try:
        response = logs_client.get_log_events(
            logGroupName="/aws/batch/job", logStreamName=logStreamName
        )
    except logs_client.meta.client.exceptions.ResourceNotFoundException:
        return f"Log stream {logStreamName} does not exist. Please try again in a few minutes"

    logs = pd.DataFrame.from_dict(response["events"])
    logs.timestamp = logs.timestamp.transform(
        lambda x: datetime.fromtimestamp(x / 1000)
    )
    logs.drop("ingestionTime", axis=1, inplace=True)
    return logs


def get_last_batch_job_id(
    batch_environment, job_name, boto_session=boto3.session.Session()
):
    batch = boto_session.client("batch")
    jobs = []
    for queue in batch_environment.get_stack_outputs(filter="JobQueue").values():
        hits = batch.list_jobs(
            jobQueue=queue, filters=[{"name": "JOB_NAME", "values": [job_name]}]
        ).get("jobSummaryList", [])
        jobs.extend(hits) if hits != [] else next

    if jobs == []:
        return ""
    else:
        return jobs[-1].get("jobId", [])
