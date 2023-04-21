import os
import requests
from Bio import AlignIO
from Bio.PDB.PDBList import PDBList
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.PDBIO import Select
import matplotlib.pyplot as plt
import py3Dmol
import io
import numpy as np
from batchfold.utils import residue_constants
from batchfold.utils import protein
import json

def download_rcsb_pdb_file(pdb_code, output_dir, model = None, chain = None, file_format="pdb"):
    """ Download a pdb file from rcsb.org """

    pdb_code = str.upper(pdb_code)
    pdbl = PDBList()
    os.makedirs(output_dir, exist_ok=True)
    ent_filename = pdbl.retrieve_pdb_file(
        pdb_code=pdb_code, file_format=file_format, pdir=output_dir, overwrite=True
    )

    if os.path.exists(ent_filename):
        pdb_filename = os.path.join(output_dir, pdb_code + ".pdb")
        os.rename(ent_filename, pdb_filename)

        if model is not None or chain is not None:
            extract_chain(pdb_filename, model, chain)
        return pdb_filename
    else:
        raise Exception("Desired structure doesn't exist at rcsb.org")

def download_rcsb_fasta_file(pdb_code, output_dir):
    """ Download a FASTA file from rcsb.org """
    pdb_code = str.upper(pdb_code)
    os.makedirs(output_dir, exist_ok=True)
    r = requests.get(f"https://www.rcsb.org/fasta/entry/{pdb_code}")
    fasta_filename = os.path.join(output_dir, pdb_code + ".fasta")
    with open(fasta_filename, "wb") as f:
        f.write(r.content)
    return fasta_filename

def get_bfactors(pdb_file):
    p = PDBParser()
    structure = p.get_structure("X", pdb_file)
    output = {}
    for model in structure:
        chains = []
        for chain in model:
            bfactors = []
            for residue in chain:
                if "CA" in residue.child_dict:
                    # print(residue["CA"].get_bfactor())
                    bfactors.append(residue["CA"].get_bfactor())
            chains.append(bfactors)
        output[model.id] = chains    
    return output

def overwrite_b_factors(pdb_str: str, bfactors: np.ndarray) -> str:
    """Overwrites the B-factors in pdb_str with contents of bfactors array.

    Args:
      pdb_str: An input PDB string.
      bfactors: A numpy array with shape [1, n_residues, 37]. We assume that the
        B-factors are per residue; i.e. that the nonzero entries are identical in
        [0, i, :].

    Returns:
      A new PDB string with the B-factors replaced.
    """
    if bfactors.shape[-1] != residue_constants.atom_type_num:
        raise ValueError(
            f"Invalid final dimension size for bfactors: {bfactors.shape[-1]}."
        )

    parser = PDBParser(QUIET=True)
    handle = io.StringIO(pdb_str)
    structure = parser.get_structure("", handle)

    curr_resid = ("", "", "")
    idx = -1
    for atom in structure.get_atoms():
        atom_resid = atom.parent.get_id()
        if atom_resid != curr_resid:
            idx += 1
            if idx >= bfactors.shape[0]:
                raise ValueError(
                    "Index into bfactors exceeds number of residues. "
                    "B-factors shape: {shape}, idx: {idx}."
                )
        curr_resid = atom_resid
        atom.bfactor = bfactors[idx, residue_constants.atom_order["CA"]]

    new_pdb = io.StringIO()
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)
    pdb_io.save(new_pdb)
    return new_pdb.getvalue()

def plot_banded_pdb(pdb_file, show_sidechains = False, width = 800, height = 600):
    with open(pdb_file) as f:
            best_pdb = f.read()
    target_protein = protein.from_pdb_string(best_pdb)
    plddt_list = target_protein.b_factors[:,0]
    atom_mask = target_protein.atom_mask
    banded_b_factors = []
    for plddt in plddt_list:
        for idx, (min_val, max_val, _) in enumerate(residue_constants.PLDDT_BANDS):
            if plddt >= min_val and plddt <= max_val:
                banded_b_factors.append(idx)
                break

    banded_b_factors = (
            np.array(banded_b_factors)[:, None] * atom_mask
    )

    to_visualize_pdb = overwrite_b_factors(best_pdb, banded_b_factors)
    # Color the structure by per-residue pLDDT
    color_map = {i: bands[2] for i, bands in enumerate(residue_constants.PLDDT_BANDS)}
    view = py3Dmol.view(width, height)
    view.addModelsAsFrames(to_visualize_pdb)
    style = {"cartoon": {"colorscheme": {"prop": "b", "map": color_map}}}
    if show_sidechains:
        style["stick"] = {}
    view.setStyle({"model": -1}, style)
    view.zoomTo()
    view.show()
    return None

def plot_plddt_legend():
    """Plots the legend for pLDDT."""

    thresh = [
        "Very low (pLDDT < 50)",
        "Low (70 > pLDDT > 50)",
        "Confident (90 > pLDDT > 70)",
        "Very high (pLDDT > 90)",
    ]

    colors = [x[2] for x in residue_constants.PLDDT_BANDS]

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

def plot_metrics(pdb_file, pae=None):
    with open(pdb_file) as f:
        best_pdb = f.read()
    target_protein = protein.from_pdb_string(best_pdb)
    plddt_list = target_protein.b_factors[:,0]

    fig = plt.figure(figsize=(8, 3))
    ax1 = fig.add_subplot(121)
    ax1.plot(plddt_list)
    ax1.set_title("Predicted LDDT")
    ax1.set_xlabel("Residue")
    ax1.set_ylabel("pLDDT")

    if pae is not None:
        ax2 = fig.add_subplot(122)
        ax2.imshow(pae, vmin=0.0, vmax=pae.max(), cmap="Greens_r")
        ax2.set_title("Predicted Aligned Error")
        ax2.set_xlabel("Scored residue")
        ax2.set_ylabel("Aligned residue")

def get_best_alphafold_model(ranking_file):
    with open(ranking_file, "r") as f:
        j = json.load(f)
    return(j["order"][0])

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


class SelectChain(Select):
    def __init__(self, chain):
        self.chain = chain

    def accept_chain(self, chain):
        if chain.__repr__ == self.chain:
            return 1
        else:
            return 0
    
def extract_chain(input_file, model = None, chain = None):
    p = PDBParser()
    structure = p.get_structure("X", input_file)
    if model is not None:
        structure = structure[model]
    else:
        structure = structure[0]
    if chain is not None:
        structure = structure[chain]
    io = PDBIO()
    io.set_structure(structure)
    new_filename = os.path.splitext(input_file)[0] + "_" + chain + ".pdb"
    io.save(new_filename)
    return new_filename
