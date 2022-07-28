# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""
Helper functions for the AWS-Alphafold notebook.
"""
from Bio import SeqIO
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

def get_pdb_data(pdb_list):
    
    pdbl = PDBList() 
    parser = PDBParser(PERMISSIVE = True, QUIET = True) 
    writer = PDBIO()
    ppb = PPBuilder()

    for id in pdb_list:
        structure_id, chain_id = id.split("_")
        print(f"structure ID is {structure_id}")
        print(f"chain ID is {chain_id}")
        
        filename = pdbl.retrieve_pdb_file(pdb_code=structure_id, file_format="pdb", pdir="data")
        structure = parser.get_structure(structure_id, filename)
        description = structure.header["name"]
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
    """ List files in dir with a given extension"""
    absdir = os.path.abspath(dir)
    return [os.path.join(absdir,f) for f in os.listdir(absdir) if f.endswith(extension)]

# Color bands for visualizing plddt
PLDDT_BANDS = [
  (0, 50, '#FF7D45'),
  (50, 70, '#FFDB13'),
  (70, 90, '#65CBF3'),
  (90, 100, '#0053D6')
]    

def msa_plot(id, sto_path, mgnify_max_hits=501):
    """ Create a plot of MSA data."""
    msas = []
    deletion_matrices = []
    full_msa = []
    for sto_file in list_files_by_extension(".sto", dir=sto_path):
        db_name = os.path.basename(sto_file).split("_")[0]
        unsorted_results = []
        with open(sto_file, "r") as f:
            result = f.read()
            msa, deletion_matrix, target_names = parsers.parse_stockholm(result)
            zipped_results = zip(msa, deletion_matrix, target_names)
            zipped_results = [x for x in zipped_results if x[2] != id]
            unsorted_results.extend(zipped_results)
            db_msas, db_deletion_matrices, _ = zip(*unsorted_results)
            if db_msas:
                if db_name == 'mgnify':
                    db_msas = db_msas[:mgnify_max_hits]
                    db_deletion_matrices = db_deletion_matrices[:mgnify_max_hits]
                full_msa.extend(db_msas)
                msas.append(db_msas)
                deletion_matrices.append(db_deletion_matrices)
                msa_size = len(set(db_msas))
                print(f'{msa_size} Sequences Found in {db_name}')
                
    deduped_full_msa = list(dict.fromkeys(full_msa))
    total_msa_size = len(deduped_full_msa)
    print(f'\n{total_msa_size} Sequences Found in Total\n')

    import numpy as np
    aa_map = {restype: i for i, restype in enumerate('ABCDEFGHIJKLMNOPQRSTUVWXYZ-')}
    msa_arr = np.array([[aa_map[aa] for aa in seq] for seq in deduped_full_msa])
    num_alignments, num_res = msa_arr.shape

    fig = plt.figure(figsize=(12, 3))
    plt.title(f'Per-Residue Count of Non-Gap Amino Acids in the MSA for {id}')
    plt.plot(np.sum(msa_arr != aa_map['-'], axis=0), color='black')
    plt.ylabel('Non-Gap Count')
    plt.yticks(range(0, num_alignments + 1, max(1, int(num_alignments / 3))))
    # plt.show()
    return plt

def plot_plddt_legend():
    """Plots the legend for pLDDT."""

    thresh = [
                'Very low (pLDDT < 50)',
                'Low (70 > pLDDT > 50)',
                'Confident (90 > pLDDT > 70)',
                'Very high (pLDDT > 90)']

    colors = [x[2] for x in PLDDT_BANDS]

    plt.figure(figsize=(2, 2))
    for c in colors:
        plt.bar(0, 0, color=c)
    plt.legend(thresh, frameon=False, loc='center', fontsize=20)
    plt.xticks([])
    plt.yticks([])
    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    plt.title('Model Confidence', fontsize=20, pad=20)
    return plt


def pdb_plot(pdb_path, show_sidechains = True):
    """ Create a plot of PDB structures"""

    pkl_file = list_files_by_extension(extension=".pkl", dir=pdb_path)[0]
    result = np.load(pkl_file, allow_pickle=True)

    pae_outputs = (
          result['predicted_aligned_error'],
          result['max_predicted_aligned_error']
      )

    # Construct multiclass b-factors to indicate confidence bands
    # 0=very low, 1=low, 2=confident, 3=very high    
    relaxed_pdb = list_files_by_extension(extension="_relaxed.pdb", dir=pdb_path)[0]

    with open(relaxed_pdb) as f:
        best_pdb = f.read()
    
    banded_b_factors = []
    for plddt in result['plddt']:
        for idx, (min_val, max_val, _) in enumerate(PLDDT_BANDS):
            if plddt >= min_val and plddt <= max_val:
                banded_b_factors.append(idx)
                break
    banded_b_factors = np.array(banded_b_factors)[:, None] * result['final_atom_mask']
    to_visualize_pdb = utils.overwrite_b_factors(best_pdb, banded_b_factors)


    # Color the structure by per-residue pLDDT
    color_map = {i: bands[2] for i, bands in enumerate(PLDDT_BANDS)}
    view = py3Dmol.view(width=800, height=600)
    view.addModelsAsFrames(to_visualize_pdb)
    style = {'cartoon': {'colorscheme': {'prop': 'b','map': color_map}}}

    if show_sidechains:
        style['stick'] = {}
    view.setStyle({'model': -1}, style)
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
    plt.plot(result['plddt'])
    plt.title('Predicted LDDT')
    plt.xlabel('Residue')
    plt.ylabel('pLDDT')

    if num_plots == 2:
        plt.subplot(1, 2, 2)
        pae, max_pae = pae_outputs
        plt.imshow(pae, vmin=0., vmax=max_pae, cmap='Greens_r')
        plt.colorbar(fraction=0.046, pad=0.04)
        plt.title('Predicted Aligned Error')
        plt.xlabel('Scored residue')
        plt.ylabel('Aligned residue')

    return plt