__author__ = "Ryan Kibler, Sergey Ovchinnikov, Nate Bennet, Philip Leung, Adam Broerman"  # TODO others?
# most of the code is copied from krypton's colabfold https://colab.research.google.com/drive/1teIx4yDyrSm0meWvM9X9yYNLHGO7f_Zy#scrollTo=vJxiCSLc7IWD
# The initial guess stuff is from Nate Bennett with maybe some helper code from Adam Broerman
# pae code is lifted from Nate
# it contains alphafold2-multimer but don't use it
# krypton is basically lead author without knowing it

from asyncio import format_helpers
import time

time_checkpoint = time.time()
import argparse
import os
import itertools

# from Bio import SeqIO
parser = argparse.ArgumentParser()


# This hack is probably unnecessary with AF2-multimer since they've switched to jax for feature processing
# tell Jax and Tensorflow to use the same memory. This allows us to run larger structures
os.environ["TF_FORCE_UNIFIED_MEMORY"] = "1"
os.environ["XLA_PYTHON_CLIENT_MEM_FRACTION"] = "2.0"

SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))


def validate_file(parser, path):
    """
    Check for file existance and read files first so that we can fail early before loading alphafold, etc

    https://stackoverflow.com/a/11541450
    """
    if not os.path.exists(path):
        parser.error("The file %s does not exist!" % path)
    else:
        if (
            path.endswith(".pdb")
            or path.endswith(".pdb.gz")
            or path.lower().endswith(".fasta")
            or path.lower().endswith(".fa")
            or path.lower().endswith(".silent")
        ):
            return path
        else:
            parser.error(
                "Only PDB files, silent files, and FASTA files are allowed. You supplied: %s"
                % path
            )


parser.add_argument(
    "input_files",
    metavar="PATH",
    nargs="+",
    type=lambda x: validate_file(parser, x),
    help="Paths to PDB files or FASTA files to run AlphaFold2 predictions on. All chains in a PDB file will be predicted as a multichain prediction. To specify chainbreaks in FASTA format, separate sequences with '/' or ':'",
)

parser.add_argument(
    "--info_collector_path",
    type=str,
    default="/net/scratch/db/af2/data", #TODO where should this be?
    help="Path where InfoCollector will drop off data. Default = '/net/scratch/db/af2/data'.",
)

parser.add_argument(
    "--info_collector_config",
    type=str,
    default="/net/scratch/db/af2/config.json", #TODO where should this be?
    help="Path to the InfoCollector configuration file. Default = '/net/scratch/db/af2/config.json'.",
)

# could try using a type here (like input files) to assert that the value is greater than 1. Instead right now we assert below.
parser.add_argument(
    "--mock_msa_depth",
    default=1,
    help="fake the msa. Lower is faster, but potentially less accurate. Range [1,inf). AF2 default is 512. Our Default = 1.",
    type=int,
)

parser.add_argument(
    "--models",
    choices=["1", "2", "3", "4", "5", "all"],
    default="all",
    nargs="+",
    help="Deepmind provided five sets of weights/models. You can choose any combination of models to run.",
)

parser.add_argument(
    "--type",
    choices=["monomer", "monomer_ptm", "multimer", "multimer_v2"],
    default="monomer_ptm",
    #help="The flavor of alphafold weights to use. 'monomer' is the original AF2. 'ptm' is the original AF2 with an extra head that predicts pTMscore. 'multimer' is AF2-Multimer. 'multimer_v2' is updated AF2-Multimer that is supposed to yield fewer clashes. The use of multimer weights with standard AF2 probably won't work",
    help="This option does not do anything anymore. The code will always use monomer_ptm weights. If you choose anything other than that, the program will exit without running predictions."
)

parser.add_argument(
    "--version",
    choices=["monomer", "multimer"],
    default="monomer",
    #help="The version of AF2 Module to use. Both versions can predict both mulimers. When used to predict multimers, the 'monomer' version is equivalent to AF2-Gap. The 'multimer' versions are equivalent to AF2-Multimer and should not be used with the monomer weight types.",
    help="This option does not do anything anymore. The code will always use monomer code. If you choose anything other than that, the program will exit without running predictions. If you want to use the Multimer code, use ColabFold instead. "
)

parser.add_argument(
    "--nstruct",
    help="Number of independent outputs to generate PER MODEL. It will make predictions with seeds starting at 'seed_start' and increasing by one until n outputs are generated (like seed_range = range(seed_start,seed_start + nstruct)). Default=1",
    default=1,
    type=int,
)

parser.add_argument(
    "--seed_start", type=int, help="Seed to start at. Default=0", default=0
)

parser.add_argument(
    "--num_ensemble",
    type=int,
    default=1,
    help="number of times to process the input features and combine. default = 1. Deepmind used 8 for casp. Expert Option.",
)

parser.add_argument(
    "--max_recycles",
    type=int,
    default=3,
    help="max number of times to run evoformer. Default is 3. Single domain proteins need fewer runs. Multidomain or PPI may need more",
)

parser.add_argument(
    "--recycle_tol",
    type=float,
    default=0.0,
    help="Stop recycling early if CA-RMSD difference between current output and previous is < recycle_tol. Default = 0.0 (no early stopping)",
)

# #An idea in current colab fold.
# parser.add_argument(
#     "--prediction_threshold",
#     nargs=2,
#     metavar=('value','type'),
#     help="Continue recycling until the prediction is above the threshold or the num_recycles == max_recycles. Type choices are ['mean_plddt','mean_pae','rmsd_prev']",
# )

# unknown if this currently works
# parser.add_argument("--show_images", action="store_true")

parser.add_argument(
    "--output_pae",
    action="store_true",
    help="dump the PAE matrix to disk. This is useful for investigating interresidue relationships.",
)

parser.add_argument(
    "--output_summary",
    action="store_true",
    help="write a 1-line summary of each prediction to disk under output_dir named 'reports.txt'.",
)

# # unknown if this currently works
# parser.add_argument(
#     "--save_intermediates",
#     action="store_true",
#     help="save intermediate structures between recycles. This is useful for making folding movies/trajectories",
# )

parser.add_argument(
    "--amber_relax",
    action="store_true",
    help="Amber relax is unsupported and enabling this option will cause the program to exit. This option is left in for backwards compatibility.",
)

parser.add_argument(
    "--overwrite",
    action="store_true",
    help="overwrite existing files. Default is to skip predictions which would result in files that already exist. This is useful for checkpointing and makes the script more backfill friendly.",
)
parser.add_argument(
    "--initial_guess",
    nargs="?",
    const=True,
    default=False,
    help="use the initial guess from the input PDB file. This is useful for trying to focus predictions toward a known conformation. If no path is provided, the input_file must be a PDB or silent. If a path is provided, the input must be a fasta.",
)
parser.add_argument(
    "--reference_pdb",
    type=str,
    help="reference PDB to use for RMSD calculations. Coordinates (after alignment) and chain order will be updated to that of this reference, unless the input_files are PDB files",
)

parser.add_argument(
    "--simple_rmsd",
    action="store_true",
    #help="compute RMSD directly with the alphafold prediction and without trying to rearrange chain orders.",
    help="This option does not do anything anymore. The code will always use MMalign for superimposing outputs."
)

# sidechain_relax_parser = parser.add_mutually_exclusive_group(required=False)
# sidechain_relax_parser.add_argument("--amber_relax",help="run Amber relax on each output prediction")
# sidechain_relax_parser.add_argument("--rosetta_relax",help="run Rosetta relax (sidechain only) on each output prediction")

parser.add_argument(
    "--enable_dropout",
    action="store_true",
    help="Introduce structural diversity by enabling dropout",
)
parser.add_argument(
    "--pct_seq_mask",
    type=float,
    default=0.15,
    help="percent of sequence to make during inference. Default = 0.15. Setting to 0 might reduce prediction stocasticity.",
)

parser.add_argument(
    "--out_dir",
    type=str,
    default="output/",
    help="Directory to output models and data.",
)

# Pass the location of alphafold_weights as an argument
parser.add_argument(
    "--alphafold_weights",
    type=str,
    required=True,
    help="Path to the AlphaFold data directory.",
)

args = parser.parse_args()

ALPHAFOLD_DATADIR = args.alphafold_weights
assert os.path.exists(ALPHAFOLD_DATADIR), f"AlphaFold data directory does not exist: {ALPHAFOLD_DATADIR}"


#manage removed options
if args.version != parser.get_default("version") or args.type != parser.get_default("type"):
    exit("ERROR: multimer functionality is deprecated because AF2-multimer performs poorly on single sequence MSAs. Use colabfold instead to run multimer with MSAs. It is up to you to decide if it is theoretically/morally correct to use MSAs with de novo proteins. The non-ptm monomer weights have also been removed as they offer no benefit over the ptm weights.")

if args.amber_relax:
    print("ERROR: Amber relax is currently broken and I don't intend to fix it.")
    exit(1)

#adding this to keep code working later on while I figure out how to make it work
args.save_intermediates = False

assert args.mock_msa_depth > 0

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / "silent_tools"))
import silent_tools  # installed as a submodule
from dataclasses import dataclass
from typing import Union, Tuple, Dict
import numpy as np
# from matplotlib import pyplot as plt

from info_collection import InfoCollector

# plt.switch_backend("agg")

os.makedirs(args.out_dir, exist_ok=True)




#from David Juergens
def np_kabsch(A,B):
    """
    Numpy version of kabsch algorithm. Superimposes B onto A

    Parameters:
        (A,B) np.array - shape (N,3) arrays of xyz crds of points


    Returns:
        rms - rmsd between A and B
        R - rotation matrix to superimpose B onto A
        rB - the rotated B coordinates
    """
    A = np.copy(A)
    B = np.copy(B)

    def centroid(X):
        # return the mean X,Y,Z down the atoms
        return np.mean(X, axis=0, keepdims=True)

    def rmsd(V,W, eps=0):
        # First sum down atoms, then sum down xyz
        N = V.shape[-2]
        return np.sqrt(np.sum((V-W)*(V-W), axis=(-2,-1)) / N + eps)


    N, ndim = A.shape

    # move to centroid
    A = A - centroid(A)
    B = B - centroid(B)

    # computation of the covariance matrix
    C = np.matmul(A.T, B)

    # compute optimal rotation matrix using SVD
    U,S,Vt = np.linalg.svd(C)


    # ensure right handed coordinate system
    d = np.eye(3)
    d[-1,-1] = np.sign(np.linalg.det(Vt.T@U.T))

    # construct rotation matrix
    R = Vt.T@d@U.T

    # get rotated coords
    rB = B@R

    # calculate rmsd
    rms = rmsd(A,rB)

    return rms, rB, R


# test_pdb_1 = "/home/rdkibler/megahelix.pdb"
# test_pdb_2 = "/home/rdkibler/1ubq.pdb"
# test_pdb_3 = "/home/rdkibler/4ncu.pdb"

from collections import defaultdict
import numpy as np

aa_3to1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLU": "E",
    "GLN": "Q",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}

class ParsedPDB:
    """
    A simple class for storing the parsed contents of a PDB file and handling PDB manipulations.
    """
    def __init__(self):
        self.coords = np.array([])
        self.atom_details = []
        self.name = None

    def copy(self):
        new = ParsedPDB()
        new.coords = np.copy(self.coords)
        new.atom_details = [list(atom) for atom in self.atom_details]
        new.name = self.name
        return new

    def parse_pdbstr(self, pdbstr):
        self.coords = []
        self.atom_details = []
        for line in pdbstr.split("\n"):
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom = line.startswith("ATOM")
                #00000000001111111111222222222233333333334444444444555555555566666666667777777777
                #01234567890123456789012345678901234567890123456789012345678901234567890123456789
                #ATOM      4  CA  GLY B   2      -4.825  -1.541  -6.447  1.00  0.00           C
                atom_num = int(line[6:11])
                atom_name = line[12:16]
                resnum = int(line[22:26])
                resname = line[17:20]
                chain = line[21]
                
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                try:
                    occupancy = float(line[54:60])
                except ValueError:
                    occupancy = 0.0
                
                try:
                    bfactor = float(line[60:66])
                except ValueError:
                    bfactor = 0.0

                try:    
                    element = line[76:78]
                except ValueError:
                    element = "  "
                
                try:
                    charge = line[78:80]
                except ValueError:
                    charge = "  "
                
                self.coords.append([x,y,z])
                self.atom_details.append([atom, atom_num, atom_name, resnum, resname, chain, bfactor, element, occupancy, charge])

        #convert coords to a numpy array
        self.coords = np.array(self.coords)


    def parse_pdb_file(self,path):
        self.name = path.split("/")[-1].split(".pdb")[0]

        with open(path) as f:
            pdbstr = f.read()
        self.parse_pdbstr(pdbstr)

    def reorder_records(self):
        #rearrange the coords and atom_details so they are in chain-sorted order
        chains = sorted(list(set([atom[5] for atom in self.atom_details])))
        new_coords = []
        new_atom_details = []

        coords_by_chain = defaultdict(list)
        atom_details_by_chain = defaultdict(list)
        for xyz,atom in zip(self.coords,self.atom_details):
            coords_by_chain[atom[5]].append(xyz)
            atom_details_by_chain[atom[5]].append(atom)

        for chain in chains:
            new_coords.extend(coords_by_chain[chain])
            new_atom_details.extend(atom_details_by_chain[chain])

        self.coords = np.array(new_coords)
        self.atom_details = new_atom_details

        self.renumber()

    def remap_chains(self, chain_mapping):
        #chain_mapping is a dict of old_chain:new_chain
        for atom in self.atom_details:
            atom[5] = chain_mapping[atom[5]]
    
        self.reorder_records()

    def get_bfactors(self):
        #return per-residue b-factors?
        bfacs = []
        for atom in self.atom_details:
            if atom[0] == "ATOM" and atom[2] == " CA ":
                bfacs.append(atom[6])
        return np.array(bfacs)

    def get_pdbstr(self):
        atom_num = 1
        buffer = ""
        for xyz,atom in zip(self.coords,self.atom_details):
            x,y,z = xyz
            #0:is_atom
            #1:atom_num
            #2:atom_name
            #3:resnum
            #4:resname
            #5:chain
            #00000000001111111111222222222233333333334444444444555555555566666666667777777777
            #01234567890123456789012345678901234567890123456789012345678901234567890123456789
            #ATOM      4  CA  GLY B   2      -4.825  -1.541  -6.447  1.00  0.00           C
            

            atom_name = atom[2]
            resnum = atom[3]
            resname = atom[4]
            chain = atom[5]
            bfactor = atom[6]
            element = atom[7]
            occupancy = atom[8]
            charge = atom[9]
            #resnum = f"{resnum_int:>4}"
            ATOM_HETATM = "ATOM  " if atom[0] else "HETATM"

            buffer += f"{ATOM_HETATM}{atom_num:>5} {atom_name} {resname} {chain}{resnum:>4}    {x:>8.3f}{y:>8.3f}{z:>8.3f}{occupancy:>6.2f}{bfactor:>6.2f}          {element:>2}{charge:>2}\n"
            atom_num += 1
        return buffer

    def make_pdb_file(self, path):
        with open(path, 'w') as f:
            f.write(self.get_pdbstr())
            

    def get_ca_indices_bychain(self):
        #return a dict of lists of indices of CA atoms
        ca_indices_bychain = defaultdict(list)
        for i,atom in enumerate(self.atom_details):
            if atom[2] == ' CA ':
                ca_indices_bychain[atom[5]].append(i)
        return ca_indices_bychain

    def renumber(self):
        previous_resid = None
        new_resid = 0

        #renumber the residues in the order they appear
        for i,atom in enumerate(self.atom_details):
            resid = int(atom[3])
            if resid != previous_resid:
                new_resid += 1
            self.atom_details[i][3] = str(new_resid)
            previous_resid = resid

    def get_chain_range_map(self):
        #returns a dict of chain:(resnum_start:resnum_stop)
        chain_range_map = {}
        for chain,indices in self.get_ca_indices_bychain().items():
            indices = sorted(indices)
            index_start = indices[0]
            index_stop = indices[-1]
            resnum_start = int(self.atom_details[index_start][3])
            resnum_stop = int(self.atom_details[index_stop][3])
            chain_range_map[chain] = (resnum_start, resnum_stop)
        return chain_range_map
    
    def get_CA_coords(self):
        CA_indices = []
        for i,atom in enumerate(self.atom_details):
            if atom[2] == ' CA ':
                CA_indices.append(i)
        return self.coords[CA_indices]

    def rmsd_static(self, other):
        #returns the CA rmsd between two ParsedPDB objects
        #assumes that the two objects have the same number of atoms in the same order
        return np.sqrt(np.mean((self.get_CA_coords() - other.get_CA_coords())**2))

    def rmsd_kabsch(self,other):
        #uses kabsch algorithm to superimpose other onto self and returns the rmsd
        #assumes that the two objects have the same number of atoms in the same order
        #the coordinates will NOT be updated
        self_ca = self.get_CA_coords()
        other_ca = other.get_CA_coords()
        # print()
        # print()
        # print(self_ca.shape, other_ca.shape)
        # print()
        # print()
        rms, rB, R = np_kabsch(self_ca, other_ca)
        return rms
    
    def get_seq(self):
        #returns a fasta string of the sequence
        seq = ""
        prev_chain = self.atom_details[0][5]
        for atom in self.atom_details:
            if atom[5] != prev_chain:
                seq += "/"
                prev_chain = atom[5]
            if atom[2] == ' CA ':
                seq += aa_3to1[atom[4]]
        return seq
    
    def total_length(self):
        #returns the total length of the sequence
        return len(self.get_seq().replace("/",""))

def get_chain_permutations(chains: list) -> list:
    """
    Gets all permutations of the chains.
    """
    import itertools

    return list(itertools.permutations(chains))





import subprocess
import os

mmalign_exe = f"{SCRIPTDIR}/mmalign/MMalign"

#ensure it exists
if not os.path.exists(mmalign_exe):
    print(f"ERROR: {mmalign_exe} does not exist. Please run 'g++ -static -O3 -ffast-math -lm -o MMalign MMalign.cpp' in the mmalign directory")
    raise FileNotFoundError

from typing import Tuple
import tempfile


def MMalign(
    model:ParsedPDB, reference:ParsedPDB
) -> Tuple[float, float, Dict[str,str]]:
    """
    Aligns two models using MMalign. Works with single and multiple chains.
    Returns the RMSD,TMscore, and the aligned model.
    """

    import random
    import string

    #decend into a temporary directory
    with tempfile.TemporaryDirectory() as tempdir:
        #save the two models to pdb files
        model_path = f"{tempdir}/model.pdb"
        reference_path = f"{tempdir}/reference.pdb"
        model.make_pdb_file(model_path)
        reference.make_pdb_file(reference_path)

        aligned_model_path = f"{tempdir}/output.pdb"

        #run mmalign
        mmalign_proc = subprocess.Popen([mmalign_exe, model_path, reference_path,'-outfmt',"2","-o",aligned_model_path], stdout=subprocess.PIPE)
        mmalign_output = mmalign_proc.communicate()[0].decode("utf-8")

        #parse the output
        data_line = mmalign_output.split("\n")[1]
        rmsd = float(data_line.split()[4])
        tmscore = float(data_line.split()[3])
        model_chain_order_1 = data_line.split()[0].split(":")[1:]
        reference_chain_order_2 = data_line.split()[1].split(":")[1:]
        
        #might be backwards, who knows
        chain_order_mapping = dict(zip(model_chain_order_1,reference_chain_order_2))

        #update model
        model.parse_pdb_file(aligned_model_path)
        return rmsd, tmscore, chain_order_mapping

def compute_per_residue_lddt(query_path:str, reference_path:str):
    raise NotImplementedError


# parsed_pdb_test_1 = ParsedPDB()
# parsed_pdb_test_1.parse_pdb_file(test_pdb_1)
# seq1 = parsed_pdb_test_1.get_seq()
# assert seq1 == "A" * 50

# parsed_pdb_test_2 = ParsedPDB()
# parsed_pdb_test_2.parse_pdb_file(test_pdb_2)

# parsed_pdb_test_3 = ParsedPDB()
# parsed_pdb_test_3.parse_pdb_file(test_pdb_3)

# rmsd, tmscore, chain_order_mapping = MMalign(parsed_pdb_test_1, parsed_pdb_test_2)
# print(rmsd, tmscore, chain_order_mapping)

# rmsd, tmscore, chain_order_mapping = MMalign(parsed_pdb_test_1, parsed_pdb_test_3)
# print(rmsd, tmscore, chain_order_mapping)

# rmsd, tmscore, chain_order_mapping = MMalign(parsed_pdb_test_2, parsed_pdb_test_3)
# print(rmsd, tmscore, chain_order_mapping)


# rmsd, tmscore, chain_order_mapping = MMalign(parsed_pdb_test_3, parsed_pdb_test_3)
# print(rmsd, tmscore, chain_order_mapping)
# print(parsed_pdb_test_3.rmsd_static(parsed_pdb_test_3))
# print(parsed_pdb_test_3.rmsd_kabsch(parsed_pdb_test_3))





def convert_pdb_chainbreak_to_new_chain(pdbstring):
    previous_resid = 0
    chain_num = 0
    new_pdbstring = ""
    import string

    alphabet = string.ascii_uppercase + string.digits + string.ascii_lowercase
    for line in pdbstring.split("\n"):
        if line[:4] == "ATOM":
            resid = int(line[22:26])
            if resid - previous_resid > 1:
                chain_num += 1
                if chain_num >= len(alphabet):
                    raise Exception(
                        "Too many chains to convert to new chain format. "
                        "Decrease the number of chains or increase the alphabet size."
                    )
            new_pdbstring += line[:21] + f"{alphabet[chain_num]: >1}" + line[22:] + "\n"
            previous_resid = resid
        else:
            new_pdbstring += line + "\n"
    return new_pdbstring


@dataclass()
class PredictionTarget:
    name: str
    seq: str
    parsed_pdb: str = None
    input_path: str = None

    def __lt__(self, other):
        return len(self) < len(other)

    def __len__(self):
        return len(self.seq.replace("/", ""))


def parse_fasta(path):
    if path.endswith(".gz"):
        import gzip

        filehandle = gzip.open(path, "rt")
    else:
        filehandle = open(path, "rt")

    outputs = []

    seq = ""
    name = ""
    for line in filehandle:
        if line.startswith(">"):
            if len(seq) > 0:
                outputs.append(PredictionTarget(name, seq))
            name = line[1:].strip()
            seq = ""
        else:
            seq += line.strip()
    seq = seq.replace(":","/")
    if len(seq) > 0:
        # This should always be true for a well formatted fasta file
        outputs.append(PredictionTarget(name, seq))

    filehandle.close()

    return outputs


unique_name_counter = 0


# def get_unique_name():
#     global unique_name_counter
#     unique_name_counter += 1
#     return f"struct{unique_name_counter}"


def parse_pdb(path):
    name = path.split("/")[-1].split(".pdb")[0]
    parsed_pdb = ParsedPDB()
    parsed_pdb.parse_pdb_file(path)
    seq = parsed_pdb.get_seq()

    return [PredictionTarget(name, seq, parsed_pdb, path)]


def parse_silent(path):
    outputs = []
    index = silent_tools.get_silent_index(path)

    tags = index["tags"]

    structures = silent_tools.get_silent_structures(path, index, tags)

    for name, structure in zip(tags, structures):

        chain_per_res = silent_tools.get_chain_ids(structure)

        # only gonna grab C-alphas
        seq = "".join(silent_tools.get_sequence_chunks(structure))
        # atoms = silent_tools.sketch_get_cas_protein_struct(structure)
        atoms = silent_tools.sketch_get_atoms(structure, 1)
        pdbstring = silent_tools.write_pdb_atoms(
            atoms, seq, ["CA"], chain_ids=chain_per_res
        )
        parsed_pdb = ParsedPDB()
        parsed_pdb.parse_pdbstr("\n".join(pdbstring))
        parsed_pdb.name = name
        seq = parsed_pdb.get_seq() #new seq is properly formatted

        outputs.append(PredictionTarget(name, seq, parsed_pdb, path))

    return outputs


def parse_file(path):
    targets = []
    if path.endswith(".gz"):
        filename = path[:-3]
    else:
        filename = path

    if filename.endswith(".silent"):
        targets.extend(parse_silent(path))
    elif filename.endswith(".fa") or filename.endswith(".fasta"):
        targets.extend(parse_fasta(path))
    elif filename.endswith(".pdb"):
        targets.extend(parse_pdb(path))

    return targets


query_targets = []
for file in args.input_files:
    query_targets.extend(parse_file(file))


from alphafold.model import model
from alphafold.model import config
from alphafold.model import data
from alphafold.common import protein
from alphafold.data import parsers

# I don't know if this is a good idea.
if args.version == "multimer" or args.version == "multimer_v2":
    from alphafold.data import pipeline_multimer
from alphafold.data import pipeline

import colabfold as cf
from collections import defaultdict
import tqdm
import jax
from jax.lib import xla_bridge

device = xla_bridge.get_backend().platform
print("using ", device)


# if args.amber_relax:
#     from alphafold.relax import relax

#     RELAX_MAX_ITERATIONS = 0
#     RELAX_ENERGY_TOLERANCE = 2.39
#     RELAX_STIFFNESS = 10.0
#     RELAX_EXCLUDE_RESIDUES = []
#     RELAX_MAX_OUTER_ITERATIONS = 3

#     amber_relaxer = relax.AmberRelaxation(
#         max_iterations=RELAX_MAX_ITERATIONS,
#         tolerance=RELAX_ENERGY_TOLERANCE,
#         stiffness=RELAX_STIFFNESS,
#         exclude_residues=RELAX_EXCLUDE_RESIDUES,
#         max_outer_iterations=RELAX_MAX_OUTER_ITERATIONS,
#     )


longest = max([len(tgt) for tgt in query_targets])

if longest < 400 and device != "cpu":
    # catch the user's eye
    plural = "s are" if len(query_targets) > 1 else " is"
    print(
        "=======================================================================================\n"
        + f"WARNING: Your query{plural} shorter than 400 residues. This is a very small protein.\n"
        + "You may want to use the CPU to conserve GPU resources for those who need them.\n"
        + "Remember that you can launch far more jobs in parallel on CPUs than you can on GPUs...\n"
        + "See this example of how prediction time scales on CPU vs GPU: \n"
        + "https://docs.google.com/spreadsheets/d/1jTGITpIx6fJehAplUkXtePOp7me3Dpq_pPKHn68F7XY\n"
        + "======================================================================================="
    )

seed_range = list(range(args.seed_start, args.seed_start + args.nstruct))

# # initial guess and multimer are not compatible
# if args.initial_guess and args.version == "multimer":
#     print("WARNING: initial guess and multimer are not compatible. ")
#     exit(1)

# # TODO initial guess needs a pdb file if and only if args.input_file is a fasta file
# if type(args.initial_guess) == str:  # check input_file type
#     if args.input_files[0].endswith(".pdb"):
#         print("WARNING: initial guess was provided a PDB and input_file was a PDB")
#         print(
#             "No followup argument is needed for initial guess when input_file is a PDB"
#         )
#         exit(1)
#     elif args.input_files[0].endswith(".silent"):
#         print("WARNING: initial guess was provided a PDB and input_file was a .silent")
#         print(
#             "No followup argument is needed for initial guess when input_file is a silent"
#         )
#         exit(1)
#     else:
#         pass
# else:
#     pass
# if args.input_files[0].endswith(".fa"):
#     if args.initial_guess is True:
#         print("WARNING: initial guess needs a PDB if input_file was a fasta")
#         exit(1)
#     else:
#         pass
# else:
#     pass


# blatently stolen from https://github.com/sokrypton/ColabFold/blob/8e6b6bb582f40a4fea06b19fc001d3d9ca208197/colabfold/alphafold/msa.py#L15
# by konstin i think
# no worries, I plan on going and actually forking colabfold eventually.
from alphafold.model.features import FeatureDict
from alphafold.model.tf import shape_placeholders
import tensorflow as tf
from typing import Mapping, Any

NUM_RES = shape_placeholders.NUM_RES
NUM_MSA_SEQ = shape_placeholders.NUM_MSA_SEQ
NUM_EXTRA_SEQ = shape_placeholders.NUM_EXTRA_SEQ
NUM_TEMPLATES = shape_placeholders.NUM_TEMPLATES


def make_fixed_size(feat, runner, max_length):
    """pad input features"""
    cfg = runner.config

    if cfg.model.global_config.multimer_mode:
        # shape_schema = ?
        # pad_size_map = {
        #     shape_placeholders.NUM_RES: max_length,
        #     shape_placeholders.NUM_MSA_SEQ: cfg.model.embeddings_and_evoformer.num_msa,
        #     shape_placeholders.NUM_EXTRA_SEQ: cfg.model.embeddings_and_evoformer.num_extra_msa,
        #     shape_placeholders.NUM_TEMPLATES: 0,
        # }
        print("Warning: padding sequences in multimer mode is not implemented yet")
        return feat
    else:
        shape_schema = {k: [None] + v for k, v in dict(cfg.data.eval.feat).items()}
        pad_size_map = {
            shape_placeholders.NUM_RES: max_length,
            shape_placeholders.NUM_MSA_SEQ: cfg.data.eval.max_msa_clusters,
            shape_placeholders.NUM_EXTRA_SEQ: cfg.data.common.max_extra_msa,
            shape_placeholders.NUM_TEMPLATES: 0,
        }


    for k, v in feat.items():
        # Don't transfer this to the accelerator.
        if k == "extra_cluster_assignment":
            continue
        shape = list(v.shape)

        schema = shape_schema[k]
        assert len(shape) == len(schema), (
            f"Rank mismatch between shape and shape schema for {k}: "
            f"{shape} vs {schema}"
        )
        pad_size = [pad_size_map.get(s2, None) or s1 for (s1, s2) in zip(shape, schema)]
        padding = [(0, p - tf.shape(v)[i]) for i, p in enumerate(pad_size)]
        if padding:
            feat[k] = tf.pad(v, padding, name=f"pad_to_fixed_{k}")
            feat[k].set_shape(pad_size)
    return {k: np.asarray(v) for k, v in feat.items()}


#######################################################################################################################
# Adapted from code by Nate Bennett for providing initial guess for the alphafold model
import jax.numpy as jnp
from alphafold.common import residue_constants
from alphafold.data import templates
import collections


def af2_get_atom_positions(parsed_pdb) -> Tuple[np.ndarray, np.ndarray]:
    """Gets atom positions and mask."""

    lines = parsed_pdb.get_pdbstr().splitlines()

    # indices of residues observed in the structure
    idx_s = [
        int(l[22:26]) for l in lines if l[:4] == "ATOM" and l[12:16].strip() == "CA"
    ]
    num_res = len(idx_s)

    all_positions = np.zeros([num_res, residue_constants.atom_type_num, 3])
    all_positions_mask = np.zeros(
        [num_res, residue_constants.atom_type_num], dtype=np.int64
    )

    residues = collections.defaultdict(list)
    # 4 BB + up to 10 SC atoms
    xyz = np.full((len(idx_s), 14, 3), np.nan, dtype=np.float32)
    for l in lines:
        if l[:4] != "ATOM":
            continue
        resNo, atom, aa = int(l[22:26]), l[12:16], l[17:20]

        residues[resNo].append(
            (atom.strip(), aa, [float(l[30:38]), float(l[38:46]), float(l[46:54])])
        )

    for resNo in residues:
        pos = np.zeros([residue_constants.atom_type_num, 3], dtype=np.float32)
        mask = np.zeros([residue_constants.atom_type_num], dtype=np.float32)

        for atom in residues[resNo]:
            atom_name = atom[0]
            x, y, z = atom[2]
            if atom_name in residue_constants.atom_order.keys():
                pos[residue_constants.atom_order[atom_name]] = [x, y, z]
                mask[residue_constants.atom_order[atom_name]] = 1.0
            elif atom_name.upper() == "SE" and res.get_resname() == "MSE":
                # Put the coordinates of the selenium atom in the sulphur column.
                pos[residue_constants.atom_order["SD"]] = [x, y, z]
                mask[residue_constants.atom_order["SD"]] = 1.0

        idx = idx_s.index(resNo)  # This is the order they show up in the pdb
        all_positions[idx] = pos
        all_positions_mask[idx] = mask
    # _check_residue_distances(
    #     all_positions, all_positions_mask, max_ca_ca_distance) # AF2 checks this but if we want to allow massive truncations we don't want to check this

    return all_positions, all_positions_mask


def af2_all_atom(parsed_pdb,pad_to=None):
    template_seq = parsed_pdb.get_seq().replace("/", "")

    all_atom_positions, all_atom_mask = af2_get_atom_positions(parsed_pdb)

    all_atom_positions = np.split(all_atom_positions, all_atom_positions.shape[0])

    templates_all_atom_positions = []

    # Initially fill will all zero values
    pad_length = pad_to if pad_to is not None else len(template_seq)

    assert pad_length >= len(template_seq)

    for _ in range(pad_length):
        templates_all_atom_positions.append(
            jnp.zeros((residue_constants.atom_type_num, 3))
        )

    for idx, i in enumerate(template_seq):
        templates_all_atom_positions[idx] = all_atom_positions[idx][
            0
        ]  # assign target indices to template coordinates

    return jnp.array(templates_all_atom_positions)


def mk_mock_template(query_sequence):
    # mock template features
    output_templates_sequence = []
    output_confidence_scores = []
    templates_all_atom_positions = []
    templates_all_atom_masks = []

    for _ in query_sequence:
        templates_all_atom_positions.append(
            np.zeros((templates.residue_constants.atom_type_num, 3))
        )
        templates_all_atom_masks.append(
            np.zeros(templates.residue_constants.atom_type_num)
        )
        output_templates_sequence.append("-")
        output_confidence_scores.append(-1)
    output_templates_sequence = "".join(output_templates_sequence)
    templates_aatype = templates.residue_constants.sequence_to_onehot(
        output_templates_sequence, templates.residue_constants.HHBLITS_AA_TO_ID
    )

    template_features = {
        "template_all_atom_positions": np.array(templates_all_atom_positions)[None],
        "template_all_atom_masks": np.array(templates_all_atom_masks)[None],
        "template_sequence": [f"none".encode()],
        "template_aatype": np.array(templates_aatype)[None],
        "template_confidence_scores": np.array(output_confidence_scores)[None],
        "template_domain_names": [f"none".encode()],
        "template_release_date": [f"none".encode()],
    }

    return template_features


#######################################################################################################################


#### load up reference pdb, if present
if args.reference_pdb is not None:
    parsed_pdb_reference = ParsedPDB()
    parsed_pdb_reference.parse_pdb_file(args.reference_pdb)

######## initial guess logic ########

if type(args.initial_guess) == bool and args.initial_guess:
    # if initial_guess is requested but no pdb is provided, we will use each 
    #input structure as the initial guess. This means that the inputs
    #must be pdb files.
    for input_path in args.input_files:
        if not input_path.endswith(".pdb"):
            raise ValueError(
                f"initial_guess is True but no pdb was provided. {input_path} is not a pdb file. What are we supposed to initialize with?"
            )

#### load up initial guess pdb, if present
if type(args.initial_guess) == str:
    parsed_pdb_initial_guess = ParsedPDB()
    parsed_pdb_initial_guess.parse_pdb_file(args.initial_guess)

    # check that all the provided fasta sequences are the same length
    initial_guess_length = parsed_pdb_initial_guess.total_length()
    for tgt in query_targets:
        if len(tgt) != initial_guess_length:
            raise ValueError(
                f"All provided fasta sequences must be the same length as the initial guess structure in order to use initial guess with fasta inputs. Initial guess length: {initial_guess_length}, {tgt.name} length: {len(tgt)}"
            )


######################################

max_length = max([len(tgt) for tgt in query_targets])


if args.type == "multimer" or args.type == "multimer_v2":
    model_name = "model_5_multimer"
else:
    model_name = "model_5_ptm" if args.type == "monomer_ptm" else "model_5"

if "all" in args.models:
    model_names = ["model_1", "model_2", "model_3", "model_4", "model_5"]
else:
    model_names = [f"model_{model_num}" for model_num in args.models]


cfg = config.model_config(model_name)
params = data.get_model_haiku_params(model_name, data_dir=ALPHAFOLD_DATADIR)

if args.version == "multimer" or args.version == "multimer_v2":
    cfg.model.num_ensemble_eval = args.num_ensemble
    # cfg.model.embeddings_and_evoformer.num_extra_msa = args.mock_msa_depth
    cfg.model.embeddings_and_evoformer.masked_msa.replace_fraction = args.pct_seq_mask

    # templates are enabled by default, but I'm not supplying them, so disable
    cfg.model.embeddings_and_evoformer.template.enabled = False

    cfg.model.embeddings_and_evoformer.initial_guess = bool(
        args.initial_guess
    )  # new for initial guessing

else:
    cfg.data.eval.num_ensemble = args.num_ensemble
    cfg.data.eval.max_msa_clusters = args.mock_msa_depth
    cfg.data.common.max_extra_msa = args.mock_msa_depth
    cfg.data.eval.masked_msa_replace_fraction = args.pct_seq_mask
    cfg.data.common.num_recycle = args.max_recycles
    cfg.model.embeddings_and_evoformer.initial_guess = bool(
        args.initial_guess
    )  # new for initial guessing
    # do I also need to turn on template?

cfg.model.recycle_tol = args.recycle_tol
cfg.model.num_recycle = args.max_recycles


model_runner = model.RunModel(
    cfg,
    params,
    is_training=args.enable_dropout,
    return_representations=args.save_intermediates,
)


output_counter = 0
with tqdm.tqdm(total=len(query_targets)) as pbar1:
    for target in query_targets:
        pbar1.set_description(f"Input: {target.name}")

        # this is a lazy hack to avoid replacing var names.
        full_sequence = target.seq.replace("/", "")
        query_sequences = target.seq.split("/")
        name = target.name

        #############################
        # define input features
        #############################

        if not args.initial_guess:  # initial guess is False by default
            initial_guess = None
        elif type(args.initial_guess) == str:  # use the provided pdb
            initial_guess = af2_all_atom(parsed_pdb_initial_guess,pad_to=max_length)
        else:  # use the target structure
            print("Using target structure as initial guess")
            initial_guess = af2_all_atom(target.parsed_pdb,pad_to=max_length)

        num_res = len(full_sequence)
        feature_dict = {}
        msas = [parsers.Msa([full_sequence], [[0] * len(full_sequence)], [name])]

        if args.version == "multimer" or args.version == "multimer_v2":

            feature_dict = pipeline_multimer.DataPipelineFaker().process(
                query_sequences
            )
        else:
            feature_dict.update(
                pipeline.make_sequence_features(full_sequence, name, num_res)
            )
            feature_dict.update(pipeline.make_msa_features(msas))

            Ls = [len(chain_seq) for chain_seq in query_sequences]
            Ls_plot = sum([[len(seq)] for seq in query_sequences], [])
            # this introduces a bug where the plot just doesn't work for multimer version

            feature_dict["residue_index"] = cf.chain_break(
                feature_dict["residue_index"], Ls
            )

            if args.initial_guess:
                feature_dict.update(mk_mock_template(query_sequences))

        ###########################
        # run alphafold
        ###########################
        def parse_results(prediction_result, processed_feature_dict):

            # figure note... it would be nice to use prediction_result["structure_module"]["final_atom_mask"] to mask out everything in prediction_result that shouldn't be there due to padding.
            b_factors = (
                prediction_result["plddt"][:, None]
                * prediction_result["structure_module"][
                    "final_atom_mask"
                ]  # I think not needed b/c I truncated the vector earlier
            )            

            # but for now let's focus on truncating the results we most care about to the length of the target sequence
            prediction_result["plddt"] = prediction_result["plddt"][: len(target.seq)]
            if "predicted_aligned_error" in prediction_result:
                prediction_result["predicted_aligned_error"] = prediction_result[
                    "predicted_aligned_error"
                ][: len(target.seq), : len(target.seq)]

            dist_bins = jax.numpy.append(0, prediction_result["distogram"]["bin_edges"])
            dist_mtx = dist_bins[prediction_result["distogram"]["logits"].argmax(-1)]
            contact_mtx = jax.nn.softmax(prediction_result["distogram"]["logits"])[
                :, :, dist_bins < 8
            ].sum(-1)

            out = {
                "unrelaxed_protein": protein.from_prediction(
                    processed_feature_dict,
                    prediction_result,
                    b_factors=b_factors,
                    remove_leading_feature_dimension=(args.type != "multimer" and args.type != "multimer_v2"),
                ),
                "plddt": prediction_result["plddt"],
                "mean_plddt": prediction_result["plddt"].mean(),
                "dists": dist_mtx,
                "adj": contact_mtx,
            }

            if "ptm" in prediction_result:
                out.update(
                    {
                        "pae": prediction_result["predicted_aligned_error"],
                        "pTMscore": prediction_result["ptm"],
                    }
                )
            if args.type == "multimer" or args.type == "multimer_v2":
                out.update(
                    {
                        "pTMscore": prediction_result["ptm"],
                        "pae": prediction_result["predicted_aligned_error"],
                        "iptm": prediction_result["iptm"],
                    }
                )
            return out

        total = len(model_names) * len(seed_range)

        with tqdm.tqdm(total=total) as pbar2:
            outs = {}

            def report(key):
                pbar2.update(n=1)
                o = outs[key]
                out_dict = {}
                out_dict["mean_plddt"] = o["mean_plddt"]

                out_dict["recycles"] = o["recycles"]
                out_dict["tol"] = o["tol"]
                out_dict["model"] = key.split("_")[1]
                out_dict["type"] = args.type
                out_dict["seed"] = key.split("_")[-1]

                output_line = f"{name} {key} recycles:{o['recycles']} tol:{o['tol']:.2f} mean_plddt:{o['mean_plddt']:.2f}"
                if args.type == "monomer_ptm" or args.type == "multimer" or args.type == "multimer_v2":
                    output_line += f" pTMscore:{o['pTMscore']:.2f}"

                prefix = f"{name}_{key}"
                fout_name = os.path.join(args.out_dir, f"{prefix}_unrelaxed.pdb")

                output_pdbstr = protein.to_pdb(o["unrelaxed_protein"])

                output_pdbstr = convert_pdb_chainbreak_to_new_chain(output_pdbstr)
                parsed_pdb_output = ParsedPDB()
                parsed_pdb_output.parse_pdbstr(output_pdbstr)
                parsed_pdb_output.renumber()
                # output_pdbstr = renumber(output_pdbstr)

                import string

                alphabet = string.ascii_uppercase + string.digits + string.ascii_lowercase
                #chain_range_map = get_chain_range_map(output_pdbstr)
                chain_range_map = parsed_pdb_output.get_chain_range_map()

                num_chains = len(chain_range_map)

                final_chain_order = list(
                    alphabet[:num_chains]
                )  # initialize with original order, basically, for the default case where there is no refernce or input pdb file

                final_chain_order_mapping = {
                    old_chain: new_chain
                    for old_chain, new_chain in zip(alphabet, final_chain_order)
                }

                #output_pymol_name = "temp_target"
                #pymol.cmd.read_pdbstr(output_pdbstr, oname=output_pymol_name)

                if args.reference_pdb is not None:
                    # if args.simple_rmsd:
                    #     rmsd = pymol_align(
                    #         "temp_target",
                    #         reference_pdb_name, "super"
                    #     )
                    # else:
                    #     rmsd, output_pdbstr, final_chain_order = pymol_multichain_align(
                    #         "temp_target", reference_pdb_name, "super"
                    #     )  # use super here b/c sequence is not guaranteed to be very similar

                    rmsd, tmscore, final_chain_order_mapping = MMalign(parsed_pdb_output, parsed_pdb_reference)
                    parsed_pdb_output.remap_chains(final_chain_order_mapping)
                    kabsch_rmsd = parsed_pdb_output.rmsd_kabsch(parsed_pdb_reference)

                    #out_dict["mmalign_rmsd_to_reference"] = rmsd
                    out_dict["rmsd_to_reference"] = kabsch_rmsd
                    out_dict["tmscore_to_reference"] = tmscore

                    #send this back up for the info-recorder
                    if target.parsed_pdb is None:
                        #outs[key]["mmalign_rmsd_to_input"] = rmsd 
                        outs[key]["rmsd_to_input"] = kabsch_rmsd
                        outs[key]["tmscore_to_input"] = tmscore

                    #pymol.cmd.delete("temp_target")
                    output_line += f" rmsd_to_reference:{rmsd:0.2f}"
                
                # parsed_pdb_output.make_pdb_file("TEST_RAW_OUTPUT.pdb")
                # target.parsed_pdb.make_pdb_file("TEST_RAW_TARGET.pdb")
                if target.parsed_pdb is not None:
                    
                    rmsd, tmscore, final_chain_order_mapping = MMalign(parsed_pdb_output, target.parsed_pdb)
                    parsed_pdb_output.remap_chains(final_chain_order_mapping)
                    # parsed_pdb_output.make_pdb_file("TEST_REMAPPED_OUTPUT.pdb")
                    kabsch_rmsd = parsed_pdb_output.rmsd_kabsch(target.parsed_pdb)

                    #out_dict["mmalign_rmsd_to_input"] = rmsd
                    out_dict["rmsd_to_input"] = kabsch_rmsd
                    out_dict["tmscore_to_input"] = tmscore

                    #send this back up for the info-recorder
                    #outs[key]["mmalign_rmsd_to_input"] = rmsd
                    outs[key]["rmsd_to_input"] = kabsch_rmsd
                    outs[key]["tmscore_to_input"] = tmscore

                    #pymol.cmd.delete("temp_target")
                    output_line += f" rmsd_to_input:{rmsd:0.2f}"
                
                #re-extract the per-residue lddt values from the b-factor column of the pdb
                plddt_list = parsed_pdb_output.get_bfactors()
                outs[key]["plddts"] = plddt_list

                # with open(fout_name, "w") as f:
                #     f.write(output_pdbstr)
                # pymol.cmd.save(fout_name, output_pymol_name)
                # pymol.cmd.delete(output_pymol_name)
                parsed_pdb_output.make_pdb_file(fout_name)

                # TO BE IMPLEMENTED
                # #compute real per-residue lddts
                # if target.pymol_obj_name is not None:
                #     query_path = fout_name
                #     with tempfile.NamedTemporaryFile() as reference_file:
                #         pymol.cmd.save(reference_file.name, target.pymol_obj_name)
                #         reference_path = reference_file.name
                #         lddt_list = compute_per_residue_lddt(query_path,reference_path)
                #         out_dict["lddts"] = lddt_list
                #         out_dict["mean_lddt"] = np.mean(lddt_list)
                #         outs[key]['lddts'] = lddt_list #send this back up for the info-recorder. this is getting messy
                


                if args.type == "monomer_ptm":
                    # calculate mean PAE for interactions between each chain pair, taking into account the changed chain order
                    pae = o["pae"]

                    # first, truncate the matrix to the full length of the sequence (without chainbreak characters "/"). It can sometimes be too long because of padding inputs
                    sequence_length = len(target.seq.replace("/", "").replace("U", ""))
                    pae = pae[:sequence_length, :sequence_length]

                    if args.output_pae:
                        out_dict["pae"] = pae

                    interaction_paes = []
                    for chain_1, chain_2 in itertools.permutations(
                        final_chain_order, 2
                    ):
                        chain_1_range_start, chain_1_range_stop = chain_range_map[
                            chain_1
                        ]
                        chain_2_range_start, chain_2_range_stop = chain_range_map[
                            chain_2
                        ]

                        final_chain_1 = final_chain_order_mapping[chain_1]
                        final_chain_2 = final_chain_order_mapping[chain_2]
                        interaction_pae = np.mean(
                            pae[
                                chain_1_range_start:chain_1_range_stop,
                                chain_2_range_start:chain_2_range_stop,
                            ]
                        )
                        interaction_paes.append(interaction_pae)
                        out_dict[
                            f"mean_pae_interaction_{final_chain_1}{final_chain_2}"
                        ] = interaction_pae

                    # average all the interaction PAEs
                    out_dict["mean_pae_interaction"] = np.mean(interaction_paes)

                    # calculate mean intra-chain PAE per chain
                    intra_chain_paes = []
                    for chain in alphabet[:num_chains]:
                        chain_range_start, chain_range_stop = chain_range_map[chain]
                        intra_chain_pae = np.mean(
                            pae[
                                chain_range_start:chain_range_stop,
                                chain_range_start:chain_range_stop,
                            ]
                        )
                        intra_chain_paes.append(intra_chain_pae)
                        out_dict[f"mean_pae_intra_chain_{chain}"] = intra_chain_pae

                    # average all the intrachain PAEs
                    out_dict["mean_pae_intra_chain"] = np.mean(intra_chain_paes)

                    # average all the PAEs
                    out_dict["mean_pae"] = np.mean(pae)

                    out_dict["pTMscore"] = o["pTMscore"]
                elif args.type == "multimer" or args.type == "multimer_v2":
                    raise NotImplementedError
                    out_dict["ptm"] = o["pTMscore"]
                    out_dict["iptm"] = o["iptm"]

                # if args.show_images:
                #     fig = cf.plot_protein(o["unrelaxed_protein"], Ls=Ls_plot, dpi=200)
                #     plt.savefig(
                #         os.path.join(args.out_dir, f"{prefix}.png"),
                #         bbox_inches="tight",
                #     )
                #     plt.close(fig)

                # if args.amber_relax:
                #     # Relax the prediction.
                #     relaxed_pdb_str, _, _ = amber_relaxer.process(
                #         prot=o["unrelaxed_protein"]
                #     )

                #     # Save the relaxed PDB.
                #     relaxed_output_path = os.path.join(
                #         args.out_dir, f"{prefix}_relaxed.pdb"
                #     )
                #     with open(relaxed_output_path, "w") as f:
                #         f.write(relaxed_pdb_str)

                # np.savez_compressed(os.path.join(args.out_dir,f'{prefix}_prediction_results.npz'),**out_dict)

                # cast devicearray to serializable type
                for key in out_dict:
                    out_dict[key] = np.array(out_dict[key]).tolist()

                import json

                # output as nicely formatted json
                global time_checkpoint
                elapsed_time = time.time() - time_checkpoint
                output_line += f" elapsed time (s): {elapsed_time}"

                if args.output_summary:
                    with open(os.path.join(args.out_dir, "reports.txt"), "a") as f:
                        f.write(output_line + "\n")
                print(output_line)

                out_dict["elapsed_time"] = elapsed_time

                with open(
                    os.path.join(args.out_dir, f"{prefix}_prediction_results.json"),
                    "w",
                ) as f:
                    json.dump(out_dict, f, indent=2)

                time_checkpoint = time.time()

            #######################################################################

            # go through each random_seed
            for seed in seed_range:

                # prep input features
                processed_feature_dict = model_runner.process_features(
                    feature_dict, random_seed=seed
                )

                # pad input features
                # Pad sequences to the same length
                ##sequence padding
                # I'm not sure if this is compatible with multimer version or not, but I'll stick it here for now
                # model_config = model_runner.config
                # eval_cfg = model_config.data.eval
                # crop_feats = {k: [None] + v for k, v in dict(eval_cfg.feat).items()}
                # print(crop_feats)
                # feature_dict = make_fixed_size(
                #     feature_dict,
                #     crop_feats,
                #     args.mock_msa_depth,
                #     args.mock_msa_depth,
                #     max_length,
                # )
                processed_feature_dict = make_fixed_size(
                    processed_feature_dict, model_runner, max_length
                )

                # go through each model
                for num, model_name in enumerate(model_names):


                    info_recorder = InfoCollector(args.info_collector_path, args.info_collector_config)
                    info_recorder['sequence'] = target.seq
                    info_recorder['padded-length'] = max_length
                    info_recorder['seed'] = seed
                    info_recorder['source'] = target.input_path
                    if target.parsed_pdb is not None:
                        #get the pdb string of the input structure as pymol has saved it
                        info_recorder['input-pdb'] = target.parsed_pdb.name
                    info_recorder['model-num'] =  model_name
                    info_recorder['used-msa'] = False
                    info_recorder['used-initial-guess'] = bool(args.initial_guess)
                    info_recorder['used-templates'] = False
                    info_recorder['output-number'] = output_counter
                    info_recorder['af2-version'] = args.type

                    model_mod = ""
                    if args.type == "monomer_ptm":
                        model_mod = "_ptm"
                    elif args.type == "multimer":
                        model_mod = "_multimer"
                    elif args.type == "multimer_v2":
                        model_mod = "_multimer_v2"
                    model_name = model_name + model_mod
                    key = f"{model_name}_seed_{seed}"
                    pbar2.set_description(f"Running {key}")

                    # check if this prediction/seed has already been done
                    prefix = f"{name}_{key}"
                    if not args.overwrite and os.path.exists(
                        os.path.join(args.out_dir, f"{prefix}_prediction_results.json")
                    ):
                        print(f"{prefix}_prediction_results.json already exists")
                        continue

                    # replace model parameters
                    params = data.get_model_haiku_params(
                        model_name, data_dir=ALPHAFOLD_DATADIR
                    )
                    for k in model_runner.params.keys():
                        model_runner.params[k] = params[k]

                    # predict
                    if args.initial_guess:
                        prediction_result, (r, t) = cf.to(
                            model_runner.predict(
                                processed_feature_dict,
                                random_seed=seed,
                                initial_guess=initial_guess,
                            ),
                            device,
                        )  # is this ok?
                    else:
                        # a quick hack because the multimer version of the model_runner doesn't have initial_guess in its signature (is that the term?).
                        # the fix will be to update Multimer code to accept initial_guess deep down in the actual code
                        prediction_result, (r, t) = cf.to(
                            model_runner.predict(
                                processed_feature_dict, random_seed=seed
                            ),
                            device,
                        )  # is this ok?

                    # save results
                    outs[key] = parse_results(prediction_result, processed_feature_dict)
                    outs[key].update({"recycles": r, "tol": t})
                    report(key)
                    output_counter += 1

                    info_recorder['pLDDT'] = outs[key]['plddts']

                    #pLDDT is currently a list of all atoms, but we want to report the pLDDT of the Ca atoms only


                    #info_recorder['LDDT'] = outs[key]['lddts'] #to be implemented
                    info_recorder['num-recycles'] = outs[key]['recycles']
                    info_recorder['pae-matrix'] = outs[key]['pae'].tolist()
                    # if 'rmsd_to_input' in outs[key].keys():
                    #     info_recorder['rmsd'] = outs[key]['rmsd_to_input']
                    #     info_recorder['TMscore'] = outs[key]['tmscore_to_input']
                    #if "mmalign_rmsd_to_input" in outs[key].keys():
                    if "rmsd_to_input" in outs[key].keys():
                        #info_recorder['mmalign_rmsd'] = outs[key]['mmalign_rmsd_to_input']
                        info_recorder['rmsd'] = outs[key]['rmsd_to_input']
                        info_recorder['TMscore'] = outs[key]['tmscore_to_input']
                    

                    info_recorder['pTMscore'] = outs[key]['pTMscore']

                    info_recorder.report()

                    del prediction_result, params
                del processed_feature_dict

        pbar1.update(1)
