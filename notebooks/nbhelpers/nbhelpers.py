# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""
Helper functions for the AWS-Alphafold notebook.
"""
from datetime import datetime
import boto3
import uuid
import sagemaker
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
import string
from string import ascii_uppercase, ascii_lowercase
import py3Dmol
import json
import re

boto_session = boto3.session.Session()
sm_session = sagemaker.session.Session(boto_session)
region = boto_session.region_name
s3 = boto_session.client("s3", region_name=region)
batch = boto_session.client("batch", region_name=region)
cfn = boto_session.client("cloudformation", region_name=region)
logs_client = boto_session.client("logs")


def create_job_name(suffix=None):

    """
    Define a simple job identifier
    """

    if suffix == None:
        return datetime.now().strftime("%Y%m%dT%H%M%S")
    else:
        ## Ensure that the suffix conforms to the Batch requirements, (only letters,
        ## numbers, hyphens, and underscores are allowed).
        suffix = re.sub("\W", "_", suffix)
        return datetime.now().strftime("%Y%m%dT%H%M%S") + "_" + suffix


def upload_fasta_to_s3(
    sequences,
    ids,
    bucket=sm_session.default_bucket(),
    job_name=uuid.uuid4(),
    region="us-east-1",
):

    """
    Create a fasta file and upload it to S3.
    """

    file_out = "_tmp.fasta"
    with open(file_out, "a") as f_out:
        for i, seq in enumerate(sequences):
            seq_record = SeqRecord(Seq(seq), id=ids[i])
            SeqIO.write(seq_record, f_out, "fasta")

    object_key = f"{job_name}/{job_name}.fasta"
    response = s3.upload_file(file_out, bucket, object_key)
    os.remove(file_out)
    s3_uri = f"s3://{bucket}/{object_key}"
    print(f"Sequence file uploaded to {s3_uri}")
    return object_key


def list_alphafold_stacks():
    af_stacks = []
    for stack in cfn.list_stacks(
        StackStatusFilter=["CREATE_COMPLETE", "UPDATE_COMPLETE"]
    )["StackSummaries"]:
        if "alphafold-cfn-batch.yaml" in stack.get("TemplateDescription", []):
            af_stacks.append(stack)
    return af_stacks

def get_batch_resources(stack_name):
    """
    Get the resource names of the Batch resources for running Alphafold jobs.
    """

    # stack_name = af_stacks[0]["StackName"]
    stack_resources = cfn.list_stack_resources(StackName=stack_name)
    cpu_job_queue_spot = None
    for resource in stack_resources["StackResourceSummaries"]:
        if resource["LogicalResourceId"] == "GPUFoldingJobDefinition":
            gpu_job_definition = resource["PhysicalResourceId"]
        if resource["LogicalResourceId"] == "PrivateGPUJobQueue":
            gpu_job_queue = resource["PhysicalResourceId"]
        if resource["LogicalResourceId"] == "CPUFoldingJobDefinition":
            cpu_job_definition = resource["PhysicalResourceId"]
        if resource["LogicalResourceId"] == "PrivateCPUJobQueueOnDemand":
            cpu_job_queue_od = download_job_queue = resource["PhysicalResourceId"]        
        if resource["LogicalResourceId"] == "PrivateCPUJobQueueSpot":
            cpu_job_queue_spot = resource["PhysicalResourceId"]                    
        if resource["LogicalResourceId"] == "CPUDownloadJobDefinition":
            download_job_definition = resource["PhysicalResourceId"]
    return {
        "gpu_job_definition": gpu_job_definition,
        "gpu_job_queue": gpu_job_queue,
        "cpu_job_definition": cpu_job_definition,
        "cpu_job_queue_od": cpu_job_queue_od,
        "cpu_job_queue_spot": cpu_job_queue_spot,
        "download_job_definition": download_job_definition,
        "download_job_queue": download_job_queue,
    }


def get_batch_job_info(jobId):

    """
    Retrieve and format information about a batch job.
    """

    job_description = batch.describe_jobs(jobs=[jobId])

    output = {
        "jobArn": job_description["jobs"][0]["jobArn"],
        "jobName": job_description["jobs"][0]["jobName"],
        "jobId": job_description["jobs"][0]["jobId"],
        "status": job_description["jobs"][0]["status"],
        "createdAt": datetime.utcfromtimestamp(
            job_description["jobs"][0]["createdAt"] / 1000
        ).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "dependsOn": job_description["jobs"][0]["dependsOn"],
        "tags": job_description["jobs"][0]["tags"],
    }

    if output["status"] in ["STARTING", "RUNNING", "SUCCEEDED", "FAILED"]:
        output["logStreamName"] = job_description["jobs"][0]["container"][
            "logStreamName"
        ]
    return output


def get_batch_logs(logStreamName):

    """
    Retrieve and format logs for batch job.
    """

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


def download_dir(client, bucket, local="data", prefix=""):
    """Recursively download files from S3."""

    paginator = client.get_paginator("list_objects_v2")
    file_count = 0
    for result in paginator.paginate(Bucket=bucket, Delimiter="/", Prefix=prefix):
        if result.get("CommonPrefixes") is not None:
            for subdir in result.get("CommonPrefixes"):
                download_dir(client, bucket, local, subdir.get("Prefix"))
        for file in result.get("Contents", []):
            dest_pathname = os.path.join(local, file.get("Key"))
            if not os.path.exists(os.path.dirname(dest_pathname)):
                os.makedirs(os.path.dirname(dest_pathname))
            client.download_file(bucket, file.get("Key"), dest_pathname)
            file_count += 1
    print(f"{file_count} files downloaded from s3.")
    return local


def download_results(bucket, job_name, local="data"):
    """Download MSA information from S3"""
    return download_dir(s3, bucket, local, job_name)


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


def display_structure(
    pdb_path,
    color="lDDT",
    show_sidechains=False,
    show_mainchains=False,
    chains=1,
    vmin=50,
    vmax=90,
):
    """
    Display the predicted structure in a Jupyter notebook cell
    """
    if color not in ["chain", "lDDT", "rainbow"]:
        raise ValueError("Color must be 'LDDT' (default), 'chain', or 'rainbow'")

    plot_pdb(
        pdb_path,
        show_sidechains=show_sidechains,
        show_mainchains=show_mainchains,
        color=color,
        chains=chains,
        vmin=vmin,
        vmax=vmax,
    ).show()
    if color == "lDDT":
        plot_plddt_legend().show()


def submit_batch_alphafold_job(
    job_name,
    fasta_paths,
    s3_bucket,
    data_dir="/mnt/data_dir/fsx",
    output_dir="alphafold",
    bfd_database_path="/mnt/bfd_database_path/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt",
    mgnify_database_path="/mnt/mgnify_database_path/mgy_clusters_2018_12.fa",
    pdb70_database_path="/mnt/pdb70_database_path/pdb70",
    obsolete_pdbs_path="/mnt/obsolete_pdbs_path/obsolete.dat",
    template_mmcif_dir="/mnt/template_mmcif_dir/mmcif_files",
    pdb_seqres_database_path="/mnt/pdb_seqres_database_path/pdb_seqres.txt",
    small_bfd_database_path="/mnt/small_bfd_database_path/bfd-first_non_consensus_sequences.fasta",
    uniclust30_database_path="/mnt/uniclust30_database_path/uniclust30_2018_08/uniclust30_2018_08",
    uniprot_database_path="/mnt/uniprot_database_path/uniprot.fasta",
    uniref90_database_path="/mnt/uniref90_database_path/uniref90.fasta",
    max_template_date=datetime.now().strftime("%Y-%m-%d"),
    db_preset="reduced_dbs",
    model_preset="monomer",
    benchmark=False,
    use_precomputed_msas=False,
    features_paths=None,
    run_features_only=False,
    logtostderr=True,
    cpu=4,
    memory=16,
    gpu=1,
    depends_on=None,
    stack_name=None,
    use_spot_instances=False,
    run_relax=True,
    num_multimer_predictions_per_model=1
):

    if stack_name is None:
        stack_name = list_alphafold_stacks()[0]["StackName"]
    batch_resources = get_batch_resources(stack_name)

    container_overrides = {
        "command": [
            f"--fasta_paths={fasta_paths}",
            f"--uniref90_database_path={uniref90_database_path}",
            f"--mgnify_database_path={mgnify_database_path}",
            f"--data_dir={data_dir}",
            f"--template_mmcif_dir={template_mmcif_dir}",
            f"--obsolete_pdbs_path={obsolete_pdbs_path}",
            f"--output_dir={output_dir}",
            f"--max_template_date={max_template_date}",
            f"--db_preset={db_preset}",
            f"--model_preset={model_preset}",
            f"--s3_bucket={s3_bucket}",
            f"--run_relax={run_relax}",
        ],
        "resourceRequirements": [
            {"value": str(cpu), "type": "VCPU"},
            {"value": str(memory * 1000), "type": "MEMORY"},
        ],
    }

    if model_preset == "multimer":
        container_overrides["command"].append(
            f"--uniprot_database_path={uniprot_database_path}"
        )
        container_overrides["command"].append(
            f"--pdb_seqres_database_path={pdb_seqres_database_path}"
        )
        container_overrides["command"].append(
            f"--num_multimer_predictions_per_model={num_multimer_predictions_per_model}"
        )
        print("If multimer prediction failes due to Amber relaxation, re-run with run_relax=False")
    else:
        container_overrides["command"].append(
            f"--pdb70_database_path={pdb70_database_path}"
        )

    if db_preset == "reduced_dbs":
        container_overrides["command"].append(
            f"--small_bfd_database_path={small_bfd_database_path}"
        )
    else:
        container_overrides["command"].append(
            f"--uniclust30_database_path={uniclust30_database_path}"
        )
        container_overrides["command"].append(
            f"--bfd_database_path={bfd_database_path}"
        )

    if benchmark:
        container_overrides["command"].append("--benchmark")

    if use_precomputed_msas:
        container_overrides["command"].append("--use_precomputed_msas")

    if features_paths is not None:
        container_overrides["command"].append(f"--features_paths={features_paths}")

    if run_features_only:
        container_overrides["command"].append("--run_features_only")

    if logtostderr:
        container_overrides["command"].append("--logtostderr")

    if gpu > 0:
        if use_spot_instances:
            print("Spot instance queue not available for GPU jobs. Using on-demand queue instead.")
        job_definition = batch_resources["gpu_job_definition"]
        job_queue = batch_resources["gpu_job_queue"]
        container_overrides["resourceRequirements"].append(
            {"value": str(gpu), "type": "GPU"}
        )
    else:
        job_definition = batch_resources["cpu_job_definition"]
        if use_spot_instances and batch_resources["cpu_job_queue_spot"] is not None:
            job_queue = batch_resources["cpu_job_queue_spot"]
        elif use_spot_instances and batch_resources["cpu_job_queue_spot"] is None:
            print("Spot instance queue not available. Using on-demand queue instead.")
            job_queue = batch_resources["cpu_job_queue_od"]
        else:
            job_queue = batch_resources["cpu_job_queue_od"]

    print(container_overrides)
    if depends_on is None:
        response = batch.submit_job(
            jobDefinition=job_definition,
            jobName=job_name,
            jobQueue=job_queue,
            containerOverrides=container_overrides,
        )
    else:
        response = batch.submit_job(
            jobDefinition=job_definition,
            jobName=job_name,
            jobQueue=job_queue,
            containerOverrides=container_overrides,
            dependsOn=[{"jobId": depends_on, "type": "SEQUENTIAL"}],
        )

    return response

def get_run_metrics(bucket, job_name):
    timings_uri = sagemaker.s3.s3_path_join(bucket, job_name, "timings.json")
    ranking_uri = sagemaker.s3.s3_path_join(bucket, job_name, "ranking_debug.json")
    downloader = sagemaker.s3.S3Downloader()
    timing_dict = json.loads(downloader.read_file(f"s3://{timings_uri}"))
    ranking_dict = json.loads(downloader.read_file(f"s3://{ranking_uri}"))

    timing_df = pd.DataFrame.from_dict(
        timing_dict, orient="index", columns=["duration_sec"]
    )
    ranking_plddts_df = pd.DataFrame.from_dict(
        ranking_dict["plddts"], orient="index", columns=["plddts"]
    )
    order_df = pd.DataFrame.from_dict(ranking_dict["order"])
    return (timing_df, ranking_plddts_df, order_df)


def validate_input(input_sequences):
    output = []
    for sequence in input_sequences:
        sequence = sequence.upper().strip()
        if re.search("[^ARNDCQEGHILKMFPSTWYV]", sequence):
            raise ValueError(
                f"Input sequence contains invalid amino acid symbols." f"{sequence}"
            )
        output.append(sequence)

    if len(output) == 1:
        print("Using the monomer models.")
        model_preset = "monomer"
        return output, model_preset
    elif len(output) > 1:
        print("Using the multimer models.")
        model_preset = "multimer"
        return output, model_preset
    else:
        raise ValueError("Please provide at least one input sequence.")



### ---------------------------------------------
# Original Copyright 2021 Sergey Ovchinnikov https://github.com/sokrypton/ColabFold
# Modifications Copyright 2022 Amazon.com, Inc. or its affiliates. All Rights Reserved.

pymol_color_list = [
    "#33ff33",
    "#00ffff",
    "#ff33cc",
    "#ffff00",
    "#ff9999",
    "#e5e5e5",
    "#7f7fff",
    "#ff7f00",
    "#7fff7f",
    "#199999",
    "#ff007f",
    "#ffdd5e",
    "#8c3f99",
    "#b2b2b2",
    "#007fff",
    "#c4b200",
    "#8cb266",
    "#00bfbf",
    "#b27f7f",
    "#fcd1a5",
    "#ff7f7f",
    "#ffbfdd",
    "#7fffff",
    "#ffff7f",
    "#00ff7f",
    "#337fcc",
    "#d8337f",
    "#bfff3f",
    "#ff7fff",
    "#d8d8ff",
    "#3fffbf",
    "#b78c4c",
    "#339933",
    "#66b2b2",
    "#ba8c84",
    "#84bf00",
    "#b24c66",
    "#7f7f7f",
    "#3f3fa5",
    "#a5512b",
]

alphabet_list = list(ascii_uppercase + ascii_lowercase)

def plot_pdb(
    pred_output_path,
    show_sidechains=False,
    show_mainchains=False,
    color="lDDT",
    chains=None,
    Ls=None,
    vmin=50,
    vmax=90,
    color_HP=False,
    size=(800, 480),
):

    """
    Create a 3D view of a pdb structure
    Copied from https://github.com/sokrypton/ColabFold/blob/main/beta/colabfold.py
    """

    if chains is None:
        chains = 1 if Ls is None else len(Ls)

    view = py3Dmol.view(
        js="https://3dmol.org/build/3Dmol.js", width=size[0], height=size[1]
    )
    view.addModel(open(pred_output_path,'r').read(),'pdb')
    if color == "lDDT":
        view.setStyle(
            {
                "cartoon": {
                    "colorscheme": {
                        "prop": "b",
                        "gradient": "roygb",
                        "min": vmin,
                        "max": vmax,
                    }
                }
            }
        )
    elif color == "rainbow":
        view.setStyle({"cartoon": {"color": "spectrum"}})
    elif color == "chain":
        for n, chain, color in zip(range(chains), alphabet_list, pymol_color_list):
            view.setStyle({"chain": chain}, {"cartoon": {"color": color}})
    if show_sidechains:
        BB = ["C", "O", "N"]
        HP = [
            "ALA",
            "GLY",
            "VAL",
            "ILE",
            "LEU",
            "PHE",
            "MET",
            "PRO",
            "TRP",
            "CYS",
            "TYR",
        ]
        if color_HP:
            view.addStyle(
                {"and": [{"resn": HP}, {"atom": BB, "invert": True}]},
                {"stick": {"colorscheme": "yellowCarbon", "radius": 0.3}},
            )
            view.addStyle(
                {"and": [{"resn": HP, "invert": True}, {"atom": BB, "invert": True}]},
                {"stick": {"colorscheme": "whiteCarbon", "radius": 0.3}},
            )
            view.addStyle(
                {"and": [{"resn": "GLY"}, {"atom": "CA"}]},
                {"sphere": {"colorscheme": "yellowCarbon", "radius": 0.3}},
            )
            view.addStyle(
                {"and": [{"resn": "PRO"}, {"atom": ["C", "O"], "invert": True}]},
                {"stick": {"colorscheme": "yellowCarbon", "radius": 0.3}},
            )
        else:
            view.addStyle(
                {
                    "and": [
                        {"resn": ["GLY", "PRO"], "invert": True},
                        {"atom": BB, "invert": True},
                    ]
                },
                {"stick": {"colorscheme": f"WhiteCarbon", "radius": 0.3}},
            )
            view.addStyle(
                {"and": [{"resn": "GLY"}, {"atom": "CA"}]},
                {"sphere": {"colorscheme": f"WhiteCarbon", "radius": 0.3}},
            )
            view.addStyle(
                {"and": [{"resn": "PRO"}, {"atom": ["C", "O"], "invert": True}]},
                {"stick": {"colorscheme": f"WhiteCarbon", "radius": 0.3}},
            )
    if show_mainchains:
        BB = ["C", "O", "N", "CA"]
        view.addStyle(
            {"atom": BB}, {"stick": {"colorscheme": f"WhiteCarbon", "radius": 0.3}}
        )
    view.zoomTo()
    return view

def plot_plddt_legend(dpi=100):

    """
    Create 3D Plot legend
    Copied from https://github.com/sokrypton/ColabFold/blob/main/beta/colabfold.py
    """

    thresh = [
        "plDDT:",
        "Very low (<50)",
        "Low (60)",
        "OK (70)",
        "Confident (80)",
        "Very high (>90)",
    ]
    plt.figure(figsize=(1, 0.1), dpi=dpi)
    ########################################
    for c in ["#FFFFFF", "#FF0000", "#FFFF00", "#00FF00", "#00FFFF", "#0000FF"]:
        plt.bar(0, 0, color=c)
    plt.legend(
        thresh,
        frameon=False,
        loc="center",
        ncol=6,
        handletextpad=1,
        columnspacing=1,
        markerscale=0.5,
    )
    plt.axis(False)
    return plt
### ---------------------------------------------
