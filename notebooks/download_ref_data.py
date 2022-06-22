# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from nbhelpers import nbhelpers
import boto3
import argparse

batch = boto3.client("batch")


def _parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument("--batch_substack_name", type=str, default=None)
    parser.add_argument("--job_name", type=str, default="download_job")
    parser.add_argument("--script", type=str, default="all")
    parser.add_argument("--cpu", type=int, default=4)
    parser.add_argument("--memory", type=int, default=16)
    parser.add_argument("--download_dir", type=str, default="/fsx")
    parser.add_argument("--download_mode", type=str, default="reduced_dbs")

    return parser.parse_known_args()


def submit_download_data_job(
    batch_substack_name,
    job_name,
    script,
    cpu,
    memory,
    download_dir,
    download_mode,
):

    if batch_substack_name is None:
        batch_substack_name = nbhelpers.list_alphafold_stacks()[0]["StackName"]
    batch_resources = nbhelpers.get_batch_resources(batch_substack_name)
    job_definition = batch_resources["download_job_definition"]
    job_queue = batch_resources["download_job_queue"]

    container_overrides = {
        "command": [script, download_dir, download_mode],
        "resourceRequirements": [
            {"value": str(cpu), "type": "VCPU"},
            {"value": str(memory * 1000), "type": "MEMORY"},
        ],
    }

    response = batch.submit_job(
        jobDefinition=job_definition,
        jobName=job_name,
        jobQueue=job_queue,
        containerOverrides=container_overrides,
    )

    return response


if __name__ == "__main__":

    ### Command line parser
    args, _ = _parse_args()

    if args.script.upper() == "ALL":
        all_scripts = [
            "download_alphafold_params_s3.sh",
            "download_bfd_s3.sh",
            "download_mgnify_s3.sh",
            "download_pdb70_s3.sh",
            "download_pdb_mmcif_s3.sh",
            "download_pdb_seqres_s3.sh",
            "download_small_bfd_s3.sh",
            "download_uniclust30_s3.sh",
            "download_uniprot.sh",
            "download_uniref90.sh",
        ]

        responses = []
        for script in all_scripts:
            script_response = submit_download_data_job(
                batch_substack_name=args.batch_substack_name,
                job_name=args.job_name + "-" + script.replace(".sh",""),
                script=script,
                cpu=args.cpu,
                memory=args.memory,
                download_dir=args.download_dir,
                download_mode=args.download_mode,
            )
            responses.append(script_response)
        response = str(responses)

    elif args.script.upper() == "PARAMETERS_ONLY":
        parameters_only_script = "download_alphafold_params_s3.sh"
        response = submit_download_data_job(
            batch_substack_name=args.batch_substack_name,
            job_name=args.job_name,
            script=parameters_only_script,
            cpu=args.cpu,
            memory=args.memory,
            download_dir=args.download_dir,
            download_mode=args.download_mode,
        )

    else:
        response = submit_download_data_job(
            batch_substack_name=args.batch_substack_name,
            job_name=args.job_name,
            script=args.script,
            cpu=args.cpu,
            memory=args.memory,
            download_dir=args.download_dir,
            download_mode=args.download_mode,
        )
    print(response)
