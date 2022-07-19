import boto3
from batchfold.batchfold_environment import BatchFoldEnvironment
from batchfold.download_job import DownloadJob
import logging
from datetime import datetime
from __future__ import print_function
import urllib3
import json

LOGGER = logging.getLogger()
LOGGER.setLevel(logging.INFO)
http = urllib3.PoolManager()

def main(job_queue_name="GravitonSpotJobQueue"):
    """ Download all data to file system """

    boto_session = boto3.session.Session(profile_name="bloyal+proteinfolding-Admin")
    region = boto_session.region_name
    batch_environment = BatchFoldEnvironment(boto_session = boto_session, region_name=region)
    
    download_test_submission = batch_environment.submit_job(
        DownloadJob(job_name="download_test" + datetime.now().strftime("%Y%m%dT%H%M%S"), script="./scripts/download_test.sh"),
        job_queue_name=job_queue_name,
    )

    download_alphafold_params_submission = batch_environment.submit_job(
        DownloadJob(job_name="download_alphafold_params" + datetime.now().strftime("%Y%m%dT%H%M%S"), script="./scripts/download_alphafold_params.sh"),
        job_queue_name=job_queue_name,
    )

    download_bfd_submission = batch_environment.submit_job(
        DownloadJob(job_name="download_bfd" + datetime.now().strftime("%Y%m%dT%H%M%S"), script="./scripts/download_bfd.sh"),
        job_queue_name=job_queue_name,
    )

    download_mgnify_submission = batch_environment.submit_job(
        DownloadJob(job_name="download_mgnify" + datetime.now().strftime("%Y%m%dT%H%M%S"), script="./scripts/download_mgnify.sh"),
        job_queue_name=job_queue_name,
    )

    download_openfold_params_submission = batch_environment.submit_job(
        DownloadJob(job_name="download_openfold_params" + datetime.now().strftime("%Y%m%dT%H%M%S"), script="./scripts/download_openfold_params.sh"),
        job_queue_name=job_queue_name,
    )

    download_pdb70_submission = batch_environment.submit_job(
        DownloadJob(job_name="download_pdb70" + datetime.now().strftime("%Y%m%dT%H%M%S"), script="./scripts/download_pdb70.sh"),
        job_queue_name=job_queue_name,
    )

    download_pdb_mmcif_submission = batch_environment.submit_job(
        DownloadJob(job_name="download_pdb_mmcif" + datetime.now().strftime("%Y%m%dT%H%M%S"), script="./scripts/download_pdb_mmcif.sh"),
        job_queue_name=job_queue_name,
    )

    download_pdb_seqres_submission = batch_environment.submit_job(
        DownloadJob(job_name="download_pdb_seqres" + datetime.now().strftime("%Y%m%dT%H%M%S"), script="./scripts/download_pdb_seqres.sh"),
        job_queue_name=job_queue_name,
    )

    download_small_bfd_submission = batch_environment.submit_job(
        DownloadJob(job_name="download_small_bfd" + datetime.now().strftime("%Y%m%dT%H%M%S"), script="./scripts/download_small_bfd.sh"),
        job_queue_name=job_queue_name,
    )

    download_uniclust30_submission = batch_environment.submit_job(
        DownloadJob(job_name="download_uniclust30" + datetime.now().strftime("%Y%m%dT%H%M%S"), script="./scripts/download_uniclust30.sh"),
        job_queue_name=job_queue_name,
    )

    download_uniprot_submission = batch_environment.submit_job(
        DownloadJob(job_name="download_uniprot" + datetime.now().strftime("%Y%m%dT%H%M%S"), script="./scripts/download_uniprot.sh"),
        job_queue_name=job_queue_name,
    )

    download_uniref30_submission = batch_environment.submit_job(
        DownloadJob(job_name="download_uniref30" + datetime.now().strftime("%Y%m%dT%H%M%S"), script="./scripts/download_uniref30.sh"),
        job_queue_name=job_queue_name,
    )

    download_uniref90_submission = batch_environment.submit_job(
        DownloadJob(job_name="download_uniref90" + datetime.now().strftime("%Y%m%dT%H%M%S"), script="./scripts/download_uniref90.sh"),
        job_queue_name=job_queue_name,
    )

    download_colabfold_envdb_submission = batch_environment.submit_job(
        DownloadJob(job_name="download_colabfold_envdb" + datetime.now().strftime("%Y%m%dT%H%M%S"), script="./scripts/download_colabfold_envdb.sh"),
        job_queue_name=job_queue_name,
    )

    prep_mmseqs_dbs_submission = batch_environment.submit_job(
        DownloadJob(job_name="prep_mmseqs_dbs" + datetime.now().strftime("%Y%m%dT%H%M%S"), script="./scripts/prep_mmseqs_dbs.sh", memory=500, cpu=64),
        job_queue_name=job_queue_name,
        depends_on=[download_uniref30_submission, download_colabfold_envdb_submission],
    )

    response = [
        download_test_submission,
        download_alphafold_params_submission,
        download_bfd_submission,
        download_mgnify_submission,
        download_mgnify_submission,
        download_openfold_params_submission,
        download_pdb70_submission,
        download_pdb_mmcif_submission,
        download_pdb_seqres_submission,
        download_small_bfd_submission,
        download_uniclust30_submission,
        download_uniprot_submission,
        download_uniref30_submission,
        download_uniref90_submission,
        download_colabfold_envdb_submission,
        prep_mmseqs_dbs_submission,
    ]

    return(response)

def lambda_handler(event, context):
    try:
        LOGGER.info("REQUEST RECEIVED:\n %s", event)
        LOGGER.info("REQUEST RECEIVED:\n %s", context)
        if event["RequestType"] == "Create":
            LOGGER.info("CREATE!")
            
            main()
            
            send(
                event, context, "SUCCESS", {"response": "Resource creation successful!"}
            )
        elif event["RequestType"] == "Update":
            LOGGER.info("UPDATE!")
            send(
                event,
                context,
                "SUCCESS",
                {"response": "Resource update successful!"},
            )
        elif event["RequestType"] == "Delete":
            LOGGER.info("DELETE!")
            send(
                event,
                context,
                "SUCCESS",
                {"response": "Resource deletion successful!"},
            )
        else:
            LOGGER.info("FAILED!")
            send(
                event,
                context,
                "FAILED",
                {"response": "Unexpected event received from CloudFormation"},
            )
    except:
        LOGGER.info("FAILED!")
        send(
            event,
            context,
            "FAILED",
            {"response": "Exception during processing"},
        )

def send(event, context, responseStatus, responseData, physicalResourceId=None, noEcho=False, reason=None):
    responseUrl = event['ResponseURL']

    print(responseUrl)

    responseBody = {
        'Status' : responseStatus,
        'Reason' : reason or "See the details in CloudWatch Log Stream: {}".format(context.log_stream_name),
        'PhysicalResourceId' : physicalResourceId or context.log_stream_name,
        'StackId' : event['StackId'],
        'RequestId' : event['RequestId'],
        'LogicalResourceId' : event['LogicalResourceId'],
        'NoEcho' : noEcho,
        'Data' : responseData
    }

    json_responseBody = json.dumps(responseBody)

    print("Response body:")
    print(json_responseBody)

    headers = {
        'content-type' : '',
        'content-length' : str(len(json_responseBody))
    }

    try:
        response = http.request('PUT', responseUrl, headers=headers, body=json_responseBody)
        print("Status code:", response.status)


    except Exception as e:

        print("send(..) failed executing http.request(..):", e)


if __name__ == "__main__":
    main()