import boto3
import os
from datetime import datetime
from batchfold.batchfold_environment import BatchFoldEnvironment
from batchfold.nextflow_job import NextFlowJob
from batchfold.utils import utils


def main():

    # Create AWS clients
    boto_session = boto3.session.Session()
    batch_environment = BatchFoldEnvironment(boto_session=boto_session)
    s3 = boto_session.client("s3")

    S3_BUCKET = batch_environment.default_bucket
    print(f" S3 bucket name is {S3_BUCKET}") 

    random_str = datetime.now().strftime("%Y%m%d%s")
    nextflow_script = "run_rfdesign_esmfold_multiple_sequences.nf"
    asset_prefix = "assets_input"
    input_prefix = 'pd1-demo'
    rf_design_prefix = "myrfdesign_hallucination_" + random_str
    esmfold_prefix = "FinalESMFoldOutput_" + random_str

    # copy pdb structures to S3
    s3.upload_file('pd1_demo/pd1.pdb', S3_BUCKET, os.path.join(input_prefix, 'pd1.pdb'))
    s3.upload_file('pd1_demo/pdl1.pdb', S3_BUCKET, os.path.join(input_prefix, 'pdl1.pdb'))

    # copy NextFlow script to S3
    s3.upload_file(nextflow_script, S3_BUCKET, os.path.join(asset_prefix, nextflow_script))

    # copy additional scripts to s3
    utils.upload_dir(bucket = S3_BUCKET, local_path = "bin", prefix = os.path.join(asset_prefix, "bin"), boto_session=boto_session)

    # Submit job
    job_name = "NexFlowJob_" + datetime.now().strftime("%Y%m%d%s")
    nextflow_job = NextFlowJob(
        boto_session=boto_session,
        job_name = job_name,
        assets_s3_uri = os.path.join("s3://", S3_BUCKET, asset_prefix),
        nextflow_script = nextflow_script,
        params = {
            "s3_input": os.path.join("s3://", S3_BUCKET, input_prefix + '/'), # The slash is important here - you want to download the "folder" recursively
            "rf_design_output": os.path.join("s3://", S3_BUCKET, rf_design_prefix),
            "esmfold_output": os.path.join("s3://", S3_BUCKET, esmfold_prefix),
        }
    )

    nextflow_submission = batch_environment.submit_job(
        nextflow_job, job_queue_name="CPUOnDemandJobQueue"
    )

    print(nextflow_submission)


if __name__ == "__main__":
    main()
