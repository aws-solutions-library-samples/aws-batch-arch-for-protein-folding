import boto3
import sys
import time


def wait_until_job_is_done(client, job_id):
    x = 1
    while x == 1:
        response = client.describe_jobs(jobs=[job_id])
        the_status = response["jobs"][0]["status"]
        if the_status in ["SUCCEEDED", "FAILED"]:
            return ()
        else:
            time.sleep(10)  # wait a bit before checking the status again


def main():
    client = boto3.client("batch")

    jobs_file = sys.argv[1]
    jobs_list = open(jobs_file).readlines()
    jobs_list = [i.rstrip() for i in jobs_list]

    for i in jobs_list:
        wait_until_job_is_done(client, i)


if __name__ == "__main__":
    main()
