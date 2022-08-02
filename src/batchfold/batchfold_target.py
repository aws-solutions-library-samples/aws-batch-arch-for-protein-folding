from __future__ import annotations

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0


from typing import List, Dict
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaIterator
from attrs import define, field
import boto3
import re
import os

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from attrs import define, field
import boto3
import re
import os
from io import  StringIO

@define
class BatchFoldTarget:
    """Define sequence targets"""

    target_id: str
    s3_bucket: str = field(kw_only=True, default="")
    s3_base_prefix: str = field(kw_only=True)
    s3_fastas_prefix: str = "fastas"
    s3_msas_prefix: str = "msas"
    s3_predictions_prefix: str = "predictions"
    boto_session: boto3.session.Session = boto3.DEFAULT_SESSION or boto3.Session()
    sequences: Dict = field(kw_only=True)

    @s3_bucket.validator
    def check_bucket_exists(self, attribute, value):
        """Validate that s3 bucket exists"""
        self.boto_session.client("s3").head_bucket(Bucket=value)

    @s3_base_prefix.default
    def define_s3_base_prefix(self) -> str:
        return self.target_id

    @sequences.default
    def init_sequences(self) -> Dict:
        return {}

    def __attrs_post_init__(self) -> None:
        """ Check and see if data already exists in S3 for this target and if so, load it"""
        s3 = self.boto_session.client("s3")
        bucket = self.s3_bucket
        key = os.path.join(self.s3_base_prefix, self.s3_fastas_prefix, self.target_id+".fasta")
        try:
            s3.head_object(
                Bucket=bucket, 
                Key=key
                )
        except:
            return None
        else:
            fasta_str = s3.get_object(Bucket=bucket, Key=key)['Body'].read().decode("utf-8")
            with StringIO(fasta_str) as stringio:
                    seq_records = SeqIO.parse(stringio, "fasta")
                    for seq_record in seq_records:
                        self.sequences[seq_record.id] = seq_record
            return None


    def add_sequence(self, seq: str, seq_id: str = "", description: str = "") -> BatchFoldTarget:
        """Add sequence to target object."""

        seq_id = seq_id or self.target_id
        self.validate_sequence(seq)
        seq_record = SeqRecord(
            seq=Seq(seq),
            id=seq_id,
            description=description,
            annotations={"molecule_type": "protein"},
        )

        self.sequences[seq_record.id] = seq_record
        return self

    def add_fasta(self, path: str) -> BatchFoldTarget:
        """Add an existing fasta file to the target object."""
        with open(path) as handle:
            for seq_record in SeqIO.parse(handle, "fasta"):
                self.sequences[seq_record.id] = seq_record
        
        return(self.upload_fasta(path))

    def upload_fasta(self) -> BatchFoldTarget:
        """Create and upload a fasta file to s3"""
        file_out = f"{self.target_id}.fasta"
        with open(file_out, "w") as f_out:
            SeqIO.write(list(self.sequences.values()), f_out, "fasta")
        s3_fasta_key = os.path.join(self.s3_base_prefix, self.s3_fastas_prefix, file_out)
        self.boto_session.client("s3").upload_file(file_out, self.s3_bucket, s3_fasta_key)
        os.remove(file_out)
        return os.join("s3://", self.s3_bucket, s3_fasta_key)

    def get_fasta_s3_uri(self) -> str:
        """Get the s3 uri for the fasta file"""
        return os.path.join("s3://", self.s3_bucket, self.s3_base_prefix, self.s3_fastas_prefix, self.target_id + ".fasta")

    def get_msas_s3_uri(self) -> str:
        """Get the s3 uri for the msas folder"""
        return os.path.join("s3://", self.s3_bucket, self.s3_base_prefix, self.s3_msas_prefix)

    def get_predictions_s3_uri(self) -> str:
        """Get the s3 uri for the predictions folder"""
        return os.path.join("s3://", self.s3_bucket, self.s3_base_prefix, self.s3_predictions_prefix)

    def download_fastas(self, local_path: str = ".") -> str:
        """Download fasta files from s3"""

        prefix = os.path.join(self.s3_base_prefix, self.s3_fastas_prefix)
        return self._download_dir(bucket = self.s3_bucket, local_path = local_path, prefix = prefix)

    def download_msas(self, local_path: str = ".", type: str = "all") -> str:
        """Download msa files from s3"""
        
        prefix = os.path.join(self.s3_base_prefix, self.s3_msas_prefix)
        if type in ["jackhmmer", "mmseqs2"]:
            prefix += "/" + type           

        return self._download_dir(bucket = self.s3_bucket, local_path = local_path, prefix = prefix)

    def download_predictions(self, local_path: str = ".", job: str = "") -> str:
        """Download prediction files from s3"""
        
        prefix = os.path.join(self.s3_base_prefix, self.s3_predictions_prefix)
        if job != "":
            prefix = os.path.join(prefix, job)           

        return self._download_dir(bucket = self.s3_bucket, local_path = local_path, prefix = prefix)

    def download_all(self, local_path: str = ".") -> str:
        """Download all files from s3"""       

        return self._download_dir(bucket = self.s3_bucket, local_path = local_path, prefix = self.s3_base_prefix)


    def _download_dir(self, client = None, bucket = None, local_path=".", prefix=""):
        """Recursively download files from S3."""
        
        bucket = bucket or self.s3_bucket
        client = client or self.boto_session.client("s3")
        paginator = client.get_paginator("list_objects_v2")
        file_count = 0
        for result in paginator.paginate(Bucket=bucket, Delimiter="/", Prefix=prefix):
            if result.get("CommonPrefixes") is not None:
                for subdir in result.get("CommonPrefixes"):
                    self._download_dir(
                        client=client, 
                        bucket=bucket, 
                        local_path=local_path, 
                        prefix=subdir.get("Prefix")
                        )
            for file in result.get("Contents", []):
                dest_pathname = os.path.join(local_path, file.get("Key"))
                if not os.path.exists(os.path.dirname(dest_pathname)):
                    os.makedirs(os.path.dirname(dest_pathname))
                client.download_file(bucket, file.get("Key"), dest_pathname)
                file_count += 1
        print(f"{file_count} files downloaded from s3.")
        return os.path.abspath(local_path)
        
    def validate_sequence(self, sequence: str) -> bool:
        """Validate that the sequences only contains valid amino acid symbols."""
        sequence = sequence.upper().strip()
        if re.search("[^ARNDCQEGHILKMFPSTWYV]", sequence):
            raise ValueError(
                f"Input sequence contains invalid amino acid symbols." f"{sequence}"
            )
        return True

    def list_job_names(self, bucket = None, client = None, job_type=""):

        bucket = bucket or self.s3_bucket
        client = client or self.boto_session.client("s3")
        paginator = client.get_paginator("list_objects_v2")

        if job_type not in [ "", "OpenFold", "AlphaFold2", "AlphaFold"]:
            raise ValueError("Valid job types are 'jackhmmer', 'mmseqs2', 'openfold', and 'alphafold2'")

        # Create a PageIterator from the Paginator
        page_iterator = paginator.paginate(
            Bucket=self.s3_bucket,
            Prefix=self.target_id+"/predictions")
        jobs = []
        for page in page_iterator:
            for key in page["Contents"]:
                jobs.append(re.search("predictions/(.*)/", key["Key"]).group(1))

        jobs = [*set(jobs)]

        if job_type != "":
            jobs = [job for job in jobs if job_type in job]

        jobs.sort()
        return jobs