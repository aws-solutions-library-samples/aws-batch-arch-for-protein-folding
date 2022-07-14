from __future__ import annotations

from attr import define
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from attrs import define, field
import boto3
import re
import os

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
    sequences: list = []

    @s3_bucket.validator
    def check_bucket_exists(self, attribute, value):
        """Validate that s3 bucket exists"""
        self.boto_session.client("s3").head_bucket(Bucket=value)

    @s3_base_prefix.default
    def define_s3_base_prefix(self) -> None:
        return self.target_id

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
        self.sequences.append(seq_record)
        return self

    def upload_fasta(self) -> BatchFoldTarget:
        """Create and upload a fasta file to s3"""
        file_out = f"{self.target_id}.fasta"
        with open(file_out, "w") as f_out:
            SeqIO.write(self.sequences, f_out, "fasta")
        s3_fasta_key = os.path.join(self.s3_base_prefix, self.s3_fastas_prefix, file_out)
        self.boto_session.client("s3").upload_file(file_out, self.s3_bucket, s3_fasta_key)
        os.remove(file_out)
        return self

    def get_fasta_s3_uri(self) -> str:
        """Get the s3 uri for the fasta file"""
        return os.path.join("s3://", self.s3_bucket, self.s3_base_prefix, self.s3_fastas_prefix, self.target_id + ".fasta")

    def get_msas_s3_uri(self) -> str:
        """Get the s3 uri for the msas folder"""
        return os.path.join("s3://", self.s3_bucket, self.s3_base_prefix, self.s3_msas_prefix)

    def get_predictions_s3_uri(self) -> str:
        """Get the s3 uri for the predictions folder"""
        return os.path.join("s3://", self.s3_bucket, self.s3_base_prefix, self.s3_predictions_prefix)

    def download_fasta(self, filename: str) -> str:
        """Download fasta file from s3"""
        self.boto_session.resource("s3").meta.client.download_file(
            self.s3_bucket, os.path.join(self.s3_base_prefix, self.s3_fastas_prefix, self.target_id + ".fasta"), filename
        )
        return os.path.abspath(filename)

    @staticmethod
    def validate_sequence(sequence: str) -> bool:
        """Validate that the sequences only contains valid amino acid symbols."""
        sequence = sequence.upper().strip()
        if re.search("[^ARNDCQEGHILKMFPSTWYV]", sequence):
            raise ValueError(
                f"Input sequence contains invalid amino acid symbols." f"{sequence}"
            )
        return True
