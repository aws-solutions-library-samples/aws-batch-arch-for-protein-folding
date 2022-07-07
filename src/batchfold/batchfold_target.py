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

    id: str
    s3_bucket: str = field()
    s3_prefix: str = ""
    boto_session: boto3.session.Session = boto3.DEFAULT_SESSION or boto3.Session()
    sequences: list = []

    @s3_bucket.validator
    def _check_bucket_exists(self, attribute, value):
        """Validate that s3 bucket exists"""
        self.boto_session.client("s3").head_bucket(Bucket=value)

    def add_sequence(self, id: str, seq: str, description: str = "") -> SeqRecord:
        """Add sequence to target object."""

        self.validate_sequence(seq)
        seq_record = SeqRecord(
            seq=Seq(seq),
            id=id,
            description=description,
            annotations={"molecule_type": "protein"},
        )
        self.sequences.append(seq_record)
        return seq_record

    def upload_fasta(self) -> str:
        """Create and upload a fasta file to s3"""
        file_out = f"{self.id}.fasta"
        with open(file_out, "w") as f_out:
            SeqIO.write(self.sequences, f_out, "fasta")
        s3_key = os.path.join(self.s3_prefix, file_out)
        self.boto_session.client("s3").upload_file(file_out, self.s3_bucket, s3_key)
        os.remove(file_out)
        return self.get_s3_uri()

    def get_s3_uri(self) -> str:
        """Get the s3 uri for the fasta file"""
        return os.path.join("s3://", self.s3_bucket, self.s3_prefix, self.id + ".fasta")

    def download_fasta(self, filename: str) -> str:
        """Download fasta file from s3"""
        self.boto_session.resource("s3").meta.client.download_file(
            self.s3_bucket, os.path.join(self.s3_prefix, self.id + ".fasta"), filename
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
