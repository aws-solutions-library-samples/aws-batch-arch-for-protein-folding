from batchfold.utils import utils
import pytest
# import boto3
# from datetime import datetime
import os

@pytest.mark.skip(reason="In progress")
def test_download_pdb_file():

    pdb_file_name = utils.download_rcsb_pdb_file(pdb_code="7FCC", output_dir="tests/data")
    assert os.path.exists(pdb_file_name) is True

    os.remove(pdb_file_name)
    assert os.path.exists(pdb_file_name) is False

@pytest.mark.skip(reason="In progress")
def test_download_fasta_file():

    fasta_file_name = utils.download_rcsb_fasta_file(pdb_code="7FCC", output_dir="tests/data")
    assert os.path.exists(fasta_file_name) is True

    os.remove(fasta_file_name)
    assert os.path.exists(fasta_file_name) is False

@pytest.mark.skip(reason="In progress")
@pytest.fixture()
def monomer():
    pdb_file_name = utils.download_rcsb_pdb_file(pdb_code="7FCC", output_dir="tests/data")
    return(pdb_file_name)

@pytest.mark.skip(reason="In progress")
@pytest.fixture()
def multimer():
    pdb_file_name = utils.download_rcsb_pdb_file(pdb_code="4ZQK", output_dir="tests/data")
    return(pdb_file_name)

@pytest.mark.skip(reason="In progress")
def test_monomer_b_factors(monomer):
    bfactors1 = utils.get_bfactors(monomer)
    assert type(bfactors1) is dict
    assert len(bfactors1) == 1
    assert type(bfactors1[0]) is list
    assert len(bfactors1[0]) == 1
    assert type(bfactors1[0][0]) is list
    assert len(bfactors1[0][0]) == 146

@pytest.mark.skip(reason="In progress")
def test_multimer_b_factors(multimer):
    bfactors2 = utils.get_bfactors(multimer)
    assert type(bfactors2) is dict
    assert len(bfactors2) == 1
    assert type(bfactors2[0]) is list
    assert len(bfactors2[0]) == 2
    assert type(bfactors2[0][0]) is list
    assert len(bfactors2[0][0]) == 115
    assert type(bfactors2[0][1]) is list
    assert len(bfactors2[0][0]) == 115

