from batchfold.batchfold_job import BatchFoldJob

def test_job_init():
    job = BatchFoldJob(
        name="MyNewJob",
        definition="MyJobDef",
        cpu = 6,
        memory = 9,
        gpu = 3,
        command = ["echo hello"],
        depends_on="MyOldJob"
    )

    assert job.name == "MyNewJob"
    assert job.definition == "MyJobDef"
    assert job.depends_on == "MyOldJob"
    assert job.container_overrides["command"] == ["echo hello"]
    assert len(job.container_overrides["resourceRequirements"]) == 3
    assert job.container_overrides["resourceRequirements"][0]["type"] == "VCPU"
    assert job.container_overrides["resourceRequirements"][0]["value"] == "6"
    assert job.container_overrides["resourceRequirements"][1]["type"] == "MEMORY"
    assert job.container_overrides["resourceRequirements"][1]["value"] == "9000"
    assert job.container_overrides["resourceRequirements"][2]["type"] == "GPU"
    assert job.container_overrides["resourceRequirements"][2]["value"] == "3"
