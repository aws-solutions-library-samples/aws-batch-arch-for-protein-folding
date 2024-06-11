Todo:

Working Dockerfile is located on Superfold EC2 instance

Here are the changes to be made:

export USER=Something
ln -s [FileSystem] /net/scratch/db/af2/data/
ln -s [FileSystem] /net/scratch/db/af2/config.json

echo [FileSystem] > alphafold_params.pth

== /mnt/efs

Fix the	error:
<00:00, 232.73s/it]
Input: test:   0%|                                                                                                                       | 0/1 [03:52<?, ?i>
Traceback (most recent call last):
  File "run_superfold.py", line 1578, in <module>
    report(key)
  File "run_superfold.py", line 1325, in report
    target.parsed_pdb.make_pdb_file("TEST_RAW_TARGET.pdb")
AttributeError: 'NoneType' object has no attribute 'make_pdb_file'

                target.parsed_pdb.make_pdb_file("TEST_RAW_TARGET.pdb")
                if target.parsed_pdb is not None:

Looks like 1325 is accessing target.parsed_pdb before checking on 1326 if target.parsed_pdb is not None...

Test with:
run_superfold.py ~/input_test/test.fa --models 4 --out_dir ~/output_test/ --max_recycles 5​