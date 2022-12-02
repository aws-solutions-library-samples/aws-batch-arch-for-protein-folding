import torch
import esm
from timeit import default_timer as timer
import logging
import json
import argparse
import sys

logger = logging.getLogger()
logger.setLevel(logging.INFO)

formatter = logging.Formatter(
    "%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%y/%m/%d %H:%M:%S",
)

console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

logger.info("Loading model")
model = esm.pretrained.esmfold_v1()
logger.info("Configuring model for CPU inference.")
model = model.eval()
model.esm.float()
model.cpu()

def lambda_handler(event, context):
    """
    sequence
    chunk_size
    num_recycles
    """
    logger.info(f"Event is:\n{event}")

    if "chunk_size" in event:
        chunk_size = event["chunk_size"]
    else:
        chunk_size = 128

    if "num_recycles" in event:
        num_recycles = event["num_recycles"]
    else:
        num_recycles = 4

    start = timer()

    model.set_chunk_size(chunk_size)
    logger.info("Starting inference")
    output = model.infer([event["sequence"]], num_recycles=num_recycles)
    logger.info("Generating pdb string")
    pdb = model.output_to_pdb(output)

    tottime = timer() - start
    metrics = {
        "time": tottime,
        "mean_plddt": output["mean_plddt"].item(),
        "ptm": output["ptm"].item(),
    }

    logger.info(f"Metrics:\n{metrics}")
    logger.info(f"PDB:\n{pdb}")

    return {
        "statusCode": 200,
        "headers": {"Content-Type": "application/json"},
        "body": json.dumps({"metrics": metrics, "pdb": pdb}),
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('sequence', type=str)
    parser.add_argument('--chunk_size', type=int, default=128)
    parser.add_argument('--num_recycles', type=int, default=4)

    args = parser.parse_args()

    event = {
        "sequence": args.sequence,
        "chunk_size": args.chunk_size,
        "num_recycles": args.num_recycles
    }
    lambda_handler(event, {})