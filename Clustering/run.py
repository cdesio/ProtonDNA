"""Cluster DNA damage from a replay ROOT file, one voxel at a time."""
import argparse
from runClustering import runClustering

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--filename", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--sugar", required=True)
parser.add_argument("--copy", type=int, required=True, help="TrackingGrid voxel ID")
parser.add_argument("--min", type=float, default=5.0)
parser.add_argument("--max", type=float, default=37.5)
parser.add_argument("--indirect", type=float, default=0.405)
parser.add_argument("--continuous", action="store_true", help="Use continuous DNA indexing")
args = parser.parse_args()
runClustering(args.filename, args.output, args.min, args.max, args.indirect, args.sugar,
              simulationType="ps", continuous=args.continuous, part1_CopyNum=args.copy)
