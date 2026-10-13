#!/bin/bash
# Run the prescribed-automorphism MIP search (kohnert.py) over one shard of targets.
# usage: codes/fleet/run_shard.sh SHARD SECS_PER_TARGET [WORKERS]
#   e.g. codes/fleet/run_shard.sh s07 500 4
# Log: runs/fleet_SHARD.log  (one line per target; "*** FOUND" lines mark candidates,
# which kohnert.py also appends to codes/results/candidates.jsonl).
set -e
cd "$(dirname "$0")/../.."
python3 -c "import highspy, numpy" 2>/dev/null || pip install -q highspy numpy
make -s -C codes/src libmindist.so
mkdir -p runs
cd codes/scripts
exec python3 kbatch.py ../fleet/shards/"$1".txt "$2" "${3:-4}" ../../runs/fleet_"$1".log 7777
