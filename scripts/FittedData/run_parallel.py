#!/usr/bin/env python3
"""
run_bwgen_parallel.py

Parallel Python script to automate bwgen runs for various centrality and event counts.
"""
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse

parser = argparse.ArgumentParser(description="Parallel bwgen runner")
parser.add_argument(
    "-c", "--config",
    help="Configuration YAML file",
    default="../../configs/default.yaml"
)
args = parser.parse_args()

# Define centralities and corresponding event counts
centralities = [55, 45, 35, 25, 15]
events = [500000, 300000, 150000, 100000, 50000]
# centralities = [45, 35]
# events = [500000, 200000]
# centralities = [25, 35]
# events = [100000, 200000]

# Path to bwgen executable
bwgen = '../../build/bwgen'

# Worker function
def run_bwgen(cent, n):
    cmd = [bwgen, '-C', str(cent), '-n', str(n), '-c', args.config]
    print(f"Starting bwgen for centrality={cent}, events={n} with config={args.config}")
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"Completed bwgen for centrality={cent}:\n{result.stdout}")
        return cent, None
    except subprocess.CalledProcessError as e:
        import sys
        print(f"Failed bwgen for centrality={cent}:\n{e.stderr}", file=sys.stderr)
        return cent, e.stderr


def main():
    with ThreadPoolExecutor(max_workers=len(centralities)) as executor:
        futures = {executor.submit(run_bwgen, cent, n): cent for cent, n in zip(centralities, events)}
        for future in as_completed(futures):
            cent = futures[future]
            err = None
            try:
                cent, err = future.result()
            except Exception as ex:
                print(f"Centrality={cent} run failed with exception: {ex}", file=sys.stderr)
                continue

            if err:
                print(f"Error in centrality={cent}:\n{err}", file=sys.stderr)
            else:
                print(f"Centrality={cent} completed successfully.")

    print("All parallel runs completed.")

if __name__ == '__main__':
    main()
