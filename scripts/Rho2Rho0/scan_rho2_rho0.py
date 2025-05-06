#!/usr/bin/env python3
"""
scan_rho2_rho0.py

Sweep the rho2scale factor for bwgen from 0.5 to 1.5 (inclusive) using 5 evenly‑spaced
points.  By default the script iterates over the five centrality bins

    55  45  35  25  15

and automatically assigns the corresponding number of events

    500 000, 300 000, 200 000, 150 000, 50 000.

You can override the centralities (and optionally the events list) from the command line.

Usage
-----
    ./scan_rho2_rho0.py [options]

Options
-------
  -p, --program    Path to the bwgen executable               [default: ../../build/bwgen]
  -c, --centrality Comma/space‑separated list of bins         [default: 55 45 35 25 15]
  -n, --nEvents    Comma/space‑separated list of event counts (must match centrality list)
                    If omitted while using the default centralities, the default event
                    counts above are applied.
  -b, --betascale  Comma/space‑separated list of betascale values (default: 0.5 … 1.5)
  -j, --jobs       Number of parallel jobs (default: 6)
  --dry            Only print commands (do not execute)
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path
import numpy as np
import concurrent.futures


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Parameter scan for rho2scale.")
    parser.add_argument(
        "-p",
        "--program",
        default="../../build/bwgen",
        help="Path to the bwgen executable (default: %(default)s)",
    )
    parser.add_argument(
        "-c",
        "--centrality",
        type=str,
        default=None,
        help="Centrality bin(s), e.g. '55,45,35'  (default: 55 45 35 25 15)",
    )
    parser.add_argument(
        "-n",
        "--nEvents",
        type=str,
        default=None,
        help="Event count(s) matching centrality list (default depends on centrals)",
    )
    parser.add_argument(
        "-b",
        "--betascale",
        type=str,
        default=None,
        help="betascale values; comma/space separated; default same as rho2scale list",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=1,
        help="Number of parallel jobs (default: 1)",
    )
    parser.add_argument(
        "--dry",
        action="store_true",
        help="Print commands without executing them",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # Determine centrality bins and corresponding event counts
    if args.centrality is None:
        centralities = [55, 45, 35, 25, 15]
    else:
        centralities = [int(c) for c in args.centrality.replace(',', ' ').split()]

    if args.nEvents is None:
        if args.centrality is None:
            events_list = [500_000, 300_000, 200_000, 150_000, 50_000]
        else:
            events_list = [50_000] * len(centralities)
    else:
        events_list = [int(n) for n in args.nEvents.replace(',', ' ').split()]

    if len(events_list) != len(centralities):
        sys.exit("ERROR: The nEvents list must have the same length as the centrality list.")

    bwgen_path = Path(args.program).expanduser()
    if not bwgen_path.is_file():
        sys.exit(f"bwgen executable not found: {bwgen_path}")

    # 5 points from 0.5 to 1.5 inclusive
    rho_scales = np.linspace(0.5, 1.5, 5)
    if args.betascale is None:
        beta_scales = rho_scales
    else:
        beta_scales = np.array([float(b) for b in args.betascale.replace(',', ' ').split()])
    print(f"Scanning rho2scale values: {', '.join(f'{s:.2f}' for s in rho_scales)}")
    print(f"Scanning betascale values:  {', '.join(f'{b:.2f}' for b in beta_scales)}")

    tasks: list[list[str]] = []

    for cent, nev in zip(centralities, events_list):
        print(f"\n=== Centrality {cent}  |  nEvents = {nev:,} ===")

        for beta in beta_scales:
            for rho in rho_scales:
                cmd = [
                    str(bwgen_path),
                    "-C", str(cent),
                    "-R", f"{rho:.3f}",
                    "-B", f"{beta:.3f}",
                    "-n", str(nev),
                ]

                print(" ".join(cmd))
                tasks.append(cmd)

    # If dry run, we are done
    if args.dry:
        return

    def run_cmd(cmd: list[str]) -> int:
        return subprocess.run(cmd, check=False).returncode

    if args.jobs <= 1:
        # Serial execution
        for cmd in tasks:
            rc = run_cmd(cmd)
            if rc != 0:
                print(f"[ERROR] Command failed (exit {rc}) -> {' '.join(cmd)}")
                sys.exit(rc)
    else:
        # Parallel execution
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            future_to_cmd = {pool.submit(run_cmd, cmd): cmd for cmd in tasks}
            for fut in concurrent.futures.as_completed(future_to_cmd):
                rc = fut.result()
                cmd = future_to_cmd[fut]
                if rc != 0:
                    print(f"[ERROR] Command failed (exit {rc}) -> {' '.join(cmd)}")
                    # Do not exit immediately; report and continue


if __name__ == "__main__":
    main()
