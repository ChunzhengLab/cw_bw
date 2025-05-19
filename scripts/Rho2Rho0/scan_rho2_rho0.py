#!/usr/bin/env python3
# scan_rho2_rho0.py

from __future__ import annotations
import argparse
import os
import subprocess
import sys
from pathlib import Path
import numpy as np
import concurrent.futures

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Parameter scan for rho2scale and betascale.")
    parser.add_argument("-p", "--program", default="../../build/bwgen", help="Path to the bwgen executable")
    parser.add_argument("-c", "--centrality", type=str, default=None, help="Centrality bins, e.g. '55,45'")
    parser.add_argument("-n", "--nEvents", type=str, default=None, help="Event counts per bin")
    parser.add_argument("-b", "--betascale", type=str, default=None, help="Betascale values")
    parser.add_argument("-r", "--rho2scale", type=str, default=None, help="Rho2scale values")
    parser.add_argument("-j", "--jobs", type=int, default=6, help="Number of parallel jobs")
    parser.add_argument("--dry", action="store_true", help="Print commands only")
    parser.add_argument("--mode", choices=["rho_only", "beta_only", "pair", "both", "single"], default="rho_only",
                        help="Scan mode: rho_only, beta_only, pair, both, single (default: rho_only)")
    return parser.parse_args()

def main() -> None:
    args = parse_args()

    if args.centrality is None:
        centralities = [55, 45, 35, 25, 15]
    else:
        centralities = [int(c) for c in args.centrality.replace(',', ' ').split()]

    if args.nEvents is None:
        events_list = [500_000, 300_000, 200_000, 150_000, 50_000][:len(centralities)]
    else:
        events_list = [int(n) for n in args.nEvents.replace(',', ' ').split()]

    if len(events_list) != len(centralities):
        sys.exit("ERROR: The nEvents list must match centrality list length.")

    bwgen_path = Path(args.program).expanduser()
    if not bwgen_path.is_file():
        sys.exit(f"bwgen not found: {bwgen_path}")

    rho_list = np.linspace(0.9, 1.1, 5) if args.rho2scale is None else \
        np.array([float(x) for x in args.rho2scale.replace(',', ' ').split()])
    beta_list = np.linspace(0.9, 1.1, 5) if args.betascale is None else \
        np.array([float(x) for x in args.betascale.replace(',', ' ').split()])

    print(f"Mode = {args.mode}")
    print(f"rho2scale:  {', '.join(f'{x:.2f}' for x in rho_list)}")
    print(f"betascale:  {', '.join(f'{x:.2f}' for x in beta_list)}")

    tasks: list[list[str]] = []

    for cent, nev in zip(centralities, events_list):
        print(f"\n=== Centrality {cent}  |  nEvents = {nev:,} ===")

        pairs = []
        if args.mode == "rho_only":
            pairs = [(rho, 1.0) for rho in rho_list]
        elif args.mode == "beta_only":
            pairs = [(1.0, beta) for beta in beta_list]
        elif args.mode == "pair":
            if len(rho_list) != len(beta_list):
                sys.exit("In pair mode, rho2scale and betascale must be same length.")
            pairs = list(zip(rho_list, beta_list))
        elif args.mode == "both":
            pairs = [(r, b) for r in rho_list for b in beta_list]
        elif args.mode == "single":
            pairs = [(rho, 1.0) for rho in rho_list] + [(1.0, beta) for beta in beta_list]
        else:
            sys.exit(f"Unknown mode: {args.mode}")

        for rho, beta in pairs:
            cmd = [
                str(bwgen_path),
                "-C", str(cent),
                "-R", f"{rho:.3f}",
                "-B", f"{beta:.3f}",
                "-n", str(nev),
                "-c", "../../configs/default.yaml"
            ]
            print(" ".join(cmd))
            tasks.append(cmd)

    if args.dry:
        return

    def run_cmd(cmd: list[str]) -> int:
        return subprocess.run(cmd, check=False).returncode

    if args.jobs <= 1:
        for cmd in tasks:
            rc = run_cmd(cmd)
            if rc != 0:
                print(f"[ERROR] Command failed (exit {rc}) -> {' '.join(cmd)}")
                sys.exit(rc)
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            future_to_cmd = {pool.submit(run_cmd, cmd): cmd for cmd in tasks}
            for fut in concurrent.futures.as_completed(future_to_cmd):
                rc = fut.result()
                cmd = future_to_cmd[fut]
                if rc != 0:
                    print(f"[ERROR] Command failed (exit {rc}) -> {' '.join(cmd)}")

if __name__ == "__main__":
    main()
