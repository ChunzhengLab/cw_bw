#!/usr/bin/env python3
# scan_pars.py

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
    parser.add_argument(
        "-f", "--frac-lbc",
        type=str,
        default=None,
        help="FracLBC values (0 to 0.16) to scan, comma-separated"
    )
    parser.add_argument("-j", "--jobs", type=int, default=8, help="Number of parallel jobs")
    parser.add_argument("--dry", action="store_true", help="Print commands only")
    parser.add_argument(
        "--mode",
        choices=["rho_only", "beta_only", "frac_only", "all"],
        default="rho_only",
        help="Scan mode: rho_only, beta_only, frac_only, all (default: rho_only)"
    )
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

    rho_list = (
        np.array([0.0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.0])
        if args.rho2scale is None
        else np.array([float(x) for x in args.rho2scale.replace(',', ' ').split()])
    )
    beta_list = (
        np.array([0.0, 0.4, 0.8, 0.96, 1.0, 1.06, 1.1, 1.16])
        if args.betascale is None
        else np.array([float(x) for x in args.betascale.replace(',', ' ').split()])
    )

    # Prepare fracLBC scan values
    frac_list = (
        np.array([0.0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16])
        if args.frac_lbc is None
        else np.array([float(x) for x in args.frac_lbc.replace(',', ' ').split()])
    )

    print(f"Mode = {args.mode}")
    print(f"rho2scale:  {', '.join(f'{x:.2f}' for x in rho_list)}")
    print(f"betascale:  {', '.join(f'{x:.2f}' for x in beta_list)}")
    print(f"fracLBC:  {', '.join(f'{x:.2f}' for x in frac_list)}")

    tasks: list[list[str]] = []
    for cent, nev in zip(centralities, events_list):
        print(f"\n=== Centrality {cent}  |  nEvents = {nev:,} ===")
        # Base command without parameter overrides
        base_cmd = [
            str(bwgen_path),
            "-C", str(cent),
            "-n", str(nev),
            "-c", "../../configs/default.yaml"
        ]

        if args.mode == "rho_only":
            for rho in rho_list:
                cmd = base_cmd + ["-R", f"{rho:.3f}"]
                print(" ".join(cmd))
                tasks.append(cmd)

        elif args.mode == "beta_only":
            for beta in beta_list:
                cmd = base_cmd + ["-B", f"{beta:.3f}"]
                print(" ".join(cmd))
                tasks.append(cmd)

        elif args.mode == "frac_only":
            for frac in frac_list:
                cmd = base_cmd + ["-f", f"{frac:.3f}"]
                print(" ".join(cmd))
                tasks.append(cmd)

        elif args.mode == "all":
            # Run all three only modes sequentially
            for rho in rho_list:
                cmd = base_cmd + ["-R", f"{rho:.3f}"]
                print(" ".join(cmd))
                tasks.append(cmd)
            for beta in beta_list:
                cmd = base_cmd + ["-B", f"{beta:.3f}"]
                print(" ".join(cmd))
                tasks.append(cmd)
            for frac in frac_list:
                cmd = base_cmd + ["-f", f"{frac:.3f}"]
                print(" ".join(cmd))
                tasks.append(cmd)

        else:
            sys.exit(f"Unknown mode: {args.mode}")

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
                try:
                    future_to_cmd = {pool.submit(run_cmd, cmd): cmd for cmd in tasks}
                    for fut in concurrent.futures.as_completed(future_to_cmd):
                        rc = fut.result()
                        cmd = future_to_cmd[fut]
                        if rc != 0:
                            print(f"[ERROR] Command failed (exit {rc}) -> {' '.join(cmd)}")
                except KeyboardInterrupt:
                    print("\n[INFO] Interrupted by user. Cancelling remaining tasks.")
                    pool.shutdown(wait=False)
                    sys.exit(1)

if __name__ == "__main__":
    main()
