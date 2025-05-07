#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
run_banch.py

并行批量驱动 `bwgen` 的执行。
每组中心度和 fLBC 值组合调用一次 bwgen。
输出日志（默认）：../../bwgen_log/results_cent<cent>_fLBC<val>.log

使用示例：
  python run_banch.py
  python run_banch.py --flbc-list 0.03 0.05
  python run_banch.py --cent-list 15 35 55
  python run_banch.py --no-log
"""

import argparse
import itertools
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
import os

# 日志设置
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] [%(name)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger("BatchDriver")

# 默认配置
DEFAULT_CENTRALITIES = [15, 25, 35, 45, 55]
DEFAULT_FLBC_VALUES = [0.02, 0.04, 0.06, 0.08, 0.12]

def generate_commands(events: int, bwgen_exec: str, flbc_values, centralities):
    """生成所有 bwgen 命令行"""
    for cent, flbc in itertools.product(centralities, flbc_values):
        yield [
            bwgen_exec,
            "-C", str(cent),
            "-f", f"{flbc:.3f}",
            "-n", str(events),
            "-c", "../../configs/default.yaml"
        ]

def run_command(cmd, no_log):
    cmd_str = " ".join(cmd)
    logger.info(f"Running: {cmd_str}")

    if no_log:
        try:
            subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
            logger.info(f"✅ Finished: {cmd_str}")
        except subprocess.CalledProcessError:
            logger.error(f"❌ Failed: {cmd_str}")
        return

    # 否则保存日志
    cent, flbc = None, None
    for i, arg in enumerate(cmd):
        if arg == "-C" and i + 1 < len(cmd):
            cent = cmd[i + 1]
        if arg == "-f" and i + 1 < len(cmd):
            flbc = cmd[i + 1]

    if cent and flbc:
        log_name = f"results_cent{cent}_fLBC{float(flbc):.3f}.log"
    else:
        log_name = "unknown.log"

    log_dir = Path(__file__).resolve().parent.parent / "bwgen_log"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / log_name

    with open(log_path, "w") as logfile:
        try:
            subprocess.run(cmd, stdout=logfile, stderr=subprocess.STDOUT, check=True)
            logger.info(f"✅ Finished: {cmd_str}")
        except subprocess.CalledProcessError:
            logger.error(f"❌ Failed: {cmd_str} (see {log_path.name})")

def main():
    parser = argparse.ArgumentParser(description="并行批量运行 bwgen")
    parser.add_argument("--dry", action="store_true", help="仅打印命令，不执行")
    parser.add_argument("--events", type=int, default=50000, help="每次运行的事件数")
    parser.add_argument("--bwgen", default="../../build/bwgen", help="bwgen 可执行路径")
    parser.add_argument(
        "--max-workers",
        type=int,
        default=max(1, (os.cpu_count() or 1) - 1),
        help="最大并发数（默认：CPU 核心数 - 1）"
    )
    parser.add_argument("--flbc-list", type=float, nargs="+", help="自定义 fLBC 值列表，例如: 0.03 0.05 0.07")
    parser.add_argument("--cent-list", type=int, nargs="+", help="自定义中心度列表，例如: 15 35 55")
    parser.add_argument("--no-log", action="store_true", help="不保存日志，子进程输出将被丢弃")

    args = parser.parse_args()
    flbc_values = args.flbc_list if args.flbc_list else DEFAULT_FLBC_VALUES
    centralities = args.cent_list if args.cent_list else DEFAULT_CENTRALITIES

    commands = list(generate_commands(args.events, args.bwgen, flbc_values, centralities))

    if args.dry:
        for cmd in commands:
            logger.info("$ " + " ".join(cmd))
        return

    logger.info(f"🧵 Launching with up to {args.max_workers} parallel workers...")

    with ProcessPoolExecutor(max_workers=args.max_workers) as executor:
        futures = [executor.submit(run_command, cmd, args.no_log) for cmd in commands]
        for _ in as_completed(futures):
            pass

    logger.info("✅ All jobs completed.")

if __name__ == "__main__":
    main()
