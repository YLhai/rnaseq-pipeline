"""Small helpers for command execution and logging."""

from __future__ import annotations

import logging
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List, Optional, Sequence


def setup_logging(level: str = "INFO", log_file: Optional[Path] = None) -> None:
    handlers: List[logging.Handler] = [logging.StreamHandler(sys.stderr)]
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file, encoding="utf-8"))
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=handlers,
        force=True,
    )


def which(cmd: str) -> Optional[str]:
    return shutil.which(cmd)


def command_available(cmd: str) -> bool:
    return which(cmd) is not None


def _format_cmd(cmd: Sequence[str]) -> str:
    return " ".join(str(c) for c in cmd)


def run_cmd(
    cmd: Sequence[str],
    log_file: Optional[Path] = None,
    cwd: Optional[Path] = None,
    dry_run: bool = False,
    check: bool = True,
) -> int:
    """Run a command and stream stdout/stderr to a log file."""
    logging.info("RUN: %s", _format_cmd(cmd))
    if dry_run:
        return 0
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        with log_file.open("a", encoding="utf-8") as fh:
            fh.write(f"\n$ {_format_cmd(cmd)}\n")
            fh.flush()
            proc = subprocess.run(
                [str(c) for c in cmd],
                stdout=fh,
                stderr=subprocess.STDOUT,
                cwd=None if cwd is None else str(cwd),
                check=False,
            )
    else:
        proc = subprocess.run(
            [str(c) for c in cmd],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            cwd=None if cwd is None else str(cwd),
            text=True,
            check=False,
        )
        if proc.stdout:
            logging.debug(proc.stdout)
    if check and proc.returncode != 0:
        raise RuntimeError(
            f"command failed with exit code {proc.returncode}: {_format_cmd(cmd)}"
        )
    return proc.returncode


def run_pipe(
    cmds: Sequence[Sequence[str]],
    log_file: Optional[Path] = None,
    dry_run: bool = False,
) -> int:
    """Run a pipeline, e.g. hisat2 | samtools sort."""
    logging.info("PIPE: %s", " | ".join(_format_cmd(c) for c in cmds))
    if dry_run:
        return 0
    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        log_handle = log_file.open("ab")
    else:
        log_handle = subprocess.DEVNULL
    procs: List[subprocess.Popen] = []
    prev_stdout = None
    try:
        for i, cmd in enumerate(cmds):
            last = i == len(cmds) - 1
            stdout = log_handle if last else subprocess.PIPE
            proc = subprocess.Popen(
                [str(c) for c in cmd],
                stdin=prev_stdout,
                stdout=stdout,
                stderr=log_handle,
            )
            if prev_stdout is not None:
                prev_stdout.close()
            prev_stdout = proc.stdout
            procs.append(proc)
        codes = [p.wait() for p in procs]
    finally:
        if log_file is not None:
            log_handle.close()
    for code in codes:
        if code != 0:
            raise RuntimeError(f"pipeline failed with exit code {code}")
    return 0


def write_tsv(path: Path, header: Iterable[str], rows: Iterable[Iterable[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join("" if v is None else str(v) for v in row) + "\n")
