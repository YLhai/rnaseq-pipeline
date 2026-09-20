"""Command line interface."""

from __future__ import annotations

import argparse
import logging
import shutil
import sys
from importlib import resources
from pathlib import Path

from . import __version__
from .config import ConfigError, load_config
from .steps import ALL_STEPS, Pipeline
from .utils import command_available, setup_logging


def _template(name: str) -> str:
    return resources.files("rnaseq_pipeline").joinpath("templates", name).read_text(
        encoding="utf-8"
    )


def cmd_init(args: argparse.Namespace) -> int:
    outdir = Path(args.output).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    config_path = outdir / "config.yaml"
    samples_path = outdir / "samples.tsv"
    if config_path.exists() and not args.force:
        print(f"config already exists: {config_path} (use --force to overwrite)")
    else:
        config_path.write_text(_template("config.example.yaml"), encoding="utf-8")
        print(f"wrote {config_path}")
    if samples_path.exists() and not args.force:
        print(f"samples already exists: {samples_path} (use --force to overwrite)")
    else:
        samples_path.write_text(_template("samples.example.tsv"), encoding="utf-8")
        print(f"wrote {samples_path}")
    print("\nNext: edit samples.tsv and config.yaml, then run:")
    print(f"  rnaseq-pipeline check -c {config_path}")
    print(f"  rnaseq-pipeline run -c {config_path}")
    return 0


def cmd_check(args: argparse.Namespace) -> int:
    try:
        cfg = load_config(args.config)
    except ConfigError as exc:
        print(f"CONFIG ERROR:\n{exc}", file=sys.stderr)
        return 2
    print(f"project        : {cfg.project_name}")
    print(f"config         : {cfg.path}")
    print(f"samples        : {len(cfg.samples)}")
    print(f"groups         : {sorted({s.group for s in cfg.samples})}")
    print(f"outdir         : {cfg.outdir}")
    print(f"genome         : {cfg.genome_fasta}")
    print(f"GTF            : {cfg.gtf}")
    print(f"contrasts      : {', '.join(c['name'] for c in cfg.contrasts)}")
    print(f"categories     : {', '.join(cfg.categories) if cfg.categories else '(none)'}")
    print("\nExternal tools:")
    for tool in ("fastp", "hisat2", "hisat2-build", "samtools",
                 "featureCounts", "stringtie", "Rscript", "multiqc"):
        print(f"  {tool:<14} {'OK' if command_available(tool) else 'MISSING'}")
    print("\nConfiguration looks valid (use --strict to also check input files).")
    if args.strict:
        missing = []
        for label, p in (
            ("genome", cfg.genome_fasta), ("GTF", cfg.gtf),
            ("annotation", cfg.annotation_table),
        ):
            if p is not None and not p.exists():
                missing.append(f"{label}: {p}")
        for s in cfg.samples:
            if not s.fastq_1.exists():
                missing.append(f"{s.sample} R1: {s.fastq_1}")
            if s.fastq_2 is not None and not s.fastq_2.exists():
                missing.append(f"{s.sample} R2: {s.fastq_2}")
        if missing:
            print("\nMissing input files:", file=sys.stderr)
            for item in missing:
                print(f"  {item}", file=sys.stderr)
            return 3
        print("All configured input files exist.")
    return 0


def cmd_run(args: argparse.Namespace) -> int:
    try:
        cfg = load_config(args.config)
    except ConfigError as exc:
        print(f"CONFIG ERROR:\n{exc}", file=sys.stderr)
        return 2
    setup_logging(cfg.log_level, cfg.logs_dir / "pipeline.log")
    steps = [s.strip() for s in args.steps.split(",") if s.strip()]
    invalid = [s for s in steps if s not in ALL_STEPS]
    if invalid:
        print(f"unknown steps: {invalid}; valid: {ALL_STEPS}", file=sys.stderr)
        return 2
    pipeline = Pipeline(cfg, dry_run=args.dry_run, force=args.force)
    pipeline.run(steps)
    return 0


def cmd_figures(args: argparse.Namespace) -> int:
    try:
        cfg = load_config(args.config)
    except ConfigError as exc:
        print(f"CONFIG ERROR:\n{exc}", file=sys.stderr)
        return 2
    setup_logging(cfg.log_level, cfg.logs_dir / "pipeline.log")
    if args.formats:
        cfg.figure_formats = [s.strip() for s in args.formats.split(",") if s.strip()]
    pipeline = Pipeline(cfg, dry_run=args.dry_run, force=args.force)
    pipeline.run(["figures"])
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rnaseq-pipeline",
        description="Reusable bulk RNA-seq pipeline (fastp/HISAT2/featureCounts/DESeq2).",
    )
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)

    p_init = sub.add_parser("init", help="create config.yaml and samples.tsv")
    p_init.add_argument("-o", "--output", default="rnaseq_project")
    p_init.add_argument("--force", action="store_true")
    p_init.set_defaults(func=cmd_init)

    p_check = sub.add_parser("check", help="validate configuration and tools")
    p_check.add_argument("-c", "--config", default="config.yaml")
    p_check.add_argument("--strict", action="store_true",
                         help="also check that input files exist")
    p_check.set_defaults(func=cmd_check)

    p_run = sub.add_parser("run", help="run the pipeline")
    p_run.add_argument("-c", "--config", default="config.yaml")
    p_run.add_argument("--steps", default=",".join(ALL_STEPS))
    p_run.add_argument("--dry-run", action="store_true")
    p_run.add_argument("--force", action="store_true")
    p_run.set_defaults(func=cmd_run)

    p_fig = sub.add_parser("figures", help="generate publication figures from existing results")
    p_fig.add_argument("-c", "--config", default="config.yaml")
    p_fig.add_argument("--formats", default="png,pdf",
                       help="comma-separated output formats: png,pdf")
    p_fig.add_argument("--dry-run", action="store_true")
    p_fig.add_argument("--force", action="store_true")
    p_fig.set_defaults(func=cmd_figures)
    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
