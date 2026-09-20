"""Configuration loading and validation."""

from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml


class ConfigError(ValueError):
    """Raised when the configuration is invalid."""


def _as_path(base: Path, value: Any) -> Optional[Path]:
    if value in (None, ""):
        return None
    raw = str(value)
    # On POSIX (e.g. WSL), accept Windows-style absolute paths and map them
    # to /mnt/<drive>/... so one config can be shared across platforms.
    if os.name != "nt":
        match = re.match(r"^([A-Za-z]):[\\/](.*)$", raw)
        if match:
            rest = match.group(2).replace("\\", "/")
            raw = f"/mnt/{match.group(1).lower()}/{rest}"
    p = Path(raw).expanduser()
    if not p.is_absolute():
        p = base / p
    return p


def _as_list(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [v.strip() for v in value.split(",") if v.strip()]
    return [str(v).strip() for v in value if str(v).strip()]


@dataclass
class Sample:
    sample: str
    group: str
    replicate: str
    fastq_1: Path
    fastq_2: Optional[Path] = None
    extra: Dict[str, str] = field(default_factory=dict)

    @property
    def paired(self) -> bool:
        return self.fastq_2 is not None


@dataclass
class Config:
    path: Path
    base_dir: Path
    raw: Dict[str, Any]
    samples: List[Sample] = field(default_factory=list)

    # project
    project_name: str = "rnaseq_project"
    outdir: Path = Path("results")
    threads: int = 8
    jobs: int = 2
    force: bool = False
    keep_intermediate: bool = True
    log_level: str = "INFO"
    input_is_clean: bool = False

    # reference
    genome_fasta: Optional[Path] = None
    gtf: Optional[Path] = None
    hisat2_index: Optional[Path] = None
    build_index: bool = True
    annotation_format: str = "GTF"
    feature_type: str = "exon"
    gene_attribute: str = "gene_id"

    # annotation
    annotation_table: Optional[Path] = None
    annotation_delimiter: str = "\t"
    annotation_from_gtf: bool = False
    gene_id_column: str = "GeneID"
    product_column: str = "Product"
    go_column: str = "GO Terms"
    kegg_gene_table: Optional[Path] = None
    go_gene_table: Optional[Path] = None
    ko_pathway_table: Optional[Path] = None
    pathway_names: Optional[Path] = None

    # analysis
    design_variable: str = "group"
    covariates: List[str] = field(default_factory=list)
    contrasts: List[Dict[str, str]] = field(default_factory=list)
    padj_cutoff: float = 0.05
    log2fc_cutoff: float = 1.0
    run_stringtie: bool = False
    run_enrichment: bool = True
    strandness: int = 0
    figure_formats: List[str] = field(default_factory=lambda: ["png", "pdf"])
    figure_enrichment_threshold: str = "strict"
    categories: Dict[str, Dict[str, Any]] = field(default_factory=dict)

    @property
    def clean_dir(self) -> Path:
        return self.outdir / "clean"

    @property
    def qc_dir(self) -> Path:
        return self.outdir / "qc"

    @property
    def bam_dir(self) -> Path:
        return self.outdir / "bam"

    @property
    def counts_dir(self) -> Path:
        return self.outdir / "counts"

    @property
    def de_dir(self) -> Path:
        return self.outdir / "DE"

    @property
    def annotated_dir(self) -> Path:
        return self.outdir / "DE_annotated"

    @property
    def category_dir(self) -> Path:
        return self.outdir / "category_analysis"

    @property
    def enrichment_dir(self) -> Path:
        return self.outdir / "enrichment"

    @property
    def logs_dir(self) -> Path:
        return self.outdir / "logs"

    @property
    def stringtie_dir(self) -> Path:
        return self.outdir / "stringtie"

    @property
    def index_prefix(self) -> Path:
        if self.hisat2_index is not None:
            return self.hisat2_index
        return self.outdir / "reference" / "hisat2_index" / "genome"

    @property
    def featurecounts_file(self) -> Path:
        return self.counts_dir / "featureCounts_counts.txt"

    @property
    def metadata_file(self) -> Path:
        return self.counts_dir / "metadata.tsv"

    @property
    def contrast_file(self) -> Path:
        return self.counts_dir / "contrasts.tsv"

    def ensure_dirs(self) -> None:
        for p in (
            self.outdir, self.clean_dir, self.qc_dir, self.bam_dir,
            self.counts_dir, self.de_dir, self.annotated_dir,
            self.category_dir, self.enrichment_dir, self.logs_dir,
        ):
            p.mkdir(parents=True, exist_ok=True)


def load_samples(path: Path, delimiter: str = "\t") -> List[Sample]:
    if not path.exists():
        raise ConfigError(f"samples file not found: {path}")
    samples: List[Sample] = []
    with path.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if not reader.fieldnames:
            raise ConfigError(f"samples file has no header: {path}")
        cols = {c.strip().lower(): c for c in reader.fieldnames}
        required = ("sample", "group", "fastq_1")
        missing = [c for c in required if c not in cols]
        if missing:
            raise ConfigError(
                f"samples file missing required columns {missing}; "
                "expected: sample, group, replicate, fastq_1, fastq_2"
            )
        for row in reader:
            sample = (row.get(cols["sample"]) or "").strip()
            if not sample:
                continue
            fq1 = _as_path(path.parent, row.get(cols["fastq_1"]))
            fq2_col = cols.get("fastq_2")
            fq2 = _as_path(path.parent, row.get(fq2_col)) if fq2_col else None
            rep_col = cols.get("replicate")
            rep = (row.get(rep_col) or "").strip() if rep_col else ""
            if fq1 is None:
                raise ConfigError(f"sample {sample}: fastq_1 is empty")
            samples.append(
                Sample(
                    sample=sample,
                    group=(row.get(cols["group"]) or "").strip(),
                    replicate=rep,
                    fastq_1=fq1,
                    fastq_2=fq2,
                    extra={
                        col: (row.get(col) or "").strip()
                        for name, col in cols.items()
                        if name not in {"sample", "group", "replicate", "fastq_1", "fastq_2"}
                    },
                )
            )
    if not samples:
        raise ConfigError(f"no samples found in {path}")
    names = [s.sample for s in samples]
    if len(names) != len(set(names)):
        raise ConfigError("sample names must be unique")
    return samples


def load_config(path: str | Path) -> Config:
    cfg_path = Path(path).expanduser().resolve()
    if not cfg_path.exists():
        raise ConfigError(f"config file not found: {cfg_path}")
    with cfg_path.open(encoding="utf-8") as fh:
        try:
            raw = yaml.safe_load(fh) or {}
        except yaml.YAMLError as exc:
            mark = getattr(exc, "problem_mark", None)
            location = (
                f"line {mark.line + 1}, column {mark.column + 1}"
                if mark is not None else "unknown location"
            )
            problem = getattr(exc, "problem", str(exc))
            raise ConfigError(
                f"invalid YAML in {cfg_path} ({location}): {problem}\n"
                "Common causes:\n"
                "  - 'project: value' followed by indented keys; use 'project:' on its own line\n"
                "  - an unquoted value containing a colon followed by a space, e.g. 'name: my: project'\n"
                "  - a full-width Chinese colon character instead of ':'\n"
                "  - tab indentation; YAML requires spaces"
            ) from exc
    base = cfg_path.parent

    project = raw.get("project", {}) or {}
    samples_value = raw.get("samples", "samples.tsv")
    samples_path = _as_path(base, samples_value)
    samples_delimiter = str(raw.get("samples_delimiter", "\t"))
    samples = load_samples(samples_path, delimiter=samples_delimiter)

    reference = raw.get("reference", {}) or {}
    annotation = raw.get("annotation", {}) or {}
    analysis = raw.get("analysis", {}) or {}
    deg = analysis.get("deg", {}) or {}
    categories = analysis.get("categories", {}) or {}
    figures = analysis.get("figures", {}) or {}

    cfg = Config(
        path=cfg_path,
        base_dir=base,
        raw=raw,
        samples=samples,
        project_name=str(project.get("name", "rnaseq_project")),
        outdir=_as_path(base, project.get("outdir", "results")),
        threads=int(project.get("threads", 8)),
        jobs=int(project.get("jobs", 2)),
        force=bool(project.get("force", False)),
        keep_intermediate=bool(project.get("keep_intermediate", True)),
        log_level=str(project.get("log_level", "INFO")),
        input_is_clean=bool(project.get(
            "input_is_clean",
            project.get("skip_fastp", False),
        )),
        genome_fasta=_as_path(base, reference.get("genome_fasta")),
        gtf=_as_path(base, reference.get("gtf")),
        hisat2_index=_as_path(base, reference.get("hisat2_index")),
        build_index=bool(reference.get("build_index", True)),
        annotation_format=str(reference.get("annotation_format", "GTF")),
        feature_type=str(reference.get("feature_type", "exon")),
        gene_attribute=str(reference.get("gene_attribute", "gene_id")),
        annotation_table=_as_path(base, annotation.get("table")),
        annotation_delimiter=str(annotation.get("delimiter", "\t")),
        annotation_from_gtf=bool(annotation.get("from_gtf", False)),
        gene_id_column=str(annotation.get("gene_id_column", "GeneID")),
        product_column=str(annotation.get("product_column", "Product")),
        go_column=str(annotation.get("go_column", "GO Terms")),
        kegg_gene_table=_as_path(base, annotation.get("kegg_gene_table")),
        go_gene_table=_as_path(base, annotation.get("go_gene_table")),
        ko_pathway_table=_as_path(base, annotation.get("ko_pathway_table")),
        pathway_names=_as_path(base, annotation.get("pathway_names")),
        design_variable=str(analysis.get("design_variable", "group")),
        covariates=_as_list(analysis.get("covariates", [])),
        contrasts=list(analysis.get("contrasts", []) or []),
        padj_cutoff=float(deg.get("padj", 0.05)),
        log2fc_cutoff=float(deg.get("log2fc", 1.0)),
        run_stringtie=bool(analysis.get("stringtie", False)),
        run_enrichment=bool(analysis.get("enrichment", True)),
        strandness=int(analysis.get("strandness", 0)),
        figure_formats=_as_list(figures.get("formats", ["png", "pdf"])) or ["png", "pdf"],
        figure_enrichment_threshold=str(
            figures.get("enrichment_threshold", "strict")
        ),
        categories=dict(categories),
    )
    # Normalize design variables to the actual sample-table column names,
    # so that users can write either "Batch" or "batch" in the config.
    sample_cols = {c.lower(): c for c in _sample_columns(cfg)}
    cfg.design_variable = sample_cols.get(cfg.design_variable.lower(), cfg.design_variable)
    cfg.covariates = [sample_cols.get(c.lower(), c) for c in cfg.covariates]
    cfg.ensure_dirs()
    validate_config(cfg, check_files=False)
    return cfg


def validate_config(cfg: Config, check_files: bool = True) -> List[str]:
    problems: List[str] = []
    if cfg.genome_fasta is None and cfg.hisat2_index is None:
        problems.append("reference.genome_fasta or reference.hisat2_index is required")
    if cfg.gtf is None:
        problems.append("reference.gtf is required")
    if not cfg.contrasts:
        problems.append("analysis.contrasts must contain at least one contrast")
    groups = {s.group for s in cfg.samples}
    for i, c in enumerate(cfg.contrasts):
        if not c.get("name") or not c.get("numerator") or not c.get("denominator"):
            problems.append(f"analysis.contrasts[{i}] needs name/numerator/denominator")
            continue
        for key in ("numerator", "denominator"):
            if c[key] not in groups:
                problems.append(
                    f"contrast {c.get('name')}: group '{c[key]}' not found in samples.tsv"
                )
    if cfg.design_variable not in _sample_columns(cfg):
        lower_map = {c.lower(): c for c in _sample_columns(cfg)}
        if cfg.design_variable.lower() not in lower_map:
            problems.append(
                f"analysis.design_variable '{cfg.design_variable}' is not available in samples"
            )
        else:
            cfg.design_variable = lower_map[cfg.design_variable.lower()]
    for cov in cfg.covariates:
        lower_map = {c.lower(): c for c in _sample_columns(cfg)}
        if cov not in _sample_columns(cfg) and cov.lower() not in lower_map:
            problems.append(f"covariate '{cov}' is not available in samples")
    if check_files:
        for label, p in (
            ("genome_fasta", cfg.genome_fasta),
            ("gtf", cfg.gtf),
            ("annotation.table", cfg.annotation_table),
            ("annotation.kegg_gene_table", cfg.kegg_gene_table),
            ("annotation.go_gene_table", cfg.go_gene_table),
            ("annotation.ko_pathway_table", cfg.ko_pathway_table),
            ("annotation.pathway_names", cfg.pathway_names),
        ):
            if p is not None and not p.exists():
                problems.append(f"{label} not found: {p}")
        for s in cfg.samples:
            if not s.fastq_1.exists():
                problems.append(f"{s.sample}: fastq_1 not found: {s.fastq_1}")
            if s.fastq_2 is not None and not s.fastq_2.exists():
                problems.append(f"{s.sample}: fastq_2 not found: {s.fastq_2}")
        if cfg.hisat2_index is not None and not Path(str(cfg.hisat2_index) + ".1.ht2").exists():
            problems.append(
                f"hisat2 index prefix seems incomplete: {cfg.hisat2_index}.1.ht2 not found"
            )
    if problems:
        raise ConfigError("Configuration problems:\n- " + "\n- ".join(problems))
    return problems


def _sample_columns(cfg: Config) -> List[str]:
    cols = ["sample", "group", "replicate"]
    if cfg.samples:
        for key in cfg.samples[0].extra:
            if key not in cols:
                cols.append(key)
            # Allow case-insensitive matching against header names.
            if key.lower() not in [c.lower() for c in cols]:
                cols.append(key.lower())
    return cols
