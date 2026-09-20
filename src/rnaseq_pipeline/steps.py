"""Pipeline orchestration."""

from __future__ import annotations

import concurrent.futures
import csv
import json
import logging
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

from . import downstream
from .config import Config
from .utils import command_available, run_cmd, run_pipe, which, write_tsv


ALL_STEPS = [
    "fastp", "index", "align", "count", "stringtie",
    "deseq", "annotate", "categories", "enrichment", "figures", "multiqc",
]


class Pipeline:
    def __init__(self, cfg: Config, dry_run: bool = False, force: bool = False):
        self.cfg = cfg
        self.dry_run = dry_run
        self.force = force or cfg.force
        cfg.ensure_dirs()
        self.per_job_threads = max(1, cfg.threads // max(1, cfg.jobs))

    # ------------------------------------------------------------------
    # fastp
    # ------------------------------------------------------------------
    def run_fastp(self) -> None:
        if self.cfg.input_is_clean:
            logging.info("input_is_clean=true; fastp is not needed")
            return
        jobs = []
        for s in self.cfg.samples:
            out1 = self.cfg.clean_dir / f"{s.sample}_R1.fastq.gz"
            out2 = self.cfg.clean_dir / f"{s.sample}_R2.fastq.gz"
            qc_json = self.cfg.qc_dir / f"{s.sample}.json"
            complete = (
                out1.exists() and out1.stat().st_size > 0
                and qc_json.exists() and qc_json.stat().st_size > 0
                and (not s.paired or (out2.exists() and out2.stat().st_size > 0))
            )
            if complete and not self.force:
                logging.info("fastp skip (exists): %s", s.sample)
                continue
            # Remove partial outputs left by a previous failed run.
            if not self.dry_run:
                for f in (out1, out2):
                    if f.exists():
                        f.unlink()
            cmd = [
                "fastp",
                "-i", str(s.fastq_1),
                "-o", str(out1),
                "-q", "20", "-u", "30", "-n", "5",
                "-w", str(self.per_job_threads),
                "-h", str(self.cfg.qc_dir / f"{s.sample}.html"),
                "-j", str(self.cfg.qc_dir / f"{s.sample}.json"),
            ]
            if s.paired:
                cmd += [
                    "-I", str(s.fastq_2),
                    "-O", str(out2),
                    "--detect_adapter_for_pe",
                ]
            jobs.append((s.sample, cmd))
        self._run_jobs(jobs, self.cfg.logs_dir / "fastp")

    # ------------------------------------------------------------------
    # HISAT2 index
    # ------------------------------------------------------------------
    def build_index(self) -> None:
        prefix = self.cfg.index_prefix
        if Path(str(prefix) + ".1.ht2").exists() and not self.force:
            logging.info("HISAT2 index exists: %s", prefix)
            return
        if self.cfg.hisat2_index is not None and not self.cfg.build_index:
            logging.info("using existing HISAT2 index: %s", prefix)
            return
        if self.cfg.genome_fasta is None:
            raise ValueError("reference.genome_fasta is required to build a HISAT2 index")
        prefix.parent.mkdir(parents=True, exist_ok=True)
        cmd = [
            "hisat2-build", "-p", str(self.cfg.threads),
            str(self.cfg.genome_fasta), str(prefix),
        ]
        run_cmd(cmd, self.cfg.logs_dir / "hisat2_build.log", dry_run=self.dry_run)

    # ------------------------------------------------------------------
    # alignment
    # ------------------------------------------------------------------
    def run_align(self) -> None:
        jobs = []
        prefix = self.cfg.index_prefix
        for s in self.cfg.samples:
            bam = self.cfg.bam_dir / f"{s.sample}.sorted.bam"
            if bam.exists() and not self.force:
                logging.info("align skip (exists): %s", s.sample)
                continue
            hisat = [
                "hisat2", "-p", str(self.per_job_threads), "--dta",
                "-x", str(prefix),
            ]
            in1, in2 = self._input_paths(s)
            if not self.dry_run:
                if not in1.exists():
                    raise FileNotFoundError(
                        f"{s.sample}: input FASTQ not found: {in1}"
                    )
                if s.paired and (in2 is None or not in2.exists()):
                    raise FileNotFoundError(
                        f"{s.sample}: input FASTQ not found: {in2}"
                    )
            if s.paired:
                hisat += ["-1", str(in1), "-2", str(in2)]
            else:
                hisat += ["-U", str(in1)]
            sort = [
                "samtools", "sort", "-@", str(self.per_job_threads),
                "-o", str(bam),
            ]
            jobs.append((s.sample, hisat, sort, bam))
        if self.dry_run:
            for sample, hisat, sort, bam in jobs:
                run_pipe([hisat, sort], self.cfg.logs_dir / f"{sample}.hisat2.log", True)
            return
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.cfg.jobs) as pool:
            futures = {}
            for sample, hisat, sort, bam in jobs:
                futures[pool.submit(
                    run_pipe, [hisat, sort],
                    self.cfg.logs_dir / f"{sample}.hisat2.log", False,
                )] = sample
            errors = []
            for future in concurrent.futures.as_completed(futures):
                sample = futures[future]
                try:
                    future.result()
                except Exception as exc:  # noqa: BLE001
                    errors.append((sample, exc))
            if errors:
                sample, exc = errors[0]
                raise RuntimeError(
                    f"alignment '{sample}' failed: {exc}\n"
                    f"See log: {self.cfg.logs_dir / (sample + '.hisat2.log')}"
                ) from exc
        for s in self.cfg.samples:
            bam = self.cfg.bam_dir / f"{s.sample}.sorted.bam"
            run_cmd(
                ["samtools", "index", str(bam)],
                self.cfg.logs_dir / "samtools_index.log",
                dry_run=self.dry_run,
            )
            if not self.cfg.keep_intermediate:
                for f in (
                    self.cfg.clean_dir / f"{s.sample}_R1.fastq.gz",
                    self.cfg.clean_dir / f"{s.sample}_R2.fastq.gz",
                ):
                    if f.exists():
                        f.unlink()

    # ------------------------------------------------------------------
    # featureCounts
    # ------------------------------------------------------------------
    def run_featurecounts(self) -> None:
        bams = [self.cfg.bam_dir / f"{s.sample}.sorted.bam" for s in self.cfg.samples]
        missing = [str(b) for b in bams if not b.exists()]
        if missing and not self.dry_run:
            raise FileNotFoundError(f"missing BAM files: {missing[:3]}")
        if self.cfg.featurecounts_file.exists() and not self.force:
            logging.info("featureCounts skip (exists): %s", self.cfg.featurecounts_file)
            return
        if self.cfg.gtf is None:
            raise ValueError("reference.gtf is required")
        paired = {s.paired for s in self.cfg.samples}
        if len(paired) > 1:
            raise ValueError("featureCounts requires all samples to be paired-end or all single-end")
        cmd = [
            "featureCounts", "-T", str(self.cfg.threads),
            "-F", str(self.cfg.annotation_format),
            "-t", str(self.cfg.feature_type),
            "-g", str(self.cfg.gene_attribute),
            "-s", str(self.cfg.strandness),
        ]
        if paired == {True}:
            cmd += ["-p"]
        cmd += [
            "-a", str(self.cfg.gtf),
            "-o", str(self.cfg.featurecounts_file),
        ] + [str(b) for b in bams]
        run_cmd(cmd, self.cfg.logs_dir / "featurecounts.log", dry_run=self.dry_run)
        if not self.dry_run:
            self._normalize_featurecounts_header()
        self._write_metadata()

    # ------------------------------------------------------------------
    # StringTie (optional)
    # ------------------------------------------------------------------
    def run_stringtie(self) -> None:
        if not self.cfg.run_stringtie:
            logging.info("StringTie disabled (analysis.stringtie=false)")
            return
        if self.cfg.gtf is None:
            raise ValueError("reference.gtf is required")
        jobs = []
        for s in self.cfg.samples:
            bam = self.cfg.bam_dir / f"{s.sample}.sorted.bam"
            outdir = self.cfg.stringtie_dir / s.sample
            outdir.mkdir(parents=True, exist_ok=True)
            gtf_out = outdir / f"{s.sample}.gtf"
            if gtf_out.exists() and not self.force:
                continue
            cmd = [
                "stringtie", "-p", str(self.per_job_threads), "-e", "-B",
                "-G", str(self.cfg.gtf), "-o", str(gtf_out), str(bam),
            ]
            jobs.append((s.sample, cmd))
        self._run_jobs(jobs, self.cfg.logs_dir / "stringtie")

    # ------------------------------------------------------------------
    # DESeq2
    # ------------------------------------------------------------------
    def run_deseq2(self) -> None:
        if not self.cfg.featurecounts_file.exists() and not self.dry_run:
            raise FileNotFoundError(self.cfg.featurecounts_file)
        if not self.dry_run:
            self._write_metadata()
        rscript = which("Rscript")
        if rscript is None:
            raise RuntimeError("Rscript not found in PATH")
        r_script = Path(__file__).resolve().parent / "R" / "run_deseq2.R"
        cmd = [
            rscript, str(r_script),
            "--counts", str(self.cfg.featurecounts_file),
            "--metadata", str(self.cfg.metadata_file),
            "--contrasts", str(self.cfg.contrast_file),
            "--outdir", str(self.cfg.de_dir),
            "--design-var", self.cfg.design_variable,
            "--covariates", ",".join(self.cfg.covariates),
            "--padj", str(self.cfg.padj_cutoff),
            "--log2fc", str(self.cfg.log2fc_cutoff),
        ]
        run_cmd(cmd, self.cfg.logs_dir / "deseq2.log", dry_run=self.dry_run)

    # ------------------------------------------------------------------
    # annotation / categories / enrichment
    # ------------------------------------------------------------------
    def run_annotation(self) -> None:
        self._prepare_annotation()
        downstream.annotate_degs(self.cfg, dry_run=self.dry_run)

    def run_categories(self) -> None:
        self._prepare_annotation()
        downstream.extract_categories(self.cfg, dry_run=self.dry_run)

    def run_enrichment(self) -> None:
        if not self.cfg.run_enrichment:
            logging.info("enrichment disabled (analysis.enrichment=false)")
            return
        rscript = which("Rscript")
        if rscript is None:
            raise RuntimeError("Rscript not found in PATH")
        r_script = Path(__file__).resolve().parent / "R" / "run_enrichment.R"
        self._prepare_annotation()
        go_map_path = self.cfg.go_gene_table
        if go_map_path is None and self.cfg.annotation_table is not None:
            if self.dry_run:
                go_map_path = self.cfg.counts_dir / "go_gene_mapping.tsv"
            else:
                go_map_path = self._prepare_go_mapping()
        payload = {
            "deg_files": {
                c["name"]: str(self.cfg.de_dir / f"{c['name']}_DEGs.csv")
                for c in self.cfg.contrasts
            },
            "background": str(self.cfg.counts_dir / "background_genes.txt"),
            "go_map": str(go_map_path) if go_map_path else "",
            "kegg_gene": str(self.cfg.kegg_gene_table) if self.cfg.kegg_gene_table else "",
            "ko_pathway": str(self.cfg.ko_pathway_table) if self.cfg.ko_pathway_table else "",
            "pathway_names": str(self.cfg.pathway_names) if self.cfg.pathway_names else "",
            "outdir": str(self.cfg.enrichment_dir),
            "padj_method": "BH",
            "pvalue_cutoff": 1.0,
            "qvalue_cutoff": 1.0,
            "min_gene_set_size": 3,
            "max_gene_set_size": 500,
        }
        if not go_map_path and not self.cfg.kegg_gene_table:
            logging.info("enrichment skip: no go_gene_table or kegg_gene_table configured")
            return
        if not self.dry_run:
            background = self._read_background_genes()
            (self.cfg.counts_dir / "background_genes.txt").write_text(
                "\n".join(background) + "\n", encoding="utf-8"
            )
            payload_file = self.cfg.counts_dir / "enrichment_config.json"
            payload_file.write_text(json.dumps(payload, ensure_ascii=False, indent=2),
                                    encoding="utf-8")
        else:
            payload_file = self.cfg.counts_dir / "enrichment_config.json"
        run_cmd(
            [rscript, str(r_script), "--config", str(payload_file)],
            self.cfg.logs_dir / "enrichment.log",
            dry_run=self.dry_run,
        )

    # ------------------------------------------------------------------
    # MultiQC
    # ------------------------------------------------------------------
    def run_multiqc(self) -> None:
        if not command_available("multiqc"):
            logging.info("MultiQC not found; skipping")
            return
        outdir = self.cfg.qc_dir / "multiqc"
        cmd = ["multiqc", "--force", "-o", str(outdir), str(self.cfg.outdir)]
        run_cmd(cmd, self.cfg.logs_dir / "multiqc.log", dry_run=self.dry_run)

    # ------------------------------------------------------------------
    # publication figures
    # ------------------------------------------------------------------
    def run_figures(self) -> None:
        rscript = which("Rscript")
        if rscript is None:
            raise RuntimeError("Rscript not found in PATH")
        r_script = Path(__file__).resolve().parent / "R" / "make_figures.R"
        figures_dir = self.cfg.outdir / "figures"
        payload = {
            "de_dir": str(self.cfg.de_dir),
            "category_dir": str(self.cfg.category_dir),
            "enrichment_dir": str(self.cfg.enrichment_dir),
            "counts_file": str(self.cfg.featurecounts_file),
            "metadata_file": str(self.cfg.metadata_file),
            "outdir": str(figures_dir),
            "contrasts": self.cfg.contrasts,
            "padj": self.cfg.padj_cutoff,
            "log2fc": self.cfg.log2fc_cutoff,
            "formats": self.cfg.figure_formats,
            "enrichment_threshold": self.cfg.figure_enrichment_threshold,
        }
        payload_file = self.cfg.counts_dir / "figures_config.json"
        if not self.dry_run:
            payload_file.write_text(
                json.dumps(payload, ensure_ascii=False, indent=2),
                encoding="utf-8",
            )
        run_cmd(
            [rscript, str(r_script), "--config", str(payload_file)],
            self.cfg.logs_dir / "figures.log",
            dry_run=self.dry_run,
        )

    # ------------------------------------------------------------------
    # run / summary
    # ------------------------------------------------------------------
    def run(self, steps: Optional[Iterable[str]] = None) -> None:
        steps = list(steps or ALL_STEPS)
        if self.cfg.input_is_clean and "fastp" in steps:
            logging.info("fastp skipped because project.input_is_clean=true")
            steps = [s for s in steps if s != "fastp"]
        dispatch = {
            "fastp": self.run_fastp,
            "index": self.build_index,
            "align": self.run_align,
            "count": self.run_featurecounts,
            "stringtie": self.run_stringtie,
            "deseq": self.run_deseq2,
            "annotate": self.run_annotation,
            "categories": self.run_categories,
            "enrichment": self.run_enrichment,
            "figures": self.run_figures,
            "multiqc": self.run_multiqc,
        }
        for step in steps:
            if step not in ALL_STEPS:
                raise ValueError(f"unknown step: {step}")
            logging.info("=== step: %s ===", step)
            dispatch[step]()
        self.write_summary()

    def write_summary(self) -> None:
        path = self.cfg.outdir / "pipeline_summary.md"
        lines = [
            f"# {self.cfg.project_name}",
            "",
            "## Samples",
            "",
            "| sample | group | replicate | R1 | R2 |",
            "|---|---|---|---|---|",
        ]
        for s in self.cfg.samples:
            lines.append(
                f"| {s.sample} | {s.group} | {s.replicate} | {s.fastq_1} | "
                f"{s.fastq_2 or ''} |"
            )
        deg_summary = self.cfg.de_dir / "DEG_summary.csv"
        if deg_summary.exists():
            lines += ["", "## DEG summary", "", deg_summary.read_text(encoding="utf-8")]
        lines += [
            "",
            "## Main outputs",
            "",
            f"- counts: `{self.cfg.featurecounts_file}`",
            f"- DE results: `{self.cfg.de_dir}`",
            f"- annotated DEGs: `{self.cfg.annotated_dir}`",
            f"- category tables: `{self.cfg.category_dir}`",
            f"- enrichment: `{self.cfg.enrichment_dir}`",
            "",
        ]
        if not self.dry_run:
            path.write_text("\n".join(lines), encoding="utf-8")
            logging.info("summary written: %s", path)

    # ------------------------------------------------------------------
    # internal helpers
    # ------------------------------------------------------------------
    def _run_jobs(self, jobs, log_dir: Path) -> None:
        if self.dry_run:
            for name, cmd in jobs:
                run_cmd(cmd, log_dir / f"{name}.log", dry_run=True)
            return
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.cfg.jobs) as pool:
            futures = {
                pool.submit(run_cmd, cmd, log_dir / f"{name}.log", None, False): name
                for name, cmd in jobs
            }
            errors = []
            for future in concurrent.futures.as_completed(futures):
                name = futures[future]
                try:
                    future.result()
                except Exception as exc:  # noqa: BLE001
                    errors.append((name, exc))
            if errors:
                name, exc = errors[0]
                raise RuntimeError(
                    f"job '{name}' failed: {exc}\n"
                    f"See log: {log_dir / (name + '.log')}"
                ) from exc

    def _write_metadata(self) -> None:
        extra_cols: List[str] = []
        for s in self.cfg.samples:
            for key in s.extra:
                if key not in extra_cols:
                    extra_cols.append(key)
        header = ["sample", "group", "replicate"] + extra_cols
        rows = [
            [s.sample, s.group, s.replicate] + [s.extra.get(c, "") for c in extra_cols]
            for s in self.cfg.samples
        ]
        write_tsv(self.cfg.metadata_file, header, rows)
        write_tsv(
            self.cfg.contrast_file,
            ["name", "numerator", "denominator"],
            [[c["name"], c["numerator"], c["denominator"]] for c in self.cfg.contrasts],
        )

    def _normalize_featurecounts_header(self) -> None:
        """Replace BAM-path column names with plain sample names."""
        path = self.cfg.featurecounts_file
        if not path.exists():
            return
        lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
        if len(lines) < 2:
            return
        header = lines[1].rstrip("\n").split("\t")
        n_samples = len(self.cfg.samples)
        if len(header) < 6 + n_samples:
            logging.warning(
                "featureCounts header has %d columns; expected at least %d",
                len(header), 6 + n_samples,
            )
            return
        header[6:6 + n_samples] = [s.sample for s in self.cfg.samples]
        lines[1] = "\t".join(header) + "\n"
        path.write_text("".join(lines), encoding="utf-8")
        logging.info("normalized featureCounts sample columns: %s",
                     ", ".join(s.sample for s in self.cfg.samples))

    def _read_background_genes(self) -> List[str]:
        genes: List[str] = []
        with self.cfg.featurecounts_file.open(encoding="utf-8") as fh:
            first = fh.readline()
            if first.startswith("#"):
                fh.readline()
            for line in fh:
                gene = line.split("\t", 1)[0].strip()
                if gene:
                    genes.append(gene)
        return genes

    def _input_paths(self, sample):
        """Return the FASTQ paths used by alignment."""
        if self.cfg.input_is_clean:
            return sample.fastq_1, sample.fastq_2
        return (
            self.cfg.clean_dir / f"{sample.sample}_R1.fastq.gz",
            self.cfg.clean_dir / f"{sample.sample}_R2.fastq.gz",
        )

    def _prepare_go_mapping(self) -> Optional[Path]:
        """Build a two-column gene->GO mapping from the annotation table."""
        if self.cfg.annotation_table is None or not self.cfg.annotation_table.exists():
            return None
        out = self.cfg.counts_dir / "go_gene_mapping.tsv"
        rows = []
        with self.cfg.annotation_table.open(newline="", encoding="utf-8-sig") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                gene = (row.get(self.cfg.gene_id_column) or "").strip()
                if not gene:
                    continue
                seen = set()
                for go in re.findall(r"GO:\d+", row.get(self.cfg.go_column, "") or ""):
                    if go not in seen:
                        rows.append((gene, go))
                        seen.add(go)
        if not rows:
            return None
        with out.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(["gene", "go"])
            writer.writerows(rows)
        logging.info("generated GO mapping from annotation table: %s", out)
        return out

    def _prepare_annotation(self) -> Optional[Path]:
        """Use an existing annotation table, or derive one from GTF/GFF."""
        if self.cfg.annotation_table is not None and self.cfg.annotation_table.exists():
            return self.cfg.annotation_table
        if not self.cfg.annotation_from_gtf or self.cfg.gtf is None:
            return None
        out = self.cfg.outdir / "annotation" / "annotation_from_gtf.tsv"
        if out.exists() and not self.force:
            self.cfg.annotation_table = out
            return out
        genes = {}
        with self.cfg.gtf.open(encoding="utf-8") as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                attrs = {}
                for item in parts[8].strip().strip(";").split(";"):
                    item = item.strip()
                    if not item:
                        continue
                    if "=" in item:
                        key, value = item.split("=", 1)
                    else:
                        bits = item.split(" ", 1)
                        if len(bits) != 2:
                            continue
                        key, value = bits
                    attrs[key.strip()] = value.strip().strip('"')
                gene = (
                    attrs.get(self.cfg.gene_attribute)
                    or attrs.get("gene_id")
                    or attrs.get("ID")
                    or attrs.get("Parent")
                )
                if not gene:
                    continue
                product = (
                    attrs.get("product")
                    or attrs.get("gene_name")
                    or attrs.get("Name")
                    or attrs.get("description")
                    or attrs.get("gene")
                    or ""
                )
                existing = genes.get(gene, "")
                if product and (not existing or existing == "hypothetical protein"):
                    genes[gene] = product
                else:
                    genes.setdefault(gene, "")
        if not genes:
            logging.warning("could not derive annotation from %s", self.cfg.gtf)
            return None
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow([self.cfg.gene_id_column, self.cfg.product_column,
                             self.cfg.go_column])
            for gene in sorted(genes):
                writer.writerow([gene, genes[gene] or "hypothetical protein", ""])
        self.cfg.annotation_table = out
        logging.info("derived annotation from GTF: %s", out)
        return out
