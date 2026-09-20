"""Annotation of DEGs and custom functional gene-set classification."""

from __future__ import annotations

import csv
import logging
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .config import Config


def _open_csv(path: Path, delimiter: str = ","):
    return path.open(newline="", encoding="utf-8-sig")


def _read_table(path: Path, delimiter: str = "\t") -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


def _split_values(value: str) -> List[str]:
    if not value:
        return []
    out = []
    for part in re.split(r"[;,]+", value):
        part = part.strip()
        if part and part not in out:
            out.append(part)
    return out


def _lookup(row: Dict[str, str], name: str) -> str:
    """Lookup that ignores case, spaces and underscores in column names."""
    if name in row:
        return row.get(name, "") or ""
    target = re.sub(r"[^a-z0-9]", "", str(name).lower())
    for key, value in row.items():
        if re.sub(r"[^a-z0-9]", "", str(key).lower()) == target:
            return value or ""
    return ""


def _read_gene_stats(cfg: Config) -> Dict[str, Dict[str, Dict[str, object]]]:
    """gene -> contrast -> stats."""
    stats: Dict[str, Dict[str, Dict[str, object]]] = defaultdict(dict)
    for contrast in cfg.contrasts:
        name = contrast["name"]
        path = cfg.de_dir / f"{name}_all_results.csv"
        if not path.exists():
            continue
        with _open_csv(path) as fh:
            for row in csv.DictReader(fh):
                gene = (row.get("gene") or "").strip()
                if not gene:
                    continue
                try:
                    lfc = float(row.get("log2FoldChange") or "nan")
                except ValueError:
                    lfc = float("nan")
                try:
                    padj = float(row.get("padj") or "nan")
                except ValueError:
                    padj = float("nan")
                is_deg = (
                    not math.isnan(padj) and not math.isnan(lfc)
                    and padj < cfg.padj_cutoff and abs(lfc) >= cfg.log2fc_cutoff
                )
                stats[gene][name] = {
                    "log2FoldChange": None if math.isnan(lfc) else lfc,
                    "padj": None if math.isnan(padj) else padj,
                    "DEG": is_deg,
                }
    return stats


def _load_annotation_lookup(cfg: Config) -> Dict[str, Dict[str, str]]:
    lookup: Dict[str, Dict[str, str]] = {}
    if cfg.annotation_table is None or not cfg.annotation_table.exists():
        return lookup
    with cfg.annotation_table.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter=cfg.annotation_delimiter)
        for row in reader:
            gene = (_lookup(row, cfg.gene_id_column) or "").strip()
            if not gene:
                continue
            ent = lookup.setdefault(gene, {})
            for key, value in row.items():
                if value is None:
                    continue
                value = value.strip()
                if not value:
                    continue
                if key in ent and value not in ent[key].split(";"):
                    ent[key] = f"{ent[key]};{value}"
                else:
                    ent.setdefault(key, value)
    return lookup


def _load_gene_go(cfg: Config, lookup: Dict[str, Dict[str, str]]) -> Dict[str, List[str]]:
    gene_go: Dict[str, List[str]] = defaultdict(list)
    for gene, ent in lookup.items():
        text = _lookup(ent, cfg.go_column)
        for go in re.findall(r"GO:\d+", text or ""):
            if go not in gene_go[gene]:
                gene_go[gene].append(go)
    if cfg.go_gene_table and cfg.go_gene_table.exists():
        with cfg.go_gene_table.open(newline="", encoding="utf-8-sig") as fh:
            reader = csv.reader(fh, delimiter="\t")
            for row in reader:
                if len(row) < 2:
                    continue
                gene, go = row[0].strip(), row[1].strip()
                if gene and go and go not in gene_go[gene]:
                    gene_go[gene].append(go)
    return gene_go


def _load_gene_kegg(cfg: Config) -> Tuple[Dict[str, List[str]], Dict[str, str]]:
    gene_paths: Dict[str, List[str]] = defaultdict(list)
    if not cfg.kegg_gene_table or not cfg.kegg_gene_table.exists():
        return gene_paths, {}
    ko_path: Dict[str, List[str]] = defaultdict(list)
    if cfg.ko_pathway_table and cfg.ko_pathway_table.exists():
        with cfg.ko_pathway_table.open(newline="", encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                ko = parts[0].strip().removeprefix("ko:")
                path = parts[1].strip().removeprefix("path:ko")
                if path and path not in ko_path[ko]:
                    ko_path[ko].append(path)
    path_names: Dict[str, str] = {}
    if cfg.pathway_names and cfg.pathway_names.exists():
        with cfg.pathway_names.open(newline="", encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t", 1)
                if len(parts) == 2:
                    path_names[parts[0].strip().removeprefix("ko")] = parts[1].strip()
    with cfg.kegg_gene_table.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene = (_lookup(row, "GeneID") or _lookup(row, "gene") or "").strip()
            ko = (_lookup(row, "KO") or "").strip().removeprefix("ko:")
            if not gene or not ko:
                continue
            for path in ko_path.get(ko, []):
                if path not in gene_paths[gene]:
                    gene_paths[gene].append(path)
    return gene_paths, path_names


def annotate_degs(cfg: Config, dry_run: bool = False) -> None:
    if cfg.annotation_table is None:
        logging.info("annotation.table not configured; skipping DEG annotation")
        return
    lookup = _load_annotation_lookup(cfg)
    if not lookup:
        logging.warning("annotation table is empty or unreadable: %s", cfg.annotation_table)
        return
    cfg.annotated_dir.mkdir(parents=True, exist_ok=True)
    for contrast in cfg.contrasts:
        name = contrast["name"]
        src = cfg.de_dir / f"{name}_DEGs.csv"
        if not src.exists():
            logging.warning("DEG file not found: %s", src)
            continue
        dst = cfg.annotated_dir / f"{name}_DEGs_annotated.csv"
        if dry_run:
            logging.info("DRY RUN: annotate %s -> %s", src, dst)
            continue
        with _open_csv(src) as fh:
            reader = csv.DictReader(fh)
            base_fields = list(reader.fieldnames or [])
            all_extra = {
                key for ent in lookup.values() for key in ent.keys()
            }
            extra = [c for c in sorted(all_extra) if c not in base_fields]
            with dst.open("w", newline="", encoding="utf-8") as out:
                writer = csv.DictWriter(out, fieldnames=base_fields + extra)
                writer.writeheader()
                for row in reader:
                    ent = lookup.get((row.get("gene") or "").strip(), {})
                    for key in extra:
                        row.setdefault(key, ent.get(key, ""))
                    writer.writerow(row)
        logging.info("wrote %s", dst)
    if not dry_run:
        lookup_path = cfg.annotated_dir / "gene_annotation_lookup.csv"
        fields = sorted({k for ent in lookup.values() for k in ent})
        with lookup_path.open("w", newline="", encoding="utf-8") as out:
            writer = csv.DictWriter(out, fieldnames=["gene"] + fields)
            writer.writeheader()
            for gene in sorted(lookup):
                writer.writerow({"gene": gene, **lookup[gene]})


def _category_matches(
    category: Dict[str, object],
    product: str,
    go_text: str,
    go_ids: Sequence[str],
    pathways: Sequence[str],
    path_names: Dict[str, str],
) -> List[str]:
    evidence: List[str] = []
    go_set = {str(x).removeprefix("GO:").lstrip(":") for x in category.get("go", []) or []}
    go_set = {f"GO:{g}" if not g.startswith("GO:") else g for g in go_set}
    for go in go_ids:
        if go in go_set:
            evidence.append(f"GO:{go}")
    kegg_set = {str(x).removeprefix("ko:") for x in category.get("kegg", []) or []}
    for path in pathways:
        if path in kegg_set or f"ko{path}" in kegg_set:
            name = path_names.get(path, "")
            evidence.append(f"KEGG:ko{path}" + (f" ({name})" if name else ""))
    text = f"{product or ''} {go_text or ''}".lower()
    for kw in category.get("keywords", []) or []:
        if str(kw).lower() in text and f"annotation:{kw}" not in evidence:
            evidence.append(f"annotation:{kw}")
    return evidence


def extract_categories(cfg: Config, dry_run: bool = False) -> None:
    if not cfg.categories:
        logging.info("no analysis.categories configured; skipping category analysis")
        return
    lookup = _load_annotation_lookup(cfg)
    gene_go = _load_gene_go(cfg, lookup)
    gene_paths, path_names = _load_gene_kegg(cfg)
    stats = _read_gene_stats(cfg)
    cfg.category_dir.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    for key, category in cfg.categories.items():
        label = str(category.get("label", key))
        candidates: Dict[str, List[str]] = {}
        for gene, ent in lookup.items():
            product = _lookup(ent, cfg.product_column)
            evidence = _category_matches(
                category, product, _lookup(ent, cfg.go_column),
                gene_go.get(gene, []), gene_paths.get(gene, []),
                path_names,
            )
            if evidence:
                candidates[gene] = evidence
        if not candidates:
            logging.info("category %s: no candidate genes", label)
            continue

        matrix_rows = []
        deg_rows = []
        for gene in sorted(candidates):
            ent = lookup.get(gene, {})
            row = {
                "category": label,
                "gene": gene,
                "product": _lookup(ent, cfg.product_column),
                "reasons": "; ".join(candidates[gene]),
            }
            any_deg = False
            for contrast in cfg.contrasts:
                name = contrast["name"]
                st = stats.get(gene, {}).get(name, {})
                lfc = st.get("log2FoldChange")
                padj = st.get("padj")
                is_deg = bool(st.get("DEG"))
                any_deg = any_deg or is_deg
                row[f"{name}_log2FC"] = "" if lfc is None else round(float(lfc), 3)
                row[f"{name}_padj"] = "" if padj is None else f"{float(padj):.5g}"
                row[f"{name}_DEG"] = "yes" if is_deg else "no"
            matrix_rows.append(row)
            if any_deg:
                for contrast in cfg.contrasts:
                    name = contrast["name"]
                    st = stats.get(gene, {}).get(name, {})
                    if not st.get("DEG"):
                        continue
                    lfc = float(st["log2FoldChange"])
                    deg_rows.append({
                        "category": label,
                        "comparison": name,
                        "gene": gene,
                        "regulation": "up" if lfc > 0 else "down",
                        "log2FoldChange": round(lfc, 3),
                        "padj": f"{float(st['padj']):.5g}",
                        "product": _lookup(ent, cfg.product_column),
                        "reasons": "; ".join(candidates[gene]),
                    })
                    summary_rows.append({
                        "category": label, "comparison": name,
                        "direction": "up" if lfc > 0 else "down",
                    })

        if dry_run:
            logging.info("DRY RUN: category %s -> %s", label, cfg.category_dir)
            continue

        matrix_fields = ["category", "gene", "product", "reasons"]
        for contrast in cfg.contrasts:
            name = contrast["name"]
            matrix_fields += [f"{name}_log2FC", f"{name}_padj", f"{name}_DEG"]
        with (cfg.category_dir / f"combined_{key}_matrix.csv").open(
            "w", newline="", encoding="utf-8"
        ) as out:
            writer = csv.DictWriter(out, fieldnames=matrix_fields)
            writer.writeheader()
            writer.writerows(matrix_rows)

        deg_fields = [
            "category", "comparison", "gene", "regulation",
            "log2FoldChange", "padj", "product", "reasons",
        ]
        with (cfg.category_dir / f"{key}_DEGs.csv").open(
            "w", newline="", encoding="utf-8"
        ) as out:
            writer = csv.DictWriter(out, fieldnames=deg_fields)
            writer.writeheader()
            writer.writerows(deg_rows)

    if not dry_run:
        counts: Dict[Tuple[str, str, str], int] = defaultdict(int)
        for row in summary_rows:
            counts[(row["category"], row["comparison"], row["direction"])] += 1
        with (cfg.category_dir / "category_summary.csv").open(
            "w", newline="", encoding="utf-8"
        ) as out:
            writer = csv.writer(out)
            writer.writerow(["category", "comparison", "up", "down", "total"])
            for key, category in cfg.categories.items():
                label = str(category.get("label", key))
                for contrast in cfg.contrasts:
                    name = contrast["name"]
                    up = counts.get((label, name, "up"), 0)
                    down = counts.get((label, name, "down"), 0)
                    writer.writerow([label, name, up, down, up + down])
        logging.info("category analysis written to %s", cfg.category_dir)
