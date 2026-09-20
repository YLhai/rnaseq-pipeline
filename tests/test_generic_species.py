import csv
import sys
import textwrap
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rnaseq_pipeline import downstream
from rnaseq_pipeline.config import load_config
from rnaseq_pipeline.steps import Pipeline


class GenericSpeciesTest(unittest.TestCase):
    def test_input_is_clean_uses_sample_paths(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "genome.fa").write_text(">chr1\nACGT\n", encoding="utf-8")
            (root / "genes.gtf").write_text(
                'chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id "g1";\n',
                encoding="utf-8",
            )
            (root / "samples.tsv").write_text(
                "sample\tgroup\treplicate\tfastq_1\tfastq_2\n"
                "c1\tcontrol\t1\tclean/c1_R1.fq.gz\tclean/c1_R2.fq.gz\n"
                "t1\ttreat\t1\tclean/t1_R1.fq.gz\tclean/t1_R2.fq.gz\n",
                encoding="utf-8",
            )
            (root / "config.yaml").write_text(
                textwrap.dedent(
                    """
                    project:
                      name: clean_input
                      outdir: results
                      input_is_clean: true
                    samples: samples.tsv
                    reference:
                      genome_fasta: genome.fa
                      gtf: genes.gtf
                    analysis:
                      contrasts:
                        - name: treat_vs_control
                          numerator: treat
                          denominator: control
                    """
                ),
                encoding="utf-8",
            )
            cfg = load_config(root / "config.yaml")
            self.assertTrue(cfg.input_is_clean)
            in1, in2 = Pipeline(cfg)._input_paths(cfg.samples[0])
            self.assertEqual(in1, cfg.samples[0].fastq_1)
            self.assertEqual(in2, cfg.samples[0].fastq_2)

    def test_featurecounts_header_normalization(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "genome.fa").write_text(">chr1\nACGT\n", encoding="utf-8")
            (root / "genes.gtf").write_text(
                'chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id "g1";\n',
                encoding="utf-8",
            )
            (root / "samples.tsv").write_text(
                "sample\tgroup\treplicate\tfastq_1\tfastq_2\n"
                "c1\tcontrol\t1\tc1_R1.fq.gz\tc1_R2.fq.gz\n"
                "t1\ttreat\t1\tt1_R1.fq.gz\tt1_R2.fq.gz\n",
                encoding="utf-8",
            )
            (root / "config.yaml").write_text(
                textwrap.dedent(
                    """
                    project:
                      name: counts_header
                      outdir: results
                    samples: samples.tsv
                    reference:
                      genome_fasta: genome.fa
                      gtf: genes.gtf
                    analysis:
                      contrasts:
                        - name: treat_vs_control
                          numerator: treat
                          denominator: control
                    """
                ),
                encoding="utf-8",
            )
            cfg = load_config(root / "config.yaml")
            cfg.ensure_dirs()
            cfg.featurecounts_file.write_text(
                "# Program:featureCounts\n"
                "Geneid\tChr\tStart\tEnd\tStrand\tLength\t"
                "/tmp/results/bam/c1.sorted.bam\t/tmp/results/bam/t1.sorted.bam\n"
                "g1\tchr1\t1\t4\t+\t4\t10\t20\n",
                encoding="utf-8",
            )
            Pipeline(cfg)._normalize_featurecounts_header()
            header = cfg.featurecounts_file.read_text(
                encoding="utf-8"
            ).splitlines()[1].split("\t")
            self.assertEqual(header[6:], ["c1", "t1"])

    def test_annotation_from_gtf(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "genome.fa").write_text(">chr1\nACGT\n", encoding="utf-8")
            (root / "genes.gtf").write_text(
                'chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id "g1"; gene_name "CAT1";\n'
                'chr1\ttest\texon\t5\t8\t.\t+\t.\tgene_id "g2"; product "SOD";\n',
                encoding="utf-8",
            )
            (root / "samples.tsv").write_text(
                "sample\tgroup\treplicate\tfastq_1\tfastq_2\n"
                "c1\tcontrol\t1\tc1_R1.fq.gz\tc1_R2.fq.gz\n"
                "t1\ttreat\t1\tt1_R1.fq.gz\tt1_R2.fq.gz\n",
                encoding="utf-8",
            )
            (root / "config.yaml").write_text(
                textwrap.dedent(
                    """
                    project:
                      name: gtf_only
                      outdir: results
                    samples: samples.tsv
                    reference:
                      genome_fasta: genome.fa
                      gtf: genes.gtf
                    annotation:
                      table: ""
                      from_gtf: true
                    analysis:
                      contrasts:
                        - name: treat_vs_control
                          numerator: treat
                          denominator: control
                    """
                ),
                encoding="utf-8",
            )
            cfg = load_config(root / "config.yaml")
            path = Pipeline(cfg)._prepare_annotation()
            self.assertIsNotNone(path)
            text = path.read_text(encoding="utf-8")
            self.assertIn("CAT1", text)
            self.assertIn("SOD", text)

    def test_csv_annotation_and_extra_sample_columns(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "genome.fa").write_text(">chr1\nACGT\n", encoding="utf-8")
            (root / "genes.gtf").write_text(
                'chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id "g1";\n',
                encoding="utf-8",
            )
            (root / "samples.tsv").write_text(
                "Sample\tGroup\tReplicate\tBatch\tfastq_1\tfastq_2\n"
                "c1\tcontrol\t1\tB1\tc1_R1.fq.gz\tc1_R2.fq.gz\n"
                "c2\tcontrol\t2\tB1\tc2_R1.fq.gz\tc2_R2.fq.gz\n"
                "t1\ttreat\t1\tB2\tt1_R1.fq.gz\tt1_R2.fq.gz\n"
                "t2\ttreat\t2\tB2\tt2_R1.fq.gz\tt2_R2.fq.gz\n",
                encoding="utf-8",
            )
            (root / "annotation.csv").write_text(
                "gene_id,product,go_terms\n"
                'g1,catalase,"GO:0004096; GO:0006979"\n'
                "g2,hypothetical protein,\n",
                encoding="utf-8",
            )
            (root / "config.yaml").write_text(
                textwrap.dedent(
                    """
                    project:
                      name: generic
                      outdir: results
                    samples: samples.tsv
                    reference:
                      genome_fasta: genome.fa
                      gtf: genes.gtf
                      annotation_format: GTF
                      feature_type: exon
                      gene_attribute: gene_id
                    annotation:
                      table: annotation.csv
                      delimiter: ","
                      gene_id_column: GeneID
                      product_column: Product
                      go_column: GO Terms
                    analysis:
                      design_variable: batch
                      contrasts:
                        - name: treat_vs_control
                          numerator: treat
                          denominator: control
                      categories:
                        antioxidant:
                          label: Antioxidant
                          go: ["GO:0004096"]
                          keywords: ["catalase"]
                    """
                ),
                encoding="utf-8",
            )

            cfg = load_config(root / "config.yaml")
            self.assertEqual(cfg.design_variable, "Batch")
            self.assertEqual(cfg.annotation_delimiter, ",")

            # Minimal DESeq2-like outputs for one contrast.
            de_dir = cfg.de_dir
            de_dir.mkdir(parents=True, exist_ok=True)
            (de_dir / "treat_vs_control_all_results.csv").write_text(
                "baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,gene\n"
                "100,1.5,0.2,7,1e-10,1e-8,g1\n"
                "50,0.1,0.2,0.5,0.6,0.8,g2\n",
                encoding="utf-8",
            )
            (de_dir / "treat_vs_control_DEGs.csv").write_text(
                "baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,gene\n"
                "100,1.5,0.2,7,1e-10,1e-8,g1\n",
                encoding="utf-8",
            )
            cfg.counts_dir.mkdir(parents=True, exist_ok=True)
            cfg.featurecounts_file.write_text(
                "# command\n"
                "Geneid\tChr\tStart\tEnd\tStrand\tLength\tc1\tc2\tt1\tt2\n"
                "g1\tchr1\t1\t4\t+\t4\t10\t12\t50\t48\n"
                "g2\tchr1\t1\t4\t+\t4\t5\t6\t6\t5\n",
                encoding="utf-8",
            )

            downstream.annotate_degs(cfg)
            annotated = cfg.annotated_dir / "treat_vs_control_DEGs_annotated.csv"
            self.assertTrue(annotated.exists())
            with annotated.open(encoding="utf-8") as fh:
                row = next(csv.DictReader(fh))
            self.assertEqual(row["product"], "catalase")

            downstream.extract_categories(cfg)
            summary = cfg.category_dir / "category_summary.csv"
            self.assertTrue(summary.exists())
            text = summary.read_text(encoding="utf-8")
            self.assertIn("Antioxidant,treat_vs_control,1,0,1", text)


if __name__ == "__main__":
    unittest.main()
