import sys
import textwrap
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from rnaseq_pipeline.config import ConfigError, load_config
from rnaseq_pipeline.downstream import _category_matches


class ConfigTest(unittest.TestCase):
    def _make_project(self, root: Path) -> Path:
        (root / "genome.fa").write_text(">chr1\nACGT\n", encoding="utf-8")
        (root / "genes.gtf").write_text(
            'chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id "g1";\n',
            encoding="utf-8",
        )
        (root / "samples.tsv").write_text(
            "sample\tgroup\treplicate\tfastq_1\tfastq_2\n"
            "c1\tcontrol\t1\tr1.fq.gz\tr2.fq.gz\n"
            "c2\tcontrol\t2\tr3.fq.gz\tr4.fq.gz\n"
            "t1\ttreat\t1\tr5.fq.gz\tr6.fq.gz\n"
            "t2\ttreat\t2\tr7.fq.gz\tr8.fq.gz\n",
            encoding="utf-8",
        )
        (root / "config.yaml").write_text(
            textwrap.dedent(
                """
                project:
                  name: test
                  outdir: results
                  threads: 2
                  jobs: 1
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
        return root / "config.yaml"

    def test_load_config(self):
        with TemporaryDirectory() as tmp:
            path = self._make_project(Path(tmp))
            cfg = load_config(path)
            self.assertEqual(cfg.project_name, "test")
            self.assertEqual(len(cfg.samples), 4)
            self.assertEqual(cfg.contrasts[0]["name"], "treat_vs_control")
            self.assertEqual(cfg.threads, 2)

    def test_invalid_group(self):
        with TemporaryDirectory() as tmp:
            path = self._make_project(Path(tmp))
            text = path.read_text(encoding="utf-8").replace("numerator: treat", "numerator: missing")
            path.write_text(text, encoding="utf-8")
            with self.assertRaises(ConfigError):
                load_config(path)


class CategoryTest(unittest.TestCase):
    def test_keyword_and_go_match(self):
        evidence = _category_matches(
            {"go": ["GO:0004096"], "keywords": ["catalase"], "kegg": []},
            "catalase",
            "",
            ["GO:0004096"],
            [],
            {},
        )
        self.assertTrue(any("GO:0004096" in e for e in evidence))
        self.assertTrue(any("catalase" in e for e in evidence))

    def test_kegg_match(self):
        evidence = _category_matches(
            {"go": [], "keywords": [], "kegg": ["ko00900"]},
            "HMG-CoA synthase",
            "",
            [],
            ["00900"],
            {"00900": "Terpenoid backbone biosynthesis"},
        )
        self.assertEqual(len(evidence), 1)
        self.assertIn("ko00900", evidence[0])


if __name__ == "__main__":
    unittest.main()
