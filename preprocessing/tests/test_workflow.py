import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SOURCE = Path(__file__).resolve().parents[1]


class WorkflowTest(unittest.TestCase):
    def setUp(self):
        runtime = SOURCE / ".runtime"
        runtime.mkdir(exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=runtime)
        self.root = Path(self.temp.name)
        self.preprocessing = self.root / "preprocessing"
        shutil.copytree(SOURCE / "scripts", self.preprocessing / "scripts")
        shutil.copy2(SOURCE / "pixi.sh", self.preprocessing / "pixi.sh")
        (self.preprocessing / "pixi.toml").write_text("manifest\n")
        (self.preprocessing / "pixi.lock").write_text("lock\n")
        for directory in (
            "alignment/mapping",
            "library_qc/fastqc",
            "library_qc/multiqc",
            "index/star-index",
            "quantification/Output",
            "raw/fastq",
            "bin",
        ):
            (self.root / directory).mkdir(parents=True)
        self.calls = self.root / "calls"
        pixi = self.root / "bin/pixi"
        pixi.write_text(
            '#!/bin/bash\nprintf "pixi:%s|%s|%s|%s|%s\\n" "$*" "$HOME" '
            '"$PIXI_CACHE_DIR" "$RAYON_NUM_THREADS" "$PWD" >> "$CALLS"\n'
        )
        pixi.chmod(0o755)
        sbatch = self.root / "bin/sbatch"
        sbatch.write_text(
            '#!/bin/bash\nprintf "sbatch:%s|%s|%s|%s\\n" '
            '"$*" "$PWD" "$PROJECT_ROOT" "$PIXI_BIN" >> "$CALLS"\necho 12345\n'
        )
        sbatch.chmod(0o755)
        self.env = os.environ.copy()
        self.env.update(
            PATH=f"{self.root / 'bin'}:{self.env['PATH']}",
            PIXI_BIN=str(pixi),
            CALLS=str(self.calls),
            PROJECT_ROOT=str(self.root),
            SLURM_CPUS_PER_TASK="8",
            SLURM_ARRAY_TASK_ID="0",
        )

    def tearDown(self):
        self.temp.cleanup()

    def run_script(self, path, *, check=True):
        return subprocess.run(
            ["bash", str(self.preprocessing / "scripts" / path)],
            cwd=self.root,
            env=self.env,
            text=True,
            capture_output=True,
            check=check,
        )

    def install(self):
        (self.preprocessing / ".pixi/envs/default").mkdir(parents=True)
        self.run_script("slurm-install.sh")

    def create_qc_inputs(self):
        (self.root / "raw/fastq/sample.fastq.gz").touch()
        index_dir = self.root / "index/star-index"
        (index_dir / "genomeParameters.txt").write_text("versionGenome\t2.7.11b\n")
        (index_dir / "chrNameLength.txt").write_text("chr1\t100\n")
        (index_dir / "Log.out").write_text(
            "STAR version=2.7.11b\n"
            "Genome sequence total length = 100\n"
            "Processed GTF: 3 transcripts, 4 exons (non-collapsed)\n"
            "Total junctions: 2\n"
        )
        sample_dir = self.root / "alignment/mapping/sample"
        sample_dir.mkdir(parents=True)
        (sample_dir / "Log.final.out").write_text("STAR alignment summary\n")
        (sample_dir / "ReadsPerGene.out.tab").write_text("N_unmapped\t1\t1\t1\n")
        (sample_dir / "SJ.out.tab").write_text(
            "chr1\t1\t2\t1\t0\t1\t4\t0\t10\n"
            "chr1\t3\t4\t1\t0\t0\t5\t0\t10\n"
        )
        (self.root / "library_qc/fastqc/sample_fastqc.zip").write_bytes(b"fastqc")
        for mode, assigned in (("unique", 8), ("multiple", 9)):
            (self.root / f"quantification/Output/featurecounts_summary_{mode}.tsv").write_text(
                f"Status\tsample\nAssigned\t{assigned}\nUnassigned_NoFeatures\t{10 - assigned}\n"
            )

    def test_installer_and_stale_lock(self):
        self.run_script("submit-install.sh")
        calls = self.calls.read_text()
        self.assertIn("slurm-install.sh", calls)
        self.assertIn(f"|{self.preprocessing}|{self.root}|", calls)
        self.install()
        calls = self.calls.read_text()
        self.assertIn("install --manifest-path", calls)
        self.assertIn(
            "--locked --run-post-link-scripts --concurrent-downloads 8 "
            "--concurrent-solves 1",
            calls,
        )
        self.assertIn(str(self.preprocessing / ".runtime/cache/pixi"), calls)
        self.run_script("reference/submit-build-index.sh")
        (self.preprocessing / "pixi.lock").write_text("changed\n")
        result = self.run_script("reference/submit-build-index.sh", check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("missing or out of date", result.stderr)

    def test_submitters_and_offline_jobs(self):
        self.install()
        result = self.run_script("qc/submit-fastqc.sh", check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no FASTQs", result.stderr)
        (self.root / "raw/fastq/sample.fastq.gz").touch()
        self.run_script("qc/submit-fastqc.sh")
        self.create_qc_inputs()
        for path in (
            "qc/submit-multiqc.sh",
            "alignment/submit-align-fastq.sh",
            "counting/submit-feature-counting.sh",
            "alignment/submit-package-outputs.sh",
            "reference/submit-build-index.sh",
            "qc/slurm-fastqc.sh",
            "qc/slurm-multiqc.sh",
            "alignment/slurm-align-fastq.sh",
            "counting/slurm-feature-counting.sh",
        ):
            self.run_script(path)
        self.run_script("qc/slurm-multiqc.sh")
        calls = self.calls.read_text()
        self.assertIn("--array=0-0%20", calls)
        self.assertIn("slurm-multiqc.sh", calls)
        self.assertIn("--array=0-0%5", calls)
        self.assertIn("prepare-multiqc-inputs.py", calls)
        self.assertIn("multiqc --file-list", calls)
        self.assertIn(f"|{self.root / 'library_qc'}|{self.root}|", calls)
        self.assertIn(f"|{self.root / 'alignment'}|{self.root}|", calls)
        for tool in ("multiqc", "python", "STAR", "samtools", "Rscript"):
            self.assertIn(f"--as-is --offline {tool}", calls)
        self.assertIn(
            f"--quantMode GeneCounts --readFilesIn {self.root / 'raw/fastq/sample.fastq.gz'} "
            "--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate "
            f"--outFileNamePrefix {self.root / 'alignment/mapping/sample'}/",
            calls,
        )

    def test_prepare_multiqc_inputs_preserves_outputs_and_summarizes_both_modes(self):
        helper = SOURCE / "scripts/qc/prepare-multiqc-inputs.py"
        result = subprocess.run(
            [sys.executable, str(helper), "--project-root", str(self.root), "--check-only"],
            text=True,
            capture_output=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("No FASTQs", result.stderr)

        self.create_qc_inputs()
        output_dir = self.root / "multiqc-inputs"
        subprocess.run(
            [sys.executable, str(helper), "--project-root", str(self.root), "--check-only"],
            check=True,
        )
        subprocess.run(
            [
                sys.executable,
                str(helper),
                "--project-root",
                str(self.root),
                "--output-dir",
                str(output_dir),
            ],
            check=True,
        )

        report_data = json.loads((output_dir / "featurecounts_modes_mqc.yaml").read_text())
        metrics = report_data["data"]["sample"]
        self.assertEqual(metrics["featurecounts_unique_assigned_pct"], 80.0)
        self.assertEqual(metrics["featurecounts_multiple_assigned_pct"], 90.0)
        self.assertEqual(metrics["star_junctions"], 2)
        self.assertEqual(metrics["star_annotated_junctions_pct"], 50.0)
        self.assertEqual(metrics["star_junction_unique_read_support"], 9)
        self.assertEqual(metrics["star_junction_multimapping_read_support"], 0)
        self.assertEqual(
            (output_dir / "featurecounts_unique.summary").read_text(),
            (self.root / "quantification/Output/featurecounts_summary_unique.tsv").read_text(),
        )
        self.assertIn(
            "Status\tsample [multimappers]",
            (output_dir / "featurecounts_multiple.summary").read_text(),
        )
        self.assertEqual(
            (self.root / "quantification/Output/featurecounts_summary_multiple.tsv").read_text().splitlines()[0],
            "Status\tsample",
        )
        self.assertIn("STAR reference index", (output_dir / "star_index_mqc.html").read_text())


if __name__ == "__main__":
    unittest.main()
