import os
import shutil
import subprocess
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
        for directory in ("alignment", "library_qc", "index", "raw/fastq", "bin"):
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
        result = self.run_script("qc/submit-quality-control.sh", check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no FASTQs", result.stderr)
        (self.root / "raw/fastq/sample.fastq.gz").touch()
        for path in (
            "qc/submit-quality-control.sh",
            "alignment/submit-align-fastq.sh",
            "counting/submit-feature-counting.sh",
            "alignment/submit-package-outputs.sh",
            "reference/submit-build-index.sh",
            "qc/slurm-quality-control.sh",
            "alignment/slurm-align-fastq.sh",
            "counting/slurm-feature-counting.sh",
        ):
            self.run_script(path)
        self.env["STEP"] = "multiqc"
        self.run_script("qc/slurm-quality-control.sh")
        calls = self.calls.read_text()
        self.assertIn("--array=0-0%20", calls)
        self.assertIn("--dependency=afterany:12345", calls)
        self.assertIn("--array=0-0%5", calls)
        self.assertIn(f"|{self.root / 'library_qc'}|{self.root}|", calls)
        self.assertIn(f"|{self.root / 'alignment'}|{self.root}|", calls)
        for tool in ("fastqc", "multiqc", "STAR", "samtools", "Rscript"):
            self.assertIn(f"--as-is --offline {tool}", calls)


if __name__ == "__main__":
    unittest.main()
