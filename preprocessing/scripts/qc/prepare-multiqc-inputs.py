#!/usr/bin/env python3
"""Prepare small, temporary MultiQC inputs without changing pipeline outputs."""

import argparse
import csv
import html
import json
import re
import shutil
import sys
from pathlib import Path


def read_featurecounts_summary(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or reader.fieldnames[0] != "Status":
            raise ValueError(f"Unexpected featureCounts summary header in {path}")
        rows = list(reader)

    samples = reader.fieldnames[1:]
    statuses = {row["Status"]: row for row in rows}
    if not samples or "Assigned" not in statuses:
        raise ValueError(f"Missing samples or Assigned row in {path}")

    assigned = {}
    for sample in samples:
        total = sum(int(row[sample]) for row in rows)
        assigned_count = int(statuses["Assigned"][sample])
        assigned[sample] = 100.0 * assigned_count / total if total else 0.0
    return samples, assigned


def read_star_junction_stats(mapping_dir, samples):
    stats = {}
    for path in mapping_dir.glob("*/SJ.out.tab"):
        sample = path.parent.name
        if sample not in samples:
            continue
        junctions = annotated = unique_support = multimapping_support = 0
        with path.open() as handle:
            for line_number, line in enumerate(handle, 1):
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    raise ValueError(f"Malformed STAR junction row at {path}:{line_number}")
                junctions += 1
                annotated += int(fields[5] == "1")
                unique_support += int(fields[6])
                multimapping_support += int(fields[7])
        stats[sample] = {
            "star_junctions": junctions,
            "star_annotated_junctions_pct": 100.0 * annotated / junctions if junctions else 0.0,
            "star_junction_unique_read_support": unique_support,
            "star_junction_multimapping_read_support": multimapping_support,
        }
    return stats


def read_index_parameters(index_dir):
    params = {}
    with (index_dir / "genomeParameters.txt").open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t", 1)
            if len(fields) == 2 and not fields[0].startswith("#"):
                params[fields[0]] = fields[1]

    chromosome_lengths = []
    with (index_dir / "chrNameLength.txt").open() as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                raise ValueError(f"Malformed chromosome-length row at {index_dir / 'chrNameLength.txt'}:{line_number}")
            chromosome_lengths.append((fields[0], int(fields[1])))

    log_text = (index_dir / "Log.out").read_text(errors="replace")
    log_version = re.search(r"^STAR version=(.+)$", log_text, re.MULTILINE)
    genome_length = re.search(r"Genome sequence total length =\s*(\d+)", log_text)
    annotation = re.search(
        r"^\s*(\d+) transcripts.*^\s*(\d+) exons \(non-collapsed\).*^Total junctions:\s*(\d+)",
        log_text,
        re.MULTILINE | re.DOTALL,
    )

    return {
        "index_version": params.get("versionGenome", "not recorded"),
        "build_log_version": log_version.group(1).strip() if log_version else "not recorded",
        "fasta": params.get("genomeFastaFiles", "not recorded").strip(),
        "annotation": params.get("sjdbGTFfile", "not recorded").strip(),
        "sjdb_overhang": params.get("sjdbOverhang", "not recorded"),
        "chromosomes": len(chromosome_lengths),
        "chromosome_bases": sum(length for _, length in chromosome_lengths),
        "genome_bases": int(genome_length.group(1)) if genome_length else None,
        "annotation_counts": annotation.groups() if annotation else None,
    }


def validate_pipeline_inputs(project_root):
    fastq_dir = project_root / "raw" / "fastq"
    fastqc_dir = project_root / "library_qc" / "fastqc"
    mapping_dir = project_root / "alignment" / "mapping"
    count_dir = project_root / "quantification" / "Output"
    index_dir = project_root / "index" / "star-index"

    fastqs = sorted(fastq_dir.glob("*.fastq.gz"))
    if not fastqs:
        raise ValueError(f"No FASTQs found in {fastq_dir}")
    samples = [path.name[: -len(".fastq.gz")] for path in fastqs]

    required_index_files = ("genomeParameters.txt", "Log.out", "chrNameLength.txt")
    for filename in required_index_files:
        if not (index_dir / filename).is_file():
            raise ValueError(f"STAR index is incomplete; missing {index_dir / filename}")

    for sample in samples:
        required_reports = (
            fastqc_dir / f"{sample}_fastqc.zip",
            mapping_dir / sample / "Log.final.out",
            mapping_dir / sample / "ReadsPerGene.out.tab",
            mapping_dir / sample / "SJ.out.tab",
        )
        for path in required_reports:
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f"Quality reporting is not ready; missing or empty {path}")

    unique_file = count_dir / "featurecounts_summary_unique.tsv"
    multiple_file = count_dir / "featurecounts_summary_multiple.tsv"
    for path in (unique_file, multiple_file):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Feature counting is incomplete; missing or empty {path}")
    unique_samples, _ = read_featurecounts_summary(unique_file)
    multiple_samples, _ = read_featurecounts_summary(multiple_file)
    if unique_samples != samples or multiple_samples != samples:
        raise ValueError(
            "FeatureCounts summary samples do not match the FASTQ samples; "
            "wait for feature counting to finish and check for missing samples"
        )


def write_multiqc_inputs(project_root, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    count_dir = project_root / "quantification" / "Output"
    unique_file = count_dir / "featurecounts_summary_unique.tsv"
    multiple_file = count_dir / "featurecounts_summary_multiple.tsv"
    samples, unique_assigned = read_featurecounts_summary(unique_file)
    multiple_samples, multiple_assigned = read_featurecounts_summary(multiple_file)
    if samples != multiple_samples:
        raise ValueError("Unique and multimapper featureCounts summaries have different sample columns")

    shutil.copyfile(unique_file, output_dir / "featurecounts_unique.summary")
    with multiple_file.open(newline="") as source, (
        output_dir / "featurecounts_multiple.summary"
    ).open("w", newline="") as target:
        rows = list(csv.reader(source, delimiter="\t"))
        rows[0][1:] = [f"{sample} [multimappers]" for sample in rows[0][1:]]
        csv.writer(target, delimiter="\t", lineterminator="\n").writerows(rows)

    junction_stats = read_star_junction_stats(project_root / "alignment" / "mapping", set(samples))
    data = {}
    for sample in samples:
        data[sample] = {
            "featurecounts_unique_assigned_pct": unique_assigned[sample],
            "featurecounts_multiple_assigned_pct": multiple_assigned[sample],
        }
        data[sample].update(junction_stats.get(sample, {}))

    general_stats = {
        "id": "pipeline_qc_general_stats",
        "section_name": "Pipeline quality metrics",
        "description": (
            "Comparison of featureCounts assignment rates with unique reads only and with "
            "multimappers included, plus STAR splice-junction support."
        ),
        "plot_type": "generalstats",
        "headers": [
            {
                "featurecounts_unique_assigned_pct": {
                    "title": "FC assigned (unique)",
                    "description": "Percentage of featureCounts records assigned to genes, unique reads only",
                    "min": 0,
                    "max": 100,
                    "suffix": "%",
                    "scale": "RdYlGn",
                }
            },
            {
                "featurecounts_multiple_assigned_pct": {
                    "title": "FC assigned (multi)",
                    "description": "Percentage of featureCounts records assigned to genes, multimappers included",
                    "min": 0,
                    "max": 100,
                    "suffix": "%",
                    "scale": "RdYlGn",
                }
            },
            {
                "star_junctions": {
                    "title": "STAR junctions",
                    "description": "Distinct splice junctions detected",
                    "min": 0,
                }
            },
            {
                "star_annotated_junctions_pct": {
                    "title": "Annotated junctions",
                    "description": "Percentage of detected splice junctions annotated in the reference",
                    "min": 0,
                    "max": 100,
                    "suffix": "%",
                    "scale": "RdYlGn",
                }
            },
            {
                "star_junction_unique_read_support": {
                    "title": "Junction read support",
                    "description": "Unique-mapping reads supporting splice junctions",
                    "min": 0,
                }
            },
            {
                "star_junction_multimapping_read_support": {
                    "title": "Multi-map junction support",
                    "description": "Multimapping reads supporting splice junctions",
                    "min": 0,
                }
            },
        ],
        "data": data,
    }
    (output_dir / "featurecounts_modes_mqc.yaml").write_text(json.dumps(general_stats, indent=2))

    index = read_index_parameters(project_root / "index" / "star-index")
    annotation_counts = index["annotation_counts"]
    details = [
        ("Index metadata version", index["index_version"]),
        ("STAR version in build log", index["build_log_version"]),
        ("FASTA", index["fasta"]),
        ("Annotation GTF", index["annotation"]),
        ("Splice-junction overhang", f"{index['sjdb_overhang']} bp"),
        ("Reference chromosomes", str(index["chromosomes"])),
        ("Bases across chromosomes", f"{index['chromosome_bases']:,}"),
    ]
    if index["genome_bases"] is not None:
        details.append(("Genome length recorded in build log", f"{index['genome_bases']:,} bp"))
    if annotation_counts:
        details.extend(
            [
                ("GTF transcripts", f"{int(annotation_counts[0]):,}"),
                ("GTF exons", f"{int(annotation_counts[1]):,}"),
                ("GTF splice junctions", f"{int(annotation_counts[2]):,}"),
            ]
        )
    rows = "\n".join(
        f"<tr><th>{html.escape(label)}</th><td>{html.escape(value)}</td></tr>"
        for label, value in details
    )
    index_html = f"""<!--
id: "star_index_metadata"
section_name: "STAR reference index"
description: "Build parameters and annotation summary read from the generated STAR index."
-->
<table class="table table-striped">
<tbody>
{rows}
</tbody>
</table>
"""
    (output_dir / "star_index_mqc.html").write_text(index_html)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    try:
        validate_pipeline_inputs(args.project_root)
        if args.check_only:
            print("MultiQC inputs are ready.")
            return
        if args.output_dir is None:
            parser.error("--output-dir is required unless --check-only is used")
        write_multiqc_inputs(args.project_root, args.output_dir)
    except (OSError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        raise SystemExit(1) from error


if __name__ == "__main__":
    main()
