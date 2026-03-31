#!/usr/bin/env python3
"""
Batch run deepTools computeMatrix across multiple bigwig files and bed files.

Usage:
    # Command line arguments
    python batch_compute_matrix.py -b file1.bw file2.bw -r regions1.bed regions2.bed -o output_dir

    # JSON config file
    python batch_compute_matrix.py --config config.json

Example JSON config:
{
    "bigwigs": ["sample1.bw", "sample2.bw"],
    "regions": ["peaks.bed", "genes.bed"],  // or "regions_dir": "path/to/beds/" or ["dir1/", "dir2/"]
    "output_dir": "output/",
    "mode": "reference-point",
    "reference_point": "TSS",
    "upstream": 3000,
    "downstream": 3000,
    "bin_size": 50,
    "processors": 4
}
"""

import argparse
import json
import shlex
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


def load_config(config_path: str) -> dict:
    """Load configuration from a JSON file."""
    config_file = Path(config_path)
    if not config_file.exists():
        print(f"Error: Config file not found: {config_path}", file=sys.stderr)
        sys.exit(1)

    with open(config_file) as f:
        config = json.load(f)

    # Resolve relative paths in config relative to the config file's directory
    config_dir = config_file.parent

    if "bigwigs" in config:
        config["bigwigs"] = [
            str(config_dir / bw) if not Path(bw).is_absolute() else bw
            for bw in config["bigwigs"]
        ]

    if "regions_dir" in config:
        regions_dirs = config["regions_dir"]
        if isinstance(regions_dirs, str):
            regions_dirs = [regions_dirs]
        resolved_dirs = [
            str(config_dir / d) if not Path(d).is_absolute() else d
            for d in regions_dirs
        ]
        bed_files = []
        bed_subdirs = []
        for regions_dir in resolved_dirs:
            dir_beds = sorted(
                f for f in Path(regions_dir).iterdir()
                if f.suffix in (".bed", ".bedGraph", ".bedgraph")
                or f.name.endswith(".bed.gz") or f.name.endswith(".bedGraph.gz") or f.name.endswith(".bedgraph.gz")
            )
            if not dir_beds:
                print(f"Error: No .bed/.bedGraph files found in directory: {regions_dir}", file=sys.stderr)
                sys.exit(1)
            print(f"Found {len(dir_beds)} .bed/.bedGraph files in {regions_dir}")
            bed_files.extend(dir_beds)
            bed_subdirs.extend([Path(regions_dir).name] * len(dir_beds))
        config["regions"] = [str(f) for f in bed_files]
        config["bed_subdirs"] = bed_subdirs

    if "regions" in config:
        config["regions"] = [
            str(config_dir / bed) if not Path(bed).is_absolute() else bed
            for bed in config["regions"]
        ]

    if "output_dir" in config and not Path(config["output_dir"]).is_absolute():
        config["output_dir"] = str(config_dir / config["output_dir"])

    return config


CONDA_ENV = "deeptools"


def _run_single_bed(
    bed_file: str,
    bigwig_file: str,
    output_dir: str,
    mode: str,
    reference_point: str,
    upstream: int,
    downstream: int,
    bin_size: int,
    body_length: int,
    processors: int,
    skip_zeros: bool,
    missing_data_as_zero: bool,
    extra_args: list[str] | None,
    output_subdir: str = "",
) -> str:
    """Run computeMatrix for a single bed file and single bigwig. Returns a status message."""
    output_path = Path(output_dir) / output_subdir if output_subdir else Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    bed_name = Path(bed_file).stem
    bw_name = Path(bigwig_file).stem
    output_name = output_path / f"{bw_name}_{bed_name}_matrix.gz"
    output_matrix = f"{output_name}.tab"
    output_sorted_regions = output_path / f"{bw_name}_{bed_name}_sorted_regions.bed"

    cmd = [
        "computeMatrix",
        mode,
        "-S", bigwig_file,
        "-R", bed_file,
        "-o", str(output_name),
        "--outFileNameMatrix", str(output_matrix),
        "--outFileSortedRegions", str(output_sorted_regions),
        "-b", str(upstream),
        "-a", str(downstream),
        "--binSize", str(bin_size),
        "-p", str(processors),
    ]

    if mode == "reference-point":
        cmd.extend(["--referencePoint", reference_point])
    elif mode == "scale-regions":
        cmd.extend(["--regionBodyLength", str(body_length)])

    if skip_zeros:
        cmd.append("--skipZeros")

    if missing_data_as_zero:
        cmd.append("--missingDataAsZero")

    if extra_args:
        cmd.extend(extra_args)

    cmd_str = shlex.join(cmd)
    shell_cmd = f"source $(conda info --base)/etc/profile.d/conda.sh && conda activate {CONDA_ENV} && {cmd_str}"

    result = subprocess.run(
        shell_cmd,
        shell=True,
        executable="/bin/bash",
        check=True,
        capture_output=True,
        text=True,
    )

    msg = f"Successfully created: {output_matrix}"
    if result.stdout:
        msg += f"\n{result.stdout}"
    return msg


def run_compute_matrix(
    bigwig_files: list[str],
    bed_files: list[str],
    output_dir: str,
    mode: str = "reference-point",
    reference_point: str = "TSS",
    upstream: int = 3000,
    downstream: int = 3000,
    bin_size: int = 50,
    body_length: int = 5000,
    processors: int = 4,
    skip_zeros: bool = False,
    missing_data_as_zero: bool = False,
    extra_args: list[str] = None,
    parallel_beds: int = 1,
    bed_subdirs: list[str] | None = None,
    force: bool = False,
):
    """
    Run computeMatrix for each bed file with all bigwig files.

    Args:
        bigwig_files: List of paths to bigwig files
        bed_files: List of paths to bed/region files
        output_dir: Directory to save output matrix files
        mode: 'reference-point' or 'scale-regions'
        reference_point: Reference point for reference-point mode (TSS, TES, center)
        upstream: Distance upstream of reference point/start
        downstream: Distance downstream of reference point/end
        bin_size: Bin size in bp
        body_length: Body length for scale-regions mode
        processors: Number of processors to use per job
        skip_zeros: Skip regions with zero coverage
        missing_data_as_zero: Treat missing data as zero
        extra_args: Additional arguments to pass to computeMatrix
        parallel_beds: Number of bed files to process in parallel
        bed_subdirs: Optional per-bed output subdirectory names (parallel to bed_files)
        force: Rerun even if output files already exist
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    subdirs = bed_subdirs if bed_subdirs else [""] * len(bed_files)

    common_kwargs = dict(
        output_dir=output_dir,
        mode=mode,
        reference_point=reference_point,
        upstream=upstream,
        downstream=downstream,
        bin_size=bin_size,
        body_length=body_length,
        processors=processors,
        skip_zeros=skip_zeros,
        missing_data_as_zero=missing_data_as_zero,
        extra_args=extra_args,
    )

    # Skip (bed, bw) pairs whose outputs already exist (unless force=True)
    jobs_to_run = []
    total_jobs = len(bed_files) * len(bigwig_files)
    for bed, subdir in zip(bed_files, subdirs):
        bed_name = Path(bed).stem
        for bw in bigwig_files:
            bw_name = Path(bw).stem
            output_matrix = (
                output_path / subdir / f"{bw_name}_{bed_name}_matrix.gz"
                if subdir else
                output_path / f"{bw_name}_{bed_name}_matrix.gz"
            )
            if output_matrix.exists() and not force:
                print(f"Skipping {bw_name}_{bed_name} (output already exists: {output_matrix})")
            else:
                jobs_to_run.append((bed, bw, subdir))

    if not jobs_to_run:
        print("All jobs already processed, nothing to do.")
        return

    print(f"Processing {len(jobs_to_run)}/{total_jobs} (bed, bigwig) pairs ({parallel_beds} in parallel)")

    with ProcessPoolExecutor(max_workers=parallel_beds) as executor:
        futures = {
            executor.submit(_run_single_bed, bed_file=bed, bigwig_file=bw, output_subdir=subdir, **common_kwargs): (bed, bw)
            for bed, bw, subdir in jobs_to_run
        }

        failed = []
        for future in as_completed(futures):
            bed, bw = futures[future]
            try:
                msg = future.result()
                print(msg)
            except subprocess.CalledProcessError as e:
                print(f"Error processing {bw} x {bed}:", file=sys.stderr)
                print(e.stderr, file=sys.stderr)
                failed.append((bed, bw))

    if failed:
        print(f"\n{len(failed)} job(s) failed:", file=sys.stderr)
        for bed, bw in failed:
            print(f"  {bw} x {bed}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Batch run deepTools computeMatrix across bigwig and bed files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Reference-point mode (default)
  python batch_compute_matrix.py -b *.bw -r regions.bed -o output/

  # Scale-regions mode
  python batch_compute_matrix.py -b *.bw -r genes.bed -o output/ --mode scale-regions

  # Custom parameters
  python batch_compute_matrix.py -b sample1.bw sample2.bw -r peaks.bed -o output/ \\
      --upstream 5000 --downstream 5000 --bin-size 25 -p 8

  # Using a JSON config file
  python batch_compute_matrix.py --config config.json
        """,
    )

    parser.add_argument(
        "-c", "--config",
        help="JSON config file (overrides other arguments if provided)",
    )

    parser.add_argument(
        "-b", "--bigwigs",
        nargs="+",
        help="Bigwig files to process",
    )
    parser.add_argument(
        "-r", "--regions",
        nargs="+",
        help="BED files with regions of interest",
    )
    parser.add_argument(
        "--regions-dir",
        nargs="+",
        help="Directory (or directories) containing .bed files to use as regions (alternative to -r)",
    )
    parser.add_argument(
        "-o", "--output-dir",
        help="Output directory for matrix files",
    )
    parser.add_argument(
        "--mode",
        choices=["reference-point", "scale-regions"],
        default="reference-point",
        help="computeMatrix mode (default: reference-point)",
    )
    parser.add_argument(
        "--reference-point",
        choices=["TSS", "TES", "center"],
        default="TSS",
        help="Reference point for reference-point mode (default: TSS)",
    )
    parser.add_argument(
        "--upstream", "-u",
        type=int,
        default=3000,
        help="Distance upstream of reference point (default: 3000)",
    )
    parser.add_argument(
        "--downstream", "-d",
        type=int,
        default=3000,
        help="Distance downstream of reference point (default: 3000)",
    )
    parser.add_argument(
        "--bin-size",
        type=int,
        default=50,
        help="Bin size in bp (default: 50)",
    )
    parser.add_argument(
        "--body-length",
        type=int,
        default=5000,
        help="Body length for scale-regions mode (default: 5000)",
    )
    parser.add_argument(
        "-p", "--processors",
        type=int,
        default=4,
        help="Number of processors (default: 4)",
    )
    parser.add_argument(
        "--skip-zeros",
        action="store_true",
        help="Skip regions with zero coverage",
    )
    parser.add_argument(
        "--missing-data-as-zero",
        action="store_true",
        help="Treat missing data as zero",
    )
    parser.add_argument(
        "--parallel-beds",
        type=int,
        default=1,
        help="Number of bed files to process in parallel (default: 1)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Rerun computeMatrix even if output files already exist",
    )
    parser.add_argument(
        "--extra-args",
        nargs=argparse.REMAINDER,
        help="Additional arguments to pass to computeMatrix",
    )

    args = parser.parse_args()

    # Load config from JSON if provided
    if args.config:
        config = load_config(args.config)
        bigwigs = config.get("bigwigs", [])
        regions = config.get("regions", [])
        regions_subdirs = config.get("bed_subdirs")
        output_dir = config.get("output_dir")
        mode = config.get("mode", "reference-point")
        reference_point = config.get("reference_point", "TSS")
        upstream = config.get("upstream", 3000)
        downstream = config.get("downstream", 3000)
        bin_size = config.get("bin_size", 50)
        body_length = config.get("body_length", 5000)
        processors = config.get("processors", 4)
        skip_zeros = config.get("skip_zeros", False)
        missing_data_as_zero = config.get("missing_data_as_zero", False)
        extra_args = config.get("extra_args")
        parallel_beds = config.get("parallel_beds", 1)
        force = config.get("force", False)
        # CLI --force overrides config
        if args.force:
            force = True
    else:
        # Use command line arguments
        if not args.bigwigs or not args.output_dir:
            parser.error("--bigwigs and --output-dir are required unless --config is provided")
        if not args.regions and not args.regions_dir:
            parser.error("--regions or --regions-dir is required unless --config is provided")

        bigwigs = args.bigwigs
        if args.regions_dir:
            bed_files = []
            regions_subdirs = []
            for regions_dir_str in args.regions_dir:
                regions_dir = Path(regions_dir_str)
                dir_beds = sorted(
                    f for f in regions_dir.iterdir()
                    if f.suffix in (".bed", ".bedGraph", ".bedgraph")
                    or f.name.endswith(".bed.gz") or f.name.endswith(".bedGraph.gz") or f.name.endswith(".bedgraph.gz")
                )
                if not dir_beds:
                    print(f"Error: No .bed/.bedGraph files found in directory: {regions_dir_str}", file=sys.stderr)
                    sys.exit(1)
                print(f"Found {len(dir_beds)} .bed/.bedGraph files in {regions_dir_str}")
                bed_files.extend(dir_beds)
                regions_subdirs.extend([regions_dir.name] * len(dir_beds))
            regions = [str(f) for f in bed_files]
        else:
            regions = args.regions
            regions_subdirs = None
        output_dir = args.output_dir
        mode = args.mode
        reference_point = args.reference_point
        upstream = args.upstream
        downstream = args.downstream
        bin_size = args.bin_size
        body_length = args.body_length
        processors = args.processors
        skip_zeros = args.skip_zeros
        missing_data_as_zero = args.missing_data_as_zero
        extra_args = args.extra_args
        parallel_beds = args.parallel_beds
        force = args.force

    # Validate required fields
    if not bigwigs:
        print("Error: No bigwig files specified", file=sys.stderr)
        sys.exit(1)
    if not regions:
        print("Error: No region files specified", file=sys.stderr)
        sys.exit(1)
    if not output_dir:
        print("Error: No output directory specified", file=sys.stderr)
        sys.exit(1)

    # Validate input files exist
    for bw in bigwigs:
        if not Path(bw).exists():
            print(f"Error: Bigwig file not found: {bw}", file=sys.stderr)
            sys.exit(1)

    for bed in regions:
        if not Path(bed).exists():
            print(f"Error: BED file not found: {bed}", file=sys.stderr)
            sys.exit(1)

    run_compute_matrix(
        bigwig_files=bigwigs,
        bed_files=regions,
        output_dir=output_dir,
        mode=mode,
        reference_point=reference_point,
        upstream=upstream,
        downstream=downstream,
        bin_size=bin_size,
        body_length=body_length,
        processors=processors,
        skip_zeros=skip_zeros,
        missing_data_as_zero=missing_data_as_zero,
        extra_args=extra_args,
        parallel_beds=parallel_beds,
        bed_subdirs=regions_subdirs,
        force=force,
    )


if __name__ == "__main__":
    main()
