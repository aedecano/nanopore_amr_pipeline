#!/usr/bin/env python3

#Usage: python3 fetch_biosampleIDs_metadata.py test_biosampleID.txt 

import sys
import shutil
import subprocess
from pathlib import Path

def main():
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print(f"Usage: {Path(sys.argv[0]).name} <samples_file> [output_dir]", file=sys.stderr)
        sys.exit(2)

    samples_file = Path(sys.argv[1])
    out_dir = Path(sys.argv[2]) if len(sys.argv) == 3 else Path("Blackwell2021_Biosample_HTMLs")

    if shutil.which("esummary") is None:
        print("Error: 'esummary' not found in PATH. Install Entrez Direct (edirect) and try again.", file=sys.stderr)
        sys.exit(1)

    if not samples_file.is_file():
        print(f"Error: samples file not found: {samples_file}", file=sys.stderr)
        sys.exit(1)

    out_dir.mkdir(parents=True, exist_ok=True)

    with samples_file.open() as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue  # skip blank lines / comments

            html_path = Path(f"{s}.html")

            try:
                with html_path.open("w") as outfh:
                    # Run esummary and stream stdout directly to the file
                    subprocess.run(
                        ["esummary", "-db", "biosample", "-id", s],
                        stdout=outfh,
                        stderr=subprocess.PIPE,
                        text=True,
                        check=True,
                    )
            except subprocess.CalledProcessError as e:
                # Keep going for the next sample, but report the failure
                sys.stderr.write(
                    f"[WARN] esummary failed for {s} (exit {e.returncode}). "
                    f"Stderr:\n{e.stderr}\n"
                )
                continue

            # Copy the generated HTML into the output directory
            try:
                shutil.copy2(html_path, out_dir / html_path.name)
            except Exception as e:
                sys.stderr.write(f"[WARN] Copy failed for {html_path} -> {out_dir}: {e}\n")

if __name__ == "__main__":
    main()