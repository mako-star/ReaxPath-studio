from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from reaxnet import __version__
from reaxnet.analysis import AnalyzerConfig
from reaxnet.parser import (
    BondsParseError,
    parse_element_mapping,
)
from reaxnet.streaming import analyze_file_streaming


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="reaxnet",
        description=(
            "Analyze a LAMMPS fix reax/c/bonds trajectory into molecular species, "
            "reaction events, canonical identifiers, and a reaction network."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("--bonds", required=True, type=Path, help="ReaxFF bonds trajectory")
    parser.add_argument(
        "--elements",
        required=True,
        help='LAMMPS atom-type mapping, for example "1:C,2:H,3:O"',
    )
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory for output files")

    parser.add_argument("--expected-interval", type=int, help="Expected saved timestep interval; inferred from first two frames if omitted")
    parser.add_argument(
        "--timestep-scale",
        type=float,
        default=1.0,
        help="Picoseconds per LAMMPS timestep (default: 1.0)",
    )
    parser.add_argument(
        "--recrossing-window",
        type=float,
        default=0.0,
        help="Filter two-step and multistep state returns within this many ps; 0 disables it",
    )


    parser.add_argument(
        "--no-identifiers",
        action="store_true",
        help="Skip optional RDKit SMILES/InChI generation",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        element_mapping = parse_element_mapping(args.elements)
        config = AnalyzerConfig(
            timestep_scale=args.timestep_scale,
            expected_interval=args.expected_interval,
            recrossing_window=args.recrossing_window,
            generate_identifiers=not args.no_identifiers,
        )
        result = analyze_file_streaming(
            args.bonds,
            element_mapping,
            args.output_dir,
            config,
        )
    except (BondsParseError, OSError, ValueError) as exc:
        parser.exit(2, f"reaxnet: error: {exc}\n")

    print(
        f"Analyzed {result.summary['frame_count']} frames in streaming mode, "
        f"found {result.summary['species_count']} species, and accepted "
        f"{result.summary['accepted_event_count']}/{result.summary['detected_event_count']} "
        "reaction events."
    )
    for name, path in result.paths.items():
        print(f"{name}: {path.resolve()}")
    return 0
