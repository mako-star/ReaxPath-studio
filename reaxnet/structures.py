"""CLI for connectivity -> structure candidates -> species registration."""
import json

from reaxnet.connectivity import build_connectivity_parser, connectivity_config_from_args, run_connectivity
from reaxnet.parser import parse_element_mapping
from reaxnet.structure import StructureConfig


def main(argv=None):
    parser = build_connectivity_parser()
    parser.description = 'CHO structure inference and species registration from ReaxFF connectivity'
    parser.add_argument('--charges', required=True, help='Explicit allowed fragment total charges, e.g. 0 or -1,0,1')
    parser.add_argument('--max-radicals', required=True, type=int, help='Maximum radical electrons per fragment; 0 means closed shell')
    parser.add_argument('--allow-charge-separation', action='store_true', help='Allow opposite formal charges in the same fragment')
    parser.add_argument('--max-search-nodes', type=int, default=50000)
    parser.add_argument('--max-candidates', type=int, default=128)
    parser.add_argument('--max-atoms', type=int, default=128)
    parser.add_argument('--max-variable-edges', type=int, default=128)
    parser.add_argument('--cache-size', type=int, default=256)
    args = parser.parse_args(argv)
    try:
        structure_config = StructureConfig(tuple(int(q) for q in args.charges.split(',')), args.max_radicals,
                                           args.allow_charge_separation, args.max_search_nodes,
                                           args.max_candidates, args.max_atoms, args.max_variable_edges,
                                           args.cache_size)
        result = run_connectivity(args.bonds, parse_element_mapping(args.elements), args.output_dir,
                                  connectivity_config_from_args(args), max_frames=args.max_frames,
                                  structure_config=structure_config)
    except (ValueError, OSError) as exc:
        parser.exit(2, f'structures: error: {exc}\n')
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
