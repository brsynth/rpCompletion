#!/usr/bin/env python

from logging      import error as logging_error
from brs_libs     import rpCache
from rpcompletion import rp2ToSBML, build_args_parser


def _cli():
    parser = build_args_parser()
    args  = parser.parse_args()

    args.pubchem_search = args.pubchem_search.lower() in ['true', 't']

    cache = rpCache(db=args.store_mode)

    try:
        rp2ToSBML(
            cache,
            args.rp2_pathways,
            args.rp2paths_compounds,
            args.rp2paths_pathways,
            args.outdir,
            int(args.upper_flux_bound),
            int(args.lower_flux_bound),
            int(args.max_subpaths_filter),
            args.pathway_id,
            args.compartment_id,
            args.species_group_id,
            args.sink_species_group_id,
            args.pubchem_search)
    except ValueError as e:
        logging_error(str(e))



if __name__ == '__main__':
    _cli()
