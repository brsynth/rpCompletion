#!/usr/bin/env python

from os import path, mkdir
from logging import error as logging_error
from sys import argv

from rpcompletion import rpCache, rpCompletion, build_args_parser


def gen_cache():
    rpCache.generate_cache('./cache-3.2')
    exit(0)

def _cli():
    parser = build_args_parser()
    args  = parser.parse_args()

    args.pubchem_search = args.pubchem_search.lower() in ['true', 't']

    rpcompletion = rpCompletion(db=args.store_mode)

    try:
        rpsbml_paths = rpcompletion.rp2ToSBML(
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
    if '--gen_cache' in argv[1:]:
        gen_cache()
    else:
        _cli()
