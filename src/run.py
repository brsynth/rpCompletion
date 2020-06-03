from sys import argv as sys_argv
from argparse import ArgumentParser as argparse_ArgumentParser
from os import path as os_path
from os import mkdir as os_mkdir

from sys import path as sys_path
sys_path.insert(0, '/home/rpCache')
from rpCache import rpCache
from rpCache import add_arguments as rpCache_add_arguments
from rpCompletion import rpCompletion

def add_arguments(parser):
    parser.add_argument('-rp2paths_compounds', type=str)
    parser.add_argument('-rp2_pathways', type=str)
    parser.add_argument('-rp2paths_pathways', type=str)
    parser.add_argument('-upper_flux_bound', type=int, default=999999)
    parser.add_argument('-lower_flux_bound', type=int, default=0)
    parser.add_argument('-max_subpaths_filter', type=int, default=10)
    parser.add_argument('-pathway_id', type=str, default='rp_pathway')
    parser.add_argument('-compartment_id', type=str, default='MNXC3')
    parser.add_argument('-species_group_id', type=str, default='central_species')
    parser.add_argument('-pubchem_search', type=str, default='False')
    parser.add_argument('-output', type=str)
    return parser

def build_parser():
    parser = argparse_ArgumentParser('Python wrapper to parse RP2 to generate rpSBML collection of unique and complete (cofactors) pathways')
    parser = rpCache_add_arguments(parser)
    parser = add_arguments(parser)

    return parser

def entrypoint(args=sys_argv[1:]):
    parser = build_parser()

    params = parser.parse_args(args)

    if params.max_subpaths_filter<0:
        logging.error('Max number of subpaths cannot be less than 0: '+str(params.max_subpaths_filter))
        exit(1)
    if params.pubchem_search=='True' or params.pubchem_search=='T' or params.pubchem_search=='true' or params.pubchem_search=='t':
        params.pubchem_search = True
    else:
        params.pubchem_search = False

    if not os_path.exists(params.output):
        os_mkdir(params.output)
    rpcompletion = rpCompletion(params.store_mode, params.print)
    rpsbml_paths = rpcompletion.rp2ToSBML(
                             params.rp2paths_compounds,
                             params.rp2_pathways,
                             params.rp2paths_pathways,
                             params.output,
                             int(params.upper_flux_bound),
                             int(params.lower_flux_bound),
                             int(params.max_subpaths_filter),
                             params.pathway_id,
                             params.compartment_id,
                             params.species_group_id,
                             params.pubchem_search)




##
#
#
if __name__ == "__main__":

    entrypoint(sys_argv[1:])
