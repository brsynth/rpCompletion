import argparse

from rpCompletion import rpCompletion

def add_arguments(parser):
    parser.add_argument('-rp2paths_compounds', type=str)
    parser.add_argument('-rp2_pathways', type=str)
    parser.add_argument('-rp2paths_pathways', type=str)
    parser.add_argument('-upper_flux_bound', type=int, default=999999)
    parser.add_argument('-lower_flux_bound', type=int, default=0)
    parser.add_argument('-maxRuleIds', type=int, default=2)
    parser.add_argument('-pathway_id', type=str, default='rp_pathway')
    parser.add_argument('-compartment_id', type=str, default='MNXC3')
    parser.add_argument('-species_group_id', type=str, default='central_species')
    parser.add_argument('-pubchem_search', type=str, default='False')
    parser.add_argument('-output', type=str)
    return parser

def build_parser():
    parser = argparse.ArgumentParser('Python wrapper to parse RP2 to generate rpSBML collection')
    parser = rpCache_add_arguments(parser)
    parser = add_arguments(parser)

    return parser

def entrypoint(args=sys.argv[1:]):
    parser = build_parser()

    params = parser.parse_args(args)

    if params.maxRuleIds<0:
        logging.error('Max rule ID cannot be less than 0: '+str(params.maxRuleIds))
        exit(1)
    if params.pubchem_search=='True' or params.pubchem_search=='T' or params.pubchem_search=='true' or params.pubchem_search=='t':
        params.pubchem_search = True
    else:
        params.pubchem_search = False
    # elif params.pubchem_search=='False' or params.pubchem_search=='F' or params.pubchem_search=='false' or params.pubchem_search=='f':
    #     params.pubchem_search = False
    # else:
    #     logging.error('Cannot interpret pubchem_search input: '+str(params.pubchem_search))
    #     exit(1)

    if not os_path.exists(params.output):
        os_mkdir(params.output)
    rpreader = rpReader(params.store_mode, params.print)
    rpsbml_paths = rpreader.rp2ToSBML(
                             params.rp2paths_compounds,
                             params.rp2_pathways,
                             params.rp2paths_pathways,
                             params.output,
                             int(params.upper_flux_bound),
                             int(params.lower_flux_bound),
                             int(params.maxRuleIds),
                             params.pathway_id,
                             params.compartment_id,
                             params.species_group_id,
                             params.pubchem_search)



    # with tarfile.open(fileobj=args.output, mode='w:xz') as tf:
    #     tf.add(outdir, arcname=os_path.basename(args.output))
    # shutil_rmtree(outdir)
        # for rpsbml_name in rpsbml_paths:
        #     data = libsbml.writeSBMLToString(rpsbml_paths[rpsbml_name].document).encode('utf-8')
        #     fiOut = io.BytesIO(data)
        #     info = tarfile.TarInfo(name=rpsbml_name)
        #     info.size = len(data)
        #     tf.addfile(tarinfo=info, fileobj=fiOut)

##
#
#
if __name__ == "__main__":

    entrypoint(sys.argv[1:])

    # if not isOK:
    #     logging.error('Function returned an error')
    # ########IMPORTANT######
    # outputTar_bytes.seek(0)
    # #######################
    # with open(outputTar, 'wb') as f:
    #     shutil.copyfileobj(outputTar_bytes, f, length=131072)
