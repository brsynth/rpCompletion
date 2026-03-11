from os import (
    path as os_path,
    mkdir as os_mkdir
)
from logging import (
    StreamHandler,
    Logger,
    getLogger,
)
from colored import fg, attr
from rr_cache import rrCache
from rplibs import build_args_parser
from rpcompletion import rp_completion
from brs_utils import init
from rpextractsink._version import __version__


def _cli():
    parser = build_args_parser(
        prog='rpcompletion',
        description='Parse RP2 pathways to generate rpSBML collection of unique and complete (cofactors) pathways'
    )
    args  = parser.parse_args()

    logger = init(parser, args, __version__)

    logger.debug('Parameters')
    logger.debug('   |--> rp2_metnet: '+str(args.rp2_metnet))
    logger.debug('   |--> sink: '+str(args.sink))
    logger.debug('   |--> rp2paths_compounds: '+str(args.rp2paths_compounds))
    logger.debug('   |--> rp2paths_pathways: '+str(args.rp2paths_pathways))
    logger.debug('   |--> outdir: '+str(args.outdir))
    logger.debug('   |--> upper_flux_bound: '+str(args.upper_flux_bound))
    logger.debug('   |--> lower_flux_bound: '+str(args.lower_flux_bound))
    logger.debug('   |--> maxsubpaths: '+str(args.maxsubpaths))


    check_args(
        args.maxsubpaths,
        args.outdir,
        logger
    )

    if args.cache_dir is None:
        cache = rrCache(
            cspace=args.cspace,
            interactive=False,
            logger=logger
        )
    else:
        cache = rrCache(
            cspace=args.cspace,
            install_dir=args.cache_dir,
            interactive=False,
            logger=logger
        )

    pathways = rp_completion(
        rp2_metnet=args.rp2_metnet,
        sink=args.sink,
        rp2paths_compounds=args.rp2paths_compounds,
        rp2paths_pathways=args.rp2paths_pathways,
        cache=cache,
        upper_flux_bound=int(args.upper_flux_bound),
        lower_flux_bound=int(args.lower_flux_bound),
        maxsubpaths=args.maxsubpaths,
        cofile=args.cofactors,
        forward=args.forward,
        logger=logger
    )

    # WRITE OUT
    if not os_path.exists(args.outdir):
        os_mkdir(args.outdir)
    # Write out selected pathways
    local_cache = {}
    for pathway in pathways:
        pathway.to_rpSBML(cache=cache, local_cache=local_cache).write_to_file(
            os_path.join(
                args.outdir,
                pathway.get_id()
            ) + '.xml'
        )

    StreamHandler.terminator = ""
    logger.info(
        '{color}{typo}Results are stored in {rst}'.format(
            color=fg('white'),
            typo=attr('bold'),
            rst=attr('reset')
        )
    )
    StreamHandler.terminator = "\n"
    logger.info(
        '{color}{outdir}\n'.format(
            color=fg('grey_70'),
            outdir=args.outdir
        )
    )


def check_args(
    max_subpaths_filter: int,
    outdir: str,
    logger: Logger = getLogger(__name__)
):
    logger.debug('Checking arguments...')
    logger.debug('   |--> max_subpaths_filter: '+str(max_subpaths_filter))
    logger.debug('   |--> outdir: '+str(outdir))

    if max_subpaths_filter < 0:
        raise ValueError('Max number of subpaths cannot be less than 0: '+str(max_subpaths_filter))

    if os_path.exists(outdir) and os_path.isfile(outdir):
        logger.error('Outdir name '+outdir+' already exists and is actually file. Stopping the process...')
        exit(-1)


if __name__ == '__main__':
    _cli()
