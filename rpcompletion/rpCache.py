from os import path as os_path
from os import mkdir as os_mkdir
from os import remove as os_rm
from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi, MolToInchiKey, AddHs
from csv import DictReader as csv_DictReader
from csv import reader as csv_reader
from logging import getLogger as logging_getLogger
from json import dump as json_dump
from json import load as json_load
from gzip import open as gzip_open
from re import findall as re_findall
from tarfile import open as tarfile_open
from shutil import move as shutil_move
from shutil import rmtree as shutil_rmtree
import sys
import time
from itertools import chain as itertools_chain
from brs_utils import print_OK, print_FAILED, download, file_length
from requests import exceptions as r_exceptions
from tarfile import open as tf_open
from redis import StrictRedis
from credisdict import CRedisDict, wait_for_redis
import redis_server
# from subprocess import Popen,PIPE
from argparse import ArgumentParser as argparse_ArgParser
from hashlib import sha512
from pathlib import Path
from colored import attr as c_attr


#######################################################
################### rpCache  ##########################
#######################################################


def add_arguments(parser):
    parser.add_argument('-sm', '--store_mode', type=str, default='file',
                        help='data storage mode: file or db')
    parser.add_argument('--gen_cache', action='store_true',
                        help='generate the cache and exits')
    parser.add_argument('-p', '--print', type=bool, default=False,
                        help='print additional informations')
    return parser

def build_parser():
    return add_arguments(argparse_ArgParser('Python script to pre-compute data'))

# def entrypoint(args=sys.argv[1:]):
#     parser = build_parser()
#
#     params = parser.parse_args(args)
#
#     rpcache = rpCache(params.store_mode, params.print)
#
# ##
# #
# #
# if __name__ == "__main__":
#     entrypoint()

## Class to generate the cache
#
# Contains all the functions that parse different files, used to calculate the thermodynamics and the FBA of the
#the other steps. These should be called only when the files have changes
class rpCache:

    logger = logging_getLogger(__name__)
    logger.info('Started instance of rpCache')

    # _input_cache_url = 'ftp://ftp.vital-it.ch/databases/metanetx/MNXref/3.2/'
    _cache_url       = 'https://gitlab.com/breakthewall/rpcache-data/-/raw/master/'

    # static attribues
    _convertMNXM = {'MNXM162231': 'MNXM6',
                    'MNXM84': 'MNXM15',
                    'MNXM96410': 'MNXM14',
                    'MNXM114062': 'MNXM3',
                    'MNXM145523': 'MNXM57',
                    'MNXM57425': 'MNXM9',
                    'MNXM137': 'MNXM588022'}

    # name: sha512sum
    _input_cache_files = {
            'chem_xref.tsv.gz':    'e558110990dcc75af943863790dc55360fd2d40ecb17d02335377671e80f0ab3738fd556acb340e03e48dd1afdec3eece1e92df1e18bc24e7445f24f778a10da',
            'reac_xref.tsv.gz':    '48b991cf4a9c2ca573d395cf35c378881ed79e87772827647bfab2f6345499698664e07195ec10b342fc0164304dbd2363cccff1a1182225e6afebce3c16448b',
            'compounds.tsv.gz':    '719716bb880257bd014e045c03eb8dd12e2bbeba3aa52e38e9632ce605817b9dc09530e81fadd25542c0a439bdb81e1dfbd3a38f35b30b061845d1a880dbfe01',
            'chem_prop.tsv.gz':    'f2d220d1f0425e5e47f01e7deccfa46b60094d43b9f62b191ffb0fab8c00ef79e87c3b71d10bdcd26020608094f24884f51b3ebc3d7d3c9a6d594c6eaa324c66',
            'retrorules_rr02_flat_all.tsv.gz':   '890bdd24042c0192b5538964d775feefcb6cff9ad5f35690bfbfc5ae09334dd19df6828cdfc7f57a2018e090571517122b99d8760128052af898c638ae667e24',
            'comp_xref.tsv.gz':    '913a827f3645fda1699676ae6c32b9d7a8debae97ce7b0c386d8447f4eee5aa721d31bfb856d4092b3d5e987a8f19a6fe4bd28ddf1c5df5f85e71c3625bd1d81',
            'rxn_recipes.tsv.gz':  'dc0624f5ed7ab0b691d9a6ba02571a5cf334cfdb3109e78c98708e31574c46aeac2a97e9433788d80490ff80337679ccfd706cbb8e71a11cdc6122573bb69b0f'
            }

    _attributes = [
            'deprecatedMNXM_mnxm',
            'deprecatedMNXR_mnxr',
            'mnxm_strc',
            'chemXref',
            'chebi_mnxm',
            'rr_reactions',
            'inchikey_mnxm',
            'compXref',
            'name_compXref',
            'full_reactions'
    ]

    # name: sha512sum
    _cache_files = {
            _attributes[0]+'.gz': '698a3e83cf4f9206ea2644c9c35a9af53957838baaae6efb245d02b6b8d0ea8b25c75008e562b99ba3e0189e50ee47655376f2d0635f6206e0015f91f0e4bad8',
            _attributes[1]+'.gz': '51554c6f6ae99c6755da7496208b3feec30547bc4cf3007d9fd30f46fa4c0cc73bad5aeb743dca07e32711c4346504296bee776d135fb18e96c891a0086fc87e',
            _attributes[2]+'.gz': '0021ef63165d75ee6b8c209ccf14b8a1b8b7b263b4077f544729c47b5525f66511c3fa578fd2089201abb61693085b9912639e62f7b7481d06ad1f38bfc2dd8e',
            _attributes[3]+'.gz': '7d559cc7389c0cb2bd10f92e6e845bb5724be64d1624adc4e447111fc63599bb69396cd0cc3066a6bb19910c00e266c97e21b1254d9a6dc9da3a8b033603fcff',
            _attributes[4]+'.gz': '587d6c5206ee94e63af6d9eaf49fd5e2ca417308b3ece8a7f47e916c42376e2c8635a031ce26dc815cd7330f2323054a44d23951e416a9a29c5a9a2ab51e8953',
            _attributes[5]+'.gz': '8783aaa65a281c4a7ab3a82a6dc99620418ed2be4a739f46db8ee304fcb3536a78fed5a955e1c373a20c3e7d3673793157c792b4429ecb5c68ddaddb1a0f7de7',
            _attributes[6]+'.gz': '8007480fc607caf41f0f9a93beb66c7caa66c37a3d01a809f6b94bc0df469cec72091e8cc0fbabb3bd8775e9776b928ecda2779fc545c7e4b9e71c504f9510ce',
            _attributes[7]+'.gz': 'afc2ad3d31366a8f7fe1604fa49c190ade6d46bc8915f30bd20fdfdfc663c979bb10ca55ad10cadec6002a17add46639c70e7adf89cb66c57ed004fd3e4f0051',
            _attributes[8]+'.gz': '81c673fe1940e25a6a9722fd74b16bc30e1590db0c40810f541ad4ffba7ae04c01268b929d4bf944e84095a0c2a1d0079d1861bc1df3e8308fbb6b35e0aaf107',
            _attributes[9]+'.gz': '599e4de4935d2ba649c0b526d8aeef6f0e3bf0ed9ee20adad65cb86b078ac139e4cc9758945c2bb6da1c6840867239c5415cb5bceeb80164798ff627aac0a985'
            }



    _ext = '.json.gz'



    ## Cache constructor
    #
    # @param self The object pointer
    # @param inputPath The path to the folder that contains all the input/output files required
    # @param db Mode of storing objects ('file' or 'redis')
    def __init__(self, db='file', print_infos=False):


        self.store_mode = db
        rpCache._db_timeout = 10

        self.dirname = os_path.dirname(os_path.abspath( __file__ ))#+"/.."
        # input_cache
        self._input_cache_dir = self.dirname+'/input_cache/'
        # cache
        self._cache_dir = self.dirname+'/cache/'

        if self.store_mode!='file':
            self.redis = StrictRedis(host=self.store_mode, port=6379, db=0, decode_responses=True)
            if not wait_for_redis(self.redis, self._db_timeout):
                # rpCache.logger.critical("Database "+self.store_mode+" is not reachable")
                # rpCache.logger.info("Trying local redis...")
                # self.redis = StrictRedis(host='localhost', port=6379, db=0, decode_responses=True)
                # if not wait_for_redis(self.redis, self._db_timeout):
                #     rpCache.logger.critical("Database on localhost is not reachable")
                #     rpCache.logger.info("Start local redis...")
                #     p1 = Popen([redis_server.REDIS_SERVER_PATH], stdout=PIPE)
                #     self.redis = StrictRedis(host='localhost', port=6379, db=0, decode_responses=True)
                #     if not wait_for_redis(self.redis, self._db_timeout):
                #         rpCache.logger.critical("Database on localhost is not reachable")
                #         exit()
                rpCache.logger.critical("Database "+self.store_mode+" is not reachable")
                exit()
            self.deprecatedMNXM_mnxm = CRedisDict('deprecatedMNXM_mnxm', self.redis)
            self.deprecatedMNXR_mnxr = CRedisDict('deprecatedMNXR_mnxr', self.redis)
            self.mnxm_strc = CRedisDict('mnxm_strc', self.redis)
            self.chemXref = CRedisDict('chemXref', self.redis)
            self.rr_reactions = CRedisDict('rr_reactions', self.redis)
            self.chebi_mnxm = CRedisDict('chebi_mnxm', self.redis)
            # rpReader attributes
            self.inchikey_mnxm = CRedisDict('inchikey_mnxm', self.redis)
            self.compXref = CRedisDict('compXref', self.redis)
            self.name_compXref = CRedisDict('name_compXref', self.redis)
            # rpCofactors attributes
            self.full_reactions = CRedisDict('full_reactions', self.redis)
        else:
            self.deprecatedMNXM_mnxm = None
            self.deprecatedMNXR_mnxr = None
            self.mnxm_strc = None
            self.chemXref = None
            self.rr_reactions = None
            self.chebi_mnxm = None
            # rpReader attributes
            self.inchikey_mnxm = None
            self.compXref = None
            self.name_compXref = None
            # rpCofactors attributes
            self.full_reactions = None



        self.print = print_infos

        try:
            if self.store_mode=='file':
                self._check_or_load_cache_in_memory(self._cache_dir)
            else:
                self._check_or_load_cache_in_db(self._cache_dir)
        except FileNotFoundError:
            print_FAILED()
            try:
                rpCache._check_or_download_cache_to_disk(self._cache_dir)
                if self.store_mode=='file':
                    self._check_or_load_cache_in_memory(self._cache_dir)
                else:
                    self._check_or_load_cache_in_db(self._cache_dir)
            except (r_exceptions.RequestException,
                    r_exceptions.InvalidSchema,
                    r_exceptions.ConnectionError):
                print_FAILED()
                rpCache.generate_cache(self._cache_dir)
                if self.store_mode=='file':
                    self._check_or_load_cache_in_memory(self._cache_dir)
                else:
                    self._check_or_load_cache_in_db(self._cache_dir)


    #####################################################
    ################# ERROR functions ###################
    #####################################################

    ## Error function for the convertion of structures
    #
    class Error(Exception):
        pass


    ## Error function for the convertion of structures
    #
    class DepictionError(Error):
        def __init__(self, message):
            #self.expression = expression
            self.message = message

    #url = 'https://www.metanetx.org/cgi-bin/mnxget/mnxref/'
    #url = 'ftp://ftp.vital-it.ch/databases/metanetx/MNXref/3.2/'

    @staticmethod
    def _check_or_download_cache_to_disk(cache_dir):
        for attr in rpCache._attributes:
            filename = attr+rpCache._ext
            if os_path.isfile(cache_dir+filename) and sha512(Path(cache_dir+filename).read_bytes()).hexdigest()==rpCache._cache_files[attr]:
                print(filename+" already downloaded ", end = '', flush=True)
                print_OK()
            else:
                filename = attr+rpCache._ext
                print("Downloading "+filename+"...", end = '', flush=True)
                start_time = time.time()
                if not os_path.isdir(cache_dir):
                    os_mkdir(cache_dir)
                download(rpCache._cache_url+filename, cache_dir+filename)
                rpCache._cache_files[attr] = True
                end_time = time.time()
                print_OK(end_time-start_time)


    def _check_or_load_cache_in_memory(self, cache_dir):
        for attribute in rpCache._attributes:
            if not getattr(self, attribute):
                filename = attribute+rpCache._ext
                print("Loading "+filename+"...", end = '', flush=True)
                data = self._load_cache_from_file(cache_dir+filename)
                print_OK()
                setattr(self, attribute, data)
            else:
                print(attribute+" already loaded in memory...", end = '', flush=True)
                print_OK()

    def _check_or_load_cache_in_db(self, cache_dir):
        for attribute in rpCache._attributes:
            if not CRedisDict.exists(self.redis, attribute):
                filename = attribute+rpCache._ext
                print("Loading "+filename+"...", end = '', flush=True)
                data = self._load_cache_from_file(cache_dir+filename)
                print_OK()
                self._store_cache_to_db(attribute, data)
            else:
                print(attribute+" already loaded in db...", end = '', flush=True)
                print_OK()


    @staticmethod
    def generate_cache(outdir):

        if not os_path.isdir(outdir):
            os_mkdir(outdir)

        url = rpCache._cache_url

        # FETCH INPUT_CACHE FILES
        input_dir = 'input-'+os_path.basename(os_path.normpath(outdir))+'/'
        for file in rpCache._input_cache_files.keys():
            rpCache._download_input_cache(url, file, input_dir)

        # GENERATE CACHE FILES AND STORE THEM TO DISK
        attribute = 'deprecatedMNXM_mnxm'
        print(c_attr('bold')+attribute+c_attr('reset'))
        deprecatedMNXM_mnxm = None
        f_deprecatedMNXM_mnxm = outdir+attribute+rpCache._ext
        if not os_path.isfile(f_deprecatedMNXM_mnxm):
            print("   Generating data...", end = '', flush=True)
            deprecatedMNXM_mnxm = rpCache._m_deprecatedMNXM_mnxm(input_dir+'chem_xref.tsv.gz')
            print_OK()
            print("   Writing data to file...", end = '', flush=True)
            rpCache._store_cache_to_file(deprecatedMNXM_mnxm, f_deprecatedMNXM_mnxm)
            print_OK()
        else:
            print("   Cache file already exists", end = '', flush=True)
            print_OK()

        attribute = 'mnxm_strc'
        print(c_attr('bold')+attribute+c_attr('reset'))
        mnxm_strc = None
        f_mnxm_strc = outdir+attribute+rpCache._ext
        if not os_path.isfile(f_mnxm_strc):
            if not deprecatedMNXM_mnxm:
                print("   Loading input data from file...", end = '', flush=True)
                deprecatedMNXM_mnxm = rpCache._load_cache_from_file(f_deprecatedMNXM_mnxm)
                print_OK()
            print("   Generating data...", end = '', flush=True)
            mnxm_strc = rpCache._m_mnxm_strc(input_dir+'/compounds.tsv.gz', input_dir+'chem_prop.tsv.gz', deprecatedMNXM_mnxm)
            print_OK()
            print("   Writing data to file...", end = '', flush=True)
            rpCache._store_cache_to_file(mnxm_strc, f_mnxm_strc)
            print_OK()
        else:
            print("   Cache file already exists", end = '', flush=True)
            print_OK()

        attribute = 'inchikey_mnxm'
        print(c_attr('bold')+attribute+c_attr('reset'))
        inchikey_mnxm = None
        f_inchikey_mnxm = outdir+attribute+rpCache._ext
        if not os_path.isfile(f_inchikey_mnxm):
            if not mnxm_strc:
                print("   Loading input data from file...", end = '', flush=True)
                mnxm_strc = rpCache._load_cache_from_file(f_mnxm_strc)
                print_OK()
            print("   Generating data...", end = '', flush=True)
            inchikey_mnxm = rpCache._m_inchikey_mnxm(mnxm_strc)
            print_OK()
            del mnxm_strc
            print("   Writing data to file...", end = '', flush=True)
            rpCache._store_cache_to_file(inchikey_mnxm, f_inchikey_mnxm)
            print_OK()
        else:
            print("   Cache file already exists", end = '', flush=True)
            print_OK()

        attribute = 'chemXref'
        print(c_attr('bold')+attribute+c_attr('reset'))
        chemXref = None
        f_chemXref = outdir+attribute+rpCache._ext
        if not os_path.isfile(f_chemXref):
            if not deprecatedMNXM_mnxm:
                print("   Loading input data from file...", end = '', flush=True)
                deprecatedMNXM_mnxm = rpCache._load_cache_from_file(f_deprecatedMNXM_mnxm)
                print_OK()
            print("   Generating data...", end = '', flush=True)
            chemXref = rpCache._m_chemXref(input_dir+'chem_xref.tsv.gz', deprecatedMNXM_mnxm)
            print_OK()
            print("   Writing data to file...", end = '', flush=True)
            rpCache._store_cache_to_file(chemXref, f_chemXref)
            print_OK()
        else:
            print("   Cache file already exists", end = '', flush=True)
            print_OK()

        attribute = 'chebi_mnxm'
        print(c_attr('bold')+attribute+c_attr('reset'))
        chebi_mnxm = None
        f_chebi_mnxm = outdir+attribute+rpCache._ext
        if not os_path.isfile(f_chebi_mnxm):
            print("   Generating data...", end = '', flush=True)
            chebi_mnxm = rpCache._m_chebi_mnxm(chemXref)
            print_OK()
            del chemXref
            print("   Writing data to file...", end = '', flush=True)
            rpCache._store_cache_to_file(chebi_mnxm, f_chebi_mnxm)
            del chebi_mnxm
            print_OK()
        else:
            print("   Cache file already exists", end = '', flush=True)
            print_OK()

        attribute = 'deprecatedMNXR_mnxr'
        print(c_attr('bold')+attribute+c_attr('reset'))
        deprecatedMNXR_mnxr = None
        f_deprecatedMNXR_mnxr = outdir+attribute+rpCache._ext
        if not os_path.isfile(f_deprecatedMNXR_mnxr):
            print("   Generating data...", end = '', flush=True)
            deprecatedMNXR_mnxr = rpCache._m_deprecatedMNXR_mnxr(input_dir+'reac_xref.tsv.gz')
            print_OK()
            print("   Writing data to file...", end = '', flush=True)
            rpCache._store_cache_to_file(deprecatedMNXR_mnxr, f_deprecatedMNXR_mnxr)
            print_OK()
        else:
            print("   Cache file already exists", end = '', flush=True)
            print_OK()

        attribute = 'rr_reactions'
        print(c_attr('bold')+attribute+c_attr('reset'))
        rr_reactions = None
        f_rr_reactions = outdir+attribute+rpCache._ext
        if not os_path.isfile(f_rr_reactions):
            if not deprecatedMNXM_mnxm:
                print("   Loading input data from file...", end = '', flush=True)
                deprecatedMNXM_mnxm = rpCache._load_cache_from_file(f_deprecatedMNXM_mnxm)
                print_OK()
            if not deprecatedMNXR_mnxr:
                print("   Loading input data from file...", end = '', flush=True)
                deprecatedMNXR_mnxr = rpCache._load_cache_from_file(f_deprecatedMNXR_mnxr)
                print_OK()
            print("   Generating data...", end = '', flush=True)
            rr_reactions = rpCache._m_rr_reactions(input_dir+'retrorules_rr02_flat_all.tsv.gz', deprecatedMNXM_mnxm, deprecatedMNXR_mnxr)
            print_OK()
            del deprecatedMNXR_mnxr
            print("   Writing data to file...", end = '', flush=True)
            rpCache._store_cache_to_file(rr_reactions, f_rr_reactions)
            print_OK()
            del rr_reactions
        else:
            print("   Cache file already exists", end = '', flush=True)
            print_OK()

        attribute = 'compXref, name_compXref'
        print(c_attr('bold')+attribute+c_attr('reset'))
        compXref = name_compXref = None
        f_compXref = outdir+'compXref'+rpCache._ext
        f_name_compXref = outdir+'name_compXref'+rpCache._ext
        if not os_path.isfile(f_compXref) or not os_path.isfile(f_name_compXref):
            print("   Generating data...", end = '', flush=True)
            compXref,name_compXref = rpCache._m_compXref(input_dir+'comp_xref.tsv.gz')
            print_OK()
            print("   Writing data to file...", end = '', flush=True)
            rpCache._store_cache_to_file(compXref, f_compXref)
            print_OK()
            del compXref
            print("   Writing data to file...", end = '', flush=True)
            rpCache._store_cache_to_file(name_compXref, f_name_compXref)
            print_OK()
            del name_compXref
        else:
            print("   Cache files already exist", end = '', flush=True)
            print_OK()

        attribute = 'full_reactions'
        print(c_attr('bold')+attribute+c_attr('reset'))
        full_reactions = None
        f_full_reactions = outdir+attribute+rpCache._ext
        if not os_path.isfile(f_full_reactions):
            print("   Generating data...", end = '', flush=True)
            if not deprecatedMNXM_mnxm:
                print("   Loading input data from file...", end = '', flush=True)
                deprecatedMNXM_mnxm = rpCache._load_cache_from_file(f_deprecatedMNXM_mnxm)
                print_OK()
            if not deprecatedMNXR_mnxr:
                print("   Loading input data from file...", end = '', flush=True)
                deprecatedMNXR_mnxr = rpCache._load_cache_from_file(f_deprecatedMNXR_mnxr)
                print_OK()
            full_reactions = rpCache._m_full_reactions(input_dir+'rxn_recipes.tsv.gz', deprecatedMNXM_mnxm, deprecatedMNXR_mnxr)
            print_OK()
            print("   Writing data to file...", end = '', flush=True)
            rpCache._store_cache_to_file(full_reactions, f_full_reactions)
            print_OK()
            del full_reactions
        else:
            print("   Cache file already exists", end = '', flush=True)
            print_OK()


    @staticmethod
    def _download_input_cache(url, file, outdir):
        if not os_path.isdir(outdir):
            os_mkdir(outdir)
        filename = outdir+file
        if not os_path.isfile(filename):
            print("Downloading "+file+"...", end = '', flush=True)
            start_time = time.time()
            rpCache.__download_input_cache(url, file, outdir)
            end_time = time.time()
            print_OK(end_time-start_time)
        else:
            print(filename+" already downloaded ", end = '', flush=True)
            print_OK()

    @staticmethod
    def __download_input_cache(url, file, outdir):

        if not os_path.isdir(outdir):
            os_mkdir(outdir)


        # 3xCommon + rpReader
        if file in ['reac_xref.tsv.gz', 'chem_xref.tsv.gz', 'chem_prop.tsv.gz', 'comp_xref.tsv.gz']:
            download(url+'metanetx/'+file, outdir+file)

        #TODO: need to add this file to the git or another location
        if file in ['compounds.tsv.gz', 'rxn_recipes.tsv.gz']:
            download(url+'rr02_more_data/'+file,
                     outdir+file)
            # tar = tarfile_open(outdir+'/rr02_more_data.tar.gz', 'r:gz')
            # tar.extractall(outdir)
            # tar.close()
            # shutil_move(outdir+'/rr02_more_data/compounds.tsv',
            #             outdir+'/rr_compounds.tsv')
            # shutil_move(outdir+'/rr02_more_data/rxn_recipes.tsv',
            #             outdir)
            # os_rm(outdir+'/rr02_more_data.tar.gz')
            # shutil_rmtree(outdir+'/rr02_more_data')

        if file=='retrorules_rr02_flat_all.tsv.gz':
            download(url+'retrorules_rr02_rp3_hs/'+file,
                     outdir+file)
            # download('https://retrorules.org/dl/preparsed/rr02/rp3/hs',
            #          outdir+'/retrorules_rr02_rp3_hs.tar.gz')
            # tar = tarfile_open(outdir+'/retrorules_rr02_rp3_hs.tar.gz', 'r:gz')
            # tar.extractall(outdir)
            # tar.close()
            # shutil_move(outdir+'/retrorules_rr02_rp3_hs/retrorules_rr02_flat_all.tsv', outdir+'/rules_rall.tsv')
            # os_rm(outdir+'/retrorules_rr02_rp3_hs.tar.gz')
            # shutil_rmtree(outdir+'/retrorules_rr02_rp3_hs')




    ##########################################################
    ################## Private Functions #####################
    ##########################################################



    ## Method to load data from file
    #
    #  Load data from file
    #
    #  @param self Object pointer
    #  @param filename File to fetch data from
    #  @return file content
    @staticmethod
    def _load_cache_from_file(filename):
        if filename.endswith('.gz') or filename.endswith('.zip'):
            fp = gzip_open(filename, 'rt', encoding='ascii')
        else:
            fp = open(filename, 'r')
        return json_load(fp)

    ## Method to store data into file
    #
    # Store data into file as json (to store dictionnary structure)
    #
    #  @param self Object pointer
    #  @param data Data to write into file
    #  @param filename File to write data into
    @staticmethod
    def _store_cache_to_file(data, filename):
        if filename.endswith('.gz') or filename.endswith('.zip'):
            fp = gzip_open(filename, 'wt', encoding='ascii')
        else:
            fp = open(filename, 'w')
        json_dump(data, fp)

    ## Method to store data into redis database
    #
    #  Assign a CRedisDict object to the attribute to copy data into the database
    #
    #  @param self Object pointer
    #  @param attr_name Attribute name (database key)
    #  @param data Content of the attribute
    def _store_cache_to_db(self, attr_name, data):
        print("Storing "+attr_name+" to db...", end = '', flush=True)
        setattr(rpCache, attr_name, CRedisDict(attr_name, self.redis, data))
        print_OK()



    ## Function to create a dictionnary of old to new chemical id's
    #
    #  Generate a one-to-one dictionnary of old id's to new ones. Private function
    #
    # TODO: check other things about the mnxm emtry like if it has the right structure etc...
    @staticmethod
    def _checkMNXMdeprecated(mnxm, deprecatedMNXM_mnxm):
        try:
            return deprecatedMNXM_mnxm[mnxm]
        except (KeyError, TypeError):
            return mnxm


    ## Function to create a dictionnary of old to new reaction id's
    #
    # TODO: check other things about the mnxm emtry like if it has the right structure etc...
    @staticmethod
    def _checkMNXRdeprecated(mnxr, deprecatedMNXR_mnxr):
        try:
            return deprecatedMNXR_mnxr[mnxr]
        except (KeyError, TypeError):
            return mnxr


    #################################################################
    ################## Public functions #############################
    #################################################################


    ## Function to parse the chem_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of old to new MetanetX identifiers to make sure that we always use the freshest id's.
    # This can include more than one old id per new one and thus returns a dictionnary. Private function
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return Dictionnary of identifiers
    #TODO: save the self.deprecatedMNXM_mnxm to be used in case there rp_paths uses an old version of MNX
    @staticmethod
    def _deprecatedMNX(xref_path):
        deprecatedMNX_mnx = {}
        with gzip_open(xref_path, 'rt') as f:
            c = csv_reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnx = row[0].split(':')
                    if mnx[0]=='deprecated':
                        deprecatedMNX_mnx[mnx[1]] = row[1]
        return deprecatedMNX_mnx

    @staticmethod
    def _m_deprecatedMNXM_mnxm(chem_xref_path):
        deprecatedMNXM_mnxm = {}
        deprecatedMNXM_mnxm = rpCache._deprecatedMNX(chem_xref_path)
        deprecatedMNXM_mnxm.update(rpCache._convertMNXM)
        deprecatedMNXM_mnxm['MNXM01'] = 'MNXM1'
        return deprecatedMNXM_mnxm

    ## Function to parse the reac_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of old to new MetanetX identifiers to make sure that we always use the freshest id's.
    # This can include more than one old id per new one and thus returns a dictionnary. Private function
    #
    #  @param self Object pointer
    #  @param reac_xref_path Input file path
    #  @return Dictionnary of identifiers
    @staticmethod
    def _m_deprecatedMNXR_mnxr(reac_xref_path):
        return rpCache._deprecatedMNX(reac_xref_path)


    ## Convert chemical depiction to others type of depictions
    #
    # Usage example:
    # - convert_depiction(idepic='CCO', otype={'inchi', 'smiles', 'inchikey'})
    # - convert_depiction(idepic='InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', itype='inchi', otype={'inchi', 'smiles', 'inchikey'})
    #  @param self The object pointer
    #  @param idepic String depiction to be converted, str
    #  @param itype type of depiction provided as input, str
    #  @param otype types of depiction to be generated, {"", "", ..}
    #  @return odepic generated depictions, {"otype1": "odepic1", ..}
    @staticmethod
    def _convert_depiction(idepic, itype='smiles', otype={'inchikey'}):
        # Import (if needed)
        if itype == 'smiles':
            rdmol = MolFromSmiles(idepic, sanitize=True)
        elif itype == 'inchi':
            rdmol = MolFromInchi(idepic, sanitize=True)
        else:
            raise NotImplementedError('"{}" is not a valid input type'.format(itype))
        if rdmol is None:  # Check imprt
            raise rpCache.DepictionError('Import error from depiction "{}" of type "{}"'.format(idepic, itype))
        # Export
        odepic = dict()
        for item in otype:
            if item == 'smiles':
                odepic[item] = MolToSmiles(rdmol)  # MolToSmiles is tricky, one mays want to check the possible options..
            elif item == 'inchi':
                odepic[item] = MolToInchi(rdmol)
            elif item == 'inchikey':
                odepic[item] = MolToInchiKey(rdmol)
            else:
                raise NotImplementedError('"{}" is not a valid output type'.format(otype))
        return odepic


    ## Function to parse the chemp_prop.tsv file from MetanetX and compounds.tsv from RetroRules. Uses the InchIkey as key to the dictionnary
    #
    #  Generate a dictionnary gaving the formula, smiles, inchi and inchikey for the components
    #
    #  @param self Object pointer
    #  @param chem_prop_path Input file path
    #  @return mnxm_strc Dictionnary of formula, smiles, inchi and inchikey
    @staticmethod
    def _m_mnxm_strc(rr_compounds_path, chem_prop_path, deprecatedMNXM_mnxm):
        mnxm_strc = {}
        for row in csv_DictReader(gzip_open(rr_compounds_path, 'rt'), delimiter='\t'):
            tmp = {'formula':  None,
                    'smiles': None,
                    'inchi': row['inchi'],
                    'inchikey': None,
                    'mnxm': rpCache._checkMNXMdeprecated(row['cid'], deprecatedMNXM_mnxm),
                    'name': None}
            try:
                resConv = rpCache._convert_depiction(idepic=tmp['inchi'], itype='inchi', otype={'smiles','inchikey'})
                for i in resConv:
                    tmp[i] = resConv[i]
            except rpCache.DepictionError as e:
                rpCache.logger.warning('Could not convert some of the structures: '+str(tmp))
                rpCache.logger.warning(e)
            mnxm_strc[tmp['mnxm']] = tmp
        with gzip_open(chem_prop_path, 'rt') as f:
            c = csv_reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnxm = rpCache._checkMNXMdeprecated(row[0], deprecatedMNXM_mnxm)
                    tmp = {'formula':  row[2],
                            'smiles': row[6],
                            'inchi': row[5],
                            'inchikey': row[8],
                            'mnxm': mnxm,
                            'name': row[1]}
                    for i in tmp:
                        if tmp[i]=='' or tmp[i]=='NA':
                            tmp[i] = None
                    if mnxm in mnxm_strc:
                        mnxm_strc[mnxm]['formula'] = row[2]
                        mnxm_strc[mnxm]['name'] = row[1]
                        if not mnxm_strc[mnxm]['smiles'] and tmp['smiles']:
                            mnxm_strc[mnxm]['smiles'] = tmp['smiles']
                        if not mnxm_strc[mnxm]['inchikey'] and tmp['inchikey']:
                            mnxm_strc[mnxm]['inchikey'] = tmp['inchikey']
                    else:
                        #check to see if the inchikey is valid or not
                        otype = set({})
                        if not tmp['inchikey']:
                            otype.add('inchikey')
                        if not tmp['smiles']:
                            otype.add('smiles')
                        if not tmp['inchi']:
                            otype.add('inchi')
                        itype = ''
                        if tmp['inchi']:
                            itype = 'inchi'
                        elif tmp['smiles']:
                            itype = 'smiles'
                        else:
                            rpCache.logger.warning('No valid entry for the convert_depiction function')
                            continue
                        try:
                            resConv = rpCache._convert_depiction(idepic=tmp[itype], itype=itype, otype=otype)
                            for i in resConv:
                                tmp[i] = resConv[i]
                        except rpCache.DepictionError as e:
                            rpCache.logger.warning('Could not convert some of the structures: '+str(tmp))
                            rpCache.logger.warning(e)
                        mnxm_strc[tmp['mnxm']] = tmp
        return mnxm_strc


    ## Function to parse the chem_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of all cross references for a given chemical id (MNX) to other database id's
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return a The dictionnary of identifiers
    #TODO: save the self.deprecatedMNXM_mnxm to be used in case there rp_paths uses an old version of MNX
    @staticmethod
    def _m_chemXref(chem_xref_path, deprecatedMNXM_mnxm):
        chemXref = {}
        with gzip_open(chem_xref_path, 'rt') as f:
            c = csv_reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnx = rpCache._checkMNXMdeprecated(row[1], deprecatedMNXM_mnxm)
                    if len(row[0].split(':'))==1:
                        dbName = 'mnx'
                        dbId = row[0]
                    else:
                        dbName = row[0].split(':')[0]
                        dbId = ''.join(row[0].split(':')[1:])
                        if dbName=='deprecated':
                            dbName = 'mnx'
                    #mnx
                    if not mnx in chemXref:
                        chemXref[mnx] = {}
                    if not dbName in chemXref[mnx]:
                        chemXref[mnx][dbName] = []
                    if not dbId in chemXref[mnx][dbName]:
                        chemXref[mnx][dbName].append(dbId)
                    ### DB ###
                    if not dbName in chemXref:
                        chemXref[dbName] = {}
                    if not dbId in chemXref[dbName]:
                        chemXref[dbName][dbId] = mnx
        return chemXref


    ## Function to parse the chem_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of all cross references for a given chemical id (MNX) to other database id's
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return a The dictionnary of identifiers
    #TODO: save the self.deprecatedMNXM_mnxm to be used in case there rp_paths uses an old version of MNX
#    def _m_chebi_mnxm(self, chemXref):
    @staticmethod
    def _m_chebi_mnxm(chemXref):
        chebi_mnxm = {}
        for mnxm in chemXref:
            if 'chebi' in chemXref[mnxm]:
                for c in chemXref[mnxm]['chebi']:
                    chebi_mnxm[c] = mnxm
        return chebi_mnxm


    ## Function to parse the rules_rall.tsv from RetroRules
    #
    #  Extract from the reactions rules the ruleID, the reactionID, the direction of the rule directed to the origin reaction
    #
    #  @param self The object pointer.
    #  @param path The input file path.
    #  @return rule Dictionnary describing each reaction rule
    @staticmethod
    def _m_rr_reactions(rules_rall_path, deprecatedMNXM_mnxm, deprecatedMNXR_mnxr):
        rr_reactions = {}
        try:
            #with gzip_open(rules_rall_path, 'r') as f:
            #    reader = csv.reader(f, delimiter = '\t')
            #    next(reader)
            #    rule = {}
            #    for row in reader:
            for row in csv_DictReader(gzip_open(rules_rall_path, 'rt'), delimiter='\t'):
                #NOTE: as of now all the rules are generated using MNX
                #but it may be that other db are used, we are handling this case
                #WARNING: can have multiple products so need to seperate them
                products = {}
                for i in row['Product_IDs'].split('.'):
                    mnxm = rpCache._checkMNXMdeprecated(i, deprecatedMNXM_mnxm)
                    if not mnxm in products:
                        products[mnxm] = 1
                    else:
                        products[mnxm] += 1
                try:
                    #WARNING: one reaction rule can have multiple reactions associated with them
                    #To change when you can set subpaths from the mutliple numbers of
                    #we assume that the reaction rule has multiple unique reactions associated
                    if row['# Rule_ID'] not in rr_reactions:
                        rr_reactions[row['# Rule_ID']] = {}
                    if row['# Rule_ID'] in rr_reactions[row['# Rule_ID']]:
                        rpCache.logger.warning('There is already reaction '+str(row['# Rule_ID'])+' in reaction rule '+str(row['# Rule_ID']))
                    rr_reactions[row['# Rule_ID']][row['Reaction_ID']] = {
                        'rule_id': row['# Rule_ID'],
                        'rule_score': float(row['Score_normalized']),
                        'reac_id': rpCache._checkMNXRdeprecated(row['Reaction_ID'], deprecatedMNXR_mnxr),
                        'subs_id': rpCache._checkMNXMdeprecated(row['Substrate_ID'], deprecatedMNXM_mnxm),
                        'rel_direction': int(row['Rule_relative_direction']),
                        'left': {rpCache._checkMNXMdeprecated(row['Substrate_ID'], deprecatedMNXM_mnxm): 1},
                        'right': products}
                except ValueError:
                    rpCache.logger.error('Problem converting rel_direction: '+str(row['Rule_relative_direction']))
                    rpCache.logger.error('Problem converting rule_score: '+str(row['Score_normalized']))
            return rr_reactions
        except FileNotFoundError as e:
                rpCache.logger.error('Could not read the rules_rall file ('+str(rules_rall_path)+')')
                return {}


    @staticmethod
    def _m_inchikey_mnxm(mnxm_strc):
        inchikey_mnxm = {}
        for mnxm in mnxm_strc:
            inchikey = mnxm_strc[mnxm]['inchikey']
            if not inchikey: inchikey = 'NO_INCHIKEY'
            if not inchikey in inchikey_mnxm:
                inchikey_mnxm[inchikey] = []
            inchikey_mnxm[inchikey].append(mnxm)
        return inchikey_mnxm

    # rpReader
    ## Function to parse the compXref.tsv file of MetanetX
    #
    #  Generate a dictionnary of compartments id's (MNX) to other database id's
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return a The dictionnary of identifiers
    #TODO: save the self.deprecatedMNXM_mnxm to be used in case there rp_paths uses an old version of MNX
    @staticmethod
    def _m_compXref(compXref_path):
        compXref = {}
        name_compXref = {}
        try:
            with gzip_open(compXref_path, 'rt') as f:
                c = csv_reader(f, delimiter='\t')
                #not_recognised = []
                for row in c:
                    #cid = row[0].split(':')
                    if not row[0][0]=='#':
                        #collect the info
                        mnxc = row[1]
                        if len(row[0].split(':'))==1:
                            dbName = 'mnx'
                            dbCompId = row[0]
                        else:
                            dbName = row[0].split(':')[0]
                            dbCompId = ''.join(row[0].split(':')[1:])
                            dbCompId = dbCompId.lower()
                        if dbName=='deprecated':
                            dbName = 'mnx'
                        #create the dicts
                        if not mnxc in compXref:
                            compXref[mnxc] = {}
                        if not dbName in compXref[mnxc]:
                            compXref[mnxc][dbName] = []
                        if not dbCompId in compXref[mnxc][dbName]:
                            compXref[mnxc][dbName].append(dbCompId)
                        #create the reverse dict
                        if not dbCompId in name_compXref:
                            name_compXref[dbCompId] = mnxc
        except FileNotFoundError:
            rpCache.logger.error('compXref file not found')
            return {}
        return compXref,name_compXref


    ## Generate complete reactions from the rxn_recipes.tsv from RetroRules
    #
    #  These are the compplete reactions from which the reaction rules are generated from. This is used to
    # reconstruct the full reactions from monocomponent reactions
    #
    #  @param self The pointer object
    #  @param rxn_recipes_path Path to the recipes file
    #  @return Boolean that determines the success or failure of the function
    @staticmethod
    def _m_full_reactions(rxn_recipes_path, deprecatedMNXM_mnxm, deprecatedMNXR_mnxr):
        #### for character matching that are returned
        DEFAULT_STOICHIO_RESCUE = {"4n": 4, "3n": 3, "2n": 2, 'n': 1,
                           '(n)': 1, '(N)': 1, '(2n)': 2, '(x)': 1,
                           'N': 1, 'm': 1, 'q': 1,
                           '0.01': 1, '0.1': 1, '0.5': 1, '1.5': 1,
                           '0.02': 1, '0.2': 1,
                           '(n-1)': 0, '(n-2)': -1}
        reaction = {}
        try:
            for row in csv_DictReader(gzip_open(rxn_recipes_path, 'rt'), delimiter='\t'):
                tmp = {} # makes sure that if theres an error its not added
                #parse the reaction equation
                if not len(row['Equation'].split('='))==2:
                    rpCache.logger.warning('There should never be more or less than a left and right of an equation')
                    rpCache.logger.warnin(row['Equation'])
                    continue
                ######### LEFT ######
                #### MNX id
                tmp['left'] = {}
                # if row['#Reaction_ID']=="MNXR141948":
                #     print(row)
                #     exit()
                for spe in re_findall(r'(\(n-1\)|\d+|4n|3n|2n|n|\(n\)|\(N\)|\(2n\)|\(x\)|N|m|q|\(n\-2\)|\d+\.\d+) ([\w\d]+)@\w+', row['Equation'].split('=')[0]):
                    #1) try to rescue if its one of the values
                    try:
                        tmp['left'][rpCache._checkMNXMdeprecated(spe[1], deprecatedMNXM_mnxm)] = DEFAULT_STOICHIO_RESCUE[spe[0]]
                    except KeyError:
                        #2) try to convert to int if its not
                        try:
                            tmp['left'][rpCache._checkMNXMdeprecated(spe[1], deprecatedMNXM_mnxm)] = int(spe[0])
                        except ValueError:
                            rpCache.logger.warning('Cannot convert '+str(spe[0]))
                            continue
                ####### RIGHT #####
                ####  MNX id
                tmp['right'] = {}
                for spe in re_findall(r'(\(n-1\)|\d+|4n|3n|2n|n|\(n\)|\(N\)|\(2n\)|\(x\)|N|m|q|\(n\-2\)|\d+\.\d+) ([\w\d]+)@\w+', row['Equation'].split('=')[1]):
                    #1) try to rescue if its one of the values
                    try:
                        tmp['right'][rpCache._checkMNXMdeprecated(spe[1], deprecatedMNXM_mnxm)] = DEFAULT_STOICHIO_RESCUE[spe[0]]
                    except KeyError:
                        #2) try to convert to int if its not
                        try:
                            tmp['right'][rpCache._checkMNXMdeprecated(spe[1], deprecatedMNXM_mnxm)] = int(spe[0])
                        except ValueError:
                            rpCache.logger.warning('Cannot convert '+str(spe[0]))
                            continue
                ####### DIRECTION ######
                try:
                    tmp['direction'] = int(row['Direction'])
                except ValueError:
                    rpCache.logger.error('Cannot convert '+str(row['Direction'])+' to int')
                    continue
                ### add the others
                tmp['main_left'] = row['Main_left'].split(',')
                tmp['main_right'] = row['Main_right'].split(',')
                reaction[rpCache._checkMNXRdeprecated(row['#Reaction_ID'], deprecatedMNXR_mnxr)] = tmp
            return reaction
        except FileNotFoundError:
            rpCache.logger.error('Cannot find file: '+str(rxn_recipes_path))
            return False
