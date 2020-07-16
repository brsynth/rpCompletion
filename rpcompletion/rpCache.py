from os import path as os_path
from os import mkdir as os_mkdir
from os import remove as os_rm
from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi, MolToInchiKey, AddHs
from csv import DictReader as csv_DictReader
from csv import reader as csv_reader
from logging import getLogger as logging_getLogger
from pickle import load as pickle_load
from pickle import loads as pickle_loads
from pickle import dumps as pickle_dumps
from gzip import open as gzip_open
from urllib.request import urlretrieve as urllib_request_urlretrieve
from re import findall as re_findall
from tarfile import open as tarfile_open
from shutil import move as shutil_move
from shutil import rmtree as shutil_rmtree
import sys
import time
from itertools import chain as itertools_chain
from brs_utils import print_OK, print_FAILED, download_and_extract_gz
from tarfile import open as tf_open
from redis import StrictRedis
from credisdict import CRedisDict, wait_for_redis
import redis_server
from subprocess import run as proc_run
from subprocess import Popen,PIPE

#######################################################
################### rpCache  ##########################
#######################################################


def add_arguments(parser):
    parser.add_argument('-sm', '--store_mode', type=str, default='file',
                        help='data storage mode: file or db')
    parser.add_argument('-p', '--print', type=bool, default=False,
                        help='print additional informations')
    return parser

def build_parser():
    return add_arguments(argparse_ArgumentParser('Python script to pre-compute data'))

def entrypoint(args=sys.argv[1:]):
    parser = build_parser()

    params = parser.parse_args(args)

    rpcache = rpCache(params.store_mode, params.print)

##
#
#
if __name__ == "__main__":
    entrypoint()

## Class to generate the cache
#
# Contains all the functions that parse different files, used to calculate the thermodynamics and the FBA of the
#the other steps. These should be called only when the files have changes
class rpCache:



    ## Cache constructor
    #
    # @param self The object pointer
    # @param inputPath The path to the folder that contains all the input/output files required
    # @param db Mode of storing objects ('file' or 'redis')
    def __init__(self, db='file', print_infos=False):

        #given by Thomas
        self.logger = logging_getLogger(__name__)
        self.logger.info('Started instance of rpCache')

        self.store_mode = db
        self._db_timeout = 15

        self.dirname = os_path.dirname(os_path.abspath( __file__ ))#+"/.."
        # input_cache
        self.input_cache_dir = self.dirname+'/input_cache/'

        # cache
        self.cache_dir = self.dirname+'/cache/'
        self._cache_url = 'https://github.com/brsynth/rpCache-data/raw/master/'
        if not os_path.isdir(self.cache_dir):
            os_mkdir(self.cache_dir)

        if self.store_mode!='file':
            self.redis = StrictRedis(host=self.store_mode, port=6379, db=0, decode_responses=True)
            if not wait_for_redis(self.redis, self._db_timeout):
                self.logger.critical("Database "+self.store_mode+" is not reachable")
                self.logger.info("Trying local redis...")
                self.redis = StrictRedis(host='localhost', port=6379, db=0, decode_responses=True)
                if not wait_for_redis(self.redis, self._db_timeout):
                    self.logger.critical("Database on localhost is not reachable")
                    self.logger.info("Start local redis...")
                    p1 = Popen([redis_server.REDIS_SERVER_PATH], stdout=PIPE)
                    self.redis = StrictRedis(host='localhost', port=6379, db=0, decode_responses=True)
                    if not wait_for_redis(self.redis, self._db_timeout):
                        self.logger.critical("Database on localhost is not reachable")
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

        # Common attribues
        self.convertMNXM = {'MNXM162231': 'MNXM6',
                            'MNXM84': 'MNXM15',
                            'MNXM96410': 'MNXM14',
                            'MNXM114062': 'MNXM3',
                            'MNXM145523': 'MNXM57',
                            'MNXM57425': 'MNXM9',
                            'MNXM137': 'MNXM588022'}


        if not self.loadCache():
            raise ValueError


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


    ##########################################################
    ################## Private Functions #####################
    ##########################################################

    ## Private function to fetch the required data, parse them and generate the pickle
    #
    #  Opens the previously generated cache to the object memory
    #
    # @param The oject pointer
    # @return Boolean detemining the success of the function or not
    def loadCache(self, fetchInputFiles=False):

        attr_names = {
        #   KEY                  : [attribute(s) name(s) list, args list to the function]
            'deprecatedMNXM_mnxm': [['deprecatedMNXM_mnxm'], ['chem_xref.tsv']],
            'deprecatedMNXR_mnxr': [['deprecatedMNXR_mnxr'], ['reac_xref.tsv']],
            'mnxm_strc': [['mnxm_strc'], ['rr_compounds.tsv', 'chem_prop.tsv']],
            'chemXref': [['chemXref'], ['chem_xref.tsv']],
            'chebi_mnxm': [['chebi_mnxm'], []],
            'rr_reactions': [['rr_reactions'], ['rules_rall.tsv']],
            'inchikey_mnxm': [['inchikey_mnxm'], []],
            'compXref': [['compXref', 'name_compXref'], ['comp_xref.tsv']],
            'full_reactions': [['full_reactions'], ['rxn_recipes.tsv']]
        }

        start_time = end_time = 0
        # For each attribute name
        for attr_name in attr_names:
            # For each real attribute
            for attr in attr_names[attr_name][0]:
                # If cache is not  loaded
                if not self.cache_loaded(attr):
                    # LOAD CACHE
                    self.check_and_fetch_cache(attr)
                    data = self.load_cache_from_file(self.cache_dir+attr+'.pickle')
                    self.store_cache_to_db(attr, data)
                    # LOAD INPUT_CACHE
                    # self.check_and_fetch_input_cache(attr)
                    # data = self.gen_cache(attributes[0], [self.input_cache_dir+input_file for input_file in attributes[1]])
                    # self.populate_cache(attr_names[attr_name], data)
                    # start_time = time.time()
                    # end_time = time.time()
                    # print_OK(end_time-start_time)

                elif not self.store_mode=='file':

                    print(" ".join(attr_names[attr_name][0])+" already in db ", end = '', flush=True)
                    print_OK()

            # Load cache from file
            if self.store_mode=='file':
                for i in range(len(attr_names[attr_name][0])):
                    _attr_name = attr_names[attr_name][0][i]
                    filename = self.cache_dir+_attr_name+'.pickle'
                    data = self.load_cache_from_file(filename)
                    setattr(self, _attr_name, data)

        return True

    def check_and_fetch_cache(self, attr):
        filename = self.cache_dir+attr
        if not os_path.exists(filename+'.pickle'):
            print("Downloading "+attr+"...", end = '', flush=True)
            start_time = time.time()
            self.fetch_cache(attr)
            end_time = time.time()
#                            print(" (%.2fs)" % (end_time - start_time))
            print_OK(end_time-start_time)
        else:
            print(filename+" already downloaded ", end = '', flush=True)
            print_OK()

    def check_and_fetch_input_cache(self, attr):
        filename = self.input_cache_dir+attr
        if not os_path.isfile(filename):
            print("Downloading "+attr+"...", end = '', flush=True)
            start_time = time.time()
            self.fetch_input_cache(attr)
            end_time = time.time()
#                            print(" (%.2fs)" % (end_time - start_time))
            print_OK(end_time-start_time)
        else:
            print(filename+" already downloaded ", end = '', flush=True)
            print_OK()

    def cache_loaded(self, attr):
        if self.store_mode=='file':
            return os_path.isfile(self.cache_dir+attr+'.pickle')
        else:
            return CRedisDict.exists(self.redis, attr)


    def fetch_cache(self, attr):
        download_and_extract_gz(self._cache_url+attr+'.pickle.tar.gz', self.cache_dir)


    def fetch_input_cache(self, file):

        if not os_path.isdir(self.input_cache_dir):
            os_mkdir(self.input_cache_dir)

        #url = 'https://www.metanetx.org/cgi-bin/mnxget/mnxref/'
        url = 'ftp://ftp.vital-it.ch/databases/metanetx/MNXref/3.2/'

        # 3xCommon + rpReader
        if file in ['reac_xref.tsv', 'chem_xref.tsv', 'chem_prop.tsv', 'comp_xref.tsv']:
            urllib_request_urlretrieve(url+file, self.input_cache_dir+'/'+file)

        #TODO: need to add this file to the git or another location
        if file in ['rr_compounds.tsv', 'rxn_recipes.tsv']:
            urllib_request_urlretrieve('https://retrorules.org/dl/this/is/not/a/secret/path/rr02',
                                       self.input_cache_dir+'/rr02_more_data.tar.gz')
            tar = tarfile_open(self.input_cache_dir+'/rr02_more_data.tar.gz', 'r:gz')
            tar.extractall(self.input_cache_dir)
            tar.close()
            shutil_move(self.input_cache_dir+'/rr02_more_data/compounds.tsv',
                        self.input_cache_dir+'/rr_compounds.tsv')
            shutil_move(self.input_cache_dir+'/rr02_more_data/rxn_recipes.tsv',
                        self.input_cache_dir)
            os_rm(self.input_cache_dir+'rr02_more_data.tar.gz')
            shutil_rmtree(self.input_cache_dir+'rr02_more_data')

        if file=='rules_rall.tsv':
            urllib_request_urlretrieve('https://retrorules.org/dl/preparsed/rr02/rp3/hs',
                                       self.input_cache_dir+'/retrorules_rr02_rp3_hs.tar.gz')
            tar = tarfile_open(self.input_cache_dir+'/retrorules_rr02_rp3_hs.tar.gz', 'r:gz')
            tar.extractall(self.input_cache_dir)
            tar.close()
            shutil_move(self.input_cache_dir+'/retrorules_rr02_rp3_hs/retrorules_rr02_flat_all.tsv', self.input_cache_dir+'/rules_rall.tsv')
            os_rm(self.input_cache_dir+'/retrorules_rr02_rp3_hs.tar.gz')
            shutil_rmtree(self.input_cache_dir+'/retrorules_rr02_rp3_hs')




    def populate_cache(self, attributes, data):

        for i in range(len(data)):
            _attr_name = attributes[0][i]
            method = getattr(self, 'store_cache_to_'+self.store_mode)
            method(_attr_name, data[i])




    ## Method to generate data to be cached
    #
    #  Call method to generate data according to attribute name(s)
    #
    #  @param self Object pointer
    #  @param attr_names Attribute names. List in case of multiple results from generate method.
    #  @param args Arguments list to pass to the generation method.
    #  @return results Data generated by the specified method
    def gen_cache(self, attr_names, args):
        try:
            results = []
            print("Generating "+" ".join(attr_names)+"...", end = '', flush=True)
            # Choose method according to attribute name
            method = getattr(self, '_m_'+attr_names[0])
            # Apply method and expand 'args' list as arguments
            # Put results in a list
            results = [method(*args)]
            if type(results[0]) is tuple:
                results = list(itertools_chain(results[0]))
            print_OK()
            return results
        except:
            print_FAILED()
            raise


    ## Method to load data from file
    #
    #  Load data from file
    #
    #  @param self Object pointer
    #  @param filename File to fetch data from
    #  @param gzip File is compressed or not
    #  @return file content
    def load_cache_from_file(self, filename, gzip=False):
        print("Loading "+filename+"...", end = '', flush=True)
        if gzip:
            print_OK()
            return pickle_load(gzip_open(filename, 'rb'))
        else:
            print_OK()
            return pickle_load(open(filename, 'rb'))

    ## Method to store data into file
    #
    # Store data into file as pickles (to store dictionnary structure)
    #
    #  @param self Object pointer
    #  @param attr_name Attribute name (filename)
    #  @param data Content of the attribute
    #  @param gzip Compress file or not
    def store_cache_to_file(self, _attr_name, data, gzip=False):
        print("Storing "+_attr_name+" to file...", end = '', flush=True)
        filename = self.cache_dir+'/'+_attr_name+'.pickle'
        pickle_obj = pickle_dumps(data)
        if gzip:
            filename += '.gz'
            with gzip_open(filename, "wb") as f:
            	f.write(pickle_obj)
        else:
            with open(filename, "wb") as f:
             	f.write(pickle_obj)
        print_OK()

    ## Method to store data into redis database
    #
    #  Assign a CRedisDict object to the attribute to copy data into the database
    #
    #  @param self Object pointer
    #  @param attr_name Attribute name (database key)
    #  @param data Content of the attribute
    def store_cache_to_db(self, attr_name, data):
        print("Storing "+attr_name+" to db...", end = '', flush=True)
        setattr(self, attr_name, CRedisDict(attr_name, self.redis, data))
        print_OK()



    ## Function to create a dictionnary of old to new chemical id's
    #
    #  Generate a one-to-one dictionnary of old id's to new ones. Private function
    #
    # TODO: check other things about the mnxm emtry like if it has the right structure etc...
    def _checkMNXMdeprecated(self, mnxm):
        try:
            return self.deprecatedMNXM_mnxm[mnxm]
        except (KeyError, TypeError):
            return mnxm


    ## Function to create a dictionnary of old to new reaction id's
    #
    # TODO: check other things about the mnxm emtry like if it has the right structure etc...
    def _checkMNXRdeprecated(self, mnxr):
        try:
            return self.deprecatedMNXR_mnxr[mnxr]
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
    def _deprecatedMNX(self, xref_path):
        deprecatedMNX_mnx = {}
        with open(xref_path) as f:
            c = csv_reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnx = row[0].split(':')
                    if mnx[0]=='deprecated':
                        deprecatedMNX_mnx[mnx[1]] = row[1]
        return deprecatedMNX_mnx

    def _m_deprecatedMNXM_mnxm(self, chem_xref_path):
        deprecatedMNXM_mnxm = {}
        deprecatedMNXM_mnxm = self._deprecatedMNX(chem_xref_path)
        deprecatedMNXM_mnxm.update(self.convertMNXM)
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
    def _m_deprecatedMNXR_mnxr(self, reac_xref_path):
        return self._deprecatedMNX(reac_xref_path)


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
    def _convert_depiction(self, idepic, itype='smiles', otype={'inchikey'}):
        # Import (if needed)
        if itype == 'smiles':
            rdmol = MolFromSmiles(idepic, sanitize=True)
        elif itype == 'inchi':
            rdmol = MolFromInchi(idepic, sanitize=True)
        else:
            raise NotImplementedError('"{}" is not a valid input type'.format(itype))
        if rdmol is None:  # Check imprt
            raise self.DepictionError('Import error from depiction "{}" of type "{}"'.format(idepic, itype))
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
    def _m_mnxm_strc(self, rr_compounds_path, chem_prop_path):
        mnxm_strc = {}
        for row in csv_DictReader(open(rr_compounds_path), delimiter='\t'):
            tmp = {'formula':  None,
                    'smiles': None,
                    'inchi': row['inchi'],
                    'inchikey': None,
                    'mnxm': self._checkMNXMdeprecated(row['cid']),
                    'name': None}
            try:
                resConv = self._convert_depiction(idepic=tmp['inchi'], itype='inchi', otype={'smiles','inchikey'})
                for i in resConv:
                    tmp[i] = resConv[i]
            except self.DepictionError as e:
                self.logger.warning('Could not convert some of the structures: '+str(tmp))
                self.logger.warning(e)
            mnxm_strc[tmp['mnxm']] = tmp
        with open(chem_prop_path) as f:
            c = csv_reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnxm = self._checkMNXMdeprecated(row[0])
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
                            self.logger.warning('No valid entry for the convert_depiction function')
                            continue
                        try:
                            resConv = self._convert_depiction(idepic=tmp[itype], itype=itype, otype=otype)
                            for i in resConv:
                                tmp[i] = resConv[i]
                        except self.DepictionError as e:
                            self.logger.warning('Could not convert some of the structures: '+str(tmp))
                            self.logger.warning(e)
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
    def _m_chemXref(self, chem_xref_path):
        chemXref = {}
        with open(chem_xref_path) as f:
            c = csv_reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnx = self._checkMNXMdeprecated(row[1])
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
    def _m_chebi_mnxm(self):
        chebi_mnxm = {}
        for mnxm in self.chemXref:
            if 'chebi' in self.chemXref[mnxm]:
                for c in self.chemXref[mnxm]['chebi']:
                    chebi_mnxm[c] = mnxm
        return chebi_mnxm


    ## Function to parse the rules_rall.tsv from RetroRules
    #
    #  Extract from the reactions rules the ruleID, the reactionID, the direction of the rule directed to the origin reaction
    #
    #  @param self The object pointer.
    #  @param path The input file path.
    #  @return rule Dictionnary describing each reaction rule
    def _m_rr_reactions(self, rules_rall_path):
        rr_reactions = {}
        try:
            #with open(rules_rall_path, 'r') as f:
            #    reader = csv.reader(f, delimiter = '\t')
            #    next(reader)
            #    rule = {}
            #    for row in reader:
            for row in csv_DictReader(open(rules_rall_path), delimiter='\t'):
                #NOTE: as of now all the rules are generated using MNX
                #but it may be that other db are used, we are handling this case
                #WARNING: can have multiple products so need to seperate them
                products = {}
                for i in row['Product_IDs'].split('.'):
                    mnxm = self._checkMNXMdeprecated(i)
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
                        self.logger.warning('There is already reaction '+str(row['# Rule_ID'])+' in reaction rule '+str(row['# Rule_ID']))
                    rr_reactions[row['# Rule_ID']][row['Reaction_ID']] = {'rule_id': row['# Rule_ID'], 'rule_score': float(row['Score_normalized']), 'reac_id': self._checkMNXRdeprecated(row['Reaction_ID']), 'subs_id': self._checkMNXMdeprecated(row['Substrate_ID']), 'rel_direction': int(row['Rule_relative_direction']), 'left': {self._checkMNXMdeprecated(row['Substrate_ID']): 1}, 'right': products}
                except ValueError:
                    self.logger.error('Problem converting rel_direction: '+str(row['Rule_relative_direction']))
                    self.logger.error('Problem converting rule_score: '+str(row['Score_normalized']))
            return rr_reactions
        except FileNotFoundError as e:
                self.logger.error('Could not read the rules_rall file ('+str(rules_rall_path)+')')
                return {}


    def _m_inchikey_mnxm(self):
        inchikey_mnxm = {}
        for mnxm in self.mnxm_strc:
            inchikey = self.mnxm_strc[mnxm]['inchikey']
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
    def _m_compXref(self, compXref_path):
        compXref = {}
        name_compXref = {}
        try:
            with open(compXref_path) as f:
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
            self.logger.error('compXref file not found')
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
    def _m_full_reactions(self, rxn_recipes_path):
        #### for character matching that are returned
        DEFAULT_STOICHIO_RESCUE = {"4n": 4, "3n": 3, "2n": 2, 'n': 1,
                           '(n)': 1, '(N)': 1, '(2n)': 2, '(x)': 1,
                           'N': 1, 'm': 1, 'q': 1,
                           '0.01': 1, '0.1': 1, '0.5': 1, '1.5': 1,
                           '0.02': 1, '0.2': 1,
                           '(n-1)': 0, '(n-2)': -1}
        reaction = {}
        try:
            for row in csv_DictReader(open(rxn_recipes_path), delimiter='\t'):
                tmp = {} # makes sure that if theres an error its not added
                #parse the reaction equation
                if not len(row['Equation'].split('='))==2:
                    self.logger.warning('There should never be more or less than a left and right of an equation')
                    self.logger.warnin(row['Equation'])
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
                        tmp['left'][self._checkMNXMdeprecated(spe[1])] = DEFAULT_STOICHIO_RESCUE[spe[0]]
                    except KeyError:
                        #2) try to convert to int if its not
                        try:
                            tmp['left'][self._checkMNXMdeprecated(spe[1])] = int(spe[0])
                        except ValueError:
                            self.logger.warning('Cannot convert '+str(spe[0]))
                            continue
                ####### RIGHT #####
                ####  MNX id
                tmp['right'] = {}
                for spe in re_findall(r'(\(n-1\)|\d+|4n|3n|2n|n|\(n\)|\(N\)|\(2n\)|\(x\)|N|m|q|\(n\-2\)|\d+\.\d+) ([\w\d]+)@\w+', row['Equation'].split('=')[1]):
                    #1) try to rescue if its one of the values
                    try:
                        tmp['right'][self._checkMNXMdeprecated(spe[1])] = DEFAULT_STOICHIO_RESCUE[spe[0]]
                    except KeyError:
                        #2) try to convert to int if its not
                        try:
                            tmp['right'][self._checkMNXMdeprecated(spe[1])] = int(spe[0])
                        except ValueError:
                            self.logger.warning('Cannot convert '+str(spe[0]))
                            continue
                ####### DIRECTION ######
                try:
                    tmp['direction'] = int(row['Direction'])
                except ValueError:
                    self.logger.error('Cannot convert '+str(row['Direction'])+' to int')
                    continue
                ### add the others
                tmp['main_left'] = row['Main_left'].split(',')
                tmp['main_right'] = row['Main_right'].split(',')
                reaction[self._checkMNXRdeprecated(row['#Reaction_ID'])] = tmp
            return reaction
        except FileNotFoundError:
            self.logger.error('Cannot find file: '+str(path))
            return False
