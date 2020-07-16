import csv
import logging
import requests

from itertools import product as itertools_product
from time import time as time_time
from time import sleep as time_sleep
from bisect import insort as bisect_insort
from argparse import ArgumentParser as argparse_ArgumentParser
from os import path as os_path
from os import mkdir as os_mkdir
from .rpCofactors import rpCofactors, add_arguments
from brs_utils import rpSBML
from io import StringIO
from copy import deepcopy

#import rpCofactors

## @package rpCompletion
#
# Collection of functions that convert the outputs from various sources to the SBML format (rpSBML) for further analyses

def build_args_parser():
    parser = argparse_ArgumentParser('Python wrapper to parse RP2 to generate rpSBML collection of unique and complete (cofactors) pathways')
    parser = _add_arguments(parser)

    return parser

def _add_arguments(parser):
    parser = add_arguments(parser)
    parser.add_argument('rp2_pathways', type=str)
    parser.add_argument('rp2paths_compounds', type=str)
    parser.add_argument('rp2paths_pathways', type=str)
    parser.add_argument('outdir', type=str)
    parser.add_argument('--upper_flux_bound', type=int, default=999999)
    parser.add_argument('--lower_flux_bound', type=int, default=0)
    parser.add_argument('--max_subpaths_filter', type=int, default=10)
    parser.add_argument('--pathway_id', type=str, default='rp_pathway')
    parser.add_argument('--compartment_id', type=str, default='MNXC3')
    parser.add_argument('--species_group_id', type=str, default='central_species')
    parser.add_argument('--sink_species_group_id', type=str, default='rp_sink_species')
    parser.add_argument('--pubchem_search', type=str, default='False')
    return parser

class Species:
    def __init__(self, inchi, inchikey, smiles, xref):
        self.inchi = inchi
        self.inchikey = inchikey
        self.smiles = smiles
        self.xref = xref

class SBML_Item:
    def __init__(self, score, index, rpsbml_obj):
        self.score = score
        self.index = index
        self.rpsbml_obj = rpsbml_obj

    def __eq__(self, sbml_item):
        return self.score == sbml_item.score
    def __lt__(self, sbml_item):
        return self.score < sbml_item.score
    def __gt__(self, sbml_item):
        return self.score > sbml_item.score


## Class to read all the input files
#
# Contains all the functions that read the cache files and input files to reconstruct the heterologous pathways
class rpCompletion(rpCofactors):
    ## InputReader constructor
    #
    #  @param self The object pointer
    def __init__(self, db='file'):
        super().__init__(db)
        # self.rpcofactors = rpCofactors(db, self.print)
        # self.logger = logging.getLogger(__name__)
        self.logger.info('Starting instance of rpCompletion')

        self.pubchem_species = {}
        #####################
        #self.pubchem_sec_count = 0
        #self.pubchem_sec_start = 0.0
        self.pubchem_min_count = 0
        self.pubchem_min_start = 0.0


    #######################################################################
    ############################# PRIVATE FUNCTIONS #######################
    #######################################################################

    def _pubChemLimit(self):
        '''
        if self.pubchem_sec_start==0.0:
            self.pubchem_sec_start = time_time()
        '''
        if self.pubchem_min_start==0.0:
            self.pubchem_min_start = time_time()
        #self.pubchem_sec_count += 1
        self.pubchem_min_count += 1
        '''
        #### requests per second ####
        if self.pubchem_sec_count>=5 and time.time()-self.pubchem_sec_start<=1.0:
            time.sleep(1.0)
            self.pubchem_sec_start = time.time()
            self.pubchem_sec_count = 0
        elif time.time()-self.pubchem_sec_start>1.0:
            self.pubchem_sec_start = time.time()
            self.pubchem_sec_count = 0
        '''
        #### requests per minute ####
        if self.pubchem_min_count>=500 and time_time()-self.pubchem_min_start<=60.0:
            logging.warning('Reached 500 requests per minute for pubchem... waiting a minute')
            time_sleep(60.0)
            self.pubchem_min_start = time_time()
            self.pubchem_min_count = 0
        elif time_time()-self.pubchem_min_start>60.0:
            self.pubchem_min_start = time_time()
            self.pubchem_min_count = 0

    ## Try to retreive the xref from an inchi structure using pubchem
    #
    #
    '''
    No more than 5 requests per second.
    No more than 400 requests per minute.
    No longer than 300 second running time per minute.
    Requests exceeding limits are rejected (HTTP 503 error)
    '''
    def _pubchemStrctSearch(self, strct, itype='inchi'):
        self._pubChemLimit()
        try:
            # print()
            # print("REQUEST")
            # print()
            r = requests.post('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/'+str(itype)+'/xrefs/SBURL/JSON', data={itype: strct})
        except requests.exceptions.ConnectionError as e:
            self.logger.warning('Overloading PubChem, waiting 5 seconds and trying again')
            time_sleep(5)
            try:
                r = requests.post('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/'+str(itype)+'/xrefs/SBURL/JSON', data={itype: strct})
            except requests.exceptions.ConnectionError as e:
                self.logger.warning(e)
                return {}
        res_list = r.json()['InformationList']['Information']
        xref = {}
        if len(res_list)==1:
            name_r = requests.get('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/'+str(res_list[0]['CID'])+'/property/IUPACName,InChI,InChIKey,CanonicalSMILES/JSON')
            name_r_list = name_r.json()
            name = name_r_list['PropertyTable']['Properties'][0]['IUPACName']
            inchi = name_r_list['PropertyTable']['Properties'][0]['InChI']
            inchikey = name_r_list['PropertyTable']['Properties'][0]['InChIKey']
            smiles = name_r_list['PropertyTable']['Properties'][0]['CanonicalSMILES']
            xref['pubchem'] = [str(res_list[0]['CID'])]
            for url in res_list[0]['SBURL']:
                if 'https://biocyc.org/compound?orgid=META&id=' in url:
                    if 'biocyc' not in xref:
                        xref['biocyc'] = []
                    xref['biocyc'].append(url.replace('https://biocyc.org/compound?orgid=META&id=', ''))
                if 'http://www.hmdb.ca/metabolites/' in url:
                    if 'hmdb' not in xref:
                        xref['hmdb'] = []
                    xref['hmdb'].append(url.replace('http://www.hmdb.ca/metabolites/', ''))
                if 'http://www.genome.jp/dbget-bin/www_bget?cpd:' in url:
                    if 'kegg_c' not in xref:
                        xref['kegg_c'] = []
                    xref['kegg_c'].append(url.replace('http://www.genome.jp/dbget-bin/www_bget?cpd:', ''))
                if 'http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:' in url:
                    if 'chebi' not in xref:
                        xref['chebi'] = []
                    xref['chebi'].append(url.replace('http://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:', ''))
        elif len(res_list)==0:
            self.logger.warning('Could not find results for: '+str(strct))
            return {}
        else:
            self.logger.warning('There are more than one result for '+str(strct)+'... Ignoring')
            return {}
        return {'name': name, 'inchi': inchi, 'inchikey': inchikey, 'smiles': smiles, 'xref': xref}



    ###############################################################
    ############################ RP2paths entry functions #########
    ###############################################################

    ###############################################################
    ############################ RP2paths entry functions #########
    ###############################################################

    ## Function to group all the functions for parsing RP2 output to SBML files
    #
    # Takes RP2paths's compounds.txt and out_paths.csv and RetroPaths's *_scope.csv files and generates SBML
    #
    # @param compounds string path to RP2paths out_paths file
    # @param scope string path to RetroPaths2's scope file output
    # @param outPaths string path to RP2paths out_paths file
    # @param max_subpaths_filter int The maximal number of subpaths per path
    # @param compartment_id string The ID of the SBML's model compartment where to add the reactions to
    # @outFolder folder where to write files
    # @return Boolean The success or failure of the function
    def rp2ToSBML(self,
                  rp2_pathways,
                  rp2paths_compounds,
                  rp2paths_pathways,
                  outdir,
                  upper_flux_bound=999999,
                  lower_flux_bound=0,
                  max_subpaths_filter=10,
                  pathway_id='rp_pathway',
                  compartment_id='MNXC3',
                  species_group_id='central_species',
                  sink_species_group_id='rp_sink_species',
                  pubchem_search=False):

        if max_subpaths_filter<0:
            raise ValueError('Max number of subpaths cannot be less than 0: '+str(max_subpaths_filter))

        if not os_path.exists(outdir):
            os_mkdir(outdir)

        rp_strc = self._compounds(rp2paths_compounds)
        rp_transformation, sink_molecules = self._transformation(rp2_pathways)
        return self.Write_rp2pathsToSBML(rp_strc,
                                         rp_transformation,
                                         sink_molecules,
                                         rp2paths_pathways,
                                         outdir,
                                         upper_flux_bound,
                                         lower_flux_bound,
                                         max_subpaths_filter,
                                         pathway_id,
                                         compartment_id,
                                         species_group_id,
                                         sink_species_group_id,
                                         pubchem_search)

    ## Function to parse the compounds.txt file
    #
    #  Extract the smile and the structure of each compounds of RP2Path output
    #  Method to parse all the RP output compounds.
    #
    #  @param self Object pointer
    #  @param path The compounds.txt file path
    #  @return rp_compounds Dictionnary of smile and structure for each compound
    def _compounds(self, path):
        #self.rp_strc = {}
        rp_strc = {}
        try:
            if isinstance(path, bytes):
                reader = csv.reader(StringIO(path.decode('utf-8')), delimiter='\t')
            else:
                reader = csv.reader(open(path, 'r', encoding='utf-8'), delimiter='\t')
            next(reader)
            for row in reader:
                rp_strc[row[0]] = {'smiles': row[1]}  #, 'structure':row[1].replace('[','').replace(']','')
                try:
                    rp_strc[row[0]]['inchi'] = self.mnxm_strc[row[0]]['inchi']
                except KeyError:
                    #try to generate them yourself by converting them directly
                    try:
                        resConv = self._convert_depiction(idepic=row[1], itype='smiles', otype={'inchi'})
                        rp_strc[row[0]]['inchi'] = resConv['inchi']
                    except NotImplementedError as e:
                        self.logger.warning('Could not convert the following SMILES to InChI: '+str(row[1]))
                try:
                    rp_strc[row[0]]['inchikey'] = self.mnxm_strc[row[0]]['inchikey']
                    #try to generate them yourself by converting them directly
                    #TODO: consider using the inchi writing instead of the SMILES notation to find the inchikey
                except KeyError:
                    try:
                        resConv = self._convert_depiction(idepic=row[1], itype='smiles', otype={'inchikey'})
                        rp_strc[row[0]]['inchikey'] = resConv['inchikey']
                    except NotImplementedError as e:
                        self.logger.warning('Could not convert the following SMILES to InChI key: '+str(row[1]))
        except (TypeError, FileNotFoundError) as e:
            self.logger.error('Could not read the compounds file ('+str(path)+')')
            raise RuntimeError
        return rp_strc


    ## Function to parse the scope.csv file
    #
    #  Extract the reaction rules from the retroPath2.0 output using the scope.csv file
    #
    #  @param self Object pointer
    #  @param path The scope.csv file path
    def _transformation(self, path):
        rp_transformation = {}
        sink_molecules = []
        #### we might pass binary in the REST version
        reader = None
        if isinstance(path, bytes):
            reader = csv.reader(StringIO(path.decode('utf-8')), delimiter=',')
        else:
            try:
                reader = csv.reader(open(path, 'r'), delimiter=',')
            except FileNotFoundError:
                self.logger.error('Could not read the compounds file: '+str(path))
                return {}
        next(reader)
        for row in reader:
            if not row[1] in rp_transformation:
                rp_transformation[row[1]] = {}
                rp_transformation[row[1]]['rule'] = row[2]
                rp_transformation[row[1]]['ec'] = [i.replace(' ', '') for i in row[11][1:-1].split(',') if not i.replace(' ', '')=='NOEC']
            if row[7]=='1':
                for i in row[8].replace(']', '').replace('[', '').replace(' ', '').split(','):
                    sink_molecules.append(i)

        # self.logger.info(rp_transformation)
        # self.logger.info(sink_molecules)
        return rp_transformation, list(set(sink_molecules))


    def _read_paths(self, rp2paths_pathways):

        #### we might pass binary in the REST version
        if isinstance(rp2paths_pathways, bytes):
            reader = csv.reader(StringIO(rp2paths_pathways.decode('utf-8')))
        else:
            reader = csv.reader(open(rp2paths_pathways, 'r'))
        next(reader)
        current_path_id = 0
        path_step = 1

        rp_paths = {}

        for row in reader:
            try:
                #Remove all illegal characters in SBML ids
                row[3] = row[3].replace("'", "").replace('-', '_').replace('+', '')
                if not int(row[0])==current_path_id:
                    path_step = 1
                else:
                    path_step += 1
                #important to leave them in order
                current_path_id = int(row[0])
            except ValueError:
                self.logger.error('Cannot convert path_id to int ('+str(row[0])+')')
                #return {}
                return False
            #################################
            ruleIds = row[2].split(',')
            if ruleIds==None:
                self.logger.warning('The rulesIds is None')
                #pass # or continue
                continue
            ###WARNING: This is the part where we select some rules over others
            # we do it by sorting the list according to their score and taking the topx
            tmp_rr_reactions = {}
            for r_id in ruleIds:
                for rea_id in self.rr_reactions[r_id]:
                    tmp_rr_reactions[str(r_id)+'__'+str(rea_id)] = self.rr_reactions[r_id][rea_id]
            # if len(ruleIds)>int(maxRuleIds):
            #     self.logger.warning('There are too many rules, limiting the number to random top '+str(maxRuleIds))
            #     try:
            #         ruleIds = [y for y,_ in sorted([(i, tmp_rr_reactions[i]['rule_score']) for i in tmp_rr_reactions])][:int(maxRuleIds)]
            #     except KeyError:
            #         self.logger.warning('Could not select topX due inconsistencies between rules ids and rr_reactions... selecting random instead')
            #         ruleIds = random.sample(tmp_rr_reactions, int(maxRuleIds))
            # else:
            ruleIds = tmp_rr_reactions
            sub_path_step = 1
            for singleRule in ruleIds:
                tmpReac = {'rule_id': singleRule.split('__')[0],
                           'rule_ori_reac': {'mnxr': singleRule.split('__')[1]},
                           'rule_score': self.rr_reactions[singleRule.split('__')[0]][singleRule.split('__')[1]]['rule_score'],
                           'right': {},
                           'left': {},
                           'path_id': int(row[0]),
                           'step': path_step,
                           'transformation_id': row[1][:-2]}
                ############ LEFT ##############
                for l in row[3].split(':'):
                    tmp_l = l.split('.')
                    try:
                        #tmpReac['left'].append({'stoichio': int(tmp_l[0]), 'name': tmp_l[1]})
                        mnxm = '' #TODO: change this
                        if tmp_l[1] in self.deprecatedMNXM_mnxm:
                            mnxm = self.deprecatedMNXM_mnxm[tmp_l[1]]
                        else:
                            mnxm = tmp_l[1]
                        tmpReac['left'][mnxm] = int(tmp_l[0])
                    except ValueError:
                        self.logger.error('Cannot convert tmp_l[0] to int ('+str(tmp_l[0])+')')
                        #return {}
                        return False
                ############## RIGHT ###########
                for r in row[4].split(':'):
                    tmp_r = r.split('.')
                    try:
                        #tmpReac['right'].append({'stoichio': int(tmp_r[0]), 'name': tmp_r[1]})
                        mnxm = '' #TODO change this
                        if tmp_r[1] in self.deprecatedMNXM_mnxm:
                            mnxm = self.deprecatedMNXM_mnxm[tmp_r[1]]  #+':'+self.rr_reactions[tmpReac['rule_id']]['left']
                        else:
                            mnxm = tmp_r[1]  #+':'+self.rr_reactions[tmpReac['rule_id']]['left']
                        tmpReac['right'][mnxm] = int(tmp_r[0])
                    except ValueError:
                        self.logger.error('Cannot convert tmp_r[0] to int ('+str(tmp_r[0])+')')
                        return {}
                #################################
                if not int(row[0]) in rp_paths:
                    rp_paths[int(row[0])] = {}
                if not int(path_step) in rp_paths[int(row[0])]:
                    rp_paths[int(row[0])][int(path_step)] = {}
                rp_paths[int(row[0])][int(path_step)][int(sub_path_step)] = tmpReac
                sub_path_step += 1

        return rp_paths

    def _unique_species(self, meta, rp_strc, pubchem_search):

        try:
            chemName = self.mnxm_strc[meta]['name']
        except KeyError:
            chemName = None

        #compile as much info as you can

        #xref
        try: xref = self.chemXref[meta]
        except KeyError: xref = {}
        spe = Species(None, None, None, xref)

        ###### Try to recover the structures ####
        pubchem = Species(None, None, None, {})

        #inchi
        try:
            spe.inchi = rp_strc[meta]['inchi']
            if not spe.xref and pubchem_search:
                try:
                    # print()
                    # print("*************")
                    # print("pubchem_species READ (INCHI)", spe.inchi)
                    # print("*************")
                    # print()
                    pubchem.inchi    = self.pubchem_species[spe.inchi]['inchi']
                    pubchem.inchikey = self.pubchem_species[spe.inchi]['inchikey']
                    pubchem.smiles   = self.pubchem_species[spe.inchi]['smiles']
                    pubchem.xref     = self.pubchem_species[spe.inchi]['xref']
                except KeyError:
                    # print()
                    # print("*************")
                    # print("pubchem_species SEARCH (INCHI)", spe.inchi)
                    # print("*************")
                    # print()
                    pubres = self._pubchemStrctSearch(spe.inchi, 'inchi')
                    if not chemName:
                        chemName = pubres['name']
                    if 'chebi' in pubres['xref']:
                        try:
                            spe.xref = self.chemXref[self.chebi_mnxm[pubres['xref']['chebi'][0]]]
                        except KeyError:
                            pass
                    # pubchem.fill_missing(pubres)
                    if not pubchem.inchi:
                        pubchem.inchi = pubres['inchi']
                    if not pubchem.xref:
                        pubchem.xref = pubres['xref']
                    if not pubchem.inchikey:
                        pubchem.inchikey = pubres['inchikey']
                    if not pubchem.smiles:
                        pubchem.smiles = pubres['smiles']
        except KeyError:
            pass
        #inchikey
        try:
            spe.inchikey = rp_strc[meta]['inchikey']
            if not spe.xref and pubchem_search:
                # print("*************")
                # print("pubchem_species SEARCH (INCHIKEY)", spe.inchi)
                # print("*************")
                # print()
                pubres = self._pubchemStrctSearch(spe.inchikey, 'inchikey')
                print(pubres)
                if not chemName:
                    chemName = pubres['name']
                if 'chebi' in pubres['xref']:
                    try:
                        spe.xref = self.chemXref[self.chebi_mnxm[pubres['xref']['chebi'][0]]]
                    except KeyError:
                        pass
                if not pubchem.xref:
                    pubchem.xref = pubres['xref']
                if not pubchem.inchi:
                    pubchem.inchi = pubres['inchi']
                if not pubchem.smiles:
                    pubchem.smiles = pubres['smiles']
        except KeyError:
            pass
        #smiles
        try:
            spe.smiles = rp_strc[meta]['smiles']
            if not spe.xref and pubchem_search:
                # print()
                # print("*************")
                # print("pubchem_species SEARCH (SMILES)", spe.inchi)
                # print("*************")
                # print()
                pubres = self._pubchemStrctSearch(spe.smiles, 'smiles')
                print(pubres)
                if not chemName:
                    chemName = pubres['name']
                if 'chebi' in pubres['xref']:
                    try:
                        spe.xref = self.chemXref[self.chebi_mnxm[pubres['xref']['chebi'][0]]]
                    except KeyError:
                        pass
                if not pubchem.xref:
                    pubchem.xref = pubres['xref']
                if not pubchem.inchi:
                    pubchem.inchi = pubres['inchi']
                if not pubchem.inchikey:
                    pubchem.inchikey = pubres['inchikey']
        except KeyError:
            pass

        if not spe.inchi:
            spe.inchi = pubchem.inchi
        if not spe.inchikey:
            spe.inchikey = pubchem.inchikey
        if not spe.smiles:
            spe.smiles = pubchem.smiles
        if not spe.xref:
            spe.xref = pubchem.xref
        if pubchem.inchi:
            # print()
            # print("*************")
            # print("pubchem_species WRITTEN")
            # print("*************")
            # print()
            self.pubchem_species[pubchem.inchi] = {'inchi': pubchem.inchi, 'smiles': pubchem.smiles, 'inchikey': pubchem.inchikey, 'xref': pubchem.xref}


        return (chemName, spe)



    #TODO: make sure that you account for the fact that each reaction may have multiple associated reactions
    ## Function to parse the out_paths.csv file
    #
    #  Reading the RP2path output and extract all the information for each pathway
    #  RP2path Metabolic pathways from out_paths.csv
    #  create all the different values for heterologous paths from the RP2path out_paths.csv file
    #  Note that path_step are in reverse order here
    #
    #  @param self Object pointer
    #  @param path The out_path.csv file path
    #  @max_subpaths_filter maximal numer of subpaths per paths
    #  @outFolder folder where to write files
    #  @return Boolean The success or failure of the function
    def Write_rp2pathsToSBML(self,
                             rp_strc, rp_transformation,
                             sink_molecules,
                             rp2paths_pathways,
                             outFolder,
                             upper_flux_bound=999999,
                             lower_flux_bound=0,
                             max_subpaths_filter=10,
                             pathway_id='rp_pathway',
                             compartment_id='MNXC3',
                             species_group_id='central_species',
                             sink_species_group_id='rp_sink_species',
                             pubchem_search=False):

        rp_paths = self._read_paths(rp2paths_pathways)
        sink_species = []

        # for each line or rp2paths_pathways:
        #     generate comb
        #     for each combinant:
        #         rank
        #         process
        #         add cofactors
        #         dedup

        #### pathToSBML ####
        try:
            mnxc = self.name_compXref[compartment_id]
        except KeyError:
            self.logger.error('Could not Xref compartment_id ('+str(compartment_id)+')')
            return False


        for pathNum in rp_paths:

            #first level is the list of lists of sub_steps
            #second is itertools all possible combinations using product
            altPathNum = 1
            # topX subpaths of the current rp2path pathway
            local_rpsbml_items = []

            for comb_path in list(itertools_product(*[[(i,y) for y in rp_paths[pathNum][i]] for i in rp_paths[pathNum]])):
                steps = []
                for i, y in comb_path:
                    steps.append(rp_paths[pathNum][i][y])
                path_id = steps[0]['path_id']
                rpsbml = rpSBML('rp_'+str(path_id)+'_'+str(altPathNum))

                #1) Create a generic Model, ie the structure and unit definitions that we will use the most
                ##### TODO: give the user more control over a generic model creation:
                #   -> special attention to the compartment
                rpsbml.genericModel(
                        'RetroPath_Pathway_'+str(path_id)+'_'+str(altPathNum),
                        'RP_model_'+str(path_id)+'_'+str(altPathNum),
                        self.compXref[mnxc],
                        compartment_id,
                        upper_flux_bound,
                        lower_flux_bound)

                #2) Create the pathway (groups)
                rpsbml.createPathway(pathway_id)
                rpsbml.createPathway(species_group_id)
                rpsbml.createPathway(sink_species_group_id)

                #3) Find all unique species and add them to the model
                all_meta = set([i for step in steps for lr in ['left', 'right'] for i in step[lr]])
                for meta in all_meta:
                    (chemName, spe) = self._unique_species(meta, rp_strc, pubchem_search)
                    if chemName:
                        chemName = chemName.replace("'", "")
                    # self.logger.info('Creating species: '+str(chemName)+' ('+str(meta)+')')
                    #pass the information to create the species
                    if meta in sink_molecules:
                        # self.logger.info('Species is sink: '+str(sink_species_group_id))
                        rpsbml.createSpecies(meta,
                                             compartment_id,
                                             chemName,
                                             spe.xref,
                                             spe.inchi,
                                             spe.inchikey,
                                             spe.smiles,
                                             species_group_id,
                                             sink_species_group_id)
                    else:
                        rpsbml.createSpecies(meta,
                                             compartment_id,
                                             chemName,
                                             spe.xref,
                                             spe.inchi,
                                             spe.inchikey,
                                             spe.smiles,
                                             species_group_id)

                #4) Add the complete reactions and their annotations
                for step in steps:
                    #add the substep to the model
                    step['sub_step'] = altPathNum
                    rpsbml.createReaction(
                            'RP'+str(step['step']), # parameter 'name' of the reaction deleted : 'RetroPath_Reaction_'+str(step['step']),
                            upper_flux_bound, lower_flux_bound, step, compartment_id,
                            rp_transformation[step['transformation_id']]['rule'],
                            {'ec': rp_transformation[step['transformation_id']]['ec']},
                            pathway_id)

                #5) Adding the consumption of the target
                targetStep = {
                        'rule_id': None,
                        'left': {[i for i in all_meta if i[:6]=='TARGET'][0]: 1},
                        'right': [],
                        'step': None,
                        'sub_step': None,
                        'path_id': None,
                        'transformation_id': None,
                        'rule_score': None,
                        'rule_ori_reac': None}
                rpsbml.createReaction('RP1_sink',
                                      upper_flux_bound, lower_flux_bound,
                                      targetStep,
                                      compartment_id)

                #6) Adding the cofactors
                self.addCofactors(rpsbml)

                #7) Filtering
                sbml_item = SBML_Item(rpsbml.getScore(),
                                      'rp_'+str(path_id)+'_'+str(altPathNum),
                                      rpsbml)
                unique = True

                # For each subpath already in local_sbml_paths
                for item in local_rpsbml_items:
                    # Compare with the new built pathway
                    if sbml_item.rpsbml_obj==item.rpsbml_obj:
                        unique = False
                        # print(rpsbml.outPathsDict())
                        # If its score is better, then replace the one already in place
                        if sbml_item.score > item.score:
                            # Remove the same pathway with worse score from the list
                            local_rpsbml_items.remove(item)
                            # Insert at the good place the new pathway (with better score)
                            bisect_insort(local_rpsbml_items, sbml_item)
                        # Leave the loop
                        break

                # If the pathway currently built is not already in local_sbml_paths
                # The new built pathway is unique
                if unique:
                    # print("INSERT UNIQUE", sbml_item.score)
                    # If its score is better, then replace the one already in place
                    # Insert at the good place the new pathway (with better score)
                    bisect_insort(local_rpsbml_items, sbml_item)

                # Keep only topX
                local_rpsbml_items = local_rpsbml_items[-max_subpaths_filter:]

                # print(sbml_item.index, [i.score for i in local_rpsbml_items])
                    # print(item.rpsbml_obj.outPathsDict())
                    # print()

                altPathNum += 1

            # for i in range(len(local_sbml_paths)):
            #     for j in range(i+1, len(local_sbml_paths)):
            #         if local_sbml_paths[i].rpsbml_obj == local_sbml_paths[j].rpsbml_obj:
            #             print("NOT UNIQUE !!!")
            #             exit()

            # Write results to files
            for rpsbml_item in local_rpsbml_items:
                rpsbml_item.rpsbml_obj.writeSBML(outFolder)

            # for item in local_rpsbml_paths:
            #     sbml_paths[item.index] = item.rpsbml_obj
            # sbml_paths += local_sbml_paths


        return True


    #######################################################################
    ############################# JSON input ##############################
    #######################################################################


    #WARNING: we are not setting any restrictions on the number of steps allowed here and instead we
    #take the rule with the highest dimension. Also assume that there is a single rule at a maximal
    #dimension
    ## Function to generate an SBLM model from a JSON file
    #
    #  Read the json files of a folder describing pathways and generate an SBML file for each
    #  TODO: remove the default MNXC3 compartment ID
    #  TODO: change the ID of all species to take a normal string and not sepcial caracters
    #  WARNING: We are only using a single rule (technically with the highest diameter)
    #
    #  @param self Object pointer
    # @param colJson Dictionnary of
    #  @return rpsbml.document the SBML document
    #TODO: update this to include _hdd parsing
    # def jsonToSBML(self,
    #                collJson,
    #                upper_flux_bound=999999,
    #                lower_flux_bound=0,
    #                pathway_id='rp_pathway',
    #                compartment_id='MNXC3',
    #                species_group_id='central_species',
    #                pubchem_search=False):
    #     #global parameters used for all parameters
    #     pathNum = 1
    #     rp_paths = {}
    #     reac_smiles = {}
    #     reac_ecs = {}
    #     species_list = {}
    #     reactions_list = {}
    #     source_cid = {}
    #     source_stochio = {}
    #     cid_inchikey = {}
    #     sink_species = {}
    #     ############################################
    #     ############## gather the data #############
    #     ############################################
    #     for json_dict in collJson:
    #         ########### construct rp_paths ########
    #         reac_smiles[pathNum] = {}
    #         reac_ecs[pathNum] = {}
    #         species_list[pathNum] = {}
    #         reactions_list[pathNum] = {}
    #         cid_inchikey[pathNum] = {}
    #         sink_species[pathNum] = {}
    #         stochio = {}
    #         inchikey_cid = {}
    #         source_species = []
    #         skip_pathway = False
    #         for node in collJson[json_dict]['elements']['nodes']:
    #             ##### compounds ####
    #             if node['data']['type']=='compound':
    #                 cid_inchikey[pathNum][node['data']['id'].replace('-', '')] = node['data']['id']
    #                 species_list[pathNum][node['data']['id'].replace('-', '')] = {'inchi': node['data']['InChI'],
    #                                 'inchikey': node['data']['id'],
    #                                 'smiles': node['data']['SMILES']}
    #                 if int(node['data']['isSource'])==1:
    #                     #TODO: there should always be only one source, to check
    #                     source_species.append(node['data']['id'].replace('-', ''))
    #                     source_cid[pathNum] = node['data']['id'].replace('-', '')
    #                 if int(node['data']['inSink'])==1:
    #                     sink_species[pathNum][node['data']['id'].replace('-', '')] = node['data']['id']
    #             ###### reactions ######
    #             elif node['data']['type']=='reaction':
    #                 #NOTE: pick the rule with the highest diameter
    #                 r_id = sorted(node['data']['Rule ID'], key=lambda x: int(x.split('-')[-2]), reverse=True)[0]
    #                 reactions_list[pathNum][node['data']['id']] = {'rule_id': r_id,
    #                     'rule_ori_reac': None,
    #                     'right': {},
    #                     'left': {},
    #                     'path_id': pathNum,
    #                     'step': None,
    #                     'sub_step': None,
    #                     'transformation_id': None,
    #                     'rule_score': node['data']['Score']}
    #                 reac_smiles[pathNum][r_id] = node['data']['Reaction SMILES']
    #                 reac_ecs[pathNum][r_id] = list(filter(None, [i for i in node['data']['EC number']]))
    #                 stochio[node['data']['id']] = {}
    #                 for i in node['data']['Stoechiometry']:
    #                     stochio[node['data']['id']][i.replace('-', '')] = node['data']['Stoechiometry'][i]
    #         ############ make main pathway ###########
    #         main_reac = {}
    #         for reaction_node in collJson[json_dict]['elements']['edges']:
    #             if not len(reaction_node['data']['source'].split('-'))==3:
    #                 if not reaction_node['data']['source'] in reactions_list[pathNum]:
    #                     self.logger.error('The following reaction was not found in the JSON elements: '+str(reaction_node['data']['source']))
    #                     skip_pathway = True
    #                     break
    #                 else:
    #                     rid = reaction_node['data']['source']
    #                     try:
    #                         cid = inchikey_cid[reaction_node['data']['target'].replace('-', '')]
    #                     except KeyError:
    #                         cid = reaction_node['data']['target'].replace('-', '')
    #                     try:
    #                         reactions_list[pathNum][rid]['left'][cid] = stochio[rid][cid]
    #                     except KeyError:
    #                         reactions_list[pathNum][rid]['left'][cid] = 1.0
    #                         self.logger.warning('The cid ('+str(cid)+') has not been detected by stochio. Setting to 1.0')
    #             if not len(reaction_node['data']['target'].split('-'))==3:
    #                 if not reaction_node['data']['target'] in reactions_list[pathNum]:
    #                     self.logger.error('The following reaction was not found in the JSON elements: '+str(reaction_node['data']['source']))
    #                     skip_pathway = True
    #                     break
    #                 else:
    #                     rid = reaction_node['data']['target']
    #                     try:
    #                         cid = inchikey_cid[reaction_node['data']['source'].replace('-', '')]
    #                     except KeyError:
    #                         cid = reaction_node['data']['source'].replace('-', '')
    #                     try:
    #                         reactions_list[pathNum][rid]['right'][cid] = stochio[rid][cid]
    #                     except KeyError:
    #                         reactions_list[pathNum][rid]['right'][cid] = 1.0
    #                         self.logger.warning('The cid ('+str(cid)+') has not been detected by stochio. Setting to 1.0')
    #         ################# calculate the steps associated with the reactions_list ######
    #         #find the source in the LAST reaction in the pathway
    #         #NOTE: this assumes that the source is contained in a single final reaction and nowhere else
    #         last_rid = None
    #         step_num = 1
    #         found_rid = []
    #         toFind_rid = list(reactions_list[pathNum].keys())
    #         for rid in reactions_list[pathNum]:
    #             if all([True if i in source_species else False for i in reactions_list[pathNum][rid]['right']]):
    #                 reactions_list[pathNum][rid]['step'] = step_num
    #                 #step_num -= 1
    #                 step_num += 1
    #                 last_rid = rid
    #                 try:
    #                     source_stochio[pathNum] = stochio[rid][source_cid[pathNum]]
    #                 except KeyError:
    #                     source_stochio[pathNum] = 1.0
    #                 found_rid.append(rid)
    #                 toFind_rid.remove(rid)
    #                 break
    #         for rid in toFind_rid[:]:
    #             if all([True if i in reactions_list[pathNum][last_rid]['left'] else False for i in reactions_list[pathNum][rid]['right']]):
    #                 reactions_list[pathNum][rid]['step'] = step_num
    #                 #step_num -= 1
    #                 step_num += 1
    #                 last_rid = rid
    #                 found_rid.append(rid)
    #                 toFind_rid.remove(rid)
    #         if not toFind_rid==[]:
    #             self.logger.error('There are reactions unaccounted for: '+str(toFind_rid))
    #             skip_pathway = True
    #             break
    #         ############# find all the alternative reactions associated with a reaction rule ###
    #         if not skip_pathway:
    #             rp_paths[pathNum] = {}
    #             for rid in reactions_list[pathNum]:
    #                 rp_paths[pathNum][reactions_list[pathNum][rid]['step']] = {}
    #                 sub_step = 1
    #                 for reac_id in self.rr_reactions[reactions_list[pathNum][rid]['rule_id']]:
    #                     tmpReac = deepcopy(reactions_list[pathNum][rid])
    #                     tmpReac['mnxr'] = reac_id
    #                     tmpReac['sub_step'] = sub_step
    #                     rp_paths[pathNum][reactions_list[pathNum][rid]['step']][sub_step] = tmpReac
    #                     sub_step += 1
    #         else:
    #             self.logger.warning('Skipping pathway '+str(pathNum))
    #         pathNum += 1
    #     ######################################
    #     ########### create the SBML's ########
    #     ######################################
    #     try:
    #         mnxc = self.name_compXref[compartment_id]
    #     except KeyError:
    #         self.logger.error('Could not Xref compartment_id ('+str(compartment_id)+')')
    #         return False
    #     sbml_paths = {}
    #     for pathNum in rp_paths:
    #         #first level is the list of lists of sub_steps
    #         #second is itertools all possible combinations using product
    #         altPathNum = 1
    #         for comb_path in list(itertools_product(*[[(i,y) for y in rp_paths[pathNum][i]] for i in rp_paths[pathNum]])):
    #             steps = []
    #             for i, y in comb_path:
    #                 steps.append(rp_paths[pathNum][i][y])
    #             path_id = steps[0]['path_id']
    #             rpsbml = rpSBML.rpSBML('rp_'+str(path_id)+'_'+str(altPathNum))
    #             #1) create a generic Model, ie the structure and unit definitions that we will use the most
    #             ##### TODO: give the user more control over a generic model creation:
    #             #   -> special attention to the compartment
    #             rpsbml.genericModel('RetroPath_Pathway_'+str(path_id)+'_'+str(altPathNum),
    #                     'RP_model_'+str(path_id)+'_'+str(altPathNum),
    #                     self.compXref[mnxc],
    #                     compartment_id,
    #                     upper_flux_bound,
    #                     lower_flux_bound)
    #             #2) create the pathway (groups)
    #             rpsbml.createPathway(pathway_id)
    #             rpsbml.createPathway(species_group_id)
    #             #3) find all the unique species and add them to the model
    #             ###################################################
    #             ############## SPECIES ############################
    #             ###################################################
    #             meta_to_cid = {}
    #             for meta in species_list[pathNum]:
    #                 #### beofre adding it to the model check to see if you can recover some MNXM from inchikey
    #                 #NOTE: only for the sink species do we try to convert to MNXM
    #                 if meta in list(sink_species[pathNum].keys()):
    #                     try:
    #                         #take the smallest MNX, usually the best TODO: review this
    #                         cid = sorted(self.inchikey_mnxm[sink_species[pathNum][meta]], key=lambda x: int(x[4:]))[0]
    #                         meta_to_cid[meta] = cid
    #                     except KeyError:
    #                         self.logger.error('Cannot find sink compound: '+str(meta))
    #                         return False
    #                 else:
    #                     cid = meta
    #                 # retreive the name of the molecule
    #                 #here we want to gather the info from rpCompletion's rp_strc and mnxm_strc
    #                 try:
    #                     chemName = self.mnxm_strc[meta]['name']
    #                 except KeyError:
    #                     chemName = None
    #                 #compile as much info as you can
    #                 #xref
    #                 try:
    #                     spe_xref = self.chemXref[meta]
    #                 except KeyError:
    #                     spe_xref = {}
    #                 ###### Try to recover the structures ####
    #                 spe_smiles = None
    #                 spe_inchi = None
    #                 spe_inchikey = None
    #                 pubchem_smiles = None
    #                 pubchem_inchi = None
    #                 pubchem_inchikey = None
    #                 pubchem_xref = {}
    #                 #inchi
    #                 try:
    #                     spe_inchi = rp_strc[meta]['inchi']
    #                     if not spe_xref and pubchem_search:
    #                         try:
    #                             pubchem_inchi = self.pubchem_species[spe_inchi]['inchi']
    #                             pubchem_inchikey = self.pubchem_species[spe_inchi]['inchikey']
    #                             pubchem_smiles = self.pubchem_species[spe_inchi]['smiles']
    #                             pubchem_xref = self.pubchem_species[spe_inchi]['xref']
    #                         except KeyError:
    #                             pubres = self._pubchemStrctSearch(spe_inchi, 'inchi')
    #                             if not chemName:
    #                                 chemName = pubres['name']
    #                             if 'chebi' in pubres['xref']:
    #                                 try:
    #                                     #WARNING: taking the first one. Better to take the smallest?
    #                                     spe_xref = self.chemXref[self.chebi_mnxm[pubres['xref']['chebi'][0]]]
    #                                 except KeyError:
    #                                     pass
    #                             if not pubchem_xref:
    #                                 pubchem_xref = pubres['xref']
    #                             if not pubchem_inchikey:
    #                                 pubchem_inchikey = pubres['inchikey']
    #                             if not pubchem_smiles:
    #                                 pubchem_smiles = pubres['smiles']
    #                 except KeyError:
    #                     pass
    #                 #inchikey
    #                 try:
    #                     spe_inchikey = rp_strc[meta]['inchikey']
    #                     if not spe_xref and pubchem_search:
    #                         pubres = self._pubchemStrctSearch(spe_inchikey, 'inchikey')
    #                         if not chemName:
    #                             chemName = pubres['name']
    #                         if 'chebi' in pubres['xref']:
    #                             try:
    #                                 spe_xref = self.chemXref[self.chebi_mnxm[pubres['xref']['chebi'][0]]]
    #                             except KeyError:
    #                                 pass
    #                         if not pubchem_xref:
    #                             pubchem_xref = pubres['xref']
    #                         if not pubchem_inchi:
    #                             pubchem_inchi = pubres['inchi']
    #                         if not pubchem_smiles:
    #                             pubchem_smiles = pubres['smiles']
    #                 except KeyError:
    #                     pass
    #                 #smiles
    #                 try:
    #                     spe_smiles = rp_strc[meta]['smiles']
    #                     if not spe_xref and pubchem_search:
    #                         pubres = self._pubchemStrctSearch(spe_smiles, 'smiles')
    #                         if not chemName:
    #                             chemName = pubres['name']
    #                         if 'chebi' in pubres['xref']:
    #                             try:
    #                                 spe_xref = self.chemXref[self.chebi_mnxm[pubres['xref']['chebi'][0]]]
    #                             except KeyError:
    #                                 pass
    #                         if not pubchem_xref:
    #                             pubchem_xref = pubres['xref']
    #                         if not pubchem_inchi:
    #                             pubchem_inchi = pubres['inchi']
    #                         if not pubchem_inchikey:
    #                             pubchem_inchikey = pubres['inchikey']
    #                 except KeyError:
    #                     pass
    #                 if not spe_inchi:
    #                     spe_inchi = pubchem_inchi
    #                 if not spe_inchikey:
    #                     spe_inchikey = pubchem_inchikey
    #                 if not spe_smiles:
    #                     spe_smiles = pubchem_smiles
    #                 if not spe_xref:
    #                     spe_xref = pubchem_xref
    #                 #pass the information to create the species
    #                 rpsbml.createSpecies(meta,
    #                         compartment_id,
    #                         chemName,
    #                         spe_xref,
    #                         spe_inchi,
    #                         spe_inchikey,
    #                         spe_smiles,
    #                         species_group_id)
    #             #4) add the reactions and their annotations
    #             for step in steps:
    #                 #add the substep to the model
    #                 step['sub_step'] = altPathNum
    #                 #### try to replace the sink compounds inchikeys with mnxm
    #                 for direc in ['right', 'left']:
    #                     step_mnxm = {}
    #                     for meta in step[direc]:
    #                         try:
    #                             step_mnxm[meta_to_cid[meta]] = step[direc][meta]
    #                         except KeyError:
    #                             step_mnxm[meta] = step[direc][meta]
    #                     step[direc] = step_mnxm
    #                 rpsbml.createReaction('RP'+str(step['step']),
    #                         upper_flux_bound,
    #                         lower_flux_bound,
    #                         step,
    #                         compartment_id,
    #                         reac_smiles[pathNum][step['rule_id']],
    #                         {'ec': reac_ecs[pathNum][step['rule_id']]},
    #                         pathway_id)
    #             #5) adding the consumption of the target
    #             targetStep = {'rule_id': None,
    #                     'left': {source_cid[pathNum]: source_stochio[pathNum]},
    #                     'right': {},
    #                     'step': None,
    #                     'sub_step': None,
    #                     'path_id': None,
    #                     'transformation_id': None,
    #                     'rule_score': None,
    #                     'rule_ori_reac': None}
    #             rpsbml.createReaction('RP1_sink',
    #                     upper_flux_bound,
    #                     lower_flux_bound,
    #                     targetStep,
    #                     compartment_id)
    #             #6) Optional?? Add the flux objectives. Could be in another place, TBD
    #             rpsbml.createFluxObj('rpFBA_obj', 'RP1_sink', 1, True)
    #             sbml_paths['rp_'+str(step['path_id'])+'_'+str(altPathNum)] = rpsbml
    #             altPathNum += 1


    #############################################################################################
    ############################### TSV data tsv ################################################
    #############################################################################################


    ## Function to parse the TSV of measured heterologous pathways to SBML
    #
    # Given the TSV of measured pathways, parse them to a dictionnary, readable to next be parsed
    # to SBML
    #
    # @param self object pointer
    # @param inFile The input JSON file
    # @param mnxHeader Reorganise the results around the target MNX products
    # @return Dictionnary of SBML
    def _parseTSV(self, inFile, remove_inchi_4p=False, mnxHeader=False):
        data = {}
        try:
            for row in csv.DictReader(open(inFile), delimiter='\t'):
                ######## path_id ######
                try:
                    pathID = int(row['pathway_ID'])
                except ValueError:
                    self.logger.error('Cannot convert pathway ID: '+str(row['pathway_ID']))
                    continue
                if not pathID in data:
                    data[pathID] = {}
                    data[pathID]['isValid'] = True
                    data[pathID]['steps'] = {}
                ####### target #########
                if not 'target' in data[pathID]:
                    data[pathID]['target'] = {}
                    data[pathID]['target']['name'] = row['target_name']
                    if remove_inchi_4p:
                        data[pathID]['target']['inchi'] = '/'.join([row['target_structure'].split('/')[i] for i in range(len(row['target_structure'].split('/'))) if i<4])
                    else:
                        data[pathID]['target']['inchi'] = row['target_structure']
                ####### step #########
                try:
                    stepID = int(row['step'])
                except ValueError:
                    self.logger.error('Cannot convert step ID: '+str(row['step']))
                    data[pathID]['isValid'] = False
                    continue
                if stepID==0:
                    continue
                elif stepID==1:
                    data[pathID]['organism'] = row['organism'].replace(' ', '')
                    data[pathID]['reference'] = row['reference'].replace(' ', '')
                data[pathID]['steps'][stepID] = {}
                ##### substrates #########
                data[pathID]['steps'][stepID]['substrates'] = []
                lenDBref = len(row['substrate_dbref'].split(';'))
                for i in row['substrate_dbref'].split(';'):
                    if i=='':
                        lenDBref -= 1
                lenStrc = len(row['substrate_structure'].split('_'))
                for i in row['substrate_structure'].split('_'):
                    if i=='':
                        lenStrc -= 1
                lenSub = len(row['substrate_name'].split(';'))
                for i in row['substrate_name'].split(';'):
                    if i=='':
                        lenSub -= 1
                if lenSub==lenStrc==lenSub:
                    for name, inchi, dbrefs in zip(row['substrate_name'].split(';'),
                            row['substrate_structure'].split('_'),
                            row['substrate_dbref'].split(';')):
                        tmp = {}
                        if remove_inchi_4p:
                            tmp['inchi'] = '/'.join([inchi.split('/')[i] for i in range(len(inchi.split('/'))) if i<4])
                        else:
                            tmp['inchi'] = inchi.replace(' ', '')
                        tmp['name'] = name
                        tmp['dbref'] = {}
                        for dbref in dbrefs.split('|'):
                            if len(dbref.split(':'))==2:
                                db_name = dbref.split(':')[0].replace(' ', '').lower()
                                db_cid = dbref.split(':')[1].replace(' ', '')
                                if not db_name in tmp['dbref']:
                                    tmp['dbref'][db_name] = []
                                tmp['dbref'][db_name].append(db_cid)
                            else:
                                self.logger.warning('Ignoring the folowing product dbref ('+str(name)+'): '+str(dbref))
                                data[pathID]['isValid'] = False
                        data[pathID]['steps'][stepID]['substrates'].append(tmp)
                else:
                    self.logger.warning('Not equal length between substrate names, their structure or dbref ('+str(name)+'): '+str(row['substrate_name'])+' <--> '+str(row['substrate_structure'])+' <--> '+str(row['substrate_dbref']))
                    data[pathID]['isValid'] = False
                    continue
                ##### products #########
                data[pathID]['steps'][stepID]['products'] = []
                lenDBref = len(row['product_dbref'].split(';'))
                for i in row['product_dbref'].split(';'):
                    if i=='':
                        lenDBref -= 1
                lenStrc = len(row['product_structure'].split('_'))
                for i in row['product_structure'].split('_'):
                    if i=='':
                        lenStrc -= 1
                lenSub = len(row['product_name'].split(';'))
                for i in row['product_name'].split(';'):
                    if i=='':
                        lenSub -= 1
                if lenSub==lenStrc==lenDBref:
                    for name, inchi, dbrefs in zip(row['product_name'].split(';'),
                            row['product_structure'].split('_'),
                            row['product_dbref'].split(';')):
                        tmp = {}
                        if remove_inchi_4p:
                            tmp['inchi'] = '/'.join([inchi.split('/')[i] for i in range(len(inchi.split('/'))) if i<4])
                        else:
                            tmp['inchi'] = inchi.replace(' ', '')
                        tmp['name'] = name
                        tmp['dbref'] = {}
                        for dbref in dbrefs.split('|'):
                            if len(dbref.split(':'))==2:
                                db_name = dbref.split(':')[0].replace(' ', '').lower()
                                db_cid = dbref.split(':')[1].replace(' ', '')
                                if not db_name in tmp['dbref']:
                                    tmp['dbref'][db_name] = []
                                tmp['dbref'][db_name].append(db_cid)
                            else:
                                data[pathID]['isValid'] = False
                                self.logger.warning('Ignoring the folowing product dbref ('+str(name)+'): '+str(dbref))
                        data[pathID]['steps'][stepID]['products'].append(tmp)
                else:
                    self.logger.warning('Not equal length between substrate names, their structure or dbref ('+str(name)+'): '+str(row['product_name'])+' <--> '+str(row['product_structure'])+' <--> '+str(row['product_dbref']))
                    data[pathID]['isValid'] = False
                if not row['uniprot']=='':
                    data[pathID]['steps'][stepID]['uniprot'] = row['uniprot'].replace(' ', '').split(';')
                if not row['EC_number']=='':
                    data[pathID]['steps'][stepID]['ec_numbers'] = [i.replace(' ', '') for i in row['EC_number'].split(';')]
                data[pathID]['steps'][stepID]['enzyme_id'] = [i.replace(' ', '') for i in row['enzyme_identifier'].split(';')]
                data[pathID]['steps'][stepID]['enzyme_name'] = row['enzyme_name'].split(';')
        except FileNotFoundError:
            self.logger.error('Cannot open the file: '+str(inFile))
        #now loop through all of them and remove the invalid paths
        toRet = deepcopy(data)
        for path_id in data.keys():
            if toRet[path_id]['isValid']==False:
                del toRet[path_id]
            else:
                del toRet[path_id]['isValid']
        #reorganise the results around the target products mnx
        if not mnxHeader:
            return toRet
        else:
            toRetTwo = {}
            for path_id in toRet:
                try:
                    final_pro_mnx = toRet[path_id]['steps'][max(toRet[path_id]['steps'])]['products'][0]['dbref']['mnx'][0]
                except KeyError:
                    self.logger.error('The species '+str(toRet[path_id]['steps'][max(toRet[path_id]['steps'])]['products'][0]['name'])+' does not contain a mnx database reference... skipping whole pathway number '+str(path_id))
                    #continue
                if not final_pro_mnx in toRetTwo:
                    toRetTwo[final_pro_mnx] = {}
                toRetTwo[final_pro_mnx][path_id] = toRet[path_id]
            return toRetTwo


    ## Parse the validation TSV to SBML
    #
    # Parse the TSV file to SBML format and adds them to the self.sbml_paths
    #
    # @param self Object pointer
    # @param inFile Input file
    # @param compartment_id compartment of the
    def TSVtoSBML(self,
                  inFile,
                  tmpOutputFolder=None,
                  upper_flux_bound=99999,
                  lower_flux_bound=0,
                  compartment_id='MNXC3',
                  pathway_id='rp_pathway',
                  species_group_id='central_species'):
        data = self._parseTSV(inFile)
        sbml_paths = {}
        header_name = inFile.split('/')[-1].replace('.tsv', '').replace('.csv', '')
        #TODO: need to exit at this loop
        for path_id in data:
            try:
                mnxc = self.name_compXref[compartment_id]
            except KeyError:
                self.logger.error('Could not Xref compartment_id ('+str(compartment_id)+')')
                return False
            rpsbml = rpSBML.rpSBML(header_name+'_'+str(path_id))
            #1) create a generic Model, ie the structure and unit definitions that we will use the most
            ##### TODO: give the user more control over a generic model creation:
            #   -> special attention to the compartment
            rpsbml.genericModel(header_name+'_Path'+str(path_id),
                                header_name+'_Path'+str(path_id),
                                self.compXref[mnxc],
                                compartment_id,
                                upper_flux_bound,
                                lower_flux_bound)
            #2) create the pathway (groups)
            rpsbml.createPathway(pathway_id)
            rpsbml.createPathway(species_group_id)
            #3) find all the unique species and add them to the model
            allChem = []
            for stepNum in data[path_id]['steps']:
                #because of the nature of the input we need to remove duplicates
                for i in data[path_id]['steps'][stepNum]['substrates']+data[path_id]['steps'][stepNum]['products']:
                    if not i in allChem:
                        allChem.append(i)
            #add them to the SBML
            for chem in allChem:
                #PROBLEM: as it stands one expects the meta to be MNX
                if 'mnx' in chem['dbref']:
                    #must list the different models
                    meta = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                else:
                    #TODO: add the species with other types of xref in annotation
                    self.logger.warning('Some species are not referenced by a MNX id and will be ignored')
                    #try CHEBI
                    try:
                        meta = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
                        meta = 'CHEBI_'+str(meta)
                    except KeyError:
                        #TODO: need to find a better way
                        self.logger.warning('Cannot determine MNX or CHEBI entry, using random')
                        tmpDB_name = list(chem['dbref'].keys())[0]
                        meta = chem['dbref'][list(chem['dbref'].keys())[0]][0]
                        meta = str(tmpDB_name)+'_'+str(meta)
                    #break
                #try to conver the inchi into the other structures
                smiles = None
                inchikey = None
                try:
                    resConv = self._convert_depiction(idepic=chem['inchi'], itype='inchi', otype={'smiles','inchikey'})
                    smiles = resConv['smiles']
                    inchikey = resConv['inchikey']
                except NotImplementedError as e:
                    self.logger.warning('Could not convert the following InChI: '+str(chem['inchi']))
                #create a new species
                #here we want to gather the info from rpCompletion's rp_strc and mnxm_strc
                try:
                    chemName = self.mnxm_strc[meta]['name']
                except KeyError:
                    chemName = meta
                #compile as much info as you can
                #xref
                try:
                    #TODO: add the xref from the document
                    spe_xref = self.chemXref[meta]
                except KeyError:
                    #spe_xref = {}
                    spe_xref = chem['dbref']
                #inchi
                try:
                    spe_inchi = self.mnxm_strc[meta]['inchi']
                except KeyError:
                    spe_inchi = chem['inchi']
                #inchikey
                try:
                    spe_inchikey = self.mnxm_strc[meta]['inchikey']
                except KeyError:
                    spe_inchikey =  resConv['inchikey']
                #smiles
                try:
                    spe_smiles = self.mnxm_strc[meta]['smiles']
                except KeyError:
                    spe_smiles = resConv['smiles']
                #pass the information to create the species
                rpsbml.createSpecies(meta,
                        compartment_id,
                        chemName,
                        spe_xref,
                        spe_inchi,
                        spe_inchikey,
                        spe_smiles,
                        species_group_id)
            #4) add the complete reactions and their annotations
            #create a new group for the measured pathway
            #need to convert the validation to step for reactions
            for stepNum in data[path_id]['steps']:
                toSend = {'left': {}, 'right': {}, 'rule_id': None, 'rule_ori_reac': None, 'rule_score': None, 'path_id': path_id, 'step': stepNum, 'sub_step': None}
                for chem in data[path_id]['steps'][stepNum]['substrates']:
                    if 'mnx' in chem['dbref']:
                        meta = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                        #try CHEBI
                    else:
                        self.logger.warning('Not all the species to have a MNX ID')
                        #break
                        try:
                            meta = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
                            meta = 'CHEBI_'+str(meta)
                        except KeyError:
                            #TODO: need to find a better way
                            self.logger.warning('Cannot determine MNX or CHEBI entry, using random')
                            tmpDB_name = list(chem['dbref'].keys())[0]
                            meta = chem['dbref'][list(chem['dbref'].keys())[0]][0]
                            meta = str(tmpDB_name)+'_'+str(meta)
                    toSend['left'][meta] = 1
                for chem in data[path_id]['steps'][stepNum]['products']:
                    if 'mnx' in chem['dbref']:
                        meta = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                        #try CHEBI
                    else:
                        self.logger.warning('Need all the species to have a MNX ID')
                        try:
                            meta = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
                            meta = 'CHEBI_'+str(meta)
                        except KeyError:
                            #TODO: need to find a better way
                            self.logger.warning('Cannot determine MNX or CHEBI entry, using random')
                            tmpDB_name = list(chem['dbref'].keys())[0]
                            meta = chem['dbref'][list(chem['dbref'].keys())[0]][0]
                            meta = str(tmpDB_name)+'_'+str(meta)
                    toSend['right'][meta] = 1
                        #break
                #if all are full add it
                reac_xref = {}
                if 'ec_numbers' in data[path_id]['steps'][stepNum]:
                    reac_xref['ec'] = data[path_id]['steps'][stepNum]['ec_numbers']
                if 'uniprot' in data[path_id]['steps'][stepNum]:
                    reac_xref['uniprot'] = data[path_id]['steps'][stepNum]['uniprot']
                rpsbml.createReaction(header_name+'_Step'+str(stepNum),
                                      upper_flux_bound,
                                      lower_flux_bound,
                                      toSend,
                                      compartment_id,
                                      None,
                                      reac_xref,
                                      pathway_id)
                if stepNum==1:
                    #adding the consumption of the target
                    targetStep = {'rule_id': None,
                            'left': {},
                            'right': {},
                            'step': None,
                            'sub_step': None,
                            'path_id': None,
                            'transformation_id': None,
                            'rule_score': None,
                            'rule_ori_reac': None}
                    for chem in data[path_id]['steps'][stepNum]['products']:
                        try:
                            #smallest MNX
                            meta = sorted(chem['dbref']['mnx'], key=lambda x : int(x.replace('MNXM', '')))[0]
                        except KeyError:
                            #try CHEBI
                            try:
                                meta = sorted(chem['dbref']['chebi'], key=lambda x : int(x))[0]
                                meta = 'CHEBI_'+str(meta)
                            except KeyError:
                                self.logger.warning('Cannot determine MNX or CHEBI entry, using random')
                                tmpDB_name = list(chem['dbref'].keys())[0]
                                meta = chem['dbref'][list(chem['dbref'].keys())[0]][0]
                                meta = str(tmpDB_name)+'_'+str(meta)
                        targetStep['left'][meta] = 1
                    rpsbml.createReaction(header_name+'_Step1_sink',
                                          upper_flux_bound,
                                          lower_flux_bound,
                                          targetStep,
                                          compartment_id)
                    rpsbml.createFluxObj('rpFBA_obj', header_name+'_Step1_sink', 1, True)
            if tmpOutputFolder:
                rpsbml.writeSBML(tmpOutputFolder)
            else:
                sbml_paths[header_name+'_Path'+str(path_id)] = rpsbml
        if tmpOutputFolder:
            return {}
        else:
            return sbml_paths
