import logging

from sys import path as sys_path
from sys import argv as sys_argv
from os import path as os_path
from os import mkdir as os_mkdir
from glob import glob as glob_glob
from copy import deepcopy as copy_deepcopy
from argparse import ArgumentParser as argparse_ArgParser

from brs_utils import rpSBML
from .rpCache import rpCache, add_arguments


def _add_arguments(parser):
    parser = add_arguments(parser)
    parser.add_argument('-input_file', type=str)
    parser.add_argument('-output_dir', type=str, default='')
    parser.add_argument('-pathway_id', type=str, default='rp_pathway')
    parser.add_argument('-compartment_id', type=str, default='MNXC3')
    return parser

def build_parser():
    parser = argparse_ArgParser('Add the missing cofactors to the monocomponent reactions to the SBML outputs of rpReader')
    parser = _add_arguments(parser)
    return parser

def entrypoint(args=sys_argv[1:]):
    parser = build_parser()

    params = parser.parse_args(args)

    rpcofactors = rpCofactors(params.store_mode, params.print)
    rpcofactors.run(params.input_file,
                    params.output_dir,
                    params.pathway_id,
                    params.compartment_id)



if __name__ == "__main__":

    entrypoint(sys_argv[1:])

## Class to add the cofactors to a monocomponent reaction to construct the full reaction
#
#
class rpCofactors(rpCache):
    ## Init method
    # Here we want to seperate what is the use input and what is parsed by the cache to make sure that
    # everything is not hadled by a single
    #
    # @param rpReader input reader object with the parsed user input and cache files required
    #DEPRECATED def __init__(self, rpReader, userXrefDbName=None):
    def __init__(self, db='file', print_infos=False):
        super().__init__(db, print_infos)
        self.logger = logging.getLogger(__name__)
        self.logger.info('Started instance of rpCofactors')

    ################################################################
    ######################### PUBLIC FUNCTIONS #####################
    ################################################################


    ## Given a dictionnary describing a monocomponent reaction, add the cofactors by comparing it with the original reaction
    #
    # @param step Dictionnary describing the reaction
    # @param reac_side String 'right' or 'left' describing the direction of the monocomponent reaction compared with the original reaction
    # @param rr_reac Dictionnary describing the monocomponent reaction from RetroRules
    # @param f_reac Dictionnary describing the full original reaction
    # @param pathway_cmp_mnxm Dictionnary used to retreive the public ID of the intermediate compounds. Resets for each individual pathway
    #
    def completeReac(self, step, rr_reac, full_reac, mono_side, rr_string, pathway_cmp_mnxm):
        if mono_side:
            ## add the unknown species to pathway_cmp_mnxm for the next steps
            rr_mono_cmp = list(rr_reac.keys())
            step_mono_cmp = list(step.keys())
            if (len(rr_mono_cmp)==1 and len(step_mono_cmp)==1):
                #this is purposely overwitten since the main cmp between reactions can change
                pathway_cmp_mnxm[step_mono_cmp[0]] = rr_mono_cmp[0]
            else:
                self.logger.warning('There should be only one compound on the left for monocomponent reaction: rr_mono_cmp: '+str(rr_mono_cmp)+' step_mono_cmp: '+str(step_mono_cmp))
                return False
        ## add the side species
        for toAdd in full_reac.keys()-rr_reac.keys():
            step.update({toAdd: full_reac[toAdd]})
            ### update the reaction rule string
            try:
                smi = self.mnxm_strc[toAdd]['smiles']
                if not smi==None:
                    for sto_add in range(int(full_reac[toAdd])):
                        rr_string += '.'+str(smi)
            except KeyError:
                self.logger.warning('Cannot find smiles structure for '+str(toAdd))
        ## Update the the stochio
        for step_spe in step:
            if step_spe in full_reac:
                if not step[step_spe]==full_reac[step_spe]:
                    stochio_diff = full_reac[step_spe]-step[step_spe]
                    step[step_spe] = full_reac[step_spe]
                    if stochio_diff<0:
                        self.logger.warning('full_reac stochio should never be smaller than step')
                        continue
                    for i in range(stochio_diff):
                        ### update the reaction rule string
                        try:
                            smi = self.mnxm_strc[step_spe]['smiles']
                            if not smi==None:
                                rr_string += '.'+str(smi)
                        except KeyError:
                            self.logger.warning('Cannot find smiles structure for '+str(toAdd))
            elif step_spe in pathway_cmp_mnxm:
                if pathway_cmp_mnxm[step_spe] in full_reac:
                    if not step[step_spe]==full_reac[pathway_cmp_mnxm[step_spe]]:
                        step[step_spe] = full_reac[pathway_cmp_mnxm[step_spe]]
            #Its fine if the stochio is not updated, better than ignoring a whole pathway
                #else:
                #    self.logger.warning('Cannot find '+str(step_spe)+' in full reaction')
                #    return False
            #else:
            #    self.logger.warning('Cannot find '+str(step_spe)+' in pathway_cmp_mnxm')
            #    return False
        return True, rr_string


    ## Add the cofactors to monocomponent reactions
    #
    # @param step Step in a pathway
    # @param pathway_cmp_mnxm Dictionnary of intermediate compounds with their public ID's
    # @return Boolean determine if the step is to be added
    def addCofactors_step(self, step, pathway_cmp_mnxm):
        reac_smiles_left = step['reaction_rule'].split('>>')[0]
        reac_smiles_right = step['reaction_rule'].split('>>')[1]
        if self.rr_reactions[step['rule_id']][step['rule_ori_reac']['mnxr']]['rel_direction']==-1:
            isSuccess, reac_smiles_left = self.completeReac(step['right'],
                    self.rr_reactions[step['rule_id']][step['rule_ori_reac']['mnxr']]['left'],
                    self.full_reactions[self._checkMNXRdeprecated(step['rule_ori_reac']['mnxr'], self.deprecatedMNXR_mnxr)]['right'],
                    True,
                    reac_smiles_left,
                    pathway_cmp_mnxm)
            if not isSuccess:
                self.logger.error('Could not recognise reaction rule for step '+str(step))
                return False
            isSuccess, reac_smiles_right = self.completeReac(step['left'],
                    self.rr_reactions[step['rule_id']][step['rule_ori_reac']['mnxr']]['right'],
                    self.full_reactions[self._checkMNXRdeprecated(step['rule_ori_reac']['mnxr'], self.deprecatedMNXR_mnxr)]['left'],
                    False,
                    reac_smiles_right,
                    pathway_cmp_mnxm)
            if not isSuccess:
                self.logger.error('Could not recognise reaction rule for step '+str(step))
                return False
        elif self.rr_reactions[step['rule_id']][step['rule_ori_reac']['mnxr']]['rel_direction']==1:
            isSuccess, reac_smiles_left = self.completeReac(step['right'],
                    self.rr_reactions[step['rule_id']][step['rule_ori_reac']['mnxr']]['left'],
                    self.full_reactions[self._checkMNXRdeprecated(step['rule_ori_reac']['mnxr'], self.deprecatedMNXR_mnxr)]['left'],
                    True,
                    reac_smiles_left,
                    pathway_cmp_mnxm)
            if not isSuccess:
                self.logger.error('Could not recognise reaction rule for step '+str(step))
                return False
            isSuccess, reac_smiles_right = self.completeReac(step['left'],
                    self.rr_reactions[step['rule_id']][step['rule_ori_reac']['mnxr']]['right'],
                    self.full_reactions[self._checkMNXRdeprecated(step['rule_ori_reac']['mnxr'], self.deprecatedMNXR_mnxr)]['right'],
                    False,
                    reac_smiles_right,
                    pathway_cmp_mnxm)
            if not isSuccess:
                self.logger.error('Could not recognise reaction rule for step '+str(step))
                return False
        else:
            self.logger.error('Relative direction can only be 1 or -1: '+str(self.rr_reactions[step['rule_id']][step['rule_ori_reac']['mnxr']]['rel_direction']))
            return False
        step['reaction_rule'] = reac_smiles_left+'>>'+reac_smiles_right
        return True


    ## Function to reconstruct the heterologous pathway
    #
    #  Read each pathway information and RetroRules information to construct heterologous pathways and add the cofactors
    #
    #  @param self Object pointer
    #  @param rpsbml rpSBML object with a single model
    #  @return Boolean if True then you keep that model for the next step, if not then ignore it
    def addCofactors(self, rpsbml, compartment_id='MNXC3', pathway_id='rp_pathway'):
        #This keeps the IDs conversions to the pathway
        pathway_cmp_mnxm = {}
        rp_path = rpsbml.outPathsDict(pathway_id)
        ori_rp_path = copy_deepcopy(rp_path)
        #We reverse the loop to ID the intermediate CMP to their original ones
        for stepNum in sorted(list(rp_path), reverse=True):
        #for stepNum in sorted(list(rp_path)):
            if self.addCofactors_step(rp_path[stepNum], pathway_cmp_mnxm):
                ###add the new cofactors to the SBML
                #remove the original species from the monocomponent reaction
                reactants = set(set(rp_path[stepNum]['left'].keys())-set(ori_rp_path[stepNum]['left'].keys()))
                products = set(set(rp_path[stepNum]['right'].keys())-set(ori_rp_path[stepNum]['right'].keys()))
                for species in reactants|products:
                    #check to make sure that they do not yet exist and if not create a new one
                    #TODO, replace the species with an existing one if it is contained in the MIRIAM annotations
                    if not rpsbml.speciesExists(species, compartment_id):
                        xref = {}
                        inchi = None
                        inchikey = None
                        smiles = None
                        try:
                            xref = self.chemXref[species]
                        except KeyError:
                            try:
                                xref = self.chemXref[self.deprecatedMNXM_mnxm[species]]
                            except KeyError:
                                #TODO: although there should not be any
                                #intermediate species here consider
                                #removing this warning
                                self.logger.warning('Cannot find the xref for this species: '+str(species))
                                pass
                        try:
                            inchi = self.mnxm_strc[species]['inchi']
                        except KeyError:
                            try:
                                inchi = self.mnxm_strc[self.deprecatedMNXM_mnxm[species]]['inchi']
                            except KeyError:
                                self.logger.warning('Cannot find the inchi for this species: '+str(species))
                                pass
                        try:
                            inchikey = self.mnxm_strc[species]['inchikey']
                        except KeyError:
                            try:
                                inchikey = self.mnxm_strc[self.deprecatedMNXM_mnxm[species]]['inchikey']
                            except KeyError:
                                self.logger.warning('Cannot find the inchikey for this species: '+str(species))
                                pass
                        try:
                            smiles = self.mnxm_strc[species]['smiles']
                        except KeyError:
                            try:
                                smiles = self.mnxm_strc[self.deprecatedMNXM_mnxm[species]]['smiles']
                            except KeyError:
                                self.logger.warning('Cannot find the smiles for this species: '+str(species))
                                pass
                        #add the new species to rpsbml
                        try:
                            chemName = self.mnxm_strc[species]['name']
                        except KeyError:
                            chemName = None
                        rpsbml.createSpecies(species,
                                compartment_id,
                                chemName,
                                xref,
                                inchi,
                                inchikey,
                                smiles)
                #add the new species to the RP reactions
                reac = rpsbml.model.getReaction(rp_path[stepNum]['reaction_id'])
                for pro in products:
                    prod = reac.createProduct()
                    prod.setSpecies(str(pro)+'__64__'+str(compartment_id))
                    prod.setConstant(True)
                    prod.setStoichiometry(rp_path[stepNum]['right'][pro])
                for sub in reactants:
                    subs = reac.createReactant()
                    subs.setSpecies(str(sub)+'__64__'+str(compartment_id))
                    subs.setConstant(True)
                    subs.setStoichiometry(rp_path[stepNum]['left'][sub])
                #replace the reaction rule with new one
                rpsbml.addUpdateBRSynth(reac, 'smiles', rp_path[stepNum]['reaction_rule'], None, True)
            else:
                #if the cofactors cannot be found delete it from the list
                self.logger.warning('Cannot find cofactors... skipping')
                return False
        return True


    ## run using HDD 3X less than the above function
    #
    #
    def run(self, input_file, output_dir='', pathway_id='rp_pathway', compartment_id='MNXC3'):
        if not os_path.isfile(input_file):
            logging.error('Input filename does not exist.')
            return False

        modelName = os_path.basename(input_file).replace('.sbml', '').replace('.xml', '').replace('.rpsbml', '')
        rpsbml = rpSBML.rpSBML(modelName)
        rpsbml.readSBML(input_file)
        self.addCofactors(rpsbml, compartment_id, pathway_id)

        outdir = output_dir
        if outdir=='':
            outdir = os_path.dirname(input_file)
        output_file = outdir+'/'+modelName+'-completed.xml'
        if not os_path.isdir(outdir):
            os_mkdir(outdir)

        rpsbml.writeSBMLToFile(output_file)

        if os_path.getsize(output_file)==0:
            logging.error('rpCofactors has not produced any results')
            return False

        return True
