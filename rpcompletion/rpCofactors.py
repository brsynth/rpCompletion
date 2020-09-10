import logging

from sys import argv as sys_argv
from os import path as os_path
from os import mkdir as os_mkdir
from copy import deepcopy
from argparse import ArgumentParser as argparse_ArgParser

from brs_libs import rpSBML, rpCache, rpCache_add_args


def add_arguments(parser):
    parser = rpCache_add_args(parser)
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
    def __init__(self, db='file'):
        super().__init__(db)
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
    # @param pathway_cmp Dictionnary used to retreive the public ID of the intermediate compounds. Resets for each individual pathway
    #
    def completeReac(self, step, rr_reac, full_reac, mono_side, rr_string, pathway_cmp):

        if mono_side:
            ## add the unknown species to pathway_cmp for the next steps
            rr_mono_cmp = list(rr_reac.keys())
            step_mono_cmp = list(step.keys())
            if (len(rr_mono_cmp)==1 and len(step_mono_cmp)==1):
                #this is purposely overwitten since the main cmp between reactions can change
                pathway_cmp[step_mono_cmp[0]] = rr_mono_cmp[0]
            else:
                self.logger.warning('There should be only one compound on the left for monocomponent reaction: rr_mono_cmp: '+str(rr_mono_cmp)+' step_mono_cmp: '+str(step_mono_cmp))
                return False

        ## add the side species
        rr_string += self.add_side_species(step, full_reac, rr_reac)

        ## Update the stochio
        return True, self.update_stochio(step, full_reac, rr_string, pathway_cmp)


    def add_side_species(self, step, full_reac, rr_reac):
        rr_string = ''
        for toAdd in full_reac.keys()-rr_reac.keys():
            step.update({toAdd: full_reac[toAdd]})
            ### update the reaction rule string
            try:
                smi = self.cid_strc[toAdd]['smiles']
                if not smi==None:
                    for sto_add in range(int(full_reac[toAdd])):
                        rr_string += '.'+str(smi)
            except KeyError:
                self.logger.warning('Cannot find smiles structure for '+str(toAdd))
        return rr_string


    def update_stochio(self, step, full_reac, rr_string, pathway_cmp):
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
                            smi = self.cid_strc[step_spe]['smiles']
                            if not smi==None:
                                rr_string += '.'+str(smi)
                        except KeyError:
                            #@Mel toAdd -> step_spe
                            self.logger.warning('Cannot find smiles structure for '+str(step_spe))
            elif step_spe in pathway_cmp:
                if pathway_cmp[step_spe] in full_reac:
                    if not step[step_spe]==full_reac[pathway_cmp[step_spe]]:
                        step[step_spe] = full_reac[pathway_cmp[step_spe]]
            #Its fine if the stochio is not updated, better than ignoring a whole pathway
                #else:
                #    self.logger.warning('Cannot find '+str(step_spe)+' in full reaction')
                #    return False
            #else:
            #    self.logger.warning('Cannot find '+str(step_spe)+' in pathway_cmp')
            #    return False
            return rr_string


    ## Add the cofactors to monocomponent reactions
    #
    # @param step Step in a pathway
    # @param pathway_cmp Dictionnary of intermediate compounds with their public ID's
    # @return Boolean determine if the step is to be added
    def addCofactors_step(self, step, pathway_cmp):
        reac_smiles_left = step['reaction_rule'].split('>>')[0]
        reac_smiles_right = step['reaction_rule'].split('>>')[1]
        if self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['rel_direction']==-1:
            try:
                isSuccess, reac_smiles_left = self.completeReac(step['right'],
                        self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['left'],
                        self.rr_full_reactions[self._checkRIDdeprecated(step['rule_ori_reac'], self.deprecatedRID_rid)]['right'],
                        True,
                        reac_smiles_left,
                        pathway_cmp)
                if not isSuccess:
                    self.logger.warning('Could not recognise reaction rule for step (1): '+str(step))
                    return False
            except KeyError:
                self.logger.warning('Could not find the full reaction for reaction (1): '+str(step))
                return False
            try:
                isSuccess, reac_smiles_right = self.completeReac(step['left'],
                        self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['right'],
                        self.rr_full_reactions[self._checkRIDdeprecated(step['rule_ori_reac'], self.deprecatedRID_rid)]['left'],
                        False,
                        reac_smiles_right,
                        pathway_cmp)
                if not isSuccess:
                    self.logger.warning('Could not recognise reaction rule for step (2): '+str(step))
                    return False
            except KeyError:
                self.logger.warning('Could not find the full reaction for reaction (2): '+str(step))
                return False
        elif self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['rel_direction']==1:
            try:
                isSuccess, reac_smiles_left = self.completeReac(step['right'],
                        self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['left'],
                        self.rr_full_reactions[self._checkRIDdeprecated(step['rule_ori_reac'], self.deprecatedRID_rid)]['left'],
                        True,
                        reac_smiles_left,
                        pathway_cmp)
                if not isSuccess:
                    self.logger.error('Could not recognise reaction rule for step (3): '+str(step))
                    return False
            except KeyError:
                self.logger.warning('Could not find the full reaction for reaction (3): '+str(step))
                return False
            try:
                isSuccess, reac_smiles_right = self.completeReac(step['left'],
                        self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['right'],
                        self.rr_full_reactions[self._checkRIDdeprecated(step['rule_ori_reac'], self.deprecatedRID_rid)]['right'],
                        False,
                        reac_smiles_right,
                        pathway_cmp)
                if not isSuccess:
                    self.logger.error('Could not recognise reaction rule for step (4): '+str(step))
                    return False
            except KeyError:
                self.logger.warning('Could not find the full reaction for reaction (4): '+str(step))
                return False
        else:
            self.logger.error('Relative direction can only be 1 or -1: '+str(self.rr_reactions[step['rule_id']][step['rule_ori_reac']]['rel_direction']))
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
        pathway_cmp = {}
        spe_conv = {}
        rpsbml_json = rpsbml.genJSON(pathway_id)
        rp_path = rpsbml.outPathsDict(pathway_id)
        ori_rp_path = deepcopy(rp_path)
        #We reverse the loop to ID the intermediate CMP to their original ones
        for stepNum in sorted(list(rp_path), reverse=True):
        #for stepNum in sorted(list(rp_path)):
            if self.addCofactors_step(rp_path[stepNum], pathway_cmp):
                ###add the new cofactors to the SBML
                #remove the original species from the monocomponent reaction
                reactants = set(set(rp_path[stepNum]['left'].keys())-set(ori_rp_path[stepNum]['left'].keys()))
                products = set(set(rp_path[stepNum]['right'].keys())-set(ori_rp_path[stepNum]['right'].keys()))
                for species in reactants|products:
                    tmp_species = self._checkCIDdeprecated(species, self.deprecatedCID_cid)
                    #check to make sure that they do not yet exist and if not create a new one
                    #TODO, replace the species with an existing one if it is contained in the MIRIAM annotations
                    if not rpsbml.speciesExists(tmp_species, compartment_id):
                        xref = {}
                        inchi = None
                        inchikey = None
                        smiles = None
                        chem_name = None
                        ###### Try to retreive the InChI ############
                        try:
                            inchi = self.cid_strc[tmp_species]['inchi']
                        except KeyError:
                            self.logger.warning('Cannot find the inchi for this species: '+str(tmp_species))
                        try:
                            inchikey = self.cid_strc[tmp_species]['inchikey']
                            #self.logger.debug('Found the inchikey: '+str(inchikey))
                            #### TODO: find a better way to check if two species are the same ####
                            isfound = False
                            for rpsbml_species in rpsbml_json['species']:
                                #TODO add a comparison by xref as well
                                #self.logger.debug(str(rpsbml_json['species'][rpsbml_species]['brsynth']['inchikey'])+' <--> '+str(inchikey))
                                if str(rpsbml_json['species'][rpsbml_species]['brsynth']['inchikey'])==str(inchikey):
                                    spe_conv[tmp_species] = rpsbml_species
                                    self.logger.debug('The species '+str(tmp_species)+' is the same as '+str(rpsbml_species))
                                    isfound = True
                                    break
                            if isfound:
                                continue
                        except KeyError:
                            self.logger.warning('Cannot find the inchikey for this species: '+str(species))
                        ##### Try to retreive the SMILES ############
                        try:
                            smiles = self.cid_strc[tmp_species]['smiles']
                        except KeyError:
                            self.logger.warning('Cannot find the smiles for this species: '+str(species))
                        ###### Try to retreive the xref, using the inchikey if the cid fails #######
                        try:
                            xref = self.cid_xref[tmp_species]
                        except KeyError:
                            try:
                                xref = self.cid_xref[tmp_species]
                            except KeyError:
                                #if you cannot find using cid, try to retreive it using its inchikey
                                try:
                                    if inchikey:
                                        #@Joan: Can you think of a better way of doing that?
                                        # WARNING here we use MNX since as of now, there are only MNX data that is parsed correctly
                                        tmp_cids = [i for i in self.inchikey_cid[inchikey] if i[:3]=='MNX']
                                        #TODO: handle multiple matches. For now we assume that multiple MNX means that there are deprecated versions of the tool
                                        if tmp_cids:
                                            xref = self.cid_xref[self._checkCIDdeprecated(tmp_cids[0], self.deprecatedCID_cid)]
                                except KeyError:
                                    self.logger.warning('Cannot find the xref for this species: '+str(species))
                                    xref = {}
                        #### Common Name ####
                        try:
                            chem_name = self.cid_name[self._checkCIDdeprecated(tmp_species, self.deprecatedCID_cid)]
                        except KeyError:
                            #if you cannot find using cid, try to retreive it using its inchikey
                            try:
                                if inchikey:
                                    #@Joan: Same question as above
                                    tmp_cids = [i for i in self.inchikey_cid[inchikey] if i[:3]=='MNX']
                                    if tmp_cids:
                                        chem_name = self.cid_name[self._checkCIDdeprecated(tmp_cids[0], self.deprecatedCID_cid)]
                            except KeyError:
                                self.logger.warning('Cannot find the name for this species: '+str(species))
                        #### Finally create the species in the SBML file ######
                        rpsbml.createSpecies(tmp_species,
                                compartment_id,
                                chem_name,
                                xref,
                                inchi,
                                inchikey,
                                smiles)
                #add the new species to the RP reactions
                reac = rpsbml.model.getReaction(rp_path[stepNum]['reaction_id'])
                pre_reactants = [i.species for i in reac.getListOfReactants()]
                pre_products = [i.species for i in reac.getListOfProducts()]
                for pro in products:
                    if self._checkCIDdeprecated(pro, self.deprecatedCID_cid) in spe_conv:
                        toadd = spe_conv[self._checkCIDdeprecated(pro, self.deprecatedCID_cid)]
                    else:
                        toadd = str(self._checkCIDdeprecated(pro, self.deprecatedCID_cid))+'__64__'+str(compartment_id)
                    #prod.setSpecies(str(self._checkCIDdeprecated(pro))+'__64__'+str(compartment_id))
                    if toadd in pre_products:
                        continue
                    prod = reac.createProduct()
                    prod.setSpecies(toadd)
                    prod.setConstant(True)
                    prod.setStoichiometry(rp_path[stepNum]['right'][pro])
                for sub in reactants:
                    if self._checkCIDdeprecated(sub, self.deprecatedCID_cid) in spe_conv:
                        toadd = spe_conv[self._checkCIDdeprecated(sub, self.deprecatedCID_cid)]
                    else:
                        toadd = str(self._checkCIDdeprecated(sub, self.deprecatedCID_cid))+'__64__'+str(compartment_id)
                    #prod.setSpecies(str(self._checkCIDdeprecated(sub))+'__64__'+str(compartment_id))
                    if toadd in pre_reactants:
                        continue
                    subs = reac.createReactant()
                    subs.setSpecies(toadd)
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
