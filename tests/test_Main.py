"""
Created on Jul 15 2020

@author: Joan Hérisson
"""

from module   import Module
from brs_libs import rpCache
from brs_libs import rpSBML
from rpcompletion import rp2ToSBML

import tempfile
import glob
import os


class Test_Main(Module):
    __test__ = True

    def setUp(self):
        self.args.cache = rpCache('file')

    # def _postexec(self):
    #     for file, hash in self.files:
    #         Module._sort_file(file)

    def test_rp2ToSBML(self):
        with tempfile.TemporaryDirectory() as tmp_output_folder:
            rp2ToSBML(self.args.cache,
                      'data/rp_pathways.csv',
                      'data/rp2paths_compounds.tsv',
                      'data/rp2paths_pathways.csv',
                      tmp_output_folder)
            #### check the SBML generation
            #should only be a single file
            all_files_results = glob.glob(os.path.join(tmp_output_folder, '*.xml'))
            self.assertTrue(len(all_files_results)==1)
            #should have a single reaction
            rpsbml = rpSBML(name='test', inFile=all_files_results[0])
            self.assertListEqual(rpsbml.readRPpathwayIDs(), ['RP1'])
            #### check the reaction
            model = rpsbml.getModel()
            reac = model.getReaction('RP1')
            self.assertIsNotNone(reac)
            self.assertListEqual(sorted([i.species for i in reac.getListOfProducts()]), sorted(['TARGET_0000000001__64__MNXC3', 'MNXM11__64__MNXC3']))
            self.assertListEqual([i.species for i in reac.getListOfReactants()], ['MNXM100__64__MNXC3'])
            ##### check the species
            self.assertListEqual(sorted([i.getId() for i in model.getListOfSpecies()]), sorted(['TARGET_0000000001__64__MNXC3', 'MNXM100__64__MNXC3', 'MNXM11__64__MNXC3']))
            #### check the reaction rules
            self.assertDictEqual(rpsbml.readRPrules(), {'RR-02-ce58f652d300d037-16-F': '[H]C([H])=C(C([H])([H])[H])C1([H])C([H])([H])C([H])=C(C([H])([H])[H])C([H])([H])C1([H])[H].O=P(O)(O)OP(=O)(O)O>>[H]OP(=O)(O[H])OP(=O)(O[H])OC([H])([H])C([H])=C(C([H])([H])[H])C([H])([H])C([H])([H])C([H])=C(C([H])([H])[H])C([H])([H])[H]'})
