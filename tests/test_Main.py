"""
Created on Jul 15 2020

@author: Joan Hérisson
"""

from unittest import TestCase

from rpcompletion import rpCompletion, build_args_parser
from tempfile import TemporaryDirectory


# Cette classe est un groupe de tests. Son nom DOIT commencer
# par 'Test' et la classe DOIT hériter de unittest.TestCase.
class Test_Main(TestCase):

    def setUp(self):
        self.rpcompletion = rpCompletion('db')
        self.data_path = 'data'

    def test_Small(self):
        tempdir = TemporaryDirectory()
        rpsbml_paths = self.rpcompletion.rp2ToSBML(
                                 self.data_path+'/rp2_pathways.csv',
                                 self.data_path+'/rp2paths_compounds.csv',
                                 self.data_path+'/rp2paths_pathways.csv',
                                 tempdir.name)
        print(rpsbml_paths)
        tempdir.cleanup()
        self.assertEqual(True, True)
