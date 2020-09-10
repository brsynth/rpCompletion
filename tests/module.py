"""
Created on Jul 15 2020

@author: Joan Hérisson
"""

from _main    import Main
from tempfile import TemporaryDirectory
from os       import path      as os_path
from os       import stat      as os_stat
from argparse import Namespace as argparse_Namespace

tempdir = TemporaryDirectory()



class Module(Main):
    __test__ = False

    mod_name  = 'rpcompletion'
    func_name = 'rp2ToSBML'
    data_path = 'data'
    args      = argparse_Namespace()
    # Args have to be set in the order as they be passed
    args.cache              = None
    args.rp2_pathways       = data_path+'/rp2_pathways.csv'
    args.rp2paths_compounds = data_path+'/rp2paths_compounds.csv'
    args.rp2paths_pathways  = data_path+'/rp2paths_pathways.csv'
    args.outdir             = tempdir.name

    # Useless to sort files since smiles could be equivalent and not equal, then chacksum will be different
    def _check(self):
        for file, size in self.files:
            self.assertTrue(os_path.isfile(file))
            self.assertEqual(os_stat(file).st_size, size)

    files = [
    (tempdir.name+'/'+'rp_10_1_sbml.xml', 37004),
    (tempdir.name+'/'+'rp_11_1_sbml.xml', 52372),
    (tempdir.name+'/'+'rp_11_2_sbml.xml', 46612),
    (tempdir.name+'/'+'rp_12_1_sbml.xml', 47382),
    (tempdir.name+'/'+'rp_13_1_sbml.xml', 46122),
    (tempdir.name+'/'+'rp_14_1_sbml.xml', 37243),
    (tempdir.name+'/'+'rp_15_1_sbml.xml', 35632),
    (tempdir.name+'/'+'rp_16_1_sbml.xml', 46142),
    (tempdir.name+'/'+'rp_16_2_sbml.xml', 46247),
    (tempdir.name+'/'+'rp_17_1_sbml.xml', 44346),
    (tempdir.name+'/'+'rp_18_1_sbml.xml', 35589),
    (tempdir.name+'/'+'rp_19_1_sbml.xml', 35485),
    (tempdir.name+'/'+'rp_1_1_sbml.xml',  41335),
    (tempdir.name+'/'+'rp_1_2_sbml.xml',  41244),
    (tempdir.name+'/'+'rp_20_1_sbml.xml', 47233),
    (tempdir.name+'/'+'rp_21_1_sbml.xml', 47019),
    (tempdir.name+'/'+'rp_22_1_sbml.xml', 44830),
    (tempdir.name+'/'+'rp_23_1_sbml.xml', 48252),
    (tempdir.name+'/'+'rp_23_2_sbml.xml', 56646),
    (tempdir.name+'/'+'rp_24_1_sbml.xml', 47453),
    (tempdir.name+'/'+'rp_24_2_sbml.xml', 49259),
    (tempdir.name+'/'+'rp_25_1_sbml.xml', 56572),
    (tempdir.name+'/'+'rp_25_2_sbml.xml', 50000),
    (tempdir.name+'/'+'rp_26_1_sbml.xml', 54496),
    (tempdir.name+'/'+'rp_27_1_sbml.xml', 64177),
    (tempdir.name+'/'+'rp_28_1_sbml.xml', 57305),
    (tempdir.name+'/'+'rp_29_1_sbml.xml', 54528),
    (tempdir.name+'/'+'rp_2_1_sbml.xml',  27076),
    (tempdir.name+'/'+'rp_30_1_sbml.xml', 52348),
    (tempdir.name+'/'+'rp_3_1_sbml.xml',  39517),
    (tempdir.name+'/'+'rp_4_1_sbml.xml',  49009),
    (tempdir.name+'/'+'rp_5_1_sbml.xml',  50635),
    (tempdir.name+'/'+'rp_5_2_sbml.xml',  44875),
    (tempdir.name+'/'+'rp_6_1_sbml.xml',  45646),
    (tempdir.name+'/'+'rp_7_1_sbml.xml',  44385),
    (tempdir.name+'/'+'rp_8_1_sbml.xml',  43248),
    (tempdir.name+'/'+'rp_8_2_sbml.xml',  37488),
    (tempdir.name+'/'+'rp_9_1_sbml.xml',  38258)
    ]
