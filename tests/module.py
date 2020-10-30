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
    data_path = 'data/lycopene'
    args      = argparse_Namespace()
    # Args have to be set in the order as they be passed
    args.cache              = None
    args.rp2_pathways       = data_path+'/1-rp2_pathways.csv'
    args.rp2paths_compounds = data_path+'/2-rp2paths_compounds.tsv'
    args.rp2paths_pathways  = data_path+'/3-rp2paths_pathways.csv'
    args.outdir             = tempdir.name

    # Useless to sort files since smiles could be equivalent and not equal, then chacksum will be different
    def _check(self):
        for file, size in self.files:
            self.assertTrue(os_path.isfile(file))
            self.assertEqual(os_stat(file).st_size, size)

    files = [
    (args.outdir+'/'+'rp_1_11_sbml.xml',  32358),
    (args.outdir+'/'+'rp_1_1_sbml.xml',   32640),
    (args.outdir+'/'+'rp_1_6_sbml.xml',   32225),
    (args.outdir+'/'+'rp_2_12_sbml.xml',  32355),
    (args.outdir+'/'+'rp_2_22_sbml.xml',  32483),
    (args.outdir+'/'+'rp_2_2_sbml.xml',   32765),
    (args.outdir+'/'+'rp_3_10_sbml.xml',  33763),
    (args.outdir+'/'+'rp_3_131_sbml.xml', 34551),
    (args.outdir+'/'+'rp_3_132_sbml.xml', 34814),
    (args.outdir+'/'+'rp_3_140_sbml.xml', 33353),
    (args.outdir+'/'+'rp_3_1_sbml.xml',   34956),
    (args.outdir+'/'+'rp_3_261_sbml.xml', 34680),
    (args.outdir+'/'+'rp_3_262_sbml.xml', 34942),
    (args.outdir+'/'+'rp_3_270_sbml.xml', 33482),
    (args.outdir+'/'+'rp_3_2_sbml.xml',   35220),
    ]
