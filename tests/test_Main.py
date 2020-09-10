"""
Created on Jul 15 2020

@author: Joan Hérisson
"""

from module   import Module
from brs_libs import rpCache


class Test_Main(Module):
    __test__ = True

    def setUp(self):
        self.args.cache = rpCache('file')

    # def _postexec(self):
    #     for file, hash in self.files:
    #         Module._sort_file(file)
