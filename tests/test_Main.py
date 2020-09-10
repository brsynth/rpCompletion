"""
Created on Jul 15 2020

@author: Joan Hérisson
"""

from module       import Module
from rpcompletion import rpCompletion


class Test_Main(Module):
    __test__ = True

    obj = rpCompletion('file')

    # def _postexec(self):
    #     for file, hash in self.files:
    #         Module._sort_file(file)
