"""
Created on Jul 15 2020

@author: Joan Hérisson
"""

from _main import Main
from tempfile import TemporaryDirectory
from os import path as os_path
from os import stat as os_stat

tempdir = TemporaryDirectory()

hashes = [
(	tempdir.name+'/'+'rp_10_1_sbml.xml', '	ac311bdd48e303d57069536a0823a69e552e0366cbfd8a3b260b7dedcf055bac	'),
(	tempdir.name+'/'+'rp_11_1_sbml.xml', '	0eae89e7572caf599eec7a1e285f1e45408605cc8258fe064c9ac5326913b944	'),
(	tempdir.name+'/'+'rp_11_2_sbml.xml', '	a9db0bb4bb4b929442a742e95b06985efac7fc122432f28db4a9487f8b15896c	'),
(	tempdir.name+'/'+'rp_12_1_sbml.xml', '	45cb4db19e052dd97ec736aadc94d6e3b81cca8415a4aefb2ce6b5d5ae0c5166	'),
(	tempdir.name+'/'+'rp_13_1_sbml.xml', '	45cb4db19e052dd97ec736aadc94d6e3b81cca8415a4aefb2ce6b5d5ae0c5166	'),
(	tempdir.name+'/'+'rp_14_1_sbml.xml', '	100c43c49265fab898a6daecd279f0d32e4c8bcafa18a2a54c3c6e08588f9266	'),
(	tempdir.name+'/'+'rp_15_1_sbml.xml', '	f9ad246dff7186f90f544184d2510c12b1c7a32dba3e15191648cbbfd4cbea00	'),
(	tempdir.name+'/'+'rp_16_1_sbml.xml', '	39cf7ee1c946dd53bc46ab792bca761dfb5b5ea742865b3fc891050aaca4a70f	'),
(	tempdir.name+'/'+'rp_16_2_sbml.xml', '	e1c525d0890b73852db1828c1513db13de2337864f2c10dbb5d2d339077d8fad	'),
(	tempdir.name+'/'+'rp_17_1_sbml.xml', '	a62a88f1d2ad538830f2a4c7ef5c30b88dfde5f2289fdca000f555a4efa2c1fc	'),
(	tempdir.name+'/'+'rp_18_1_sbml.xml', '	0ae1d4e3316a615b5524a0a0a17ce3319c5e7024fcba948576d4a9d2201c315e	'),
(	tempdir.name+'/'+'rp_19_1_sbml.xml', '	8f5b127c50f1b168da2d3b533e12696051c658034004a683ad6e3f599ca1fe77	'),
(	tempdir.name+'/'+'rp_1_1_sbml.xml', '	142f798cfde6889c9f183784215710dfd9e11ccf9f00a15d83fdad8a5b2d52d3	'),
(	tempdir.name+'/'+'rp_1_2_sbml.xml', '	b03e8df789dd270fc82d1b8edd75f3aa1b77ae575d72bbc79401f89f3f2dfa7c	'),
(	tempdir.name+'/'+'rp_20_1_sbml.xml', '	0c1eb43ced5de2f16a9e2526e6dd24a9639a1fe170ab8c8c64efb7a2c6ef0947	'),
(	tempdir.name+'/'+'rp_21_1_sbml.xml', '	9b83c89bf8b318b844cfdc3bce5a822893b5e3935be958b44de7fd2dfef51585	'),
(	tempdir.name+'/'+'rp_22_1_sbml.xml', '	0d37b4015bf067d5a3ab46f392e7daa74ebdfdf4e81c7ed85a880d4d0e1cee24	'),
(	tempdir.name+'/'+'rp_23_1_sbml.xml', '	fed08979a62deeb2305fade558efd3c2a0eb092fd0f50ddf82dd9ab99f964f87	'),
(	tempdir.name+'/'+'rp_23_2_sbml.xml', '	a0dc303b7facec9279f25cffb5f0b339e12bef7024b21091f7b46576b597cfa4	'),
(	tempdir.name+'/'+'rp_24_1_sbml.xml', '	2e5b17239067609bcda3b690c3599bb55a1eacaef8cca51d563b82fbfb8e5c44	'),
(	tempdir.name+'/'+'rp_24_2_sbml.xml', '	f3ddc301a94614d10707015ddd06d26223faadb8b39de99fe275ec4731e382fd	'),
(	tempdir.name+'/'+'rp_25_1_sbml.xml', '	1ddcd1974e0887109f91108ee8f6e63025c9a4c4006a15a5205eccd88eceeb82	'),
(	tempdir.name+'/'+'rp_25_2_sbml.xml', '	3d9ce5448d6ca69634e907afb41652a2ba296937583784558f2d727b877590ac	'),
(	tempdir.name+'/'+'rp_26_1_sbml.xml', '	cee4dc9632721e600bcb68309e79fc40af48ce845c1332b7a86f2b642c1c40ad	'),
(	tempdir.name+'/'+'rp_27_1_sbml.xml', '	6c7255757549e06c209bf98cd1836a7a77f0eb1139401002d0492ec05963d61a	'),
(	tempdir.name+'/'+'rp_28_1_sbml.xml', '	902c513668275e3c4b9b8f00f4abe3ebfebbd935f4bb153a51df309fdbdaa691	'),
(	tempdir.name+'/'+'rp_29_1_sbml.xml', '	e6e1b8cb1c7f4a4dd7427e630170de306204d5fd2880bb4ee1dea1e1558446bf	'),
(	tempdir.name+'/'+'rp_2_1_sbml.xml', '	841ac16b98f6a6c088b580895649610e27261a2376b0d6e1873dcc8f51991ff4	'),
(	tempdir.name+'/'+'rp_30_1_sbml.xml', '	e4cb9bc50fb8aa98aaae417c99a61128ba5913777e951c756c5c3d72378c67d9	'),
(	tempdir.name+'/'+'rp_3_1_sbml.xml', '	e2ddbc66d6bd1ab392de2ea4287f923cdab2f4e2fd80c56a8d378bb5014ab5c6	'),
(	tempdir.name+'/'+'rp_4_1_sbml.xml', '	870713da2635e3f7701bd18a91bc81760bb30f3b7f8c63fd85957a7af77d2f45	'),
(	tempdir.name+'/'+'rp_5_1_sbml.xml', '	67b2f9e271f9ce050a10bd5c0271120ddcadf5d558739f0dc62e1c5d14efe94c	'),
(	tempdir.name+'/'+'rp_5_2_sbml.xml', '	e6fc85a2c066690af9027bbbd87f3e0d2e07a4c8210aa1881043f7b2f3e5b2b7	'),
(	tempdir.name+'/'+'rp_6_1_sbml.xml', '	0cadfdf08e45bbb9b56153bb454e4b038d401de5e968517e4a98ab6272646312	'),
(	tempdir.name+'/'+'rp_7_1_sbml.xml', '	b4626b41bf06f8e56d3c9acebab59e7deaf704a238e8ce9ee8088bb02dbca508	'),
(	tempdir.name+'/'+'rp_8_1_sbml.xml', '	a89da14bb19c35cc26ba16a39a9f57b21f06dea67fd1570e050d1421454b8909	'),
(	tempdir.name+'/'+'rp_8_2_sbml.xml', '	3468ea65149085536e2f959d8ddea52129dbb78d101e9b5b1c9f564cf8cdf31f	'),
(	tempdir.name+'/'+'rp_9_1_sbml.xml', '	b74fb3b6f96863d44d63606334eb54e7c60c5f9b676baa43e95c2739f9ed1c7a	')
]

class Module(Main):
    __test__ = False

    mod_name  = 'rpcompletion'
    func_name = 'rp2ToSBML'
    data_path = 'data'
    args      = [data_path+'/rp2_pathways.csv',
                 data_path+'/rp2paths_compounds.csv',
                 data_path+'/rp2paths_pathways.csv',
                 tempdir.name]

    hashes = hashes

    def _check(self):
        for file, hash in hashes:
            self.assertTrue(os_path.isfile(file))
            self.assertGreater(os_stat(file).st_size, 0)
