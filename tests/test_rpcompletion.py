"""
Created on Jul 15 2020

@author: Joan Hérisson
"""

# from rptools.rpcompletion.rpCompletion import (
#     # build_side_rxn,
#     # rp2paths_to_dict
# )
import os
from json import load as json_load
from os import path as os_path
from unittest import TestCase

from brs_utils import (
    create_logger,
)
from rplibs import rpPathway
from rr_cache import rrCache

from rpcompletion.completion import (
    __infer_separator as infer_separator,
)
from rpcompletion.completion import (
    __infer_separator_from_content as infer_separator_from_content,
)
from rpcompletion.completion import (
    __infer_separator_from_extension as infer_separator_from_extension,
)
from rpcompletion.completion import (
    __read_sink as read_sink,
)
from rpcompletion.completion import (
    rp_completion,
)


class Test_rpCompletion(TestCase):

    def setUp(self):
        self.logger = create_logger(__name__, "ERROR")
        self.cache = rrCache(cspace="mnx3.1", interactive=False, logger=self.logger)
        self.data_path = os_path.join(
            os_path.dirname(__file__), "data", "input", "lycopene"
        )
        self.output_path = os_path.join(
            os_path.dirname(__file__), "data", "output", "lycopene"
        )
        self.rp2_pathways = os_path.join(self.data_path, "1-rp2_metnet.csv")
        self.sink = os_path.join(self.data_path, "2-sink.csv")
        self.sink_ids = os_path.join(self.data_path, "refs", "sink_ids.json")
        self.rp2paths_compounds = os_path.join(
            self.data_path, "3-rp2paths_compounds.tsv"
        )
        self.rp2paths_pathways = os_path.join(self.data_path, "4-rp2paths_pathways.csv")
        test_file_pattern = "rp_002_0022"
        self.rpsbml_xml = test_file_pattern + "_sbml.xml"
        self.rpsbml_json = os_path.join(
            self.data_path, "refs", test_file_pattern + ".json"
        )

        self.ref_files = [
            "001_0001",
            "001_0006",
            "001_0011",
            "002_0001",
            "002_0011",
            "002_0021",
            "003_0001",
            "003_0131",
            "003_0261",
        ]

    def test_infer_separator_from_extension(self):
        self.assertEqual(infer_separator_from_extension("file.csv"), ",")
        self.assertEqual(infer_separator_from_extension("file.tsv"), "\t")
        self.assertEqual(infer_separator_from_extension("file.txt"), None)
        self.assertEqual(
            infer_separator_from_extension("file.csv", default=";"),
            ",",
        )
        self.assertEqual(
            infer_separator_from_extension("file.unknown", default=";"),
            ";",
        )

    def test_infer_separator_from_content(self):
        # Create temporary files with different separators
        import tempfile

        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp_csv:
            tmp_csv.write("col1,col2,col3\n1,2,3\n4,5,6")
            tmp_csv_path = tmp_csv.name

        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp_tsv:
            tmp_tsv.write("col1\tcol2\tcol3\n1\t2\t3\n4\t5\t6")
            tmp_tsv_path = tmp_tsv.name

        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp_semicolon:
            tmp_semicolon.write("col1;col2;col3\n1;2;3\n4;5;6")
            tmp_semicolon_path = tmp_semicolon.name

        self.assertEqual(infer_separator_from_content(tmp_csv_path), ",")
        self.assertEqual(infer_separator_from_content(tmp_tsv_path), "\t")
        self.assertEqual(infer_separator_from_content(tmp_semicolon_path), ";")

    def test_infer_separator(self):
        # Test with known extensions
        self.assertEqual(infer_separator("file.csv"), ",")
        self.assertEqual(infer_separator("file.tsv"), "\t")
        # Test with unknown extension but content-based inference
        import tempfile

        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as tmp_csv:
            tmp_csv.write("col1,col2,col3\n1,2,3\n4,5,6")
            tmp_csv_path = tmp_csv.name

        self.assertEqual(infer_separator(tmp_csv_path), ",")

    def test_infer_separator_with_invalid_file(self):
        # Test with a non-existent file
        self.assertEqual(infer_separator("non_existent_file.csv"), ",")

    def test_read_sink(self):
        sink_molecules = read_sink(self.sink, logger=self.logger)
        with open(self.sink_ids, "r", encoding="utf-8") as f:
            expected_sink = set(json_load(f))
        self.assertEqual(set(sink_molecules), expected_sink)

    def test_read_sink_splits_dot_ids(self):
        import tempfile

        with tempfile.NamedTemporaryFile(
            mode="w+",
            suffix=".csv",
            delete=False,
            encoding="utf-8",
        ) as tmp:
            tmp.write("Name\n")
            tmp.write("MNXM1.MNXM2\n")
            tmp.write("MNXM3\n")
            tmp_path = tmp.name

        try:
            sink_molecules = read_sink(tmp_path, logger=self.logger)
            self.assertEqual(set(sink_molecules), {"MNXM1", "MNXM2", "MNXM3"})
        finally:
            os.remove(tmp_path)

    def test_rp_completion(self):
        pathways = rp_completion(
            rp2_metnet=self.rp2_pathways,
            sink=self.sink,
            rp2paths_compounds=self.rp2paths_compounds,
            rp2paths_pathways=self.rp2paths_pathways,
            cache=self.cache,
            upper_flux_bound=999999,
            lower_flux_bound=0,
            maxsubpaths=10,
            logger=self.logger,
        )
        pathways = {pathway.get_id(): pathway for pathway in pathways}
        for pathway_id in self.ref_files:
            ref_file = os_path.join(self.output_path, f"rp_{pathway_id}.xml")
            ref_pathway = rpPathway(ref_file)
            self.assertEqual(pathways[f"rp_{pathway_id}"], ref_pathway)

    def test_rp_completion_wo_cofactors(self):
        data_path = os_path.join(
            os_path.dirname(__file__), "data", "input", "wo_cofactors"
        )
        output_path = os_path.join(
            os_path.dirname(__file__), "data", "output", "wo_cofactors"
        )
        rp2_pathways = os_path.join(data_path, "1-rp2_metnet.csv")
        sink = os_path.join(data_path, "2-sink.txt")
        rp2paths_compounds = os_path.join(data_path, "3-rp2paths_compounds.tsv")
        rp2paths_pathways = os_path.join(data_path, "4-rp2paths_pathways.csv")
        cofile = os_path.join(data_path, "cofactors_mnx.tsv")
        pathways = rp_completion(
            rp2_metnet=rp2_pathways,
            sink=sink,
            rp2paths_compounds=rp2paths_compounds,
            rp2paths_pathways=rp2paths_pathways,
            cache=self.cache,
            upper_flux_bound=999999,
            lower_flux_bound=0,
            maxsubpaths=0,
            cofile=cofile,
            logger=self.logger,
        )
        pathways = {pathway.get_id(): pathway for pathway in pathways}
        ref_files_wo_cofactors = [
            "010_0025",
            "010_0026",
            "010_0030",
            "010_0067",
            "010_0068",
            "010_0072",
            "012_0056",
            "012_0059",
            "013_0012",
            "067_0013",
        ]
        for pathway_id in ref_files_wo_cofactors:
            ref_file = os_path.join(output_path, f"rp_{pathway_id}.xml")
            ref_pathway = rpPathway(ref_file)
            self.assertEqual(pathways[f"rp_{pathway_id}"], ref_pathway)

        # print(pathways[0].get_id())
        # exit()
        # for i in range(self.files)
        # # Useless to sort files since smiles could be equivalent and not equal, then checksum will be different
        # for file in listdir(temp_d):
        #     self.assertEqual(
        #         self.files[file],
        #         Path(os_path.join(temp_d, file)).stat().st_size
        #     )
        # rpsbml = rpSBML(os_path.join(temp_d, self.rpsbml_xml))
        # # print(json_dumps(rpsbml.toDict(), indent=4))
        # # self.assertTrue(False)
        # # exit()
        # with open(self.rpsbml_json, 'r') as f:
        #     self.assertDictEqual(rpsbml.toDict(), json_load(f))
        #     # self.assertEqual(os_stat(os_path.join(temp_d, file)).st_size, size)

    # def test_update_rppaths(self):
    #     path_base_id = 2
    #     path_step = 3
    #     path_variant_idx = 1
    #     rxn = {'rule_id'           : 'RR-02-ae41ec5771136ea5-14-F',
    #            'rule_ori_reac'     : 'MNXR113128',
    #            'rule_score'        : 0.7358363677022237,
    #            'left'              : {'CMPD_0000000001': 1},
    #            'right'             : {'TARGET_0000000001': 1},
    #            'step'              : path_step,
    #            'transformation_id' : 'TRS_0_0_1'}
    #     rp_paths = update_rppaths({}, path_base_id, path_variant_idx, rxn)
    #     self.assertEqual(rp_paths, {
    #                                 path_base_id: {
    #                                     path_step: {
    #                                         path_variant_idx: rxn
    #                                     }
    #                                 }
    #                                 })

    # def test_build_side_rxn(self):
    #     cid           = 'CMPD_0000000001'
    #     index         = 1
    #     deprecatedCID = {}
    #     self.assertDictEqual(
    #         build_side_rxn(
    #             str(index)+'.'+cid,
    #             deprecatedCID
    #         ),
    #         {cid: index}
    #     )

    # def test_build_side_rxn_deprecatedCID_NoMatch(self):
    #     cid           = 'CMPD_0000000001'
    #     index         = 1
    #     deprecatedCID = {'CMPD_000000001': 'FOO'}
    #     self.assertDictEqual(
    #         build_side_rxn(
    #             str(index)+'.'+cid,
    #             deprecatedCID
    #         ),
    #         {cid: index}
    #     )

    # def test_build_side_rxn_deprecatedCID_Match(self):
    #     cid           = 'CMPD_0000000001'
    #     index         = 1
    #     deprecatedCID = {'CMPD_0000000001': 'FOO'}
    #     self.assertDictEqual(
    #         build_side_rxn(
    #             str(index)+'.'+cid,
    #             deprecatedCID
    #         ),
    #         {deprecatedCID[cid]: index}
    #     )

    # def test_rp2paths_to_dict(self):
    #     with open(os_path.join(self.data_path, 'refs', 'rp2paths_pathways.json'), 'r') as read_file:
    #         # object_hook is used to convert str keys into int keys as stored in rpCompletion functions
    #         data = json_load(read_file, object_hook=lambda d: {int(k) if k.lstrip('-').isdigit() else k: v for k, v in d.items()})
    #         self.assertDictEqual(
    #             rp2paths_to_dict(
    #                 self.rp2paths_pathways,
    #                 self.cache.get('rr_reactions'), self.cache.get('deprecatedCID_cid')
    #             ),
    #             data
    #         )
