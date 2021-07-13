from src import util
import pandas as pd
import numpy as np
import tellurium as te
import unittest

IGNORE_TEST = False
IS_PLOT = False


class TestUtil(unittest.TestCase):

    def setUp(self):	
        pass

    def testGetModel(self):
        if IGNORE_TEST:
            return
        def test(model):
            rr = te.loada(model)
            data = rr.simulate()
            self.assertGreater(len(data), 0)
        #
        test(util.LINEAR_PATHWAY_MODEL)
        test(util.WOLF_MODEL)

    def testGetData(self):
        if IGNORE_TEST:
            return
        def test(df):
            self.assertGreater(len(df), 0)
            self.assertTrue(isinstance(df, pd.DataFrame))
        #
        test(util.LINEAR_PATHWAY_DF)
        test(util.WOLF_DF)

    def testGetModule(self):
        if IGNORE_TEST:
            return
        def test(module, name):
            codeStr = util.getModule(module)
            exec(codeStr, globals())
            self.assertTrue(name in globals().keys())
        #
        test("util", "getModule")
        

if __name__ == '__main__':
    unittest.main()
