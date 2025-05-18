__pysys_title__   = r""" Simple 2D channel flow using non-equilibrium extrapolation boundary conditions """ 
#                        ==========================
__pysys_purpose__ = r""" Testing proper boundary condition implementation of non-equilibrium extrapolation bcs """ 
    
__pysys_created__ = "2025-05-04"
#__pysys_skipped_reason__   = "Skipped until Bug-1234 is fixed"

#__pysys_traceability_ids__ = "Bug-1234, UserStory-456" 
#__pysys_groups__           = "myGroup, disableCoverage, performance"
#__pysys_modes__            = lambda helper: helper.inheritedModes + [ {'mode':'MyMode', 'myModeParam':123}, ]
#__pysys_parameterized_test_modes__ = {'MyParameterizedSubtestModeA':{'myModeParam':123}, 'MyParameterizedSubtestModeB':{'myModeParam':456}, }

import pysys.basetest, pysys.mappers
from pysys.constants import *

from apes.apeshelper import ApesHelper
class PySysTest(ApesHelper, pysys.basetest.BaseTest):
    def setup(self):
        self.apes.setupMusubi()

    def execute(self):
        musrun = self.apes.runMusubi(np = 5)

    def validate(self):
        self.apes.checkMusLog()
        trackfile = 'channel_DiffAlongHeight_p00000_t2.451E+00.res'
        self.assertPathExists('tracking/'+trackfile,
                              abortOnError = True)
        self.apes.assertIsClose(trackfile, dir = 'tracking')
