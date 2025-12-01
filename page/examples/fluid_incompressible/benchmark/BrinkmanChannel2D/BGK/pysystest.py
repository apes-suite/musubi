__pysys_title__   = r""" Brinkman Channel 2D Example with BGK operator """
#                        ==========================
__pysys_purpose__ = r""" Testing Navier--Stokes--Brinkman solver in a 2D channel with BGK collision operator in Musubi. """
    
__pysys_created__ = "2025-11-26"
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
        self.copy(self.input + '/params.lua', self.output)
        self.copy(self.input + '/func.lua', self.output)
        self.apes.setupMusubi()

    def execute(self):
        musrun = self.apes.runMusubi(np = 1)

    def validate(self):
        self.apes.checkMusLog()
        trackfile = 'brinkman_F0_10_p00000_t148.810E+00.res'
        self.assertPathExists('tracking/'+trackfile,
                              abortOnError = True)
        self.apes.assertIsClose(trackfile, dir = 'tracking')
