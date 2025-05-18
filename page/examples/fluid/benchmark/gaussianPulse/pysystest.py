__pysys_title__   = r""" Gaussian pulse in Pressure """ 
#                        ==========================
__pysys_purpose__ = r""" Basic functionality test """ 
    
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
        self.apes.setupMusubi(sdrfile=None)

    def execute(self):
        musrun = self.apes.runMusubi(np = 2)

    def validate(self):
        self.apes.checkMusLog()
        self.assertPathExists('tracking/gaussianPulse_pressAlongLength_p00000_t10.001E+00.res',
                              abortOnError = True)
        self.apes.assertIsClose('gaussianPulse_pressAlongLength_p00000_t10.001E+00.res',
                                dir = 'tracking')
