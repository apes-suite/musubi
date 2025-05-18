__pysys_title__   = r""" Poisueille flow through a 2D channel with drt kernel""" 
#                        ==========================
__pysys_purpose__ = r""" Testing properly working wall setups for the drt kernel""" 
    
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
        musrun = self.apes.runMusubi(np = 2)

    def validate(self):
        self.apes.checkMusLog()
        trackfile = 'channel_pressAlongLength_p00000_t13.335E+00.res'
        self.assertPathExists('tracking/'+trackfile,
                              abortOnError = True)
        self.apes.assertIsClose(trackfile, dir = 'tracking')
