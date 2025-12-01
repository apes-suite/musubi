__pysys_title__   = r""" Semi-infinite 1D passive scalar benchmark """
#                        ==========================
__pysys_purpose__ = r""" Testing advection-diffusion-reaction of a passive scalar in a semi-infinite 1D domain in Musubi. """
    
__pysys_created__ = "2025-11-27"
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
        trackfile = 'simulation_spc1_p00000_t500.565E-03.res'
        self.assertPathExists('tracking/'+trackfile,
                              abortOnError = True)
        self.apes.assertIsClose(trackfile, dir = 'tracking')
