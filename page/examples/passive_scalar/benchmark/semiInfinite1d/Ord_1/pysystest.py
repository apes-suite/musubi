__pysys_title__   = r""" Semi-infinite 1D passive scalar benchmark with 1st order source term """
#                        ==========================
__pysys_purpose__ = r""" Testing advection-diffusion-reaction of a passive scalar in a semi-infinite 1D domain with 1st order source term in Musubi. """
    
__pysys_created__ = "2026-09-04"

import pysys.basetest
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
