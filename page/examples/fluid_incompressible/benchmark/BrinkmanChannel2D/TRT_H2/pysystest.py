__pysys_title__   = r""" Brinkman Channel 2D Example with TRT operator and H2 Brinkman force """
#                        ==========================
__pysys_purpose__ = r""" Testing Navier--Stokes--Brinkman solver in a 2D channel with TRT collision operator and H2 Brinkman force in Musubi. """
    
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
        trackfile = 'brinkman_F0_100_p00000_t148.810E+00.res'
        self.assertPathExists('tracking/'+trackfile,
                              abortOnError = True)
        self.apes.assertIsClose(trackfile, dir = 'tracking')
