#!/usr/bin/python3
import unittest
import subprocess
import pandas as pd
import numpy as np
import re
import os


PROG = ['../../build/bin/cis_overlap.exe']


class turbomole(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            os.mkdir('logs')
        except:
            pass


    def compareOverlaps(self, dir1, dir2, name, opts = []):
        call = PROG + opts
        call = call + ['-swf', 'logs/' + name]
        call = call + [dir1, dir2]
        with open('logs/'+name+'.log', 'w') as log, open('logs/'+name+'.err', 'w') as err:
            subprocess.call(call, stdout=log, stderr=err)
        self.assertTrue(os.path.isfile('logs/'+name), 'Error termination.')
        new = pd.read_csv('logs/'+name, delim_whitespace=True, header=None, skiprows=1).values
        ref = pd.read_csv('refs/'+name, delim_whitespace=True, header=None, skiprows=1).values
        self.assertTrue(np.allclose(new, ref), 'Incorrect overlap.')


    def test_activate_frozen_orbitals_ket(self):
        self.compareOverlaps('tm_h2o_freeze', 'tm_h2o_adc2', 'activate_frozen_orbitals_bra')


    def test_activate_frozen_orbitals_bra(self):
        self.compareOverlaps('tm_h2o_adc2', 'tm_h2o_freeze', 'activate_frozen_orbitals_ket')


    def test_change_basis(self):
        self.compareOverlaps('tm_h2o_adc2', 'tm_h2o_basis', 'change_basis')


    def test_change_method(self):
        self.compareOverlaps('tm_h2o_adc2', 'tm_h2o_tddft', 'change_method')


    def test_change_geom(self):
        self.compareOverlaps('tm_h2o_adc2', 'tm_h2o_geom', 'change_geom')


    def test_orthogonalize_states(self):
        self.compareOverlaps('tm_h2o_adc2', 'tm_h2o_adc2', 'orthogonalize_states', ['-os'])


    def test_orthogonalize_overlap(self):
        self.compareOverlaps('tm_h2o_adc2', 'tm_h2o_adc2', 'orthogonalize_overlap', ['-os'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
