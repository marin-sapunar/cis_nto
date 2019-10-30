#!/usr/bin/python3
''' Program tests. '''
import os
import time
import unittest
import subprocess
import numpy as np
import pandas as pd


PROG = ['../../build/bin/cis_overlap.exe']
LOGDIR = 'logs/'
REFDIR = 'refs/'


class OverlapTest(unittest.TestCase):
    ''' Test overlap program.'''
    PROG_FLAG = []


    def run_calculation(self, name, call):
        ''' Run single calculation and read results + reference results.'''
        with open(LOGDIR+name+'.log', 'w') as log, open(LOGDIR+name+'.err', 'w') as err:
            subprocess.call(PROG + call, stdout=log, stderr=err)
        self.assertTrue(os.path.isfile(LOGDIR+name), 'Error termination.')
        if not os.path.isfile(REFDIR+name):
            raise FileNotFoundError(REFDIR + name + ' reference file not found.')
        new = pd.read_csv(LOGDIR+name, delim_whitespace=True, header=None, skiprows=1).values
        ref = pd.read_csv(REFDIR+name, delim_whitespace=True, header=None, skiprows=1).values
        return ref, new


    def compare_sao(self, path1, path2, name, opts=None):
        ''' Run single calculation and compare final matrix with reference.'''
        call = self.PROG_FLAG + ['-ao', '-sao', LOGDIR + name, path1, path2]
        if opts:
            call = opts + call
        ref, new = self.run_calculation(name, call)
        self.assertTrue(np.allclose(new, ref), 'Difference in AO overlap.')


    def compare_smo(self, path1, path2, name, opts=None):
        ''' Run single calculation and compare final matrix with reference.'''
        call = self.PROG_FLAG + ['-mo', '-smo', LOGDIR + name, path1, path2]
        if opts:
            call = opts + call
        ref, new = self.run_calculation(name, call)
        self.assertTrue(np.allclose(new, ref), 'Difference in MO overlap.')


    def compare_swf(self, path1, path2, name, opts=None):
        ''' Run single calculation and compare final matrix with reference.'''
        call = self.PROG_FLAG + ['-swf', LOGDIR + name, path1, path2]
        if opts:
            call = opts + call
        ref, new = self.run_calculation(name, call)
        self.assertTrue(np.allclose(new, ref), 'Difference in WF overlap.')


class Turbomole(OverlapTest):
    ''' Turbomole interface tests. '''
    PROG_FLAG = ['-in1', 'turbomole', '-in2', 'turbomole']


    def test_ao_overlap(self):
        ''' Test AO overlap calculation with turbomole.'''
        self.compare_sao('tm_h2o_adc2', 'tm_h2o_adc2', 'tm_ao_overlap')


    def test_mo_overlap(self):
        ''' Test MO overlap calculation with turbomole.'''
        self.compare_smo('tm_h2o_adc2', 'tm_h2o_adc2', 'tm_mo_overlap')


    def test_wf_overlap(self):
        ''' Test WF overlap calculation with turbomole.'''
        self.compare_swf('tm_h2o_adc2', 'tm_h2o_adc2', 'tm_wf_overlap')


    def test_frozen_core_bra(self):
        ''' Test handling of frozen orbitals in ket wave functions.'''
        self.compare_swf('tm_h2o_freeze', 'tm_h2o_adc2', 'tm_frozen_core_bra')


    def test_frozen_core_ket(self):
        ''' Test handling of frozen orbitals in bra wave functions.'''
        self.compare_swf('tm_h2o_adc2', 'tm_h2o_freeze', 'tm_frozen_core_ket')


    def test_change_basis(self):
        ''' Test handling of wave functions with different basis sets.'''
        self.compare_swf('tm_h2o_adc2', 'tm_h2o_basis', 'tm_change_basis')


    def test_change_method(self):
        ''' Test handling of CIS wave functions from different methods.'''
        self.compare_swf('tm_h2o_adc2', 'tm_h2o_tddft', 'tm_change_method')


    def test_change_geom(self):
        ''' Test handling of wave functions at different geometries.'''
        self.compare_swf('tm_h2o_adc2', 'tm_h2o_geom', 'tm_change_geom')


    def test_orthogonalize_states(self):
        ''' Test pre-orthogonalization of wave functions.'''
        self.compare_swf('tm_h2o_adc2', 'tm_h2o_adc2', 'orthogonalize_states', ['-os'])


    def test_orthogonalize_overlap(self):
        ''' Test orthogonalization of overlap matrix.'''
        self.compare_swf('tm_h2o_adc2', 'tm_h2o_adc2', 'orthogonalize_overlap', ['-oo'])


class Molden(OverlapTest):
    ''' Molden interface tests. '''
    PROG_FLAG = ['-in1', 'molden', '-in2', 'molden']


    def test_ao_overlap(self):
        ''' Test AO overlap from molcas generated molden file.'''
        self.compare_sao('molden/h2o_molcas', 'molden/h2o_molcas', 'molcas_ao_overlap')


    def test_mo_overlap(self):
        ''' Test MO overlap from molcas generated molden file.'''
        self.compare_smo('molden/h2o_molcas', 'molden/h2o_molcas', 'molcas_mo_overlap')


class Molpro(OverlapTest):
    ''' Molpro interface tests. '''
    PROG_FLAG = ['-in1', 'molpro_output', '-in2', 'molpro_output']


    def test_ao_overlap(self):
        ''' Test AO overlap from molpro output file.'''
        self.compare_sao('molpro/o2.out', 'molpro/o2.out', 'molpro_ao_overlap')


    def test_mo_overlap(self):
        ''' Test MO overlap from molpro output file.'''
        self.compare_smo('molpro/o2.out', 'molpro/o2.out', 'molpro_mo_overlap')


def check_paths():
    ''' Check that REFDIR and PROG exist. '''
    if not os.path.isdir(REFDIR):
        raise RuntimeError(REFDIR + ' directory not found.')
    if not os.path.isfile(PROG[0]):
        raise RuntimeError(PROG[0]+ ' not found.')
    try:
        os.mkdir(LOGDIR)
    except FileExistsError: # pylint: disable=undefined-variable
        old_log_time = time.gmtime(os.path.getmtime(LOGDIR))
        old_log_time = time.strftime("%y.%m.%d._%H:%M:%S", old_log_time)
        os.rename(LOGDIR, LOGDIR[:-1] + old_log_time)
        os.mkdir(LOGDIR)


if __name__ == '__main__':
    check_paths()
    unittest.main(verbosity=2)
