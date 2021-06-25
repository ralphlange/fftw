import time
import unittest

import numpy as np
from epics import PV


class TestFFTW(unittest.TestCase):

    def test_4e_pulse(self):
        """
        Test a 4 element array with a pulse
        """
        real_is_in = False
        imag_is_in = False

        def data_callback(pvname=None, **kwargs):
            nonlocal real_is_in, imag_is_in
            if pvname.endswith('real'):
                real_is_in = True
            elif pvname.endswith('imag'):
                imag_is_in = True

        data = [0, 1, 0, 0]
        wintype = PV('A1:wintype')
        inp = PV('A1:inp-real')
        outr = PV('A1:out-real', callback=data_callback)
        outi = PV('A1:out-imag', callback=data_callback)

        while not all([real_is_in, imag_is_in]):
            time.sleep(0.001)
        real_is_in = False
        imag_is_in = False

        wintype.put('None', wait=True)
        inp.put(data, wait=True)

        result = np.fft.rfft(data)

        while not all([real_is_in, imag_is_in]):
            time.sleep(0.002)

        self.assertTrue(np.allclose(result.imag, outi.get()))
        self.assertTrue(np.allclose(result.real, outr.get()))

    def test_1ke_1sine(self):
        """
        Test a 1k array with a single sine wave
        """
        real_is_in = False
        imag_is_in = False

        def data_callback(pvname=None, **kwargs):
            nonlocal real_is_in, imag_is_in
            if pvname.endswith('real'):
                real_is_in = True
            elif pvname.endswith('imag'):
                imag_is_in = True

        data = np.cos(2 * np.pi * np.arange(1024) / 1024) + np.cos(4 * 2 * np.pi * np.arange(1024) / 1024)
        wintype = PV('A2:wintype')
        inp = PV('A2:inp-real')
        outr = PV('A2:out-real', callback=data_callback)
        outi = PV('A2:out-imag', callback=data_callback)
        while not all([real_is_in, imag_is_in]):
            time.sleep(0.001)
        real_is_in = False
        imag_is_in = False

        wintype.put('None', wait=True)
        inp.put(data, wait=True)

        result = np.fft.rfft(data)

        while not all([real_is_in, imag_is_in]):
            time.sleep(0.001)

        self.assertTrue(np.allclose(result.imag, outi.get()))
        self.assertTrue(np.allclose(result.real, outr.get()))

    def test_1ke_asub_1sine(self):
        """
        Test a 1k array with a single sine wave - using aSub
        """
        real_is_in = False
        imag_is_in = False

        def data_callback(pvname=None, **kwargs):
            nonlocal real_is_in, imag_is_in
            if pvname.endswith('real'):
                real_is_in = True
            elif pvname.endswith('imag'):
                imag_is_in = True

        data = np.cos(2 * np.pi * np.arange(1024) / 1024) + np.cos(4 * 2 * np.pi * np.arange(1024) / 1024)
        wintype = PV('A3:wintype')
        inp = PV('A3:inp-real')
        outr = PV('A3:out-real', callback=data_callback)
        outi = PV('A3:out-imag', callback=data_callback)
        while not all([real_is_in, imag_is_in]):
            time.sleep(0.001)
        real_is_in = False
        imag_is_in = False

        wintype.put('None', wait=True)
        inp.put(data, wait=True)

        result = np.fft.rfft(data)

        while not all([real_is_in, imag_is_in]):
            time.sleep(0.001)

        self.assertTrue(np.allclose(result.imag, outi.get()))
        self.assertTrue(np.allclose(result.real, outr.get()))

if __name__ == '__main__':
    unittest.main()
