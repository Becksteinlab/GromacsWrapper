# -*- coding: utf-8 -*-
# numkit.integration test cases
# Part of GromacsWrapper
# Copyright (c) Oliver Beckstein <orbeckst@gmail.com>
# Published under the Modified BSD Licence.

"""
===================================
 Test cases for numkit.integration
====================================

"""

import numkit.integration
import numpy
from numpy.testing import *

class Test_simps_error(TestCase):
    # special case with constant spacing and constant error (=1.0) so
    # that one can easily calculate the correct errors
    def setUp(self):
        self.even_x = numpy.linspace(0,10,10)           # const spacing
        self.even_dx = numpy.diff(self.even_x)[0]
        self.even_dy = 1.0 * numpy.ones(10, dtype=float) # const error
        self.odd_x = numpy.linspace(0,10,11)             # const spacing
        self.odd_dx = numpy.diff(self.odd_x)[0]
        self.odd_dy = 1.0 * numpy.ones(11, dtype=float) # const error        

    def err_analytic_odd(self):
        # error for all dy and dx equal, using only Simpson's rule:
        #  err**2 = (h/3)^2(df1^2 + (4df2)^2 + (2df3)^2 + ... + (4df_N-1)^2 + (dfN)^2)
        #  err**2 = (h/3)^2 dy^2 (10(N-1)-(2^2-1)+1)
        dy = self.odd_dy[0]
        N = self.odd_dy.shape[0]
        h = self.odd_dx
        return numpy.sqrt((h/3)**2 * dy**2 * (10*(N-1)-2))

    def err_analytic_even(self):
        # err**2 = (h/6)^2((3df1)^2 + ((3+2)df2)^2 + (8df3)^2 + (4df4)^2 + (8df5)^2 + ... + (2dfN)^2)
        # special case for constant dy and constant spacing
        dy = self.even_dy[0]
        N = self.even_dy.shape[0]
        h = self.even_dx
        #err_analytic**2 = (h/6.0)**2 * dy**2 * (2*3**2 + 12 + (8**2+4**2)*0.5*(N-2) - (4**2-2**2) + 2**2)
        return numpy.sqrt((h/6.0)**2 * dy**2 * (40.*N - 58.))        

    def test_odd_const_dx(self):
        err_dx = numkit.integration.simps_error(self.odd_dy, dx=self.odd_dx)
        err_analytic = self.err_analytic_odd()
        assert_equal(err_dx, err_analytic, err_msg="Simps error for const spacing dx")

    def test_odd_const_x(self):
        err_x =  numkit.integration.simps_error(self.odd_dy, x=self.odd_x)
        err_analytic = self.err_analytic_odd()
        assert_equal(err_x, err_analytic, err_msg="Simps error for const spacing dx, from abscissas")

    def test_odd_const_x_dx(self):
        err_dx = numkit.integration.simps_error(self.odd_dy, dx=self.odd_dx)
        err_x =  numkit.integration.simps_error(self.odd_dy, x=self.odd_x)
        assert_equal(err_x, err_dx, err_msg="Simps error for const spacing: either scalar dx or abscissas x")

    def test_even_first(self):
        err_analytic = self.err_analytic_even()
        err = numkit.integration.simps_error(self.even_dy, dx=self.even_dx, even="first")
        assert_almost_equal(err, err_analytic, 
                            err_msg="Even # points, first: Simps+Trapezoid stat. error: const dx and const dy")

    def test_even_last(self):
        err_analytic = self.err_analytic_even()
        err = numkit.integration.simps_error(self.even_dy, dx=self.even_dx, even="last")
        assert_almost_equal(err, err_analytic, 
                            err_msg="Even # points, last: Trapezoid+Simps stat. error: const dx and const dy")

    def test_even_avg(self):
        err_analytic = self.err_analytic_even()
        err = numkit.integration.simps_error(self.even_dy, dx=self.even_dx, even="avg")
        assert_almost_equal(err, err_analytic, 
                            err_msg="Even # points, avg: Simps+Trapezoid stat. error: const dx and const dy")

    def Xtest_h1h2(self):
        x = numpy.array([0,1,2,4.,5])
        y = 2*x
        dy = numpy.array([1.,1,1,1,1])

        err = numkit.integration.simps_error(dy,x=x)
        # check correctness!!
        assert_almost_equal(err, 2.7613402542968153, err_msg="Simps error for un-even spacing")



