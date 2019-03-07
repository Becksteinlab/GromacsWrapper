# -*- coding: utf-8 -*-
# GromacsWrapper
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from __future__ import division, absolute_import, print_function

import pytest
from numpy.testing import assert_equal, assert_almost_equal

from gromacs.fileformats import convert


@pytest.mark.parametrize(
    's,expected',
    [(100, 100),
     ("Jabberwock", u"Jabberwock"),
     (u"Ångström", u"Ångström"),
    ]
)
def test_to_unicode(s, expected):
    output = convert.to_unicode(s)
    assert output == expected

class TestAutoconverter(object):
    def _convert(self, s, **kwargs):
        ac = convert.Autoconverter(**kwargs)
        assert ac.active is True
        return ac.convert(s)


    @pytest.mark.parametrize(
        "s,expected", [
            ('foo bar 22  boing ---', ('foo', 'bar', 22, 'boing', None)),
            ('1 2 3 4', (1, 2, 3, 4)),
            ('1    2  3 4', (1, 2, 3, 4)),
            ('True x X yes Present', (True, True, True, True, True)),
            ('False no - None none', (False, False, False, False, False))
        ],
    )
    @pytest.mark.parametrize('sep', (True, None))
    def test_convert_default(self, s, expected, sep):
        output = self._convert(s, sep=sep)
        assert_equal(output, expected)

    @pytest.mark.parametrize(
        "s,expected", [
            ('1,2,3,4', (1, 2, 3, 4)),
            ('1 2,3,4', ('1 2', 3, 4)),
        ]
    )
    def test_convert_default_sep(self, s, expected, sep=','):
        output = self._convert(s, sep=sep)
        assert_equal(output, expected)


    @pytest.mark.parametrize(
        "s,expected", [
            ('2.71213 3.14', (2.71213, 3.14)),
            ('1000 -234 987654', (1000, -234, 987654)),
        ]
    )
    @pytest.mark.parametrize('sep', (True, None))
    def test_convert_numbers(self, s, expected, sep):
        output = self._convert(s, sep=sep)
        assert_almost_equal(output, expected)


