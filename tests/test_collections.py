# -*- coding: utf-8 -*-
# GromacsWrapper
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from __future__ import division, absolute_import, print_function

import six

import pytest

import gromacs


@pytest.fixture(scope="function",
                params=[
                    [0, "foo", None, 2.7e-2, "foo"],
                    (0, "foo", None, 2.7e-2, "foo"),
                    set([0, "foo", None, 2.7e-2]),
                    ['ant', 'boar', 'ape', 'gnu'],
                    [['ant', 'spider'], ['boar', 'ape', 'gnu']]
                ])
def things(request):
    stuff = request.param
    return stuff, gromacs.collections.Collection(stuff)

@pytest.fixture(scope="function",
                params=[
                    ['ant', 'boar', 'ape', 'gnu'],
                    [u'åmeise', u'Beißfliege', u'Ürmelchen']
                ])
def textthings(request):
    stuff = request.param
    return stuff, gromacs.collections.Collection(stuff)


class TestCollection(object):
    def test_list_like(self, things):
        seq, collection = things
        assert isinstance(collection, list)
        assert len(collection) == len(seq)

    def test_tolist(self, things):
        seq, collection = things
        lst = collection.tolist()
        assert lst == list(seq)

    def test_save_load(self, things, tmpdir):
        fname = "stuff"
        seq, collection = things
        with tmpdir.as_cwd():
            collection.save(fname)
            newcollection = gromacs.collections.Collection()
            newcollection.load(fname)
        assert newcollection.tolist() == list(seq)
        assert newcollection == collection

    @pytest.mark.parametrize('method,args',
                             [
                                 ('startswith', (u'å',)),
                                 ('upper', ()),
                                 ('capitalize', ())
                             ])
    def test_method_pass_through(self, textthings, method, args):
        seq, collection = textthings
        results = getattr(collection, method)(*args)
        assert results == [getattr(elem, method)(*args) for elem in seq]

    @pytest.mark.parametrize('attribute', ['__doc__'])
    def test_attribute_pass_through(self, textthings, attribute):
        _, collection = textthings
        results = getattr(collection, attribute)
        assert results

    def test_add(self, textthings):
        seq, collection = textthings
        bigger = collection + collection
        assert isinstance(bigger, gromacs.collections.Collection)

        seq.extend(seq)
        assert bigger.tolist() == list(seq)
