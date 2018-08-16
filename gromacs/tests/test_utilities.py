# GromacsWrapper: test_example.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from __future__ import division, absolute_import, print_function

import os.path
import pytest
from six.moves import cPickle as pickle
from six.moves import StringIO
from collections import Iterable

import numpy as np

import gromacs.utilities

@pytest.fixture
def string_buffer():
    return StringIO()

def test_which(name="cat"):
    path = gromacs.utilities.which(name)
    assert os.path.basename(path) == name


@pytest.mark.parametrize('path', ["~/whatever", "$HOME/whatever"])
def test_realpath(path):
    abspath = gromacs.utilities.realpath(path)
    assert abspath.startswith(os.path.sep)
    assert abspath.startswith(os.path.expanduser("~"))

class TestAttributeDict(object):
    def setup(self):
        self.d = gromacs.utilities.AttributeDict(foo="bar", baz="boing")

    def test_attribute_get(self):
        assert self.d.foo == "bar"

    def test_dict_get(self):
        assert self.d['foo'] == "bar"
        assert self.d.get('foo') == "bar"

    def test_attribute_set(self):
        self.d.gargl = "blaster"
        assert self.d.gargl == "blaster"

    def test_dict_set(self):
        self.d['gargl'] = "blaster"
        assert self.d['gargl'] == "blaster"

    def test_pickle(self):
        try:
            dump = pickle.dumps(self.d, pickle.HIGHEST_PROTOCOL)
        except Exception as err:
            raise AssertionError("serializing failed: {}".format(err))
        try:
            d = pickle.loads(dump)
        except Exception as err:
            raise AssertionError("de-serializing failed: {0}".format(str(err)))
        assert set(d.keys()) == set(self.d.keys())
        assert set(d.values()) == set(self.d.values())


@pytest.fixture(params=[(42, int),("42", int), ([42],int) ,
                        (2.7, float), ("2.7", float), ([2.7],float),
                        ("jabberwock", str),(["foo","bar"], str),
                        ("42 42", np.integer),("2.7 2.7",float),("foo bar",str)],
                ids=["int -> int", "str -> int","int list -> int list",
                     "float -> float", "str -> float", "float list -> float list",
                     "str -> str","str list -> str list",
                     "str -> int list","str -> float list","str -> str list"])
def conversions(request):
    value, target_type = request.param
    x = gromacs.utilities.autoconvert(value)
    return x, value, target_type

def test_autoconvert(conversions):
    x, value, target_type = conversions
    if isinstance(x, Iterable):
        x = x[0]
    assert isinstance(x, target_type), \
        "Failed to convert '{0}' to type {1}".format(value, target_type)

@pytest.fixture(params=["", "gz", "bz2"])
def openanyfilename(request):
    return "lifeofbrian.txt." + request.param

class TestOpenAny(object):
    quote = """There shall, in that time, be rumours of things going astray,
    erm, and there shall be a great confusion as to where things really are,
    and nobody will really know where lieth those little things... with the
    sort of raffia work base that has an attachment. At this time, a friend
    shall lose his friend's hammer and the young shall not know where lieth
    the things possessed by their fathers that their fathers put there only
    just the night before, about eight o'clock."""

    def test_file(self, tmpdir, openanyfilename):
        filename = openanyfilename
        name, ext = os.path.splitext(filename)
        if ext in ['.gz', '.bz2']:
            quote = self.quote.encode('utf8')
        else:
            quote = self.quote

        with tmpdir.as_cwd():
            with gromacs.utilities.openany(filename, "w") as datafile:
                datafile.write(quote)
            with gromacs.utilities.openany(filename, "r") as datafile:
                data = datafile.read()
        assert data == quote

    def test_stream_write(self, string_buffer):
        with gromacs.utilities.openany(string_buffer, "w") as datafile:
            datafile.write(self.quote)
            # note: closes stream so not easily useable with StringIO (see NamedStream
            # in MDAnalysis)
            # check inside with block
            string_buffer.seek(0)
            assert string_buffer.read() == self.quote

    def test_stream_read(self, string_buffer):
        string_buffer.write(self.quote)
        string_buffer.seek(0)
        with gromacs.utilities.openany(string_buffer, "r") as datafile:
            data = datafile.read()
            # note: closes stream so not easily useable with StringIO (see NamedStream
            # in MDAnalysis)
            # check inside with block
            assert data == self.quote


@pytest.fixture
def pdb_files(tmpdir):
    for i in [1, 15, 300]:
        tmpdir.join('myfile{}.pdb'.format(i)).write('foo\nbar\n')

    currdir = os.getcwd()
    try:
        os.chdir(str(tmpdir))
        yield
    finally:
        os.chdir(currdir)


@pytest.mark.parametrize('args', [
    ('myfile*.pdb',),
    ('myfile1.pdb', 'myfile15.pdb', 'myfile300.pdb'),
    ('myfile1.pdb', 'myfile*.pdb'),
    ('myfile*.pdb', 'myotherfiles*.pdb'),
])
def test_number_pdbs(pdb_files, args):
    gromacs.utilities.number_pdbs(*args)

    assert os.path.exists('myfile0001.pdb')
    assert os.path.exists('myfile0015.pdb')
    assert os.path.exists('myfile0300.pdb')
