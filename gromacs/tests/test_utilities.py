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


def test_cat(tmpdir):
    # assert cat.noise == miaow
    # assert not cat.noise == woof
    f1 = tmpdir.join('file1.txt')
    f1.write('foo\n')
    f2 = tmpdir.join('file2.txt')
    f2.write('bar\n')
    out = tmpdir.join('out.txt')
    gromacs.utilities.cat([str(f1), str(f2)], str(out))

    assert out.read() == 'foo\nbar\n'


def test_cat_fail(tmpdir):
    with pytest.raises(OSError):
        gromacs.utilities.cat(['not_here.txt'], str(tmpdir.join('new.txt')))


def test_cat_None_out(tmpdir):
    gromacs.utilities.cat(['some_file.txt'], None)


def test_cat_None_in(tmpdir):
    gromacs.utilities.cat(None, str(tmpdir.join('out.txt')))

    assert not os.path.exists(str(tmpdir.join('out.txt')))


@pytest.fixture
def unlink_files(tmpdir):
    tmpdir.join('hello.mdp').write('foo\n')
    tmpdir.join('#hello.mdp.1#').write('foo\n')
    tmpdir.join('#hello.mdp.2#').write('foo\n')

    tmpdir.join('out.gro').write('bar\n')
    tmpdir.join('#out.gro.1#').write('bar\n')
    tmpdir.join('#out.gro.2#').write('bar\n')

    currdir = os.getcwd()
    try:
        os.chdir(str(tmpdir))
        yield
    finally:
        os.chdir(currdir)


def test_unlink(unlink_files):
    assert os.path.exists('out.gro')

    gromacs.utilities.unlink_f('out.gro')

    assert not os.path.exists('out.gro')


def test_unlink_nonexistant(unlink_files):
    assert not os.path.exists('out.xtc')
    gromacs.utilities.unlink_f('out.xtc')


def test_unlink_gmx_backups(unlink_files):
    gromacs.utilities.unlink_gmx_backups('hello.mdp')

    assert os.path.exists('hello.mdp')
    assert not os.path.exists('#hello.mdp.1#')
    assert not os.path.exists('#hello.mdp.2#')
    assert os.path.exists('out.gro')
    assert os.path.exists('#out.gro.1#')


def test_unlink_gmx(unlink_files):
    gromacs.utilities.unlink_gmx('hello.mdp')
    assert not os.path.exists('hello.mdp')
    assert not os.path.exists('#hello.mdp.1#')
    assert not os.path.exists('#hello.mdp.2#')
    assert os.path.exists('out.gro')
    assert os.path.exists('#out.gro.1#')


@pytest.mark.parametrize('iterable,expected', [
    ('this', 'this'),
    (['this', 'that'], 'this'),
    ([1, 2, 3], 1),
    (np.arange(4), 0),
])
def test_firstof(iterable, expected):
    assert gromacs.utilities.firstof(iterable) == expected


@pytest.mark.parametrize('val,ref', [
    ('a', 'ALA'), ('A', 'ALA'),
    ('ala', 'A'), ('ALA', 'A'), ('Ala', 'A'),
    ('Q', 'GLN'), ('q', 'GLN'),
    ('GLN', 'Q'), ('gln', 'Q'), ('Gln', 'Q'),
])
def test_conv_aa_code(val, ref):
    assert gromacs.utilities.convert_aa_code(val) == ref


@pytest.mark.parametrize('val', ['ALAA', ''])
def test_conv_aa_code_VE(val):
    with pytest.raises(ValueError):
        gromacs.utilities.convert_aa_code(val)


@pytest.fixture
def fileutil():
    class MyFileUtil(gromacs.utilities.FileUtils):
        default_extension = '.test'

        def __init__(self):
            self._init_filename(filename='simple.test')

    return MyFileUtil()


@pytest.mark.parametrize('filename,ext,ref', [
    (None, None, 'simple'),
    ('other', None, 'other'),
    (None, 'pdf', 'simple.pdf'),
    ('other', 'pdf', 'other.pdf'),
])
def test_FileUtils_filename(fileutil, filename, ext, ref):
    assert fileutil.filename(filename=filename, ext=ext) == ref

def test_FileUtils_filename_VE(fileutil):
    del fileutil._filename

    with pytest.raises(ValueError):
        fileutil.filename()

@pytest.fixture
def fileutil_withfiles(fileutil, tmpdir):
    curr = os.getcwd()

    tmpdir.join('exists.txt').write('hello\n')

    try:
        os.chdir(str(tmpdir))
        yield fileutil
    finally:
        os.chdir(curr)


@pytest.mark.parametrize('fn', ['exists.txt', 'nonexistant.txt'])
def test_check_file_exists_ignore(fileutil_withfiles, fn):
    assert fileutil_withfiles.check_file_exists(fn, resolve='ignore') is False


@pytest.mark.parametrize('fn', ['exists.txt', 'nonexistant.txt'])
def test_check_file_exists_force(fileutil_withfiles, fn):
    assert fileutil_withfiles.check_file_exists(fn, force=True) is False


@pytest.mark.parametrize('fn,ref', [('exists.txt', True),
                                    ('nonexistant.txt', False)])
def test_check_file_exists_indicate(fileutil_withfiles, fn, ref):
    assert fileutil_withfiles.check_file_exists(fn, resolve='indicate') is ref


@pytest.mark.parametrize('resolve', ['warn', 'warning'])
def test_check_file_exists_warn(fileutil_withfiles, resolve):
    with pytest.warns(UserWarning):
        fileutil_withfiles.check_file_exists('exists.txt', resolve=resolve)


@pytest.mark.parametrize('resolve', ['exception', 'raise'])
def test_check_file_exists_raise(fileutil_withfiles, resolve):
    with pytest.raises(IOError):
        fileutil_withfiles.check_file_exists('exists.txt', resolve=resolve)
