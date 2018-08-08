# GromacsWrapper: test_example.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

from __future__ import division, absolute_import, print_function

import os.path
import pytest

import gromacs

# use 'ls' as command and only use common BSD/GNU options
@pytest.fixture
def command(command="ls",
                 args=('-a', '-1'),
                 options={'l': True, 'v': True}):
    Command_cls = type(command.capitalize(),
                       (gromacs.core.Command,),
                       {'command_name': command})
    return Command_cls(*args, **options)


class TestCommand(object):
    def test_run_default(self, command):
        rc, out, err = command()
        assert rc == 0
        assert not out
        assert not err

    def test_run_with_args(self, command):
        rc, out, err = command('-d', os.path.curdir, b=True, F=True )
        assert rc == 0
        assert not out
        assert not err

    def test_run_capture_stdout(self, command):
        rc, out, err = command(stdout=False)
        assert rc == 0
        assert out
        assert not err

    def test_run_capture_stderr(self, command):
        rc, out, err = command('/this_does_not_exist_Foo_Bar', stderr=False)
        assert rc > 0
        assert not out
        assert err

    @pytest.mark.parametrize('inp', ("not_used", ("not", "used")))
    def test_run_with_input(self, command, inp):
        rc, out, err = command(input=inp)
        assert rc == 0
        assert not out
        assert not err

    def test_help_short(self, command, capsys):
        command.help()
        captured = capsys.readouterr()
        assert command.command_name in captured.out

    def test_help_long(self, command, capsys):
        command.help(long=True)
        captured = capsys.readouterr()
        assert command.command_name in captured.out
        assert gromacs.core.Command.__call__.__doc__ in captured.out
