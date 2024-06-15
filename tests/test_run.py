# GromacsWrapper: test_example.py
# Copyright (c) 2009 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License 3 (or higher, your choice)
# See the file COPYING for details.

import time
import pytest
import os
import signal
import threading

from .datafiles import datafile

import gromacs.run


class Test_check_mdrun_success(object):
    @staticmethod
    def test_no_logfile():
        assert gromacs.run.check_mdrun_success("bogus_file.log") is None

    @staticmethod
    def test_success_Gromacs4():
        assert gromacs.run.check_mdrun_success(datafile("gromacs4_success.log")) is True

    @staticmethod
    def test_incomplete_Gromacs4():
        assert (
            gromacs.run.check_mdrun_success(datafile("gromacs4_incomplete.log"))
            is False
        )

    @staticmethod
    def test_success_Gromacs5():
        assert gromacs.run.check_mdrun_success(datafile("gromacs5_success.log")) is True

    @staticmethod
    def test_incomplete_Gromacs5():
        assert (
            gromacs.run.check_mdrun_success(datafile("gromacs5_incomplete.log"))
            is False
        )


# The following tests need an existing Gromacs environment. They should run
# with either Gromacs 4 or Gromacs 5


def test_MDRunner():
    try:
        mdrun = gromacs.run.MDrunner()
    except OSError:
        raise RuntimeError("This test requires a Gromacs environment.")

    rc = mdrun.run(mdrunargs={"version": True})
    assert rc == 0, "mdrun failed to run through MDrunner"


def test_MDRunner_keyboard_interrupt(monkeypatch):
    """Test that a keyboard interrupt is handled correctly."""

    # Create a mock for subprocess.Popen
    class MockPopen:
        def __init__(self, *args, **kwargs):
            pass

        def wait(self):
            time.sleep(2)  # Simulate a long-running process

        def terminate(self):
            pass

    # Use monkeypatch to replace subprocess.Popen with MockPopen
    monkeypatch.setattr("gromacs.run.subprocess.Popen", MockPopen)

    # Initialize MDrunner
    mdrunner = gromacs.run.MDrunner()

    # Function to send SIGINT after a delay
    def send_interrupt():
        time.sleep(1)  # Short delay before sending SIGINT
        os.kill(os.getpid(), signal.SIGINT)

    # Start a thread to send the SIGINT
    interrupt_thread = threading.Thread(target=send_interrupt)
    interrupt_thread.start()

    # Run the MDrunner in the main thread to handle the signal
    try:
        mdrunner.run()
    except KeyboardInterrupt:
        # Handle the KeyboardInterrupt within the test
        pass

    # Ensure the signal handler was called
    assert mdrunner.signal_handled, "The signal handler was not called as expected"

    # Ensure the interrupt thread has finished
    interrupt_thread.join()


class Test_find_gromacs_command(object):
    # Gromacs 4 or Gromacs 5 (in this order)
    commands = ["grompp", "gmx grompp"]

    def test_find(self):
        driver, name = gromacs.run.find_gromacs_command(self.commands)
        assert driver in (
            None,
            "gmx",
        ), "find_gromacs_command() did not identify a driver"
        assert name == self.commands[0], "find_gromacs_command() did not find a command"

    @staticmethod
    def test_raises_ValueError():
        with pytest.raises(OSError):
            driver, name = gromacs.run.find_gromacs_command(["./not_a_command"])


def test_get_double_or_single_prec_mdrun():
    # tests only ship with single prec mdrun
    mdrun = gromacs.run.get_double_or_single_prec_mdrun()
    assert mdrun.command_name in (
        "mdrun",
        "mdrun_d",
    ), "gromacs.run.get_double_or_single_prec_mdrun() could not find any mdrun"
