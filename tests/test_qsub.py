import os.path
import pytest

import gromacs.qsub

def test_queuing_systems(known=("Sun Gridengine", "PBS", "LoadLeveler", 'Slurm')):
    assert len(gromacs.qsub.queuing_systems) == len(known)
    for qs in gromacs.qsub.queuing_systems:
        assert qs.name in known
        assert qs.flag("-x").endswith("-x")
        assert qs.isMine("script.sge") in (True, False)
        directories = ["run1", "run2", "run3"]
        try:
            assert qs.array(directories)
        except NotImplementedError:
            pass

@pytest.mark.parametrize("scriptfile,name", [
    ("foo.sge", "Sun Gridengine"),
    ("foo.pbs", "PBS"),
    ("foo.ll", "LoadLeveler"),
    ("foo.slu", "Slurm")])
def test_detect_queuing_system(scriptfile, name):
    qs = gromacs.qsub.detect_queuing_system(scriptfile)
    assert qs.name == name


def test_generate_submit_scripts(tmpdir):
    templates = ["local.sh", "darwin.sh"]
    with tmpdir.as_cwd():
        names = gromacs.qsub.generate_submit_scripts(templates)
        for name in names:
            assert os.path.exists(name)
    basenames = [os.path.basename(name) for name in names]
    assert basenames == templates

    # todo: check utf-8...
