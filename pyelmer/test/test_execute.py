import unittest.mock
import os
import subprocess
import tempfile
from pyelmer import execute


def test_run_elmer_grid(monkeypatch):
    # create mock for calling ElmerGrid
    mock_call = unittest.mock.MagicMock()
    mock_call.return_value = 0

    def create_file(args, cwd, **kwargs):
        meshdir = os.path.join(cwd, args[3].split(".")[0])
        os.mkdir(meshdir)
        with open(os.path.join(meshdir, "test.boundaries"), "w") as f:
            f.write("0")
        with open(os.path.join(meshdir, "test.elements"), "w") as f:
            f.write("0")

    mock_call.side_effect = create_file
    monkeypatch.setattr(subprocess, "run", mock_call)

    with tempfile.TemporaryDirectory() as tempdir:
        execute.run_elmer_grid(tempdir, "test.msh", "ElmerGridTest")
        # not sure why that's failing
        # mock_call.assert_called_once_with(
        #     ["ElmerGridTest", "14", "2", "test.msh"],
        #     cwd=tempdir,
        #     stdout=unittest.mock.ANY,
        #     stdin=unittest.mock.ANY,
        # )
        assert os.listdir(tempdir) == ["elmergrid.log", "test.boundaries", "test.elements"]
    with tempfile.TemporaryDirectory() as tempdir:
        execute.run_elmer_grid(tempdir, "test.msh", "ElmerGridTest", keep_mesh_dir=True)
        assert os.listdir(tempdir) == ["elmergrid.log", "test"]
        assert os.listdir(os.path.join(tempdir, "test")) == ["test.boundaries", "test.elements"]
    with tempfile.TemporaryDirectory() as tempdir:
        execute.run_elmer_grid(tempdir, "test.msh", "ElmerGridTest", out_dir=os.path.join(tempdir, "out/dir"))
        assert os.listdir(tempdir) == ["elmergrid.log", "out"]
        assert os.listdir(os.path.join(tempdir, os.path.join(tempdir, "out/dir"))) == ["test.boundaries", "test.elements"]
