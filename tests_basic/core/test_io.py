"""Test suite for the pylinac.io module."""

import json
import os
import os.path as osp
import unittest
import zipfile
from io import BytesIO
from unittest.mock import MagicMock

import pytest
from google.api_core.exceptions import NotFound, ServiceUnavailable
from google.cloud import storage

from pylinac import Interpolation
from pylinac.core import io
from pylinac.core.io import (
    DataUnavailableError,
    LoadableMixin,
    Path,
    SNCProfiler,
    TemporaryZipDirectory,
    URLError,
    get_url,
    is_dicom,
)
from pylinac.core.profile import SingleProfile
from tests_basic.utils import get_file_from_cloud_test_repo


class TestIO(unittest.TestCase):
    def test_temp_zip_dir(self):
        """Test the TemporaryZipDirectory."""
        zfile = get_file_from_cloud_test_repo(["VMAT", "DRMLC.zip"])

        # test context manager use; shouldn't raise
        with TemporaryZipDirectory(zfile) as tmpzip:
            files = [osp.join(tmpzip, file) for file in os.listdir(tmpzip)]
            # test that they are real files
            self.assertTrue(osp.isfile(files[0]))
            # test that both images were unpacked
            self.assertEqual(len(files), 2, msg="There were not 2 files found")

    def test_get_url(self):
        """Test the URL retreiver."""
        # test webpage
        webpage_url = "http://google.com"
        get_url(webpage_url)  # shouldn't raise
        # test file
        file_url = "https://storage.googleapis.com/pylinac_demo_files/winston_lutz.zip"
        local_file = get_url(file_url)
        osp.isfile(local_file)
        # bad URL
        with self.assertRaises(URLError):
            get_url("http://asdfasdfasdfasdfasdfasdfasdfasdf.org")

    def test_is_dicom(self):
        """Test the is_dicom function."""

        test_file = get_file_from_cloud_test_repo(["VMAT", "DRGSdmlc-105-example.dcm"])
        invalid_file = test_file.replace("DR", "DR_")
        notdicom_file = osp.abspath(__file__)

        # valid file returns True
        self.assertTrue(is_dicom(test_file))

        # return false for real file but not dicom
        self.assertFalse(is_dicom(notdicom_file))

        # test invalid path
        self.assertRaises(IOError, is_dicom, invalid_file)


class TestTempZipDir(unittest.TestCase):
    def test_dir_is_deleted_normally(self):
        """Test that the directory is deleted normally."""
        zfile = get_file_from_cloud_test_repo(["VMAT", "DRMLC.zip"])
        with TemporaryZipDirectory(zfile) as tmpzip:
            self.assertTrue(osp.isdir(tmpzip))
        self.assertFalse(osp.exists(tmpzip))

    def test_dir_remains_when_delete_false(self):
        """Test that the directory is not deleted when delete=False."""
        zfile = get_file_from_cloud_test_repo(["VMAT", "DRMLC.zip"])
        with TemporaryZipDirectory(zfile, delete=False) as tmpzip:
            self.assertTrue(osp.isdir(tmpzip))
        self.assertTrue(osp.exists(tmpzip))


class TestSNCProfiler(unittest.TestCase):
    def test_loading(self):
        path = get_file_from_cloud_test_repo(["9E-GA0.prs"])
        prof = SNCProfiler(path)
        self.assertEqual(len(prof.detectors), 254)

    def test_to_profiles(self):
        path = get_file_from_cloud_test_repo(["9E-GA0.prs"])
        prof = SNCProfiler(path)
        profs = prof.to_profiles()
        self.assertEqual(len(profs), 4)
        self.assertIsInstance(profs[0], SingleProfile)

    def test_detectors(self):
        path = get_file_from_cloud_test_repo(["6XFFF.prs"])
        prof = SNCProfiler(path)
        profs = prof.to_profiles(interpolation=Interpolation.NONE)
        self.assertEqual(len(profs), 4)
        crossplane_prof = profs[0]
        self.assertTrue(crossplane_prof[30] < crossplane_prof[31])
        self.assertTrue(crossplane_prof[32] < crossplane_prof[31])
        self.assertEqual(len(crossplane_prof), 63)
        self.assertEqual(len(profs[1]), 65)


class TestLoadableMixin(unittest.TestCase):
    """Unit tests for the LoadableMixin behaviors."""

    def _make_zip(self, files: dict[str, bytes]) -> str:
        """Create a temporary zip file containing given files. Returns path."""
        import tempfile
        import zipfile

        fd, zpath = tempfile.mkstemp(suffix=".zip")
        os.close(fd)
        with zipfile.ZipFile(zpath, "w") as zf:
            for name, data in files.items():
                zf.writestr(name, data)
        return zpath

    def test_from_zip_uses_tmpdir_when_no_glob(self):
        class Dummy(LoadableMixin):
            def __init__(self, path, **kwargs):
                # record the received argument
                self.init_arg = path

        zpath = self._make_zip({"a.txt": b"1"})
        inst = Dummy.from_zip(zpath)
        # Should have received a path-like string to the extracted folder
        self.assertTrue(isinstance(inst.init_arg, str))
        self.assertTrue(len(inst.init_arg) > 0)

    def test_from_zip_collects_files_when_glob_set(self):
        class DummyList(LoadableMixin):
            ZIP_FILE_GLOB = "*.txt"

            def __init__(self, files, **kwargs):
                # files should be an iterable of Path objects
                self.files = list(files)

        zpath = self._make_zip({"b.txt": b"b", "a.txt": b"a"})
        inst = DummyList.from_zip(zpath)
        # Should be list of Path objects sorted alphabetically
        self.assertEqual(len(inst.files), 2)
        self.assertEqual([p.name for p in inst.files], ["a.txt", "b.txt"])

    def test_from_url_dispatches_to_from_zip_or_constructor(self):
        # Prepare files
        import tempfile

        fd, nonzip = tempfile.mkstemp(suffix=".bin")
        os.close(fd)
        with open(nonzip, "wb") as f:
            f.write(b"data")
        zpath = self._make_zip({"x.txt": b"x"})

        # Monkeypatch get_url in the module
        import pylinac.core.io as io_mod

        original_get_url = io_mod.get_url

        try:

            def mock_get_url(url, destination=None, progress_bar=True):
                import shutil

                shutil.copy(zpath, destination)
                return destination

            io_mod.get_url = mock_get_url

            class FromZipDummy(LoadableMixin):
                ZIP_FILE_GLOB = "*.txt"

                def __init__(self, files, **kwargs):
                    self.files = list(files)

            inst = FromZipDummy.from_url("http://example.com/fake.zip")
            self.assertEqual(len(inst.files), 1)

            # Now non-zip
            def mock_get_url(url, destination=None, progress_bar=True):
                import shutil

                shutil.copy(nonzip, destination)
                return destination

            io_mod.get_url = mock_get_url

            class FromFileDummy(LoadableMixin):
                def __init__(self, path, **kwargs):
                    self.path = path

            inst2 = FromFileDummy.from_url("http://example.com/file.bin")
            # self.assertTrue(inst2.path.endswith(".bin"))
            self.assertEqual(inst2.path.suffix, ".bin")
        finally:
            io_mod.get_url = original_get_url
            try:
                os.remove(nonzip)
            except Exception:
                pass

    def test_from_demo_uses_retrieve_demo_file_and_errors(self):
        """Cover demo error cases and delegate behavior for demo/image helpers."""

        import pylinac.core.io as io_mod

        # Case: no DEMO_FILES configured
        class NoDemo(LoadableMixin):
            DEMO_FILES = None

        with self.assertRaises(DataUnavailableError):
            NoDemo.from_demo()

        # Case: empty DEMO_FILES list
        class EmptyDemo(LoadableMixin):
            DEMO_FILES = []

        with self.assertRaises(DataUnavailableError):
            EmptyDemo.from_demo()

        # Case: empty but truthy DEMO_FILES (force the len==0 branch)
        class TruthyEmpty(list):
            def __bool__(self):
                return True

        class WeirdDemo(LoadableMixin):
            DEMO_FILES = TruthyEmpty()

        with self.assertRaises(DataUnavailableError):
            # this should hit the len(cls.DEMO_FILES) == 0 branch
            WeirdDemo.from_demo()

        # Case: delegate to from_zip or constructor via retrieve_demo_file
        import tempfile

        fd, nonzip = tempfile.mkstemp(suffix=".dat")
        os.close(fd)
        with open(nonzip, "wb") as f:
            f.write(b"x")
        zpath = self._make_zip({"y.txt": b"y"})

        original_retrieve = io_mod.retrieve_demo_file
        try:
            # demo returns zip
            io_mod.retrieve_demo_file = lambda name, force=False: Path(zpath)

            class DemoZip(LoadableMixin):
                DEMO_FILES = ["demo.zip"]
                ZIP_FILE_GLOB = "*.txt"

                def __init__(self, files, **kwargs):
                    self.files = list(files)

            inst = DemoZip.from_demo()
            self.assertEqual(len(inst.files), 1)

            # demo returns non-zip
            io_mod.retrieve_demo_file = lambda name, force=False: Path(nonzip)

            class DemoFile(LoadableMixin):
                DEMO_FILES = ["demo.dat"]

                def __init__(self, path, **kwargs):
                    self.path = path

            inst2 = DemoFile.from_demo()
            self.assertTrue(str(inst2.path).endswith(".dat"))

            # Test from_demo_image/from_demo_images delegate to from_demo
            called = {}

            # demo returns non-zip
            io_mod.retrieve_demo_file = lambda name, force=False: Path(nonzip)

            class DemoImage(LoadableMixin):
                DEMO_FILES = ["demo.dat"]

                def __init__(self, path, **kwargs):
                    called["init"] = path

            # call the convenience wrappers
            DemoImage.from_demo_image()
            self.assertIn("init", called)
            called.clear()
            DemoImage.from_demo_images()
            self.assertIn("init", called)

            # Test from_image directly
            class FromImageDummy(LoadableMixin):
                def __init__(self, image, **kwargs):
                    self.image = image

            inst_img = FromImageDummy.from_image("rawimage")
            self.assertEqual(inst_img.image, "rawimage")

        finally:
            io_mod.retrieve_demo_file = original_retrieve
            try:
                os.remove(nonzip)
            except Exception:
                pass


def _mock_gcs(monkeypatch, generation="123"):
    """Set up a mocked anonymous GCS client and blob."""
    blob = MagicMock()
    blob.generation = generation

    def reload():
        blob.generation = generation

    blob.reload.side_effect = reload

    def download_to_filename(filename, **kwargs):
        Path(filename).write_bytes(b"downloaded demo data")

    blob.download_to_filename.side_effect = download_to_filename

    bucket = MagicMock()
    bucket.blob.return_value = blob

    client = MagicMock()
    client.bucket.return_value = bucket

    monkeypatch.setattr(
        storage.Client,
        "create_anonymous_client",
        MagicMock(return_value=client),
    )

    return client, bucket, blob


def test_retrieve_demo_file_downloads_when_missing(tmp_path, monkeypatch):
    monkeypatch.setattr(io, "_get_demo_dir", lambda: tmp_path)
    _, bucket, blob = _mock_gcs(monkeypatch)

    result = io.retrieve_demo_file("example.zip")

    assert result == tmp_path / "example.zip"
    assert result.read_bytes() == b"downloaded demo data"
    assert json.loads((tmp_path / ".example.zip.pylinac-cache.json").read_text()) == {
        "generation": "123"
    }

    bucket.blob.assert_called_once_with("example.zip")
    blob.reload.assert_called_once()
    blob.download_to_filename.assert_called_once()

    # The download must be pinned to the generation we checked.
    _, kwargs = blob.download_to_filename.call_args
    assert kwargs["if_generation_match"] == "123"


def test_retrieve_demo_file_uses_cache_when_generation_matches(tmp_path, monkeypatch):
    monkeypatch.setattr(io, "_get_demo_dir", lambda: tmp_path)
    _, _, blob = _mock_gcs(monkeypatch)

    demo_path = tmp_path / "example.zip"
    demo_path.write_bytes(b"already cached")
    metadata_path = tmp_path / ".example.zip.pylinac-cache.json"
    metadata_path.write_text(json.dumps({"generation": "123"}))

    result = io.retrieve_demo_file("example.zip")

    assert result == demo_path
    assert result.read_bytes() == b"already cached"
    blob.reload.assert_called_once()
    blob.download_to_filename.assert_not_called()


def test_retrieve_demo_file_downloads_when_generation_changes(tmp_path, monkeypatch):
    monkeypatch.setattr(io, "_get_demo_dir", lambda: tmp_path)
    _, _, blob = _mock_gcs(monkeypatch, generation="456")

    demo_path = tmp_path / "example.zip"
    demo_path.write_bytes(b"old cached data")
    metadata_path = tmp_path / ".example.zip.pylinac-cache.json"
    metadata_path.write_text(json.dumps({"generation": "123"}))

    result = io.retrieve_demo_file("example.zip")

    assert result.read_bytes() == b"downloaded demo data"
    assert json.loads(metadata_path.read_text()) == {"generation": "456"}
    blob.download_to_filename.assert_called_once()


def test_retrieve_demo_file_force_downloads_current_file(tmp_path, monkeypatch):
    monkeypatch.setattr(io, "_get_demo_dir", lambda: tmp_path)
    _, _, blob = _mock_gcs(monkeypatch)

    demo_path = tmp_path / "example.zip"
    demo_path.write_bytes(b"already cached")
    metadata_path = tmp_path / ".example.zip.pylinac-cache.json"
    metadata_path.write_text(json.dumps({"generation": "123"}))

    result = io.retrieve_demo_file("example.zip", force=True)

    assert result.read_bytes() == b"downloaded demo data"
    blob.download_to_filename.assert_called_once()


def test_retrieve_demo_file_uses_cache_when_service_unavailable(tmp_path, monkeypatch):
    monkeypatch.setattr(io, "_get_demo_dir", lambda: tmp_path)
    _, _, blob = _mock_gcs(monkeypatch)

    demo_path = tmp_path / "example.zip"
    demo_path.write_bytes(b"already cached")
    blob.reload.side_effect = ServiceUnavailable("Service unavailable")

    result = io.retrieve_demo_file("example.zip")

    assert result == demo_path
    assert result.read_bytes() == b"already cached"
    blob.download_to_filename.assert_not_called()


def test_retrieve_demo_file_raises_not_found(tmp_path, monkeypatch):
    monkeypatch.setattr(io, "_get_demo_dir", lambda: tmp_path)
    _, _, blob = _mock_gcs(monkeypatch)

    blob.reload.side_effect = NotFound("File not found")

    demo_path = tmp_path / "example.zip"
    demo_path.write_bytes(b"stale cached data")

    with pytest.raises(NotFound):
        io.retrieve_demo_file("example.zip")

    assert demo_path.read_bytes() == b"stale cached data"
    blob.download_to_filename.assert_not_called()


def test_retrieve_demo_file_redownloads_with_corrupt_metadata(tmp_path, monkeypatch):
    monkeypatch.setattr(io, "_get_demo_dir", lambda: tmp_path)
    _, _, blob = _mock_gcs(monkeypatch)

    demo_path = tmp_path / "example.zip"
    demo_path.write_bytes(b"old cached data")
    metadata_path = tmp_path / ".example.zip.pylinac-cache.json"
    metadata_path.write_text("not valid json")

    result = io.retrieve_demo_file("example.zip")

    assert result.read_bytes() == b"downloaded demo data"
    assert json.loads(metadata_path.read_text()) == {"generation": "123"}
    blob.download_to_filename.assert_called_once()


def test_retrieve_demo_file_preserves_cached_file_when_download_fails(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(io, "_get_demo_dir", lambda: tmp_path)
    _, _, blob = _mock_gcs(monkeypatch, generation="456")

    demo_path = tmp_path / "example.zip"
    demo_path.write_bytes(b"old cached data")
    metadata_path = tmp_path / ".example.zip.pylinac-cache.json"
    metadata_path.write_text(json.dumps({"generation": "123"}))

    def fail_download(filename, **kwargs):
        raise ServiceUnavailable("Download failed")

    blob.download_to_filename.side_effect = fail_download

    with pytest.raises(ServiceUnavailable):
        io.retrieve_demo_file("example.zip")

    assert demo_path.read_bytes() == b"old cached data"
    assert json.loads(metadata_path.read_text()) == {"generation": "123"}


def test_from_demo_uses_gcs_cached_retrieval(monkeypatch, tmp_path):
    """Ensure LoadableMixin.from_demo delegates demo retrieval to retrieve_demo_file."""
    from pylinac.core.io import LoadableMixin

    demo_file = tmp_path / "example.zip"
    demo_file.write_bytes(b"demo data")

    calls = []

    def mock_retrieve_demo_file(name, force=False):
        calls.append((name, force))
        return demo_file

    monkeypatch.setattr(io, "retrieve_demo_file", mock_retrieve_demo_file)

    class DemoLoader(LoadableMixin):
        DEMO_FILES = ["example.zip"]

        @classmethod
        def from_zip(cls, filename, **kwargs):
            return filename

    result = DemoLoader.from_demo()

    assert result == demo_file
    assert calls == [("example.zip", False)]


def test_from_demo_forwards_force_to_retrieve_demo_file(monkeypatch, tmp_path):
    from pylinac.core.io import LoadableMixin

    demo_file = tmp_path / "example.zip"
    demo_file.write_bytes(b"demo data")

    calls = []

    def mock_retrieve_demo_file(name, force=False):
        calls.append((name, force))
        return demo_file

    monkeypatch.setattr(io, "retrieve_demo_file", mock_retrieve_demo_file)

    class DemoLoader(LoadableMixin):
        DEMO_FILES = ["example.zip"]

        @classmethod
        def from_zip(cls, filename, **kwargs):
            return filename

    DemoLoader.from_demo(force=True)

    assert calls == [("example.zip", True)]


def test_from_url_uses_anonymous_gcs_client(monkeypatch):
    """GCS URLs should use the anonymous GCS client, not get_url."""
    from pylinac.core.io import LoadableMixin

    zip_buffer = BytesIO()
    with zipfile.ZipFile(zip_buffer, "w") as archive:
        archive.writestr("example.txt", "hello")
    zip_data = zip_buffer.getvalue()

    blob = MagicMock()

    def download_to_filename(filename, **kwargs):
        Path(filename).write_bytes(zip_data)

    blob.download_to_filename.side_effect = download_to_filename

    bucket = MagicMock()
    bucket.blob.return_value = blob

    client = MagicMock()
    client.bucket.return_value = bucket

    monkeypatch.setattr(
        storage.Client,
        "create_anonymous_client",
        MagicMock(return_value=client),
    )

    def unexpected_get_url(*args, **kwargs):
        pytest.fail("get_url should not be used for GCS URLs")

    monkeypatch.setattr(io, "get_url", unexpected_get_url)

    class ZipLoader(LoadableMixin):
        @classmethod
        def from_zip(cls, filename, **kwargs):
            with zipfile.ZipFile(filename) as archive:
                return archive.read("example.txt")

    result = ZipLoader.from_url("https://storage.googleapis.com/my-bucket/example.zip")

    assert result == b"hello"
    client.bucket.assert_called_once_with("my-bucket")
    bucket.blob.assert_called_once_with("example.zip")
    blob.download_to_filename.assert_called_once()


def test_download_url_rejects_invalid_gcs_url():
    """Malformed GCS URLs should fail clearly instead of being treated as HTTP."""
    with pytest.raises(ValueError, match="Invalid GCS URL"):
        with io._download_url("https://storage.googleapis.com/"):
            pass


def test_download_url_cleans_up_temporary_directory(monkeypatch, tmp_path):
    """The temporary download should exist inside the context and be removed after."""
    downloaded_path = None
    temp_directory = None

    def mock_get_url(url, destination=None, progress_bar=True):
        Path(destination).write_bytes(b"downloaded")
        return destination

    monkeypatch.setattr(io, "get_url", mock_get_url)

    with io._download_url("https://example.com/file.bin") as path:
        downloaded_path = path
        temp_directory = path.parent

        assert path.exists()
        assert path.read_bytes() == b"downloaded"

    assert downloaded_path is not None
    assert temp_directory is not None
    assert not downloaded_path.exists()
    assert not temp_directory.exists()


def test_retrieve_demo_file_downloads_when_metadata_is_missing(tmp_path, monkeypatch):
    """An existing file without metadata cannot be verified as current."""
    monkeypatch.setattr(io, "_get_demo_dir", lambda: tmp_path)
    _, _, blob = _mock_gcs(monkeypatch)

    demo_path = tmp_path / "example.zip"
    demo_path.write_bytes(b"unverified cached data")

    metadata_path = tmp_path / ".example.zip.pylinac-cache.json"
    assert not metadata_path.exists()

    result = io.retrieve_demo_file("example.zip")

    assert result.read_bytes() == b"downloaded demo data"
    assert json.loads(metadata_path.read_text()) == {"generation": "123"}
    blob.download_to_filename.assert_called_once()


def test_get_demo_dir_uses_platform_cache(tmp_path, monkeypatch):
    def mock_user_cache_dir(appname, appauthor=False):
        assert appname == "pylinac"
        assert appauthor is False
        return str(tmp_path)

    monkeypatch.setattr(io, "user_cache_dir", mock_user_cache_dir)

    demo_dir = io._get_demo_dir()

    assert demo_dir == tmp_path / "demo_files"
    assert demo_dir.is_dir()
