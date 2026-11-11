"""I/O helper functions for pylinac."""

from __future__ import annotations

import json
import math
import os
import os.path as osp
import struct
import tempfile
import zipfile
from collections.abc import Callable
from contextlib import contextmanager
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import BinaryIO
from urllib.error import HTTPError, URLError
from urllib.parse import unquote, urlparse
from urllib.request import urlopen, urlretrieve

import numpy as np
import pydicom
from google.api_core.exceptions import (
    DeadlineExceeded,
    InternalServerError,
    ServiceUnavailable,
    TooManyRequests,
)
from google.auth.exceptions import TransportError
from google.cloud import storage
from platformdirs import user_cache_dir
from tqdm import tqdm

from .profile import SingleProfile


def _get_demo_file_generation(name: str) -> str:
    """Get the generation number of a demo file in the GCS bucket."""
    client = storage.Client.create_anonymous_client()
    bucket = client.bucket("pylinac_demo_files")
    blob = bucket.blob(name)
    blob.reload()
    return str(blob.generation)


def _get_demo_dir() -> Path:
    """Return the directory used to cache demo files."""
    demo_dir = Path(user_cache_dir("pylinac", appauthor=False)) / "demo_files"
    demo_dir.mkdir(parents=True, exist_ok=True)
    return demo_dir


def is_dicom(file: str | Path) -> bool:
    """Boolean specifying if the file is a proper DICOM file.

    This function is a pared down version of read_preamble meant for a fast return.
    The file is read for a proper preamble ('DICM'), returning True if so,
    and False otherwise. This is a conservative approach.

    Parameters
    ----------
    file : str
        The path to the file.

    See Also
    --------
    pydicom.filereader.read_preamble
    pydicom.filereader.read_partial
    """
    with open(file, "rb") as fp:
        fp.read(0x80)
        prefix = fp.read(4)
        return prefix == b"DICM"


def is_dicom_image(file: str | Path | BinaryIO) -> bool:
    """Boolean specifying if file is a proper DICOM file with a image

    Parameters
    ----------
    file : str
        The path to the file.

    See Also
    --------
    pydicom.filereader.read_preamble
    pydicom.filereader.read_partial
    """
    result = False
    try:
        img = pydicom.dcmread(file, force=True)
        if "TransferSyntaxUID" not in img.file_meta:
            img.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
        img.pixel_array
        result = True
    except (AttributeError, TypeError, KeyError, struct.error):
        pass
    return result


def retrieve_dicom_file(file: str | Path | BinaryIO) -> pydicom.FileDataset:
    """Read and return the DICOM dataset.

    Parameters
    ----------
    file : str
        The path to the file.
    """
    ds = pydicom.dcmread(file, force=True)
    if "TransferSyntaxUID" not in ds.file_meta:
        ds.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
    return ds


class TemporaryZipDirectory(TemporaryDirectory):
    """Creates a temporary directory that unpacks a ZIP archive. Shockingly useful"""

    def __init__(self, zfile: str | Path | BinaryIO, delete: bool = True):
        """
        Parameters
        ----------
        zfile : str
            String that points to a ZIP archive.
        delete : bool
            Whether to delete the temporary directory when the context manager exits.
        """
        # this is a hack. In Python 3.12, we can pass the delete argument to TemporaryDirectory directly
        self._delete = delete
        if delete:
            super().__init__()
            name = self.name
        else:
            # don't use the super() call because it will delete the directory
            # instead, make a temporary directory ourselves and set the name
            name = tempfile.mkdtemp()
            self.name = name
        zfiles = zipfile.ZipFile(zfile)
        zfiles.extractall(path=name)

    def cleanup(self) -> None:
        if self._delete:
            super().cleanup()
        else:
            pass


def retrieve_filenames(
    directory: str | Path,
    func: Callable | None = None,
    recursive: bool = True,
    **kwargs,
) -> list[str]:
    """Retrieve file names in a directory.

    Parameters
    ----------
    directory : str
        The directory to walk over recursively.
    func : function, None
        The function that validates if the file name should be kept.
        If None, no validation will be performed and all file names will be returned.
    recursive : bool
        Whether to search only the root directory.
    kwargs
        Additional arguments passed to the func parameter.
    """
    filenames = []
    if func is None:

        def func(x):
            return True

    for pdir, _, files in os.walk(directory):
        for file in files:
            filename = osp.join(pdir, file)
            if func(filename, **kwargs):
                filenames.append(filename)
        if not recursive:
            break
    return filenames


def retrieve_demo_file(name: str, force: bool = False) -> Path:
    """Retrieve a demo file, downloading it if missing or outdated.

    Files are cached locally. The Google Cloud Storage generation number
    is used to determine whether a cached file is up to date.

    Parameters
    ----------
    name : str
        The name of the demo file in the GCS bucket.
    force : bool
        Download the file even if it already exists and is up to date.

    Returns
    -------
    Path
        The local path to the demo file.
    """
    demo_dir = _get_demo_dir()
    demo_path = demo_dir / name
    metadata_path = demo_path.with_name(f".{demo_path.name}.pylinac-cache.json")
    demo_path.parent.mkdir(parents=True, exist_ok=True)

    # Connect anonymously to the public GCS bucket.
    client = storage.Client.create_anonymous_client()
    bucket = client.bucket("pylinac_demo_files")
    blob = bucket.blob(name)

    # Retrieve the remote object's metadata.
    try:
        blob.reload()
    except (
        TransportError,
        DeadlineExceeded,
        InternalServerError,
        ServiceUnavailable,
        TooManyRequests,
    ):
        # Use the cached copy if the remote service is temporarily unavailable.
        if demo_path.is_file() and not force:
            return demo_path
        raise

    remote_generation = str(blob.generation)

    # Check whether the local file matches the remote version.
    if demo_path.is_file() and metadata_path.is_file() and not force:
        try:
            metadata = json.loads(metadata_path.read_text())
            if metadata.get("generation") == remote_generation:
                return demo_path
        except (OSError, json.JSONDecodeError):
            # Missing or invalid metadata: redownload the file.
            pass

    # Download to a temporary file before replacing the cached copy.
    with TemporaryDirectory(dir=demo_path.parent) as temp_dir:
        temp_path = Path(temp_dir) / demo_path.name

        blob.download_to_filename(
            str(temp_path),
            if_generation_match=blob.generation,
        )

        # Replace the existing file only after a successful download.
        os.replace(temp_path, demo_path)

    # Save metadata only after the successful download.
    metadata_temp = metadata_path.with_suffix(metadata_path.suffix + ".tmp")
    metadata_temp.write_text(json.dumps({"generation": remote_generation}))
    os.replace(metadata_temp, metadata_path)

    return demo_path


def is_url(url: str) -> bool:
    """Determine whether a given string is a valid URL.

    Parameters
    ----------
    url : str

    Returns
    -------
    bool
    """
    try:
        with urlopen(url) as r:
            return r.status == 200
    except Exception:
        return False


def get_url(
    url: str, destination: str | Path | None = None, progress_bar: bool = True
) -> str:
    """Download a URL to a local file.

    Parameters
    ----------
    url : str
        The URL to download.
    destination : str, None
        The destination of the file. If None is given the file is saved to a temporary directory.
    progress_bar : bool
        Whether to show a command-line progress bar while downloading.

    Returns
    -------
    filename : str
        The location of the downloaded file.

    Notes
    -----
    Progress bar use/example adapted from tqdm documentation: https://github.com/tqdm/tqdm
    """

    def my_hook(t):
        last_b = [0]

        def inner(b=1, bsize=1, tsize=None):
            if tsize is not None:
                t.total = tsize
            if b > 0:
                t.update((b - last_b[0]) * bsize)
            last_b[0] = b

        return inner

    try:
        if progress_bar:
            with tqdm(
                unit="B", unit_scale=True, miniters=1, desc=url.split("/")[-1]
            ) as t:
                filename, _ = urlretrieve(
                    url, filename=destination, reporthook=my_hook(t)
                )
        else:
            filename, _ = urlretrieve(url, filename=destination)
    except (HTTPError, URLError, ValueError) as e:
        raise e
    return filename


# this is easier with pandas, but I don't want that as a dependency at this point
class SNCProfiler:
    """Load a file from a Sun Nuclear Profiler device. This accepts .prs files."""

    def __init__(
        self,
        path: str,
        gain_row: int = 20,
        detector_row: int = 106,
        bias_row: int = 107,
        calibration_row: int = 108,
        data_row: int = -1,
        data_columns: slice = slice(5, 259),
    ):
        """
        Parameters
        ----------
        path : str
            Path to the .prs file.
        detector_row
            The row that contains the detector data.
        bias_row
            The row that contains the bias data.
        calibration_row
            The row that contains the calibration data.
        data_row
            The row that contains the data.
        data_columns
            The range of columns that the data is in. Usually, there are some columns before and after the real data.
        """
        with open(path, encoding="cp437") as f:
            raw_data = f.read().splitlines()
            self.detectors = raw_data[detector_row].split("\t")[data_columns]
            self.bias = np.array(raw_data[bias_row].split("\t")[data_columns]).astype(
                float
            )
            self.calibration = np.array(
                raw_data[calibration_row].split("\t")[data_columns]
            ).astype(float)
            self.data = np.array(raw_data[data_row].split("\t")[data_columns]).astype(
                float
            )
            self.gain = float(raw_data[gain_row].split("\t")[1])
            self.timetic = float(raw_data[data_row].split("\t")[2])
            self.integrated_dose = (
                self.calibration * (self.data - self.bias * self.timetic) / self.gain
            )

    def to_profiles(
        self, n_detectors_row: int = 63, **kwargs
    ) -> tuple[SingleProfile, SingleProfile, SingleProfile, SingleProfile]:
        """Convert the SNC data to SingleProfiles. These can be analyzed directly or passed to other modules like flat/sym.
        For the horizontal/cross-plane profile, the detectors on either side of the
        central detector are missing. We adjust the x-values to reflect this.

        Parameters
        ----------
        n_detectors_row : int
            The number of detectors in a given row. Note that they Y profile includes 2 extra detectors from the other 3.
        """

        def drop_cax_sides(vals: np.ndarray) -> np.ndarray:
            x_prof_x_vals = np.arange(start=1, stop=len(vals) + 3)
            half_idx = math.ceil(len(x_prof_x_vals) / 2) - 1
            return np.delete(x_prof_x_vals, [half_idx - 1, half_idx + 1])

        # x values are shifted by 1 to reflect the physical detector number
        y_prof_vals = self.integrated_dose[n_detectors_row : 2 * n_detectors_row + 2]
        y_prof = SingleProfile(
            y_prof_vals,
            x_values=np.arange(start=1, stop=len(y_prof_vals) + 1),
            **kwargs,
        )
        # for all but the y profile, we are missing detectors to the left and right of center because the center y-detector is too wide
        x_prof_vals = self.integrated_dose[:n_detectors_row]
        x_prof = SingleProfile(
            x_prof_vals,
            x_values=drop_cax_sides(x_prof_vals),
            **kwargs,
        )
        pos_prof_values = self.integrated_dose[
            2 * n_detectors_row + 2 : 3 * n_detectors_row + 2
        ]

        pos_prof = SingleProfile(
            pos_prof_values,
            x_values=drop_cax_sides(pos_prof_values),
            **kwargs,
        )
        neg_prof_values = self.integrated_dose[
            3 * n_detectors_row + 2 : 4 * n_detectors_row + 2
        ]
        neg_prof = SingleProfile(
            neg_prof_values,
            x_values=drop_cax_sides(neg_prof_values),
            **kwargs,
        )
        return x_prof, y_prof, pos_prof, neg_prof


def _download_gcs_blob(url: str, destination: Path) -> None:
    """Download a GCS object using the Google Cloud Storage client."""
    from urllib.parse import unquote, urlparse

    parsed = urlparse(url)

    if parsed.netloc == "storage.googleapis.com":
        parts = parsed.path.lstrip("/").split("/", 1)
        if len(parts) != 2:
            raise ValueError(f"Invalid GCS URL: {url}")
        bucket_name, blob_name = parts
    else:
        raise ValueError(f"Unsupported GCS URL: {url}")

    client = storage.Client.create_anonymous_client()
    blob = client.bucket(bucket_name).blob(unquote(blob_name))
    blob.download_to_filename(str(destination))


@contextmanager
def _download_url(
    url: str,
    progress_bar: bool = True,
):
    """Download a URL to a temporary directory and yield its local path."""
    parsed = urlparse(url)
    filename = Path(parsed.path).name or "download"

    with TemporaryDirectory() as temp_dir:
        destination = Path(temp_dir) / filename

        if parsed.netloc == "storage.googleapis.com":
            parts = parsed.path.lstrip("/").split("/", 1)
            if len(parts) != 2:
                raise ValueError(f"Invalid GCS URL: {url}")

            bucket_name, blob_name = parts
            client = storage.Client.create_anonymous_client()
            blob = client.bucket(bucket_name).blob(unquote(blob_name))
            blob.download_to_filename(str(destination))
        else:
            get_url(
                url,
                destination=destination,
                progress_bar=progress_bar,
            )

        yield destination


class DataUnavailableError(RuntimeError):
    pass


class LoadableMixin:
    """Mixin providing generic loading classmethods for common patterns.

    Provides: from_zip, from_url, from_demo, from_demo_image, from_demo_images, from_image.

    Classes may set the class attributes ``DEMO_FILES`` (a list) to enable ``from_demo*`` helpers.
    If a class expects a list of files when constructed (common for multi-image analyses), set
    ``ZIP_FILE_GLOB`` to a glob pattern (e.g. "*.dcm") so the mixin will collect matching files
    and pass the list to the constructor.
    """

    DEMO_FILES: list[str] | None = None
    ZIP_FILE_GLOB: str | None = None

    @classmethod
    def from_zip(cls, zfile: str | Path | BinaryIO, **kwargs):
        """Instantiate ``cls`` from a ZIP archive by extracting to a temporary directory.

        Behavior:
        - If ``ZIP_FILE_GLOB`` is set on the class, collect sorted matches using that glob
          from the extracted directory and call ``cls(file_list, **kwargs)``.
        - Otherwise, call ``cls(tmpdir, **kwargs)`` (constructor accepts a folder path).
        """
        with TemporaryZipDirectory(zfile) as tmpdir:
            if getattr(cls, "ZIP_FILE_GLOB", None):
                # return a list of matching files (Path objects) to the constructor
                files = sorted(Path(tmpdir).glob(cls.ZIP_FILE_GLOB))
                return cls(files, **kwargs)
            return cls(tmpdir, **kwargs)

    @classmethod
    def from_url(cls, url: str, progress_bar: bool = True, **kwargs):
        """Download ``url`` to a temporary file and instantiate ``cls``."""
        with _download_url(url, progress_bar=progress_bar) as filename:
            if zipfile.is_zipfile(filename):
                return cls.from_zip(filename, **kwargs)
            return cls(filename, **kwargs)

    @classmethod
    def from_demo(
        cls,
        name: str | None = None,
        force: bool = False,
        **kwargs,
    ):
        """Instantiate ``cls`` using the repository demo file registry.

        The demo file name is taken from the provided ``name`` argument or the class
        attribute ``DEMO_FILES`` (first entry). If no name is determinable, the
        method raises DataUnavailableError to allow demo-reliant tests to skip when
        demo assets are unavailable.
        The demo file is downloaded via :func:`retrieve_demo_file` and then
        dispatched to ``from_zip`` or the class constructor based on its
        extension.
        """
        # Priority: explicit name argument > DEMO_FILES
        demo_name = name
        if demo_name is None:
            demo_files = getattr(cls, "DEMO_FILES", None)
            if demo_files is not None and len(demo_files) > 0:
                demo_name = demo_files[0]
        if demo_name is None:
            raise DataUnavailableError(
                f"No demo file available for class {cls.__name__}"
            )
        demo_path = retrieve_demo_file(demo_name, force=force)
        # Merge class-provided demo kwargs if present (e.g., sid, dpi)
        demo_defaults = getattr(cls, "DEMO_KWARGS", {}) or {}
        merged_kwargs = {**demo_defaults, **kwargs}
        if str(demo_path).lower().endswith(".zip"):
            return cls.from_zip(demo_path, **merged_kwargs)
        return cls(demo_path, **merged_kwargs)

    @classmethod
    def from_demo_image(cls, **kwargs):
        return cls.from_demo(**kwargs)

    @classmethod
    def from_demo_images(cls, **kwargs):
        return cls.from_demo(**kwargs)

    @classmethod
    def from_image(cls, image, **kwargs):
        """Instantiate ``cls`` directly from an image-like object or path."""
        return cls(image, **kwargs)
