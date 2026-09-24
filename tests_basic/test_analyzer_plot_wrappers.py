import io
import os

import matplotlib
from matplotlib import pyplot as plt

from pylinac.core.plotable import MatplotlibPlotableMixin

# Ensure headless backend for CI
matplotlib.use("Agg")


class DummyAnalyzer(MatplotlibPlotableMixin):
    """Minimal analyzer-like object implementing plot_analyzed_image using the mixin primitives."""

    def plot_analyzed_image(self, show: bool = True, **plt_kwargs):
        fig, ax = self._subplots()
        ax.plot([0, 1], [0, 1])
        ax.set_title("dummy")
        return self._finalize_plot(fig, show=show)


def test_dummy_analyzer_plot_returns_figure():
    d = DummyAnalyzer()
    fig = d.plot_analyzed_image(show=False)
    assert fig is not None
    # must be a matplotlib Figure
    from matplotlib.figure import Figure

    assert isinstance(fig, Figure)
    plt.close(fig)


def test_dummy_analyzer_save_to_bytesio_and_path(tmp_path):
    d = DummyAnalyzer()

    # save to BytesIO
    bio = io.BytesIO()
    d.save_analyzed_image(bio)
    bio.seek(0)
    data = bio.read()
    assert len(data) > 0

    # save to path
    p = tmp_path / "dummy.png"
    d.save_analyzed_image(p)
    assert p.exists()
    assert p.stat().st_size > 0

    # clean up
    try:
        os.remove(p)
    except Exception:
        pass
