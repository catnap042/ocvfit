import numpy as np
import pytest

import ocvfit as of


@pytest.mark.parametrize("name", sorted(of.LITERATURE_OCP))
def test_literature_ocp_monotone_decreasing(name):
    x = np.linspace(0, 1, 2001)
    v = of.LITERATURE_OCP[name](x)
    assert np.all(np.isfinite(v))
    assert np.all(np.diff(v) < 0)


def test_reference_values():
    # spot checks against the published expressions
    assert of.LITERATURE_OCP["lfp_afshar2017"](0.5) == pytest.approx(3.4077 - 0.020269 * 0.5, abs=1e-6)
    assert 3.0 < of.LITERATURE_OCP["nmc811_chen2020"](0.5) < 4.0
    assert of.LITERATURE_OCP["graphite_chen2020"](0.9) < 0.15


def test_unknown_name():
    with pytest.raises(KeyError):
        of.literature_ocp("unobtainium")
