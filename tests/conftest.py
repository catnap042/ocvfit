import numpy as np
import pytest

import ocvfit as of


def _fits(cell):
    model = of.CapacityModel(cell.positive, cell.negative, v_full=cell.v_full)
    return {k: of.fit_capacity(model, cell.q[k], cell.voltage[k]) for k in ("BOL", "aged")}


@pytest.fixture(scope="session")
def demo():
    """Headline demo: LFP/graphite (Afshar 2017 LFP + Chen 2020 graphite)."""
    return of.demo_cells("lfp", noise=1e-3, seed=0)


@pytest.fixture(scope="session")
def demo_fits(demo):
    return _fits(demo)


@pytest.fixture(scope="session")
def demo_nmc():
    return of.demo_cells("nmc532", noise=1e-3, seed=0)


@pytest.fixture(scope="session")
def demo_nmc_fits(demo_nmc):
    return _fits(demo_nmc)


@pytest.fixture
def rng():
    return np.random.default_rng(1234)
