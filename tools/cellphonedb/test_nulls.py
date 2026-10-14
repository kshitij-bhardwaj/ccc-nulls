import pytest

from ccc_nulls import schema
from tools.cellphonedb import adapter, nulls


@pytest.fixture(scope="module")
def interactions(small_adata):
    pytest.importorskip("liana")
    return adapter.liana_reference(small_adata, "label", n_perms=10)


@pytest.mark.parametrize("which", ["a", "b"])
def test_null_output(small_adata, interactions, which):
    fn = getattr(nulls, f"null_{which}")
    try:
        out = fn(small_adata, interactions.copy(), "label", B=20, seed=0)
    except NotImplementedError as e:
        pytest.skip(str(e))
    schema.validate_null(interactions, out, which)
