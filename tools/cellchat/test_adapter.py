import pytest

from ccc_nulls import schema
from tools.cellchat import adapter


def test_liana_reference_schema(small_adata):
    pytest.importorskip("liana")
    df = adapter.liana_reference(small_adata, "label", n_perms=10)
    assert len(df) > 0
    assert (df.score_name == adapter.SCORE_NAME).all()
    schema.validate(df)
    if adapter.HAS_NATIVE_NULL:
        assert df.p_native.notna().all()


def test_native_pipeline_schema(small_adata, tmp_path):
    try:
        inputs = adapter.prepare_inputs(small_adata, "label", tmp_path)
        native = adapter.run_native(inputs, tmp_path)
        df = adapter.native_to_schema(native)
    except NotImplementedError as e:
        pytest.skip(str(e))
    schema.validate(df, require_p_native=adapter.HAS_NATIVE_NULL)
