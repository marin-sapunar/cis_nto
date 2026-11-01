"""Test main API imports."""
def test_api_imports():
    import cis_nto
    assert hasattr(cis_nto, "compute_overlap")
    assert hasattr(cis_nto, "write_libra_format")
    assert hasattr(cis_nto, "read_libra_format")
    print(f"cis_nto version: {cis_nto.__version__}")
