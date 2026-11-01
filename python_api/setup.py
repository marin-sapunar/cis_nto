"""setup.py - pybind11 is optional."""
try:
    import pybind11
    from pybind11.setup_helpers import Pybind11Extension
    HAS_PYBIND = True
except ImportError:
    HAS_PYBIND = False
    print("INFO: pybind11 not available. C++ extension will be skipped.")

from setuptools import setup, find_packages

ext_modules = []
if HAS_PYBIND:
    ext_modules = [
        Pybind11Extension(
            "cis_nto._core",
            sources=["cis_nto/_pybind.cpp", "../src/easywfo/easywfo_bridge.cpp"],
            include_dirs=[pybind11.get_include(), "../src/easywfo", "../src/easywfo/eigen"],
            cxx_std=17,
        ),
    ]

setup(
    packages=find_packages(where="."),
    package_dir={"cis_nto": "cis_nto"},
    ext_modules=ext_modules,
)
