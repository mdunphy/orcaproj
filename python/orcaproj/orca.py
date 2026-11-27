from __future__ import annotations

import os
from ctypes import CDLL, c_char_p, c_int, c_size_t
from ctypes.util import find_library
from pathlib import Path

import numpy as np
from numpy.ctypeslib import ndpointer

_SHLIB_NAMES = ("liborca.so", "liborca.dylib", "orca.dll")


def _find_shlib() -> str:
    env = os.environ.get("ORCA_LIB")
    if env:
        return env

    here = Path(__file__).resolve()

    # Bundled with the installed package
    for name in _SHLIB_NAMES:
        candidate = here.parent / "lib" / name
        if candidate.is_file():
            return str(candidate)

    # Found through the system loader (ldconfig, LD_LIBRARY_PATH, DYLD_*)
    found = find_library("orca")
    if found:
        return found

    # Running straight out of the source tree: lib/ next to src/ and include/
    for parent in here.parents:
        for name in _SHLIB_NAMES:
            candidate = parent / "lib" / name
            if candidate.is_file():
                return str(candidate)

    raise OSError(
        "liborca not found: build it with `make` and set ORCA_LIB to the "
        "shared library, or install the package with the library bundled"
    )


_lib = CDLL(_find_shlib())

_lib.orca_init.argtypes = []
_lib.orca_init.restype = c_int

_lib.orca_set_data_dir.argtypes = [c_char_p]
_lib.orca_set_data_dir.restype = None

_initialized = False


def _bundled_data_dir() -> str | None:
    here = Path(__file__).resolve()

    # Bundled with the installed package
    bundled = here.parent / "data"
    if (bundled / "chebcheb_data.bin").is_file():
        return str(bundled)

    # Running straight out of the source tree: data/ next to src/ and include/
    for parent in here.parents:
        candidate = parent / "data"
        if (candidate / "chebcheb_data.bin").is_file():
            return str(candidate)

    return None


def set_data_dir(path: str):
    global _initialized
    if _initialized:
        raise RuntimeError("Cannot set data dir after initialization")
    _lib.orca_set_data_dir(str(path).encode())


def _ensure_init():
    global _initialized
    if not _initialized:
        if "ORCA05_DATA_PATH" not in os.environ:
            data_dir = _bundled_data_dir()
            if data_dir is not None:
                _lib.orca_set_data_dir(data_dir.encode())
        if _lib.orca_init() != 0:
            raise RuntimeError("orca_init failed")
        _initialized = True


_f64 = ndpointer(np.float64, flags="C_CONTIGUOUS")
_i32 = ndpointer(np.int32, flags="C_CONTIGUOUS")

_lib.lonlat_to_ij_vector.argtypes = [_f64, _f64, c_size_t, _f64, _f64, _i32]
_lib.lonlat_to_ij_vector.restype = c_int

_lib.ij_to_lonlat_vector.argtypes = [_f64, _f64, c_size_t, _f64, _f64, _i32]
_lib.ij_to_lonlat_vector.restype = c_int

_lib.ij_to_theta_vector.argtypes = [_f64, _f64, c_size_t, _f64, _f64, _f64, _i32]
_lib.ij_to_theta_vector.restype = c_int


def _as_f64(x):
    if isinstance(x, np.ndarray) and x.dtype == np.float64 and x.flags.c_contiguous:
        return x
    return np.asarray(x, dtype=np.float64, order="C")


def _vec2(f, a, b):
    a = _as_f64(a)
    b = _as_f64(b)

    if a.shape != b.shape:
        raise ValueError("shape mismatch")

    if a.size == 0:
        return a.copy(), a.copy(), np.empty_like(a, dtype=np.int32)

    shape = a.shape
    a = a.ravel()
    b = b.ravel()
    n = a.size

    o1 = np.empty(n, np.float64)
    o2 = np.empty(n, np.float64)
    o3 = np.empty(n, np.int32)

    _ensure_init()
    rc = f(a, b, c_size_t(n), o1, o2, o3)
    if rc != 0:
        raise RuntimeError(f"{f.__name__} failed: {rc}")

    return o1.reshape(shape), o2.reshape(shape), o3.reshape(shape)


def _vec4(f, a, b):
    a = _as_f64(a)
    b = _as_f64(b)

    if a.shape != b.shape:
        raise ValueError("shape mismatch")

    if a.size == 0:
        return a.copy(), a.copy(), a.copy(), np.empty_like(a, dtype=np.int32)

    shape = a.shape
    a = a.ravel()
    b = b.ravel()
    n = a.size

    o1 = np.empty(n, np.float64)
    o2 = np.empty(n, np.float64)
    o3 = np.empty(n, np.float64)
    o4 = np.empty(n, np.int32)

    _ensure_init()
    rc = f(a, b, c_size_t(n), o1, o2, o3, o4)
    if rc != 0:
        raise RuntimeError(f"{f.__name__} failed: {rc}")

    return (
        o1.reshape(shape),
        o2.reshape(shape),
        o3.reshape(shape),
        o4.reshape(shape),
    )


def lonlat_to_ij(lon, lat):
    return _vec2(_lib.lonlat_to_ij_vector, lon, lat)


def ij_to_lonlat(i, j):
    return _vec2(_lib.ij_to_lonlat_vector, i, j)


def ij_to_theta(i, j):
    return _vec4(_lib.ij_to_theta_vector, i, j)
