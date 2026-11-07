import struct
from pathlib import Path

import numpy as np
import pytest

from bamci.fortran_io import (encode_record, read_direct, read_float_record,
                              read_records, record_count, write_float_record,
                              write_record)


def test_record_layout():
    payload = b"\x01\x02\x03\x04"
    encoded = encode_record(payload)
    assert encoded == struct.pack(">i", 4) + payload + struct.pack(">i", 4)


def test_roundtrip_records(tmp_path):
    path = tmp_path / "ic.unf"
    values = np.arange(12, dtype=float)
    with path.open("wb") as fh:
        write_float_record(fh, values[:6])
        write_float_record(fh, values[6:])
    assert record_count(path) == 2
    np.testing.assert_allclose(read_float_record(path, 1), values[6:])
    sizes = [len(p) for p in read_records(path)]
    assert sizes == [24, 24]


def test_corrupt_marker_detected(tmp_path):
    path = tmp_path / "bad.unf"
    path.write_bytes(struct.pack(">i", 8) + b"0" * 8 + struct.pack(">i", 4))
    with pytest.raises(ValueError):
        list(read_records(path))


def test_direct_read(tmp_path):
    path = tmp_path / "grid.dat"
    arr = np.arange(6, dtype=">f4")
    path.write_bytes(arr.tobytes())
    np.testing.assert_allclose(read_direct(path, 1, 6), np.arange(6))
