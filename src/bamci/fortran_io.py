"""Leitura/escrita de arquivos Fortran UNFORMATTED binarios (big endian).

O BAM le as condicoes iniciais como arquivo sequencial UNFORMATTED com marcadores
de 4 bytes (gfortran/ifort classicos): cada registro e
    [int32 len][len bytes][int32 len]
Sem CONVERT=, o legado so funcionava em maquina big endian (SX6). Em x86 o
arquivo resultante e lixo; aqui a ordem e sempre big endian explicita.
"""

from __future__ import annotations

import struct
from pathlib import Path
from typing import Iterator, List

import numpy as np

MARKER = ">i"


def encode_record(payload: bytes) -> bytes:
    n = len(payload)
    return struct.pack(MARKER, n) + payload + struct.pack(MARKER, n)


def write_record(fh, payload: bytes) -> None:
    fh.write(encode_record(payload))


def write_float_record(fh, values, dtype: str = ">f4") -> None:
    arr = np.ascontiguousarray(np.asarray(values).ravel(), dtype=dtype)
    write_record(fh, arr.tobytes())


def write_int_record(fh, values) -> None:
    arr = np.ascontiguousarray(np.asarray(values).ravel(), dtype=">i4")
    write_record(fh, arr.tobytes())


def read_records(path: str | Path, max_records: int | None = None) -> Iterator[bytes]:
    """Itera os registros crus de um arquivo sequencial UNFORMATTED."""
    p = Path(path)
    size = p.stat().st_size
    with p.open("rb") as fh:
        pos = 0
        count = 0
        while pos < size:
            head = fh.read(4)
            if len(head) < 4:
                raise ValueError("%s: marcador truncado em %d" % (p, pos))
            (n,) = struct.unpack(MARKER, head)
            if n < 0 or pos + 8 + n > size:
                raise ValueError("%s: registro invalido em %d (len=%d)" % (p, pos, n))
            payload = fh.read(n)
            tail = fh.read(4)
            (m,) = struct.unpack(MARKER, tail)
            if m != n:
                raise ValueError("%s: marcadores inconsistentes em %d (%d != %d)" % (p, pos, n, m))
            yield payload
            pos += 8 + n
            count += 1
            if max_records is not None and count >= max_records:
                return


def read_float_record(path: str | Path, index: int, dtype: str = ">f4") -> np.ndarray:
    for i, payload in enumerate(read_records(path)):
        if i == index:
            return np.frombuffer(payload, dtype=dtype).astype(np.float64)
    raise IndexError("registro %d inexistente em %s" % (index, path))


def record_count(path: str | Path) -> int:
    return sum(1 for _ in read_records(path))


def read_direct(path: str | Path, record: int, count: int, dtype: str = ">f4") -> np.ndarray:
    """Leitura DIRECT (sem marcadores): registro de count floats."""
    p = Path(path)
    offset = (record - 1) * count * 4
    with p.open("rb") as fh:
        fh.seek(offset)
        data = fh.read(count * 4)
    if len(data) < count * 4:
        raise ValueError("%s: registro %d truncado" % (p, record))
    return np.frombuffer(data, dtype=dtype).astype(np.float64)
