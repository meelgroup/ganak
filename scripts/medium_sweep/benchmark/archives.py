from contextlib import contextmanager
from pathlib import Path
import tarfile,lzma,re
import os
base=Path(os.environ['GANAK_DATA'])
@contextmanager
def cnf_stream(year,i):
    if year in (2022,2025):
        path=base/'heur/reference'/f'mc{year%100}.tar'
        with tarfile.open(path) as ar:
            member=next(m for m in ar if m.name.endswith(f'_{i:03}.cnf.xz'))
            with ar.extractfile(member) as packed,lzma.open(packed) as cnf:yield cnf
    else:
        path=base/f'mc{year}-track1-mc_competition'/f'mc{year}_track1_{i:03}.cnf.xz'
        with lzma.open(path) as cnf:yield cnf
