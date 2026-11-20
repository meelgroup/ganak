"""Exercise native telemetry, bucket ablation, and policy output through the CLI."""
import gzip
from pathlib import Path
import random
import subprocess
import sys
import tempfile

driver = Path(__file__).with_name('policy.py')
rng = random.Random(142)
clauses = [[v*rng.choice((-1,1)) for v in rng.sample(range(1,46),3)] for _ in range(110)]
data = 'p cnf 45 110\n' + ''.join(' '.join(map(str,c))+' 0\n' for c in clauses)
with tempfile.TemporaryDirectory() as directory:
    path = Path(directory)/'sample.cnf.gz'
    path.write_bytes(gzip.compress(data.encode()))
    for flags in ([], ['--no-buckets']):
        log = subprocess.check_output([sys.executable, driver, '--kernel-timing', *flags, path], text=True)
        answers = [line for line in log.splitlines() if line.startswith('c s exact arb int')]
        assert answers == ['c s exact arb int 62028404']
        samples = [dict(kv.split('=') for kv in line.split()[3:])
                   for line in log.splitlines() if line.startswith('c frontier sample ')]
        final = [r for r in samples if r['phase']=='end']
        assert final and any(r['phase']=='begin' for r in samples)
        assert all(int(r['lookups'])>=int(r['hits']) and 'sparse_cpu_s' in r for r in final)
        calls = sum(int(r['bucket_calls']) for r in final)
        assert (calls == 0) if flags else (calls > 0)
        assert 'c medium event=finish' in log
print('Native samples, gzip policy output and both bucket settings checked.')
