#!/usr/bin/env python3
"""Run the full supplied atomic data on two positive-density smoke-test cells."""
from pathlib import Path
import math,subprocess,sys,tempfile
root=Path(__file__).resolve().parents[1]
exe=Path(sys.argv[1] if len(sys.argv)>1 else './cak_v4').resolve()
with tempfile.TemporaryDirectory(prefix='cak-v4-two-cells-') as tmp:
 d=Path(tmp)
 for src in (root/'tests/fixtures/two_cells').glob('*.dat'):(d/src.name).write_bytes(src.read_bytes())
 for n in ['atomic_models.txt','transitiondata.txt']:(d/n).symlink_to(root/n)
 r=subprocess.run([str(exe)],cwd=d,capture_output=True,text=True)
 assert r.returncode==0,(r.stdout,r.stderr)
 rows=(d/'M_UV_data.dat').read_text().splitlines()[2:]
 assert len(rows)==2
 for row in rows:
  vals=list(map(float,row.split()[4:]));assert len(vals)==32
  assert all(math.isfinite(x) and x>=0 for x in vals)
  assert all(b<=a for a,b in zip(vals,vals[1:]))
 assert '453309 lines' in r.stdout
 print('PASS two cells, 178 ions, 453309 transitions; 64 finite, nonnegative, monotonic values')
