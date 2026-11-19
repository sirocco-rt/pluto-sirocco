#!/usr/bin/env python3
"""Small independent equation and input-contract checks; no production data required."""
from pathlib import Path
import math, os, subprocess, sys, tempfile
exe=Path(sys.argv[1] if len(sys.argv)>1 else './cak_v4').resolve()
passed=[]
def case(name, energy=10, scale=1, mode=9, changes=None, args=(), valid=True):
 with tempfile.TemporaryDirectory(prefix='cak-v4-test-') as tmp:
  d=Path(tmp)
  inputs={
   'sirocco_ion_data.dat':f'nions 1\nion 0 H 1 1\nnplasma 1\n0 0 {1e8*scale:.17g}\n',
   'sirocco_pcon_data.dat':f'nplasma 1\n0 0 1e10 45 40000 {1e-15*scale:.17g} {1e8*scale:.17g} {1e8*scale:.17g} .01 .01 .01\n',
   'sirocco_spec_data.dat':'model 9\nnbands 1\nnplasma 1\n1e13 1e17\n0 0 1e13 1e17 1 -10 0 0 0\n' if mode==9 else 'model 8\nnplasma 1\n0 0 40000 1\n',
   'atomic_models.txt':f'1 1 2 13.6\n1 0 2 0\n2 {energy} 4 0\n',
   'transitiondata.txt':'1 1 1\n1 1 2 1e8\n'}
  if changes:inputs.update(changes)
  for n,t in inputs.items():
   if t is not None:(d/n).write_text(t)
  r=subprocess.run([str(exe),*args],cwd=d,capture_output=True,text=True,env=os.environ)
  if not valid:
   assert r.returncode!=0,(name,'expected failure',r.stdout,r.stderr)
   assert not (d/'M_UV_data.dat').exists(),name
   assert 'AddressSanitizer' not in r.stderr and 'runtime error:' not in r.stderr,(name,r.stderr)
   passed.append(name);return
  assert r.returncode==0,(name,r.stdout,r.stderr)
  ls=(d/'M_UV_data.dat').read_text().splitlines(); ts=list(map(float,ls[1].split()[1:])); ms=list(map(float,ls[2].split()[4:]));assert len(ts)==len(ms)==32
  assert all(math.isfinite(x) and x>=0 for x in ms)
  assert all(b<=a for a,b in zip(ms,ms[1:]))
  assert 'frequency_scope=all' in (d/'M_UV_data.meta').read_text()
  if mode==9 and not changes:
   H=6.6260755e-27;C=2.99792458e10;KB=1.380658e-16;EV=1.6021772e-12;MH=1.6733e-24
   nu=energy*EV/H;vth=math.sqrt(2*KB*40000/MH);q=math.exp(-energy*EV/(KB*40000));part=2+4*q
   nl=1e8*scale*2/part;nu_pop=1e8*scale*4*q/part;bul=6.7819570e46/nu**3*1e8
   kappa=(2*bul*nl-bul*nu_pop)*1.580764662876770e-17/(1e-15*scale)/vth
   sigma=0.66524e-24/1.672661e-24
   expected=[nu*vth/C/(1e17-1e13)*(-math.expm1(-kappa*t/sigma))/t for t in ts]
   assert all(math.isclose(a,b,rel_tol=2e-6) for a,b in zip(ms,expected)),(name,ms,expected)
  passed.append(name);return ms
base=case('all-line UV transition')
opt=case('all-line optical transition',energy=2);assert opt[0]>0
xray=case('all-line X-ray transition',energy=200);assert xray[0]>0
assert base==case('consistent density rescaling',scale=.01)
case('fixed electron opacity default')
case('blackbody',mode=8)
case('zero radiation dilution',mode=8,changes={'sirocco_spec_data.dat':'model 8\nnplasma 1\n0 0 40000 0\n'})
case('mismatched cells',changes={'sirocco_spec_data.dat':'model 8\nnplasma 1\n9 9 40000 1\n'},valid=False)
for field,val in [(4,'0'),(5,'0'),(7,'0'),(5,'nan'),(4,'-1')]:
 row='0 0 1e10 45 40000 1e-15 1e8 1e8 .01 .01 .01'.split();row[field]=val
 case(f'invalid physical field {field}={val}',changes={'sirocco_pcon_data.dat':'nplasma 1\n'+' '.join(row)+'\n'},valid=False)
case('unsupported model',changes={'sirocco_spec_data.dat':'model 123\n'},valid=False)
case('truncated physical row',changes={'sirocco_pcon_data.dat':'nplasma 1\n0 0 1e10\n'},valid=False)
case('out-of-range atomic index',changes={'transitiondata.txt':'1 1 1\n1 1 3 1e8\n'},valid=False)
case('negative ion density',changes={'sirocco_ion_data.dat':'nions 1\nion 0 H 1 1\nnplasma 1\n0 0 -1\n'},valid=False)
case('no legacy fallback',changes={'sirocco_ion_data.dat':None,'py_ion_data.dat':'nions 1\nion 0 H 1 1\nnplasma 1\n0 0 1e8\n'},valid=False)
case('invalid spectral bounds',changes={'sirocco_spec_data.dat':'model 9\nnbands 1\nnplasma 1\n1e13 1e14\n0 0 1e13 1e17 1 -10 0 0 0\n'},valid=False)
case('duplicate cells',changes={'sirocco_ion_data.dat':'nions 1\nion 0 H 1 1\nnplasma 2\n0 0 1e8\n0 0 1e8\n'},valid=False)
case('truncated atomic data',changes={'atomic_models.txt':'1 1 2 13.6\n1 0 2 0\n'},valid=False)
nu=10*1.6021772e-12/6.6260755e-27
split=f'model 9\nnbands 2\nnplasma 1\n1e13 {nu:.17g} 1e17\n0 0 1e13 {nu:.17g} 1 -10 0 0 0 {nu:.17g} 1e17 1 -10 0 0 0\n'
assert base==case('line exactly on spectral boundary',changes={'sirocco_spec_data.dat':split})
case('power law alpha minus one',changes={'sirocco_spec_data.dat':'model 9\nnbands 1\nnplasma 1\n1e13 1e17\n0 0 1e13 1e17 1 -10 -1 0 0\n'})
case('exponential continuum',changes={'sirocco_spec_data.dat':'model 9\nnbands 1\nnplasma 1\n1e13 1e17\n0 0 1e13 1e17 2 0 0 1e-10 40000\n'})
for name in passed:print('PASS',name)
print(f'{len(passed)} cases passed')
