#!/usr/bin/env python3
"""Exercise the actual coupling reader without importing SIROCCO dependencies.

Usage: python3 tests/check_pload.py
"""
import ast
from pathlib import Path
import sys
import tempfile
import numpy as np

source = Path(__file__).resolve().parents[1] / 'pluto_sirocco_sub.py'
node = next(n for n in ast.parse(source.read_text()).body
            if isinstance(n, ast.ClassDef) and n.name == 'pload')
namespace = {'np': np}
exec(compile(ast.Module(body=[node], type_ignores=[]), str(source), 'exec'), namespace)
pload = namespace['pload']

cases = 0
for dims in ((4, 3, 1), (4, 3, 2)):
    for names in (['rho', 'prs', 'T'], ['rho', 'vx1', 'tr1', 'T', 'ch']):
        for datatype, precision in (('dbl', 'f8'), ('flt', 'f4')):
            for endian, prefix in (('little', '<'), ('big', '>')):
                with tempfile.TemporaryDirectory() as tmp:
                    folder = Path(tmp)
                    grid = '# GEOMETRY: CARTESIAN\n'
                    for size in dims:
                        grid += str(size) + '\n'
                        grid += ''.join(f'{i+1} {i} {i+1}\n' for i in range(size))
                    (folder / 'grid.out').write_text(grid)
                    metadata = f'0 1.25 0.01 123 single_file {endian} '
                    (folder / (datatype + '.out')).write_text(metadata + ' '.join(names) + '\n')
                    expected = {}
                    blocks = []
                    for index, name in enumerate(names):
                        shape = list(dims)
                        values = np.arange(np.prod(shape)).reshape(shape) + 10000*(index+1)
                        expected[name] = values[:, :, 0] if shape[2] == 1 else values
                        blocks.append(values.transpose(2, 1, 0).astype(prefix + precision).tobytes())
                    binary = b''.join(blocks)
                    datafile = folder / ('data.0000.' + datatype)
                    datafile.write_bytes(binary)
                    result = pload(0, w_dir=str(folder), datatype=datatype)
                    assert result.NStep == 123 and result.SimTime == 1.25
                    for name, values in expected.items():
                        np.testing.assert_array_equal(getattr(result, name), values)
                    for broken in (binary[:-1], binary+b'\x00'):
                        datafile.write_bytes(broken)
                        try:
                            pload(0, w_dir=str(folder), datatype=datatype)
                        except ValueError:
                            pass
                        else:
                            raise AssertionError('Incorrect file length was accepted')
                    datafile.write_bytes(binary)
                    try:
                        pload(1, w_dir=str(folder), datatype=datatype)
                    except ValueError:
                        pass
                    else:
                        raise AssertionError('Missing metadata was accepted')
                cases += 1
print(f'PASS: {cases} HD cell layouts, 2D/3D, both precisions and byte orders; invalid files rejected')

