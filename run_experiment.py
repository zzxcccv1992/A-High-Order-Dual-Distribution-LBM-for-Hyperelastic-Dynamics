"""Run the 60 by 18 sinusoidal-traction DDF experiment."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(name, '1')

from eulerian_traction import run
from parameters import DEFAULT_CONFIG, load_parameters


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=DEFAULT_CONFIG)
    parser.add_argument('--out', type=Path, default=Path('outputs/traction'))
    args = parser.parse_args()
    cfg, mapping = load_parameters(args.config)
    args.out.mkdir(parents=True, exist_ok=False)
    (args.out/'parameters.json').write_text(
        json.dumps(mapping, indent=2)+'\n', encoding='utf-8')
    run(cfg, args.out)


if __name__ == '__main__':
    main()
