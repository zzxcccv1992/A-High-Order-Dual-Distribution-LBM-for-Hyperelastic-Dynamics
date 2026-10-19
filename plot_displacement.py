"""Plot global DDF displacement fields and tip motion in dimensional units."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def plot(input_dir, output_dir):
    mapping = json.loads((input_dir/'parameters.json').read_text(encoding='utf-8'))
    scales = mapping['scales']
    requested = mapping['inputs']['output']['snapshot_times_code']
    with np.load(input_dir/'snapshots.npz', allow_pickle=False) as data:
        reference = data['reference_points']
        times = data['time']
        displacement = (data['positions']-reference[None])*scales['length_m']*1e6
    indices = []
    for target in requested:
        hits = np.flatnonzero(np.isclose(times, target, rtol=0, atol=1e-10))
        if len(hits) != 1:
            raise ValueError(f'No unique saved field at normalized time {target}')
        indices.append(int(hits[0]))
    if not indices:
        raise ValueError('At least one output time is required')
    output_dir.mkdir(parents=True, exist_ok=True)
    x, y = np.moveaxis(reference*scales['length_m']*1e3, -1, 0)
    selected = displacement[indices]
    fields = dict(ux=selected[..., 0], uy=selected[..., 1],
                  magnitude=np.linalg.norm(selected, axis=-1))
    labels = dict(ux=r'$u_x$', uy=r'$u_y$', magnitude=r'$|\mathbf{u}|$')
    for name, values in fields.items():
        rows = (len(indices)+1)//2
        fig, axes = plt.subplots(rows, 2, figsize=(11, 2.55*rows),
                                 squeeze=False, constrained_layout=True)
        vmax = max(float(np.max(np.abs(values))), np.finfo(float).eps)
        vmin = 0.0 if name == 'magnitude' else -vmax
        for ax, k, value in zip(axes.flat, indices, values):
            mesh = ax.pcolormesh(x, y, value, shading='gouraud',
                                 cmap='viridis' if name == 'magnitude' else 'RdBu_r',
                                 vmin=vmin, vmax=vmax, rasterized=True)
            ax.set_aspect('equal')
            ax.set(xlabel='Reference X (mm)', ylabel='Reference Y (mm)',
                   title=f't = {times[k]*scales["time_s"]*1e6:.3f} μs')
            ax.set_xlim(float(x.min()), float(x.max()))
            ax.set_ylim(float(y.min()), float(y.max()))
        for ax in list(axes.flat)[len(indices):]:
            ax.set_visible(False)
        fig.colorbar(mesh, ax=list(axes.flat)[:len(indices)], shrink=.85,
                     label=f'{labels[name]} (μm)')
        fig.suptitle(f'DDF global displacement: {labels[name]}')
        fig.savefig(output_dir/f'global_{name}.png', dpi=180)
        plt.close(fig)

    history = np.genfromtxt(input_dir/'history.csv', delimiter=',', names=True)
    time_us = history['time']*scales['time_s']*1e6
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True,
                             constrained_layout=True)
    axes[0].plot(time_us, history['tip_x']*scales['length_m']*1e6, label=r'$u_x$')
    axes[0].plot(time_us, history['tip_y']*scales['length_m']*1e6, label=r'$u_y$')
    axes[0].set(ylabel='Tip displacement (μm)', title='Right-boundary midpoint')
    axes[0].legend()
    axes[1].plot(time_us, history['traction_x']*scales['stress_Pa']/1e6)
    axes[1].set(xlabel='Time (μs)', ylabel='Nominal traction (MPa)')
    for ax in axes:
        ax.grid(alpha=.25)
    fig.savefig(output_dir/'tip_displacement.png', dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path('outputs/traction'))
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    plot(args.input, args.out or args.input/'figures')
