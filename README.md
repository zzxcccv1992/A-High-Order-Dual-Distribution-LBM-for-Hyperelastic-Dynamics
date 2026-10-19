# Dual-distribution LBM for hyperelastic dynamics

Python implementation of a dual-distribution-function (DDF) lattice Boltzmann method with a fixed Cartesian lattice and a tracked moving boundary. The included experiment applies sinusoidal nominal traction to a rectangular, plane-strain hyperelastic body on a 60 × 18 lattice.

## Run the experiment

```sh
python -m pip install -r requirements.txt
python run_experiment.py --out outputs/traction
python plot_displacement.py --input outputs/traction
```

Use `run_experiment.py` as the experiment entry point: it loads `experiment.json` and supplies the experiment values to the solver through `parameters.py`. The solver module's standalone defaults are general loading presets. See [the parameter mapping](docs/parameters.md) for the physical, normalized, and lattice units. Choose a new output directory for each run.

## Experiment

The reference body is 60 mm long and 18 mm high. It starts undeformed and at rest, with its left edge fixed and its top and bottom edges traction-free. The right edge carries

$$
\mathbf{T}(t)=1\,\mathrm{MPa}\,\sin(2\pi f t)\,(1,0),
\qquad f=48\,782.0520868\,\mathrm{Hz}.
$$

This is nominal traction per reference area. The run lasts one loading period, approximately 20.499343 µs. The elastic inputs are λ = 115 GPa and μ = 77 GPa, with an illustrative reference density of 7850 kg/m³.

## Files and outputs

| File | Purpose |
| --- | --- |
| `eulerian_traction.py` | DDF solver and field output |
| `experiment.json` | Physical inputs and numerical settings |
| `parameters.py` | Unit conversion and solver configuration |
| `run_experiment.py` | Run the configured experiment and save its parameter mapping |
| `plot_displacement.py` | Plot global displacement components, magnitude, and tip history |

The run writes `parameters.json`, `config.json`, `history.csv`, `snapshots.npz`, `final_eulerian.npz`, and `summary.json` to the output directory. The plotting command writes `global_ux.png`, `global_uy.png`, `global_magnitude.png`, and `tip_displacement.png` under `outputs/traction/figures`. Global displacement is shown over the full material reference grid at normalized times 0.5, 1, 1.5, 2, 3, and 4. Plot coordinates are in millimetres, displacement is in micrometres, and time is in microseconds.

In `snapshots.npz`, `positions` contains current material coordinates and `reference_points` contains their reference coordinates. Their difference is displacement in normalized length units. Multiply by 0.06 m to obtain physical displacement. The 61 × 19 material sampling points are separate from the 60 × 18 cell-centred computational resolution.

## License

See [LICENSE](LICENSE).
