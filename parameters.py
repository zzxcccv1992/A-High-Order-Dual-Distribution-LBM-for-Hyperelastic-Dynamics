"""Convert dimensional experiment inputs to normalized and lattice quantities."""
from __future__ import annotations

import json
import math
from dataclasses import asdict
from pathlib import Path

from eulerian_traction import Config

DEFAULT_CONFIG = Path(__file__).resolve().with_name('experiment.json')


def load_parameters(path=DEFAULT_CONFIG):
    inputs = json.loads(Path(path).read_text(encoding='utf-8'))
    p, n, load = inputs['physical'], inputs['numerical'], inputs['loading']
    length, height = p['length_m'], p['height_m']
    rho, lam, mu = p['rho0_kg_m3'], p['lambda_Pa'], p['mu_Pa']
    cp_target = n['longitudinal_courant']
    for name, value in dict(length=length, height=height, density=rho,
                            lam=lam, mu=mu, courant=cp_target).items():
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f'{name} must be finite and positive')
    if n['resolution'] != [60, 18]:
        raise ValueError('This experiment uses a 60 by 18 lattice')
    if not math.isclose(height/length, 18/60, rel_tol=0, abs_tol=1e-12):
        raise ValueError('Geometry must match the 60 by 18 lattice aspect ratio')
    if load['frequency_rule'] != 'cp/(2*length)' or load['cycles'] != 1:
        raise ValueError('Use one sinusoidal cycle with frequency cp/(2*length)')
    if load['direction'] != [1.0, 0.0]:
        raise ValueError('This experiment uses horizontal right-boundary traction')

    cp = math.sqrt((lam+2*mu)/rho)
    cs = math.sqrt(mu/rho)
    speed = cp/cp_target
    spacing = length/60
    dt = cp_target*spacing/cp
    stress_scale = (lam+2*mu)/(cp_target*cp_target)
    time_scale = length/speed
    frequency = cp/(2*length)
    poisson = lam/(2*(lam+mu))
    # In normalized coordinates the domain length is one and h = dt = 1/60.
    # The loading frequency is cp_target/2 cycles per normalized time unit.
    cfg = Config.for_case(
        'normal', resolution=60, length=1.0, height=height/length, rho0=1.0,
        lam=lam/stress_scale, mu=mu/stress_scale,
        amplitude=p['traction_amplitude_Pa']/stress_scale,
        omega=math.pi*cp_target, end_time=2/cp_target,
        tau_f=n['tau_f'], tau_g=n['tau_g'], padding=n['padding'],
        recovery_tolerance=n['recovery_tolerance'],
        recovery_iterations=n['recovery_iterations'],
        recovery_relaxation=n['recovery_relaxation'])
    step_count = cfg.end_time*60
    if not math.isclose(step_count, round(step_count), rel_tol=0, abs_tol=1e-10):
        raise ValueError('The duration must contain an integer number of lattice steps')
    result = dict(
        inputs=inputs,
        physical=dict(**p, youngs_modulus_Pa=2*mu*(1+poisson),
                      poisson_ratio=poisson, bulk_modulus_Pa=lam+2*mu/3,
                      cp_m_s=cp, cs_m_s=cs, lattice_speed_m_s=speed,
                      spacing_m=spacing, dt_s=dt, frequency_Hz=frequency,
                      angular_frequency_rad_s=2*math.pi*frequency,
                      period_s=1/frequency, end_time_s=1/frequency),
        scales=dict(length_m=length, time_s=time_scale, velocity_m_s=speed,
                    density_kg_m3=rho, stress_Pa=stress_scale),
        solver_config=asdict(cfg),
        normalized=dict(h=1/60, dt=1/60, frequency=cp_target/2,
                        cp=cp_target, cs=cs/speed),
        lattice=dict(spacing=1.0, dt=1.0, reference_length=60, reference_height=18,
                     rho0=cfg.rho0, lam=cfg.lam, mu=cfg.mu, amplitude=cfg.amplitude,
                     cp=cp_target, cs=cs/speed,
                     frequency_cycles_per_step=cp_target/120,
                     angular_frequency_rad_per_step=cfg.omega/60,
                     steps=round(step_count)))
    return cfg, result


if __name__ == '__main__':
    _, mapping = load_parameters()
    print(json.dumps(mapping, indent=2))
