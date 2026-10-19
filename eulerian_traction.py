"""Eulerian DDF-LBM with a tracked moving traction boundary.

Both populations collide/stream on a fixed Cartesian lattice. Material markers
provide surface geometry and output, never the bulk deformation gradient.
Units: dx=dt in normalized solver units, lattice streaming speed = 1.
"""
from __future__ import annotations

import argparse
import csv
import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from matplotlib.path import Path as PolygonPath
from scipy.ndimage import distance_transform_edt

E = np.array([[0, 0], [1, 0], [0, 1], [-1, 0], [0, -1],
              [1, 1], [-1, 1], [-1, -1], [1, -1]], dtype=int)
W = np.array([4 / 9] + [1 / 9] * 4 + [1 / 36] * 4)
OPP = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6])
CS2 = 1 / 3
I2 = np.eye(2)


@dataclass
class Config:
    case: str = 'normal'
    resolution: int = 60
    length: float = 1.0
    height: float = 0.3
    rho0: float = 1.0
    mu: float = 0.225
    lam: float = 0.15
    tau_f: float = 0.85
    tau_g: float = 0.85
    amplitude: float = 0.002
    omega: float = 2 * np.pi * 0.3
    end_time: float = 25 / 3
    padding: int = 8
    recovery_tolerance: float = 2e-11
    recovery_iterations: int = 40
    recovery_relaxation: float = 1.0
    grid_shift_y: float = 0.0
    load_direction_x: float = 1.0
    load_direction_y: float = 0.0

    @classmethod
    def for_case(cls, case='normal', poisson=None, **overrides):
        """Loading presets; mixed components each have amplitude T0."""
        if case not in ('normal', 'mixed'):
            raise ValueError('Unknown loading case')
        values = dict(case=case)
        if case == 'mixed':
            values.update(mu=.2, lam=.3, load_direction_y=1.,
                          recovery_iterations=80, recovery_relaxation=.65)
            # An exterior buffer allows the bending surface to cross grid cells.
            resolution = overrides.get('resolution', 60)
            values['padding'] = max(8, int(np.ceil(.2*resolution))+2)
        values.update(overrides)
        if poisson is not None:
            if not np.isfinite(poisson) or not -1 < poisson < .5:
                raise ValueError('Poisson ratio must lie strictly between -1 and 0.5')
            mu = values.get('mu', cls.mu)
            values['lam'] = 2*mu*poisson/(1-2*poisson)
        return cls(**values)

    def nominal_traction(self, t):
        return self.amplitude*np.sin(self.omega*t)*np.array(
            [self.load_direction_x, self.load_direction_y])


def hermite(u):
    q = np.einsum('qa,...a->...q', E, u)
    usq = np.sum(u * u, axis=-1)[..., None]
    return (1 + q / CS2 + (q*q - CS2*usq) / (2*CS2**2)
            + (q**3 - 3*CS2*usq*q) / (6*CS2**3)
            + (q**4 - 6*CS2*usq*q*q + 3*CS2**2*usq*usq) / (24*CS2**4))


def f_equilibrium(rho, u):
    return rho[..., None] * W * hermite(u)


def g_equilibrium(F, u):
    return F[..., None] * (W * hermite(u))[..., None, None, :]


def omega_f(S, u):
    q = np.einsum('qa,...a->...q', E, u)
    term = (E-u[..., None, :])/CS2 + q[..., None]*E/CS2**2
    return W * np.einsum('...qa,...a->...q', term, S)


def omega_g(G, u):
    q = np.einsum('qa,...a->...q', E, u)
    return G[..., None] * (W*(1+q/CS2))[..., None, None, :]


def det2(A):
    return A[..., 0, 0]*A[..., 1, 1]-A[..., 0, 1]*A[..., 1, 0]


def solve2(A, B):
    d = det2(A)
    if np.min(d) <= 0 or not np.all(np.isfinite(d)):
        raise FloatingPointError('Singular half-source recovery')
    inv = np.empty_like(A)
    inv[..., 0, 0], inv[..., 1, 1] = A[..., 1, 1]/d, A[..., 0, 0]/d
    inv[..., 0, 1], inv[..., 1, 0] = -A[..., 0, 1]/d, -A[..., 1, 0]/d
    return inv @ B


def stress(F, mu, lam):
    J = det2(F)
    if np.min(J) <= 1e-6 or not np.all(np.isfinite(J)):
        raise FloatingPointError(f'Invalid deformation: min J={np.nanmin(J)}')
    return mu/J[..., None, None] * (F @ np.swapaxes(F, -1, -2)-I2) + (
        lam*np.log(J)/J)[..., None, None]*I2


def traction_correction(sigma, normal, traction):
    delta = traction - np.einsum('...ij,...j->...i', sigma, normal)
    return (delta[..., :, None]*normal[..., None, :]
            + normal[..., :, None]*delta[..., None, :]
            - np.sum(delta*normal, axis=-1)[..., None, None]
            * normal[..., :, None]*normal[..., None, :])


def cross2(a, b):
    return a[..., 0]*b[..., 1]-a[..., 1]*b[..., 0]


class Solver:
    def __init__(self, cfg):
        self.cfg = cfg
        self.h = cfg.length/cfg.resolution
        self.dt = self.h
        self.ny = round(cfg.height/self.h)
        if abs(self.ny*self.h-cfg.height) > 1e-10:
            raise ValueError('Height must have an integer number of lattice cells')
        if min(cfg.tau_f, cfg.tau_g) <= .5:
            raise ValueError('Relaxation times must exceed 0.5')
        pad = cfg.padding
        self.x = (np.arange(cfg.resolution+2*pad)+.5-pad)*self.h
        self.y = (np.arange(self.ny+2*pad)+.5-pad+cfg.grid_shift_y)*self.h
        self.xy = np.stack(np.meshgrid(self.x, self.y, indexing='ij'), axis=-1)
        self.shape = self.xy.shape[:2]
        self.reference = np.stack(np.meshgrid(np.linspace(0, cfg.length, 61),
                                              np.linspace(0, cfg.height, 19), indexing='ij'), axis=-1)
        self.markers = self.reference.copy()
        # Surface samples resolve each lattice cell, independent of output markers.
        bx = np.linspace(0, cfg.length, cfg.resolution+1)
        by = np.linspace(0, cfg.height, self.ny+1)
        self.front0 = np.concatenate((np.column_stack((bx[:-1], np.zeros_like(bx[:-1]))),
            np.column_stack((np.full_like(by[:-1], cfg.length), by[:-1])),
            np.column_stack((bx[:0:-1], np.full_like(bx[:0:-1], cfg.height))),
            np.column_stack((np.zeros_like(by[:0:-1]), by[:0:-1]))))
        self.front = self.front0.copy()
        self.face_ids = np.repeat([0, 1, 2, 3], [cfg.resolution, self.ny, cfg.resolution, self.ny])
        self.ds0 = np.linalg.norm(np.roll(self.front0, -1, axis=0)-self.front0, axis=1)
        self.fixed_front = self.front0[:, 0] == 0
        self.mask = np.zeros(self.shape, bool)
        self.switches = 0
        self.geometry()
        rho = np.full(self.shape, cfg.rho0)
        self.u = np.zeros(self.shape+(2,))
        self.previous_u = None
        F = np.broadcast_to(I2, self.shape+(2, 2)).copy()
        self.f = f_equilibrium(rho, self.u)
        self.g = g_equilibrium(F, self.u)
        self.step_number = 0
        self.recovery_residual = 0.
        self.face_residual = 0.
        self.constitutive_residual = 0.

    def sample(self, a, points):
        """Bilinear sampling of any trailing tensor shape on the fixed grid."""
        shape = points.shape[:-1]
        ij = (points.reshape(-1, 2)-[self.x[0], self.y[0]])/self.h
        base = np.floor(ij).astype(int)
        base[:, 0] = np.clip(base[:, 0], 0, self.shape[0]-2)
        base[:, 1] = np.clip(base[:, 1], 0, self.shape[1]-2)
        f = np.clip(ij-base, 0, 1)
        result = np.zeros((len(ij),)+a.shape[2:])
        expand = (slice(None),)+(None,)*(a.ndim-2)
        for di, dj in ((0, 0), (0, 1), (1, 0), (1, 1)):
            w = (f[:, 0] if di else 1-f[:, 0])*(f[:, 1] if dj else 1-f[:, 1])
            result += w[expand]*a[base[:, 0]+di, base[:, 1]+dj]
        return result.reshape(shape+a.shape[2:])

    def extend(self, a):
        return a[tuple(self.nearest)]

    def velocity_extension(self, u):
        result = self.extend(u)
        left = self.xy[..., 0] < 0
        mirror = self.xy[left].copy()
        mirror[:, 0] *= -1
        result[left] = -self.sample(result, mirror)
        return result

    def geometry(self):
        oldmask = self.mask
        self.mask = PolygonPath(self.front).contains_points(self.xy.reshape(-1, 2)).reshape(self.shape)
        if not np.any(self.mask):
            raise FloatingPointError('Body left computational domain')
        if np.any(self.mask[:2]) or np.any(self.mask[-2:]) or np.any(self.mask[:, :2]) or np.any(self.mask[:, -2:]):
            raise FloatingPointError('Increase background padding: material approaches outer grid')
        self.nearest = distance_transform_edt(~self.mask, return_distances=False, return_indices=True)
        self.switches += int(np.count_nonzero(self.mask ^ oldmask)) if np.any(oldmask) else 0
        seg = np.roll(self.front, -1, axis=0)-self.front
        ds = np.linalg.norm(seg, axis=-1)
        if np.min(ds) < .1*self.h:
            raise FloatingPointError('Surface segment collapsed')
        self.normals = np.column_stack((seg[:, 1], -seg[:, 0]))/ds[:, None]
        self.stretch = ds/self.ds0
        self.cuts = []
        # Separate crossing per spatial face: both faces act at corners.
        for axis in (0, 1):
            for sign in (-1, 1):
                neighbour = np.roll(self.mask, -sign, axis=axis)
                ii, jj = np.where(self.mask & ~neighbour)
                if not len(ii):
                    continue
                p = self.xy[ii, jj]
                direction = np.zeros(2)
                direction[axis] = sign*self.h
                ap = self.front[None, :, :]-p[:, None, :]
                denominator = cross2(direction, seg)
                safe = np.where(np.abs(denominator)>1e-14, denominator, np.nan)
                q = cross2(ap, seg[None, :, :])/safe
                s = cross2(ap, direction)/safe
                valid = (q >= -1e-9) & (q <= 1+1e-9) & (s >= -1e-9) & (s <= 1+1e-9)
                candidates = np.where(valid, q, np.inf)
                k = np.argmin(candidates, axis=1)
                fraction = candidates[np.arange(len(p)), k]
                if np.any(~np.isfinite(fraction)):
                    raise FloatingPointError('Cannot intersect active cut link with surface')
                xb = p+fraction[:, None]*direction
                pg = p+direction
                n = self.normals[k]
                distance = np.sum((pg-xb)*n, axis=-1)
                image = pg-2*distance[:, None]*n
                self.cuts.append(dict(axis=axis, sign=sign, ii=ii, jj=jj, segment=k,
                                      image=image, point=xb, fraction=fraction))
        return oldmask

    def H(self, u):
        ext = self.velocity_extension(u)
        grad = np.stack([.5*(np.roll(ext, -1, axis=a)-np.roll(ext, 1, axis=a)) for a in (0, 1)], axis=-1)
        div = grad[..., 0, 0]+grad[..., 1, 1]
        return grad+div[..., None, None]*I2

    def source(self, rho, F, t, diagnostics=False):
        sigma = self.extend(stress(F, self.cfg.mu, self.cfg.lam))
        erho = self.extend(rho)
        S = np.zeros(self.shape+(2,))
        residual, constitutive = 0., 0.
        for axis in (0, 1):
            for sign in (-1, 1):
                col = np.roll(sigma[..., :, axis], -sign, axis=axis).copy()
                for cut in self.cuts:
                    if cut['axis'] != axis or cut['sign'] != sign:
                        continue
                    k = cut['segment']
                    n = self.normals[k]
                    simage = self.sample(sigma, cut['image'])
                    traction = np.zeros((len(k), 2))
                    right = self.face_ids[k] == 1
                    traction[right] = self.cfg.nominal_traction(t)[None, :]/self.stretch[k[right], None]
                    correction = traction_correction(simage, n, traction)
                    free = self.face_ids[k] != 3
                    sghost = simage.copy()
                    sghost[free] += 2*correction[free]
                    col[cut['ii'], cut['jj']] = sghost[..., :, axis]
                    if diagnostics and np.any(free):
                        sb = .5*(simage+sghost)
                        rr = np.einsum('nij,nj->ni', sb, n)-traction
                        residual = max(residual, float(np.max(np.linalg.norm(rr[free], axis=-1))))
                        raw = self.sample(sigma, cut['point'])
                        rr = np.einsum('nij,nj->ni', raw, n)-traction
                        # Exclude the clamp neighbourhood; right corner effects
                        # remain in this sampled constitutive-face diagnostic.
                        p = cut['point']
                        interior_face = free & (p[:, 0] > 3*self.h)
                        if np.any(interior_face):
                            constitutive = max(constitutive, float(np.max(np.linalg.norm(rr[interior_face], axis=-1))))
                S += sign*.5*col
            S[..., axis] += CS2*.5*(np.roll(erho, -1, axis=axis)-np.roll(erho, 1, axis=axis))
        if diagnostics:
            self.face_residual = residual
            self.constitutive_residual = constitutive
        return self.extend(S)

    def recover(self, t):
        rho = self.extend(self.f.sum(axis=-1))
        if np.min(rho) <= 0 or not np.all(np.isfinite(rho)):
            raise FloatingPointError('Invalid density')
        m = self.extend(np.einsum('...q,qa->...a', self.f, E))
        M = self.extend(self.g.sum(axis=-1))
        u = self.u.copy()
        for iteration in range(self.cfg.recovery_iterations):
            H = self.H(u)
            F = self.extend(solve2(I2-.5*H, M))
            S = self.source(rho, F, t)
            unew = (m+.5*S)/rho[..., None]
            residual = float(np.max(np.abs(unew[self.mask]-u[self.mask])))
            u = u+self.cfg.recovery_relaxation*(unew-u)
            if residual < self.cfg.recovery_tolerance:
                break
        else:
            raise FloatingPointError(f'Coupled recovery did not converge: residual={residual:.3g}')
        H = self.H(u)
        F = self.extend(solve2(I2-.5*H, M))
        S = self.source(rho, F, t, diagnostics=True)
        self.recovery_residual = float(np.max(np.abs(((m+.5*S)/rho[..., None]-u)[self.mask])))
        self.u = u
        G = H @ F
        if np.max(np.linalg.norm(u[self.mask], axis=-1))/np.sqrt(CS2) > .2:
            raise FloatingPointError('Mach number exceeds 0.2')
        return rho, u, F, S, G

    def advance(self, fields, t):
        rho, u, F, S, G = fields
        cf, cg = self.cfg.tau_f, self.cfg.tau_g
        fp = self.f-(self.f-f_equilibrium(rho, u))/cf+(1-.5/cf)*omega_f(S, u)
        gp = self.g-(self.g-g_equilibrium(F, u))/cg+(1-.5/cg)*omega_g(G, u)
        fp, gp = self.extend(fp), self.extend(gp)
        fnew, gnew = np.empty_like(fp), np.empty_like(gp)
        for q, (cx, cy) in enumerate(E):
            fnew[..., q] = np.roll(fp[..., q], (cx, cy), axis=(0, 1))
            gnew[..., q] = np.roll(gp[..., q], (cx, cy), axis=(0, 1))
            # Stationary halfway bounce-back on links crossing the true left wall.
            crossing = self.mask & (self.xy[..., 0]-cx*self.h < 0)
            fnew[..., q][crossing] = fp[..., OPP[q]][crossing]
        old_u_ext = self.velocity_extension(u)
        # Spatial midpoint plus second-order extrapolation to t+dt/2. The
        # velocity is defined at fixed spatial nodes, including when the front
        # moves; this is not integration of u at fixed reference coordinates.
        mid_u = u if self.previous_u is None else 1.5*u-.5*self.previous_u
        mid_u_ext = self.velocity_extension(mid_u)
        def move(points):
            v = self.sample(old_u_ext, points)
            midpoint = points+.5*self.dt*v
            return points+self.dt*self.sample(mid_u_ext, midpoint)
        self.markers = move(self.markers)
        self.markers[0] = self.reference[0]
        self.front = move(self.front)
        self.front[self.fixed_front] = self.front0[self.fixed_front]
        oldmask = self.geometry()
        newly_active = self.mask & ~oldmask
        if np.any(newly_active):
            oldnearest = distance_transform_edt(~oldmask, return_distances=False, return_indices=True)
            fnew[newly_active] = fnew[tuple(oldnearest[:, newly_active])]
            gnew[newly_active] = gnew[tuple(oldnearest[:, newly_active])]
        self.f, self.g = self.extend(fnew), self.extend(gnew)
        self.previous_u = u.copy()
        self.step_number += 1

    def record(self, fields, t):
        rho, u, F, S, G = fields
        sampled_rho = self.sample(self.extend(rho), self.markers)
        sampled_F = self.sample(self.extend(F), self.markers)
        J = det2(F)[self.mask]
        displacement = self.markers-self.reference
        # Geometric material-cell areas with sampled density: no binary mass count.
        a, b, c, d = self.markers[:-1, :-1], self.markers[1:, :-1], self.markers[1:, 1:], self.markers[:-1, 1:]
        area = .5*(cross2(b-a, c-a)+cross2(c-a, d-a))
        density_cell = .25*(sampled_rho[:-1, :-1]+sampled_rho[1:, :-1]+sampled_rho[1:, 1:]+sampled_rho[:-1, 1:])
        nominal = self.cfg.nominal_traction(t)
        return dict(time=t, traction=self.cfg.amplitude*np.sin(self.cfg.omega*t),
                    traction_x=float(nominal[0]), traction_y=float(nominal[1]),
                    tip_x=float(displacement[-1, 9, 0]), tip_y=float(displacement[-1, 9, 1]),
                    max_displacement=float(np.max(np.linalg.norm(displacement, axis=-1))),
                    min_J=float(np.min(J)), max_mach=float(np.max(np.linalg.norm(u[self.mask], axis=-1))/np.sqrt(CS2)),
                    rhoJ_error=float(np.max(np.abs(rho[self.mask]*J/self.cfg.rho0-1))),
                    mass=float(np.sum(area*density_cell)), face_traction_residual=self.face_residual,
                    constitutive_traction_residual=self.constitutive_residual,
                    active_nodes=int(self.mask.sum()), node_switches=self.switches,
                    recovery_residual=self.recovery_residual), dict(positions=self.markers.copy(),
                    velocity=self.sample(self.velocity_extension(u), self.markers), F=sampled_F, density=sampled_rho)


def snapshot_steps(end_time, dt, interval=1/12):
    """Retain native states bracketing each requested output time."""
    if dt <= 0 or interval <= 0 or end_time < 0:
        raise ValueError('Invalid snapshot timing')
    targets = np.r_[np.arange(int(np.floor(end_time/interval+1e-10))+1)*interval, end_time]
    raw = targets/dt
    nearest = np.rint(raw)
    raw = np.where(np.isclose(raw, nearest, rtol=0., atol=1e-10), nearest, raw)
    return set(np.clip(np.r_[np.floor(raw), np.ceil(raw)].astype(int), 0, round(end_time/dt)))


def run(cfg, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    (out/'config.json').write_text(json.dumps(asdict(cfg), indent=2)+'\n')
    start = time.perf_counter()
    solver = Solver(cfg)
    steps = round(cfg.end_time/solver.dt)
    if abs(steps*solver.dt-cfg.end_time) > 1e-10:
        raise ValueError('end_time must be a multiple of dt')
    history, snapshots, times = [], [], []
    saved_steps = snapshot_steps(cfg.end_time, solver.dt)
    for step in range(steps+1):
        t = step*solver.dt
        fields = solver.recover(t)
        row, snap = solver.record(fields, t)
        history.append(row)
        if step in saved_steps:
            snapshots.append(snap)
            times.append(t)
        if step % max(1, steps//10) == 0:
            print(json.dumps(row), flush=True)
        if step < steps:
            solver.advance(fields, t)
    with (out/'history.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(history[0]))
        writer.writeheader()
        writer.writerows(history)
    np.savez_compressed(out/'snapshots.npz', time=np.array(times), reference_points=solver.reference,
                        **{key: np.stack([s[key] for s in snapshots]) for key in snapshots[0]})
    np.savez_compressed(out/'final_eulerian.npz', grid=solver.xy, mask=solver.mask,
                        rho=fields[0], velocity=fields[1], F=fields[2], front=solver.front,
                        f=solver.f, g=solver.g)
    summary = dict(config=asdict(cfg), steps=steps, elapsed_seconds=time.perf_counter()-start,
                   max_rhoJ_error=max(r['rhoJ_error'] for r in history),
                   minimum_J=min(r['min_J'] for r in history),
                   maximum_mach=max(r['max_mach'] for r in history),
                   maximum_displacement=max(r['max_displacement'] for r in history),
                   maximum_relative_mass_drift=max(abs(r['mass']/(cfg.length*cfg.height*cfg.rho0)-1) for r in history),
                   maximum_face_traction_residual=max(r['face_traction_residual'] for r in history),
                   maximum_constitutive_traction_residual=max(r['constitutive_traction_residual'] for r in history),
                   node_switches=solver.switches,
                   maximum_recovery_residual=max(r['recovery_residual'] for r in history))
    (out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary), flush=True)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=['normal', 'mixed'], default='normal')
    parser.add_argument('--resolution', type=int, default=60)
    parser.add_argument('--poisson', type=float, default=None,
                        help='Override Poisson ratio while retaining the case shear modulus')
    parser.add_argument('--amplitude', type=float, default=.002)
    parser.add_argument('--tau-f', type=float, default=.85)
    parser.add_argument('--tau-g', type=float, default=.85)
    parser.add_argument('--end-time', type=float, default=25/3)
    parser.add_argument('--grid-shift-y', type=float, default=0.)
    parser.add_argument('--padding', type=int, default=None)
    parser.add_argument('--recovery-iterations', type=int, default=None,
                        help='Maximum algebraic recovery iterations; tolerance and dynamics stay unchanged')
    parser.add_argument('--out', default=None)
    args = parser.parse_args()
    cfg = Config.for_case(args.case, poisson=args.poisson, resolution=args.resolution, amplitude=args.amplitude, tau_f=args.tau_f,
                          tau_g=args.tau_g, end_time=args.end_time, grid_shift_y=args.grid_shift_y)
    if args.padding is not None:
        cfg.padding = args.padding
    if args.recovery_iterations is not None:
        if args.recovery_iterations < 1:
            parser.error('--recovery-iterations must be positive')
        cfg.recovery_iterations = args.recovery_iterations
    default_out = (f'outputs/mixed/n{args.resolution}' if args.case == 'mixed' else f'outputs/n{args.resolution}')
    if args.poisson is not None:
        tag = f'{args.poisson:.2f}'.replace('.', '')
        default_out = f'outputs/nu{tag}/{args.case}/n{args.resolution}'
    out = args.out or default_out
    run(cfg, out)
