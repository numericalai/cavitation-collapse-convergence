"""
Analytical anchors for 3D spherical bubble collapse (choke/valve trim cavitation).

The role of this module is exactly what thermo.py did for the water-hammer case:
provide the closed-form answers the simulation is validated against, before any
CFD is trusted. Two independent references here, in increasing physical content:

  1. RAYLEIGH COLLAPSE TIME. An empty spherical cavity of radius R0 collapsing
     under constant external pressure p_inf in an incompressible liquid collapses
     completely in
         t_c = 0.914681 * R0 * sqrt(rho_L / (p_inf - p_v)).
     This is exact and parameter-free. It is the primary anchor: the simulation
     must reproduce this collapse time in the incompressible, low-driving limit.

  2. RAYLEIGH-PLESSET HISTORY. The full R(t) for a gas-filled bubble, including
     the non-condensable gas that arrests the collapse and sets the peak
     pressure. This is the reference for the whole collapse curve, not just its
     endpoint, and it predicts the minimum radius and rebound the compressible
     simulation should match until compressibility and shock formation take over
     near the instant of minimum radius.

Neither anchor requires the CFD. If the simulation disagrees with (1) in the
regime where (1) is valid, the setup is wrong and nothing downstream matters.

Fluids: water as stiffened gas (Saurel-Petitpas-Berry 2009), vapour/gas as
ideal gas. Same water parameters as the water-hammer study, deliberately, so the
two cases share a validated EOS.

Apache License 2.0.
"""

import numpy as np
from scipy.integrate import solve_ivp


# ----------------------------------------------------------------------------
# Stiffened-gas fluids (same water as the water-hammer case)
# ----------------------------------------------------------------------------

class SGFluid:
    def __init__(self, name, gamma_sg, p_inf):
        self.name = name
        self.gamma = gamma_sg
        self.p_inf = p_inf

    @property
    def mfc_gamma(self):
        return 1.0 / (self.gamma - 1.0)

    @property
    def mfc_pi_inf(self):
        return self.gamma * self.p_inf / (self.gamma - 1.0)

    def c(self, p, rho):
        return np.sqrt(self.gamma * (p + self.p_inf) / rho)


WATER = SGFluid("water", gamma_sg=4.4, p_inf=6.0e8)
GAS = SGFluid("gas", gamma_sg=1.4, p_inf=0.0)   # non-condensable / vapour

RHO_WATER = 1000.0
RHO_GAS = 1.0
P_VAPOR = 2339.0        # water vapour pressure at ~20 C [Pa]

# Rayleigh's numerical constant, from the exact quarter-period integral
#   t_c = R0 sqrt(rho/(p_inf-p_v)) * sqrt(3/2) * Beta(5/6, 1/2)/3
# Precomputed here; verify() recomputes it from the integral as a check.
RAYLEIGH_CONST = 0.914681


# ----------------------------------------------------------------------------
# 1. Rayleigh collapse time -- the primary anchor
# ----------------------------------------------------------------------------

def rayleigh_collapse_time(R0, p_inf, rho=RHO_WATER, p_v=P_VAPOR):
    """Exact total collapse time of an empty cavity under constant p_inf."""
    dp = p_inf - p_v
    if dp <= 0:
        raise ValueError(f"Driving pressure p_inf - p_v = {dp:.4g} Pa must be > 0. "
                         "An empty cavity does not collapse without net external "
                         "pressure.")
    return RAYLEIGH_CONST * R0 * np.sqrt(rho / dp)


def rayleigh_wall_speed(R, R0, p_inf, rho=RHO_WATER, p_v=P_VAPOR):
    """Interface velocity during Rayleigh collapse, from energy conservation:
        (dR/dt)^2 = (2/3)(p_inf-p_v)/rho * [(R0/R)^3 - 1].
    Diverges as R -> 0, which is why the incompressible model is only an anchor
    for the early/mid collapse, not the final instant."""
    dp = p_inf - p_v
    return -np.sqrt((2.0 / 3.0) * dp / rho * ((R0 / R) ** 3 - 1.0))


# ----------------------------------------------------------------------------
# 2. Rayleigh-Plesset -- full history with gas cushioning
# ----------------------------------------------------------------------------

def rayleigh_plesset(R0, p_inf, p_g0, kappa=1.4, rho=RHO_WATER, p_v=P_VAPOR,
                     sigma=0.0728, mu=1.0e-3, t_end=None, n=4000):
    """Integrate the Rayleigh-Plesset equation for a gas-filled bubble.

        rho[ R R'' + 1.5 R'^2 ] = p_g0 (R0/R)^{3k} + p_v - p_inf
                                  - 2 sigma / R - 4 mu R'/R

    Returns t, R, Rdot, and the minimum radius. The gas partial pressure p_g0
    at R0 is what arrests the collapse; the minimum radius and the compression
    ratio p_g0 (R0/R_min)^{3k} set the peak pressure the simulation must capture.
    """
    if t_end is None:
        t_end = 3.0 * rayleigh_collapse_time(R0, p_inf, rho, p_v)

    def rhs(t, y):
        R, V = y
        R = max(R, 1e-9 * R0)
        p_gas = p_g0 * (R0 / R) ** (3 * kappa)
        acc = (p_gas + p_v - p_inf - 2 * sigma / R - 4 * mu * V / R) / rho
        Rddot = (acc - 1.5 * V * V) / R
        return [V, Rddot]

    # Stop if the bubble rebounds past its initial radius (one full cycle).
    def rebound(t, y):
        return y[0] - R0
    rebound.direction = 1
    rebound.terminal = False

    sol = solve_ivp(rhs, (0, t_end), [R0, 0.0], method="LSODA",
                    dense_output=True, rtol=1e-10, atol=1e-12,
                    max_step=t_end / n, events=rebound)
    t = np.linspace(0, t_end, n)
    R, V = sol.sol(t)
    i_min = int(np.argmin(R))
    return {
        "t": t, "R": R, "Rdot": V,
        "R_min": float(R[i_min]),
        "t_min": float(t[i_min]),
        "compression_ratio": float(R0 / R[i_min]),
        "peak_gas_pressure": float(p_g0 * (R0 / R[i_min]) ** (3 * kappa)),
    }


def initial_gas_pressure(R0, p_inf, p_v=P_VAPOR, sigma=0.0728, equilibrium=True):
    """Gas partial pressure at R0. If the bubble starts in equilibrium,
    p_g0 = p_inf - p_v + 2 sigma / R0. For a driven collapse the bubble is
    usually started at p_g0 = p_v (or a chosen understaturation), so expose both."""
    if equilibrium:
        return p_inf - p_v + 2 * sigma / R0
    return p_v


# ----------------------------------------------------------------------------
# Near-wall scaling (Phase 3 reference, not a solver result)
# ----------------------------------------------------------------------------

def standoff_parameter(h, R0):
    """gamma = h / R0, the distance from bubble centre to wall in initial radii.
    The organizing variable for near-wall collapse: jetting onset near gamma < 3,
    peak wall loading near gamma ~ 1. Used to plan the Phase 3 sweep; the jet
    itself comes from the simulation, not from a formula."""
    return h / R0


# ----------------------------------------------------------------------------
# Verification
# ----------------------------------------------------------------------------

def verify(verbose=True):
    ok = True
    checks = []

    # Recompute the Rayleigh constant from the exact integral and compare to the
    # tabulated value. t_c = R0 sqrt(rho/dp) * INT_0^1 dx / sqrt((1/x^3)-1)
    # over the collapse, with the standard substitution giving the Beta form.
    from scipy.integrate import quad
    integrand = lambda x: 1.0 / np.sqrt((1.0 / x**3) - 1.0)
    I, _ = quad(integrand, 0, 1)
    const_computed = np.sqrt(1.5) * I
    checks.append(("Rayleigh constant vs integral", const_computed, "-",
                   RAYLEIGH_CONST * 0.999, RAYLEIGH_CONST * 1.001))

    # Water sound speed, shared with the water-hammer case.
    checks.append(("water sound speed", WATER.c(1e5, RHO_WATER), "m/s",
                   1400.0, 1700.0))

    # Collapse time for a 1 mm bubble at 100 bar -- an order-of-magnitude sanity
    # check against the microsecond timescale expected for cavitation.
    tc = rayleigh_collapse_time(1e-3, 100e5)
    checks.append(("collapse time 1mm @ 100bar", tc * 1e6, "us", 5.0, 15.0))

    # Rayleigh-Plesset must arrest: a gas-filled bubble driven to collapse has a
    # finite minimum radius, not zero. This is the DRIVEN case -- the bubble
    # starts at vapour pressure and is exposed to high p_inf, which is the actual
    # cavitation scenario. (Starting in mechanical equilibrium gives no collapse
    # at all, which is correct physics but the wrong test.)
    rp = rayleigh_plesset(1e-3, 100e5, P_VAPOR)
    checks.append(("R-P driven minimum radius", rp["R_min"] / 1e-3, "R0",
                   1e-3, 0.5))
    checks.append(("R-P compression ratio", rp["compression_ratio"], "-",
                   3.0, 1e4))

    # R-P collapse time must be close to Rayleigh (gas cushioning delays it only
    # slightly for a strongly driven collapse).
    ratio = rp["t_min"] / tc
    checks.append(("R-P t_min / Rayleigh t_c", ratio, "-", 0.95, 1.6))

    # Peak gas pressure must be far above the driving pressure -- this is the
    # damage mechanism, and it must be orders of magnitude above p_inf.
    checks.append(("R-P peak/driving pressure", rp["peak_gas_pressure"] / 100e5,
                   "-", 10.0, 1e5))

    if verbose:
        print(f"{'quantity':<32s} {'value':>12s}  {'units':<5s} {'accept':<22s} ok")
        print("-" * 82)
    for name, val, units, lo, hi in checks:
        passed = lo <= val <= hi
        ok = ok and passed
        if verbose:
            print(f"{name:<32s} {val:>12.5g}  {units:<5s} "
                  f"{f'[{lo:.4g}, {hi:.4g}]':<22s} " + ("PASS" if passed else "FAIL"))
    return ok


if __name__ == "__main__":
    ok = verify()
    print()
    print("Rayleigh collapse time across choke operating pressures (R0 = 1 mm)")
    for p_bar in (10, 50, 100, 200, 500, 1000):
        tc = rayleigh_collapse_time(1e-3, p_bar * 1e5)
        print(f"  dp = {p_bar:>5d} bar   t_c = {tc*1e6:8.3f} us")
    print()
    print("MFC fluid_pp inputs")
    for i, f in enumerate((WATER, GAS), start=1):
        print(f"  fluid_pp({i}) {f.name:<6s} gamma={f.mfc_gamma:.10g}  "
              f"pi_inf={f.mfc_pi_inf:.10g}")
    print()
    print("VERIFICATION:", "all checks passed" if ok else "FAILED")