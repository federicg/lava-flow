# 1D MATLAB Implementation of the Lava-Flow Scheme

A didactic 1D MATLAB implementation of the finite element time-stepping scheme used in the `lava-flow` C++ code. It solves the linear advection-decay equation

$$\partial_t q + a  \partial_x q = -\chi  q$$

on a uniform grid, using the two-stage implicit-explicit (IMEX) Runge-Kutta method described in the companion paper.

---

## Files

- `lava_flow_1d.m` — main script (self-contained, no toolboxes required beyond `fsolve`)

---

## Parameters

| Variable | Description | Default |
|---|---|---|
| `a` | Advection speed | `1` |
| `chi` | Decay coefficient | `0.05` |
| `num_elements` | Number of cells | `200` |
| `T` | Final simulation time | `100` |
| `cfl` | CFL number | `1.22` |
| `flag_ini` | Initial condition (1 = constant, 2 = Gaussian, 3 = windowed checkerboard) | `2` |
| `flag_tab` | Coefficient table for the IMEX scheme (0–3, see paper) | `3` |

---

## Usage

1. Open MATLAB and navigate to this folder
2. Run:
```matlab
lava_flow_1d
```
3. The script will plot the numerical vs exact solution live during the simulation, and print the final $$L^\infty$$ error.

> **Note:** `fsolve` from the Optimization Toolbox is required for the implicit stages.

---

## Scheme

The method uses a two-stage IMEX-RK approach solving the linear advection-decay equation on a uniform grid:

- **Stage 1** (`q_2`, defined at half-integer nodes $j+1/2$) — implicit solve on the cell centers (primal grid), coupling an explicit upwind flux term with an implicit decay term
- **Stage 2** (`q_3`, defined at integer nodes $j$) — implicit solve using fluxes from both the initial solution and stage 1
- **Update** (`q^{n+1}`) — combines contributions from all stages; the mass-lumped Q1 treatment relocates the degrees of freedom to the dual grid, where `q_3` and `q^{n+1}` live, leading to a staggered structure since stage 1 lives on the primal grid

Four coefficient sets are available via `flag_tab`:

| `flag_tab` | Description |
|---|---|
| `0` | Coefficients from [Gatti et al., JCP 2024](https://doi.org/10.1016/j.jcp.2024.112798) |
| `1` | Coefficients from [Gatti et al., CAMWA 2025](https://doi.org/10.1016/j.camwa.2025.02.014) |
| `2` | Table 2 in [Gatti et al., arXiv 2025](https://arxiv.org/abs/2509.09460) |
| `3` | Table 3 in [Gatti et al., arXiv 2025](https://arxiv.org/abs/2509.09460) |

The default is `flag_tab = 2`.

---

## Exact Solution

For smooth initial data, the exact solution is

$$q(x,t) = q_0(x - at) e^{-\chi t}$$

which is used to compute the $$L^\infty$$ error at each time step.
