# Vapor-Liquid Equilibrium: Bubble Point via Bisection Method + Monte Carlo Simulation

Solves for the bubble point temperature of a binary benzene-toluene mixture using
Raoult's Law and the Antoine equation, then quantifies how sensitive that result is
to real-world uncertainty in composition and pressure using a 200-trial Monte Carlo
simulation.

## Problem

Given a benzene-toluene liquid mixture at a known composition and pressure, find the
temperature at which it begins to boil (the "bubble point"). Raoult's Law gives the
total vapor pressure as:

```
P = x1 * P1_sat(T) + (1 - x1) * P2_sat(T)
```

where each component's saturation pressure follows the Antoine equation:

```
log10(P_sat) = A - B / (T + C)
```

Because `T` appears inside two separate exponential terms, this can't be solved for
algebraically — it has to be solved numerically as a root-finding problem.

## Approach

1. **Model the physics** — Antoine equation for each component's saturation pressure,
   combined via Raoult's Law into a single residual function `f(T)`.
2. **Root-find with Bisection** — repeatedly halve a bracketing interval `[350, 400] K`
   until `f(T) ≈ 0`, converging to within a specified tolerance.
3. **Verify convergence theoretically** — the number of iterations needed for a given
   error tolerance can be computed in advance (`n > log2(interval width / tolerance)`);
   this is checked against the actual number of iterations used.
4. **Quantify uncertainty with Monte Carlo** — repeat the whole solve 200 times, each
   time sampling composition (`x1`) and pressure (`P`) from realistic uniform ranges,
   to see how much the bubble point temperature varies as a result.

## Results

| Metric | Value |
|---|---|
| Nominal bubble point (x1 = 0.4, P = 1 atm) | 368.22 K (≈ 95.07 °C) |
| Bisection iterations to converge | 13 (matches theoretical minimum) |
| Monte Carlo mean T (200 trials) | 368.22 K |
| Monte Carlo std. dev. | 2.12 K |
| Monte Carlo range | 363.93 K – 372.93 K |
| Correlation: composition (x1) vs T | strongly negative |
| Correlation: pressure (P) vs T | positive |

The negative correlation between benzene fraction and bubble point makes physical
sense — benzene is the more volatile component, so more of it in the mixture lowers
the boiling point.

## Repo contents

```
├── vle_bisection_montecarlo.py   # main solver + Monte Carlo simulation
├── figures/
│   ├── convergence.png           # bisection error vs iteration (log scale)
│   └── histogram.png             # distribution of 200 simulated bubble points
├── VLE_Bisection_Report.pdf      # full write-up: methodology, results, appendix
└── README.md
```

## Running it

```bash
python3 vle_bisection_montecarlo.py
```

No dependencies beyond the Python standard library for the core solver;
`matplotlib` is used only to regenerate the figures.

## Why this project

Built as part of a portfolio piece 
demonstrating numerical root-finding, simulation-based uncertainty analysis, and 
clear technical communication — skills that carry over directly from engineering
coursework into data analysis work.

clear technical communication — skills that carry over directly from engineering
coursework into data analysis work.
