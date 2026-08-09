"""
=============================================================================
VAPOR-LIQUID EQUILIBRIUM (RAOULT'S LAW) - BUBBLE POINT VIA BISECTION METHOD
=============================================================================
Problem:
    For a binary mixture (benzene + toluene), the bubble point temperature T
    satisfies:

        P = x1 * P1_sat(T) + (1 - x1) * P2_sat(T)

    where P1_sat and P2_sat are the saturation (vapor) pressures of benzene
    and toluene, given by the Antoine equation:

        log10(Psat) = A - B / (T + C)

    Given x1 (liquid mole fraction of benzene) and P (total pressure), we
    need to find T. Since T cannot be isolated algebraically (it appears
    inside two separate exponential terms), we solve for it numerically
    using the BISECTION METHOD.

    We then wrap this solver in a MONTE CARLO SIMULATION: we run it 200
    times, each time with randomly sampled x1 and P, to see how the bubble
    point temperature varies under uncertainty in composition and pressure.
=============================================================================
"""

import math        # for log2 (used to work out how many iterations we need)
import random       # for generating random x1 and P values (Monte Carlo)
import statistics   # for computing mean/standard deviation of results


# ----------------------------------------------------------------------
# STEP 1: ANTOINE EQUATION CONSTANTS
# ----------------------------------------------------------------------
# These constants come from standard reference tables (NIST).
# They let us calculate the saturation pressure (in bar) of each pure
# component at a given temperature T (in Kelvin):
#
#       log10(Psat) = A - B / (T + C)
#
ANTOINE = {
    "benzene": {"A": 4.72583, "B": 1660.652, "C": -1.461},
    "toluene": {"A": 4.07827, "B": 1343.943, "C": -53.773},
}

# Conversion factor: 1 atm = 1.01325 bar
# (Antoine constants above give pressure in bar, but the problem gives P in atm,
# so we need this to convert between the two units.)
BAR_PER_ATM = 1.01325


# ----------------------------------------------------------------------
# STEP 2: FUNCTION TO COMPUTE SATURATION PRESSURE OF ONE COMPONENT
# ----------------------------------------------------------------------
def p_sat(T, component):
    """
    Calculates the saturation pressure (bar) of a given component
    (either 'benzene' or 'toluene') at temperature T (Kelvin),
    using the Antoine equation.
    """
    c = ANTOINE[component]                       # look up A, B, C for this component
    return 10 ** (c["A"] - c["B"] / (T + c["C"])) # Antoine equation, solved for Psat


# ----------------------------------------------------------------------
# STEP 3: THE FUNCTION WHOSE ROOT WE WANT TO FIND (f(T) = 0)
# ----------------------------------------------------------------------
def bubble_point_residual(T, x1, P_atm):
    """
    This represents Raoult's Law rearranged into "residual" form:

        f(T) = [x1 * P1_sat(T) + (1 - x1) * P2_sat(T)] - P_target

    When f(T) = 0, the calculated mixture pressure exactly equals the
    target pressure P_target -> that T is the bubble point temperature.

    Bisection will search for the value of T that makes this function zero.
    """
    P_target_bar = P_atm * BAR_PER_ATM  # convert target pressure from atm to bar

    # Total pressure predicted by Raoult's Law at temperature T
    P_calc = x1 * p_sat(T, "benzene") + (1 - x1) * p_sat(T, "toluene")

    # Residual: difference between calculated and target pressure
    return P_calc - P_target_bar


# ----------------------------------------------------------------------
# STEP 4: THE BISECTION METHOD (GENERIC ROOT FINDER)
# ----------------------------------------------------------------------
def bisection(func, a, b, tol=0.005, max_iter=100):
    """
    Finds a root of `func` within the interval [a, b] using bisection.

    How bisection works:
      1. Check that func(a) and func(b) have opposite signs (this guarantees
         a root lies somewhere between them, by the Intermediate Value Theorem).
      2. Evaluate the function at the midpoint c = (a+b)/2.
      3. Replace whichever endpoint (a or b) has the SAME sign as func(c)
         with c. This shrinks the bracket [a, b] by half each time.
      4. Repeat until the bracket width is smaller than the tolerance `tol`.

    Returns:
        (root, n_iterations) - the estimated root and how many iterations
        it took to converge. Returns (None, 0) if no sign change is found
        (meaning no root exists in [a, b]).
    """
    fa, fb = func(a), func(b)   # evaluate function at both endpoints

    # If fa and fb have the same sign, there's no guaranteed root in [a,b]
    if fa * fb > 0:
        return None, 0

    n_iter = 0
    # Keep halving the interval until it's smaller than the tolerance
    while (b - a) / 2 > tol and n_iter < max_iter:
        c = (a + b) / 2       # midpoint of current interval
        fc = func(c)          # function value at midpoint

        if fa * fc < 0:
            # Root lies between a and c -> shrink interval to [a, c]
            b, fb = c, fc
        else:
            # Root lies between c and b -> shrink interval to [c, b]
            a, fa = c, fc

        n_iter += 1

    # Final estimate: midpoint of the (now very small) remaining interval
    return (a + b) / 2, n_iter


# ----------------------------------------------------------------------
# STEP 5: THEORETICAL ITERATION COUNT (ERROR ANALYSIS)
# ----------------------------------------------------------------------
def iterations_required(a, b, tol):
    """
    Bisection guarantees the interval width halves every iteration.
    Starting width = (b - a). After n iterations, width = (b-a) / 2^n.

    We want (b-a) / 2^n < tol, so solving for n:

        n > log2( (b-a) / tol )

    This tells us, in advance, exactly how many iterations are needed
    to guarantee a given error tolerance -- without running the loop.
    """
    return math.ceil(math.log2((b - a) / tol))


# ----------------------------------------------------------------------
# STEP 6: MONTE CARLO SIMULATION
# ----------------------------------------------------------------------
def monte_carlo_bubble_point(n_trials=200, T_bounds=(350, 400),
                              x1_bounds=(0.3, 0.5), P_bounds=(0.95, 1.05),
                              tol=0.005, seed=42):
    """
    Runs the bubble-point calculation many times (n_trials), each time with
    randomly sampled inputs:

        x1 ~ Uniform(0.3, 0.5)        (liquid mole fraction of benzene)
        P  ~ Uniform(0.95, 1.05) atm  (total pressure)

    For each (x1, P) pair, we solve for T using bisection. Collecting all
    these T values lets us see how much the bubble point temperature
    varies due to uncertainty in composition and pressure.

    A fixed `seed` is used so the random numbers are reproducible --
    running the script twice gives the same results.
    """
    rng = random.Random(seed)   # seeded random number generator (reproducible)
    results = []                # will hold a dict of results per trial

    for _ in range(n_trials):
        # Randomly sample x1 and P uniformly within their given ranges
        x1 = rng.uniform(*x1_bounds)
        P = rng.uniform(*P_bounds)

        # Solve for T using bisection, given this trial's x1 and P
        T_root, n_iter = bisection(
            lambda T: bubble_point_residual(T, x1, P),  # function to root-find
            T_bounds[0], T_bounds[1],                    # search bracket [350, 400]
            tol=tol
        )

        # Store this trial's inputs and result
        results.append({"x1": x1, "P_atm": P, "T_K": T_root, "iterations": n_iter})

    return results


# ----------------------------------------------------------------------
# STEP 7: MAIN PROGRAM - RUNS EVERYTHING AND PRINTS RESULTS
# ----------------------------------------------------------------------
def main():
    # --- PART A: Nominal case (x1 = 0.4, P = 1 atm) ---------------------
    # This is the "baseline" single calculation before we add randomness.
    T_nom, n_iter_nom = bisection(
        lambda T: bubble_point_residual(T, 0.4, 1.0),  # fix x1=0.4, P=1 atm
        350, 400,     # search for T between 350K and 400K
        tol=0.005     # half-width tolerance -> guarantees error < 0.01 K
    )
    print(f"Nominal case (x1=0.4, P=1 atm): T = {T_nom:.4f} K "
          f"in {n_iter_nom} iterations")

    # --- PART B: Error analysis ------------------------------------------
    # Confirms, mathematically, how many iterations bisection NEEDS
    # to guarantee an error smaller than 0.01 K over the interval [350,400].
    n_needed = iterations_required(350, 400, tol=0.01)
    print(f"Iterations required for error < 0.01 K: {n_needed}")

    # --- PART C: Monte Carlo simulation (200 trials) ----------------------
    results = monte_carlo_bubble_point(n_trials=200, tol=0.005)

    # Filter out any trials where bisection failed to find a root (shouldn't
    # happen here, but good practice to check)
    valid = [r for r in results if r["T_K"] is not None]
    temps = [r["T_K"] for r in valid]   # extract just the T values

    # --- PART D: Summary statistics ---------------------------------------
    print(f"\nMonte Carlo simulation: {len(valid)}/{len(results)} valid trials")
    print(f"Mean T   : {statistics.mean(temps):.4f} K")   # average bubble point
    print(f"Std dev  : {statistics.pstdev(temps):.4f} K") # spread/variability
    print(f"Min / Max: {min(temps):.4f} K / {max(temps):.4f} K")  # range

    # --- PART E: Show a few sample trials for inspection -------------------
    print("\nSample trials:")
    for r in results[:5]:
        print(f"  x1={r['x1']:.4f}, P={r['P_atm']:.4f} atm "
              f"-> T={r['T_K']:.4f} K ({r['iterations']} iters)")


# ----------------------------------------------------------------------
# STEP 8: SCRIPT ENTRY POINT
# ----------------------------------------------------------------------
# This ensures main() only runs when the script is executed directly
# (not when imported as a module into another script).
if __name__ == "__main__":
    main()
