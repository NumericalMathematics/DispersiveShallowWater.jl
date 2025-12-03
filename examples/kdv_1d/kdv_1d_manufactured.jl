using OrdinaryDiffEqRosenbrock
using DispersiveShallowWater
using SummationByPartsOperators: upwind_operators, periodic_derivative_operator

###############################################################################
# Semidiscretization of the KdV equation 

equations = KdVEquation1D(gravity = 9.81, D = 1.0)

initial_condition = initial_condition_manufactured
source_terms = source_terms_manufactured
boundary_conditions = boundary_condition_periodic

# create homogeneous mesh
coordinates_min = -1.0
coordinates_max = 1.0
N = 512
mesh = Mesh1D(coordinates_min, coordinates_max, N)

# Create solver with periodic SBP operators of accuracy order 3,
# which results in a 4th order accurate semi discretizations.
# We can set the accuracy order of the upwind operators to 3 since
# we only use central versions/combinations of the upwind operators.
D1_upwind = upwind_operators(periodic_derivative_operator;
                             derivative_order = 1, accuracy_order = 3,
                             xmin = xmin(mesh), xmax = xmax(mesh),
                             N = nnodes(mesh))
solver = Solver(D1_upwind)

semi = Semidiscretization(mesh, equations, initial_condition, solver,
                          boundary_conditions = boundary_conditions,
                          source_terms = source_terms)

tspan = (0.0, 1.0)
ode = semidiscretize(semi, tspan, split_ode = Val{false}())

summary_callback = SummaryCallback()
analysis_callback = AnalysisCallback(semi; interval = 100,
                                     extra_analysis_errors = (:conservation_error,),
                                     extra_analysis_integrals = (waterheight_total,
                                                                 waterheight, entropy))
callbacks = CallbackSet(analysis_callback, summary_callback)
saveat = range(tspan..., length = 100)

# The problem is very stiff. `Tsit5()` will cause a `maxiter` exceeded error.
alg = Rodas5()
sol = solve(ode, alg, abstol = 1e-12, reltol = 1e-12,
            save_everystep = false, callback = callbacks, saveat = saveat)

""" 
using alg = KenCarp4() I get the following Benchmarks: 
For some reason IMEX seems to perform worse - doing more steps
no IMEX: (split_ode = Val{false}())
──────────────────────────────────────────────────────────────────────────────────────
           DispersiveSWE                     Time                    Allocations
                                    ───────────────────────   ────────────────────────
         Tot / % measured:               611ms /  34.2%           40.0MiB /  10.0%

Section                     ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────────────────────
rhs!                         6.62k    196ms   93.6%  29.6μs   3.94MiB   98.4%     624B
  source terms               6.62k    108ms   51.5%  16.3μs   3.94MiB   98.4%     624B
  third-order derivatives    6.62k   44.1ms   21.1%  6.66μs     0.00B    0.0%    0.00B
  hyperbolic                 6.62k   41.1ms   19.7%  6.22μs     0.00B    0.0%    0.00B
  ~rhs!~                     6.62k   2.97ms    1.4%   450ns   1.25KiB    0.0%    0.19B
analyze solution                 3   13.3ms    6.4%  4.45ms   64.8KiB    1.6%  21.6KiB
──────────────────────────────────────────────────────────────────────────────────────

IMEX with sources in nonstiff part:
──────────────────────────────────────────────────────────────────────────────────────
           DispersiveSWE                     Time                    Allocations
                                    ───────────────────────   ────────────────────────
         Tot / % measured:               2.85s /  12.3%           98.6MiB /   5.7%

Section                     ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────────────────────
rhs_split_nonstiff!          8.93k    169ms   48.3%  18.9μs   5.31MiB   93.8%     624B
  source terms               8.93k    160ms   45.6%  17.9μs   5.31MiB   93.8%     624B
  hyperbolic                 8.93k   5.67ms    1.6%   635ns     0.00B    0.0%    0.00B
  ~rhs_split_nonstiff!~      8.93k   3.44ms    1.0%   386ns      976B    0.0%    0.11B
rhs_split_stiff!             45.4k    103ms   29.4%  2.27μs      672B    0.0%    0.01B
  third-order derivatives    45.4k   95.9ms   27.4%  2.11μs     0.00B    0.0%    0.00B
  ~rhs_split_stiff!~         45.4k   7.07ms    2.0%   156ns      672B    0.0%    0.01B
analyze solution                15   78.2ms   22.3%  5.21ms    360KiB    6.2%  24.0KiB
──────────────────────────────────────────────────────────────────────────────────────

IMEX with sources in stiff part:
──────────────────────────────────────────────────────────────────────────────────────
           DispersiveSWE                     Time                    Allocations
                                    ───────────────────────   ────────────────────────
         Tot / % measured:               4.89s /  28.2%           88.7MiB /  49.0%

Section                     ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────────────────────
rhs_split_stiff!             72.2k    1.27s   92.0%  17.6μs   43.0MiB   98.9%     624B
  source terms               72.2k    1.17s   85.2%  16.3μs   43.0MiB   98.9%     624B
  third-order derivatives    72.2k   70.8ms    5.1%   981ns     0.00B    0.0%    0.00B
  ~rhs_split_stiff!~         72.2k   22.1ms    1.6%   306ns      976B    0.0%    0.01B
analyze solution                23    102ms    7.4%  4.42ms    501KiB    1.1%  21.8KiB
rhs_split_nonstiff!          13.3k   8.86ms    0.6%   665ns      672B    0.0%    0.05B
  hyperbolic                 13.3k   6.82ms    0.5%   512ns     0.00B    0.0%    0.00B
  ~rhs_split_nonstiff!~      13.3k   2.03ms    0.1%   153ns      672B    0.0%    0.05B
──────────────────────────────────────────────────────────────────────────────────────






No with fixed time step dt = 1e-2:



no IMEX
──────────────────────────────────────────────────────────────────────────────────────
           DispersiveSWE                     Time                    Allocations
                                    ───────────────────────   ────────────────────────
         Tot / % measured:               6.33s /  57.1%           0.99GiB /   1.3%

Section                     ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────────────────────
rhs!                         22.0k    3.46s   96.0%   157μs   13.1MiB   99.8%     624B
  third-order derivatives    22.0k    1.61s   44.6%  73.1μs     0.00B    0.0%    0.00B
  hyperbolic                 22.0k    1.49s   41.2%  67.5μs     0.00B    0.0%    0.00B
  source terms               22.0k    355ms    9.8%  16.1μs   13.1MiB   99.8%     624B
  ~rhs!~                     22.0k   13.5ms    0.4%   615ns   1.25KiB    0.0%    0.06B
analyze solution                 1    146ms    4.0%   146ms   21.9KiB    0.2%  21.9KiB
──────────────────────────────────────────────────────────────────────────────────────

IMEX with sources in nonstiff part:
──────────────────────────────────────────────────────────────────────────────────────
           DispersiveSWE                     Time                    Allocations
                                    ───────────────────────   ────────────────────────
         Tot / % measured:               4.34s /  39.4%           0.98GiB /   0.0%

Section                     ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────────────────────
rhs_split_stiff!             22.2k    1.70s   99.0%  76.3μs      672B    0.2%    0.03B
  third-order derivatives    22.2k    1.69s   98.7%  76.1μs     0.00B    0.0%    0.00B
  ~rhs_split_stiff!~         22.2k   5.55ms    0.3%   250ns      672B    0.2%    0.03B
rhs_split_nonstiff!            601   15.1ms    0.9%  25.0μs    367KiB   89.4%     626B
  source terms                 601   12.9ms    0.8%  21.5μs    366KiB   89.2%     624B
  hyperbolic                   601   1.41ms    0.1%  2.35μs     0.00B    0.0%    0.00B
  ~rhs_split_nonstiff!~        601    713μs    0.0%  1.19μs      976B    0.2%    1.62B
analyze solution                 1   1.86ms    0.1%  1.86ms   42.7KiB   10.4%  42.7KiB
──────────────────────────────────────────────────────────────────────────────────────

IMEX with sources in stiff part:
──────────────────────────────────────────────────────────────────────────────────────
           DispersiveSWE                     Time                    Allocations
                                    ───────────────────────   ────────────────────────
         Tot / % measured:               4.29s /  48.3%           0.99GiB /   1.3%

Section                     ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────────────────────
rhs_split_stiff!             22.2k    2.06s   99.7%  92.9μs   13.2MiB   99.8%     624B
  third-order derivatives    22.2k    1.70s   82.1%  76.4μs     0.00B    0.0%    0.00B
  source terms               22.2k    356ms   17.2%  16.0μs   13.2MiB   99.8%     624B
  ~rhs_split_stiff!~         22.2k   8.90ms    0.4%   401ns      976B    0.0%    0.04B
analyze solution                 1   3.62ms    0.2%  3.62ms   21.9KiB    0.2%  21.9KiB
rhs_split_nonstiff!            601   1.61ms    0.1%  2.68μs      672B    0.0%    1.12B
  hyperbolic                   601   1.16ms    0.1%  1.93μs     0.00B    0.0%    0.00B
  ~rhs_split_nonstiff!~        601    448μs    0.0%   746ns      672B    0.0%    1.12B
──────────────────────────────────────────────────────────────────────────────────────


Also for dt = 1e-1 the solution for IMEX already looks bad(ish),
will for split_ode = Val{false}() it looks still good.
"""
