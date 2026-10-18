# ES-QPID: Enhanced Smith Predictor-Based Quantum-Inspired PID for Antenna Azimuth Control

Simulation code for the article

> P. Bharat Kumar, K. Srinivasa Rao, S. Kannappan, G. Jyothi, R. Lakshmi,
> *Smith Predictor-Based Quantum-Inspired PID Controller with Hyperbolic Tangent Nonlinearity for Robust Azimuth Position Control of Earth Station Antennas under Variable Time Delay*,
> (article under review; the full reference will be added after publication).

The code reproduces every table and figure of the article: nominal comparison, robust design cost and ablation, ±30 % parameter errors, delay margins, step sizes, actuator limits, 200 Monte Carlo trials per delay, and the Popov / small-gain stability analysis.

## The controller in brief

The plant is an antenna azimuth drive, G(s) = 6.63 / (s² + 14 s + 6.63), behind a communication delay of L = 0.3, 0.8 or 1.5 s. The ES-QPID controller combines:

1. **Filtered Smith predictor.** The predictor is closed around the measured antenna angle. A lead-lag mismatch filter cancels the slow plant pole in the response to wind disturbances.
2. **Quantum-inspired PID law.** The pointing error is mapped to a qubit angle θₑ = π·sat(e/e_max).
3. **Two-qubit probability weights.** P_reg = cos²(θ_r/2)·cos²(θₑ/2) weights the integral, and 1 − P_reg weights the derivative. Here θ_r is driven by the set-point transition only.
4. **Smooth bounded correction.** The correction is applied through U_q·tanh(α·θ_u), and the command is limited to ±20 V.

The baselines are PID, Smith-PID, filtered Smith-PID, QPID and Smith-QPID, plus three ablations. All of them are tuned with the same robust cost and have the same actuator limit.

## Requirements

| Use | Software | Notes |
|---|---|---|
| Run everything (recommended) | **MATLAB R2016b or later** | No toolbox needed. Folder `matlab/` |
| Run everything (free) | **GNU Octave 8.x** with the `control` package | Folder `octave/`. Install the package with `pkg install -forge control` |
| Optional: global gain search, redraw manuscript figures | **Python ≥ 3.9** with `numpy`, `scipy`, `numba`, `matplotlib` | `pip install -r requirements.txt` |

Tested with MATLAB-compatible code in GNU Octave 8.4.0 (control 4.0.1), and with Python 3.13 (numpy 2.5, scipy 1.18, numba 0.68, matplotlib 3.11). A full run takes about 10–15 minutes. Most of that time is the Monte Carlo study.

## Quick start

**MATLAB**
```matlab
cd matlab
main_option2                 % or: ES_QPID_Option2_complete   (single-file version)
```

**GNU Octave**
```octave
cd octave
pkg load control
main_option2
```

Both write `results_option2.txt` (all tables), `results_option2.mat` (all data) and the folder `figures_option2/` with the PNG figures.

**Single runs** (MATLAB or Octave):
```matlab
[m, o] = sqpid_run('ES-QPID', 0.8);                    % nominal, L = 0.8 s
[m, o] = sqpid_run('ES-QPID', 0.8, 'Lerr', 1.2);       % true delay 20 % longer than the model
[m, o] = sqpid_run('PID', 1.5, 'kg', 0.8, 'gamp', 1);  % plant gain -20 %, double wind gust
plot(o.t, o.y); disp(m.IAE)
```
The run options are `R`, `T`, `Lerr`, `kg`, `delta`, `gamp`, `gsig`, `t_on`, `t_off`, `noise`, `seed`, `Usat`, `emax` and `gains` (see `sqpid_run.m`). The controllers are `PID`, `S-PID`, `FS-PID`, `QPID`, `S-QPID`, `ES-QPID`, `ES-noweight`, `ES-nofilt` and `ES-noDw`. Tuned gains exist for L = 0.3, 0.8 and 1.5 s (`sqpid_ctrl.m`).

**Python (optional)**
```bash
pip install -r requirements.txt
cd python/figures && python make_figs.py 2 ../../results out      # manuscript-style figures
cd python/tuning  && python tune_global.py eq final_eq.json        # global gain search (several hours)
```

## Repository layout

```
matlab/        MATLAB version (base MATLAB only)
  main_option2.m              runs all experiments and figures
  ES_QPID_Option2_complete.m  the same, as one self-contained file
  sqpid_sim.m                 closed-loop simulation (true delay, Smith predictor, PID / QPID laws)
  sqpid_run.m                 one experiment with options
  sqpid_ctrl.m                controller definitions and tuned gains (Table 4)
  sqpid_cost.m                robust design cost (7 scenarios)
  sqpid_experiments.m         experiments E1-E7
  sqpid_stability.m, sqpid_e7.m   Popov test, margins, small-gain delay bound
  sqpid_metrics.m, sqpid_noise.m, sqpid_zoh.m, sqpid_figures.m, tune_gains.m
  sqpid_random.mat            exact random sequences of the paper (keep next to the code)
octave/        the original GNU Octave version (uses the control package)
python/tuning  Python/numba port of the simulation loop and the global gain search
python/figures make_figs.py - draws the manuscript figures from results/results_option2.mat
results/       results_option2.txt / .mat - the results reported in the article
figures/       Fig01-Fig12 as they appear in the article
```

## Mapping to the article

| Experiment | Article |
|---|---|
| E1 nominal runs | Tables 6–7, Figs. 4, 5, 7, 8, 9 |
| Robust design cost J and ablation | Table 8, Fig. 6 |
| E2 ±30 % gain, damping, delay errors | Tables 9–10, Fig. 10 |
| E3 delay margins | Table 11 |
| E4 step sizes / e_max | Table 12, Fig. 11 |
| E5 actuator limits | Table 13 |
| E6 Monte Carlo (200 runs per delay) | Table 14, Fig. 12 |
| E7 stability analysis | Table 3, Fig. 3 |

## Reproducibility

- MATLAB and Octave produce different random numbers. The file `matlab/sqpid_random.mat` therefore stores the exact sequences used in the article: measurement noise, turbulence and the Monte Carlo draws. With this file, the MATLAB version reproduces `results/results_option2.txt` to the last printed digit. Without it, MATLAB's `rng` is used and the last digits change slightly.
- The gains in `sqpid_ctrl.m` come from the Python global search. They were rounded to 5 significant figures and are stored in `python/tuning/locked.json`. The Python port gives the same design cost as the MATLAB/Octave code to three decimals.
- The article figures were drawn with `python/figures/make_figs.py`. `sqpid_figures.m` plots the same data with a slightly different layout.

## Citation

If you use this code, please cite it (see `CITATION.cff`). The journal reference and DOI of the article will be added after publication.

## License

No license has been chosen yet. Until a license is added, the default copyright applies and others may view the code but not reuse it.
