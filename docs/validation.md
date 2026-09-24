# Harmonic validation of the mdins pDOS

**Status:** milestone M2, complete. Every number and figure here is produced by
`tests/test_euphonic.py` and `docs/make_validation_figures.py`, which share the harness
in `tests/lj_reference.py`.

## What is being tested, and why this way

The mdins pipeline turns MD velocities into a phonon density of states. Stage A reads the
trajectory, stage B forms the velocity cross-spectral density by Welch's method, stage C
normalises it against equipartition. An error anywhere in that chain — a factor of 2π, a
wrong dump interval, a mis-scaled window — produces a spectrum that still *looks*
plausible. Self-consistency checks cannot catch it, because the same wrong constant
appears on both sides.

So the test compares against something that shares no code with the pipeline. The same
Lennard-Jones crystal is described twice:

```
ASE LennardJones ──┬── finite displacement → force constants → Euphonic → phonons
                   └── NVE MD, low T       → mdins → pdos()
```

The upper route diagonalises the dynamical matrix and never integrates an equation of
motion. The lower route integrates the real dynamics and never diagonalises anything.
The only thing they have in common is the potential. If they agree, the pipeline is
doing what it claims.

Two conditions have to hold for the comparison to be meaningful, and both are the
test's responsibility rather than the code's:

- **The reference must use the q-grid commensurate with the MD supercell.** An *n*×*n*×*n*
  supercell under periodic boundaries can host exactly *n*³ wavevectors. Sampling the
  reference more finely adds modes the MD was never able to support, and the comparison
  then fails for a reason that has nothing to do with mdins.
- **The MD must be cold enough to stay harmonic.** The reference *is* the harmonic
  answer, so real anharmonicity is a genuine physical difference, not an error. The runs
  below average 9.8 K, roughly a fortieth of argon's melting point — but see the
  systematic discussed under [Results](#results), which is consistent with that not being
  quite cold enough.

A third condition was missed in the first version of this page and is the reason the
numbers below have changed: **one NVE run is not a sample of anything.** The normal-mode
energies of a harmonic crystal are constants of the motion, fixed once and for all by the
initial `MaxwellBoltzmannDistribution` draw. Every Welch segment of a single trajectory
therefore sees the *same* mode occupation, and the inter-segment spread measures phase
and leakage noise at fixed occupation rather than the error that matters. Measured here,
it understates the run-to-run scatter by a factor of about 4.5, and individual runs land
anywhere from 7.4 to 13.3 K. Every uncertainty on this page is now the scatter of that
statistic over ten independent runs, divided by √10. Nothing is propagated from within a
run, which also disposes of the question of how errors correlate between bins or between
atoms: there is no propagation step to get wrong.

Both crystals are relaxed to the Lennard-Jones equilibrium lattice constant
(5.264747 Å) first. At a lattice constant that is not a stationary point, the residual
forces give imaginary modes in the reference and a spurious low-energy feature in the
MD — a mismatch that looks exactly like a bug in the code under test.

## Case 1: single-species argon

27 atoms (3×3×3 of the primitive FCC cell), ten independent runs of 163.8 ps of NVE at
9.8 ± 2.0 K, Welch segments of 2048 frames giving 0.101 meV resolution.

![MD and harmonic pDOS for FCC argon](figures/single-species.png)

The 78 finite modes collapse into three tight clusters, because 27 wavevectors on a
high-symmetry lattice are heavily degenerate. The MD peaks land on them. The blue curve
is the ensemble mean, and the figure carries **two** bands on purpose: the dark inner one
is what a single run reports from its own Welch segments, and the pale outer one is the
actual spread over the ten runs. The gap between them is the correction described above.

The harmonic reference is a set of delta functions, so it is drawn convolved with the
estimator's own resolution kernel — a Gaussian of 0.15 meV, being the 1.44-bin Hann main
lobe of the 0.101 meV Welch segment spacing. Without that the reference is a spike
against a peak of finite height and the eye reads the difference as disagreement. The
broadening is a plotting choice and affects nothing in the test suite, which compares
moments and coarse bins.

## Case 2: mixed species, Ar and Kr

This is the case the single-species comparison cannot address. When every atom is
equivalent, an error in the per-atom projection averages away and leaves the total DOS
correct. A crystal with two species does not allow that.

The system is an ordered Ar/Kr arrangement on the conventional FCC cell (alternating
(001) layers, L1₀), repeated 2×2×2 for 32 atoms. The choice is deliberate: **ASE's
Lennard-Jones takes a single epsilon and sigma, so it cannot tell argon from krypton.**
The two species feel identical forces and sit on a geometrically identical lattice. The
mass ratio of 2.10 is the only asymmetry in the entire problem, so every difference
between the two spectra has exactly one cause, and a projection bug cannot hide behind a
difference in the forces.

![Species-projected pDOS for ordered Ar/Kr](figures/mixed-species.png)

The lighter argon carries the top of the band and krypton the bottom, with mean energies
of 5.60 and 3.92 meV across the ensemble. Dashed lines are the Euphonic reference,
projected onto species through the mode eigenvectors and broadened as above; solid lines
are the MD. The agreement is peak by peak, not merely in the envelope.

Three details in this figure are worth naming:

- **The reference is projected per species from the eigenvectors, not taken from
  `calculate_pdos`.** The three acoustic modes at Γ are rigid translation of the whole
  cell, which `remove_com_velocity` strips out of the MD. Leaving them in the reference
  would compare a spectrum that has them against one that does not.
- **The acoustic weight is not shared equally.** Γ acoustic eigenvectors go as √*m*, so
  krypton holds 4.2% of its weight there against argon's 2.0% — a ratio of 2.10, the
  mass ratio. Masking the same fraction from both species would bias the two comparisons
  in opposite directions.
- **The ratio of integrated weights between species needs the exact centre-of-mass
  correction.** Equipartition puts 3*k*<sub>B</sub>*T*/*m*<sub>i</sub> under each atom's
  spectrum, but after COM projection the correct target is
  3*k*<sub>B</sub>*T*(1/*m*<sub>i</sub> − 1/*M*). Using the uncorrected form makes the
  Ar/Kr weight ratio 3.4% wrong instead of 1.1%.

## Results

Tolerances are **measured, not tuned**. Each statistic is computed separately for ten
independent NVE runs, and the quoted 1σ is the standard error of their mean. The test
threshold is 3σ throughout. Mixed-species moments are taken over the harmonic band (see
[Known limitations](#known-limitations)); single-species moments over the full range.

| Comparison | MD (mean of 10 runs) | Euphonic | Difference | 1σ | Deviation |
|---|---|---|---|---|---|
| Ar only, ⟨E⟩ (meV) | 5.7124 ± 0.0785 | 5.5852 | +2.28% | 1.37% | +1.6σ |
| Ar only, ⟨E²⟩ (meV²) | 35.337 ± 0.858 | 33.592 | +5.20% | 2.43% | +2.0σ |
| Ar/Kr, Ar ⟨E⟩ (meV) | 5.6037 ± 0.0413 | 5.5975 | +0.11% | 0.74% | +0.1σ |
| Ar/Kr, Ar ⟨E²⟩ (meV²) | 33.286 ± 0.461 | 33.014 | +0.82% | 1.39% | +0.6σ |
| Ar/Kr, Kr ⟨E⟩ (meV) | 3.9165 ± 0.0354 | 3.8407 | +1.97% | 0.90% | +2.1σ |
| Ar/Kr, Kr ⟨E²⟩ (meV²) | 16.766 ± 0.296 | 16.102 | +4.12% | 1.76% | +2.2σ |

![Every comparison in units of its ensemble uncertainty](figures/residuals.png)

Everything passes at 3σ, but **every deviation has the same sign**, which no amount of
random scatter explains and which the previous single-seed version of this table could
not have shown. The MD spectrum sits systematically a little above the harmonic
reference: about +2% in ⟨E⟩ and, because the second moment weights by energy squared,
roughly double that in ⟨E²⟩. That ratio is what a small excess of weight at the top of
the band produces, and about 1% of the argon weight does lie above the harmonic band top
(Welch leakage and residual anharmonicity, bounded by its own test). The residual
anharmonicity at 9.8 K is the leading suspect for the rest. The honest summary is that
the agreement is good to a couple of per cent with a systematic of known sign and
unproven cause, not that the two routes agree to within noise.

A measured tolerance is only worth anything if it is tight enough to reject a wrong
answer, so the suite includes a negative control. At 3σ the single-species first-moment
test rejects a frequency-axis error above roughly 3.5% at ten seeds — much weaker than
the 1% the within-run Welch spread would have claimed, but still far smaller than a
factor of 2π, of 2, or of any plausible unit conversion. `test_the_comparison_can_fail`
asserts that a 5% stretch is rejected.

## Known limitations

**The second moment is dominated by the band tail.** About 1% of the argon weight and
0.3% of the krypton weight lies above the harmonic band top, from Welch leakage and
residual anharmonicity. Because ⟨E²⟩ weights by energy squared, that small tail carries
a disproportionate share, and the full-range mixed-species second moments disagree by
several times as much as the band-restricted ones in the table. The moments quoted above
for the mixed-species case are therefore taken over the band, cut at the top plus five
resolution widths. This is a choice of statistic, not a filter on the result — the
excluded weight is bounded by a separate test
(`test_little_weight_lies_outside_the_harmonic_band`), so it cannot grow unnoticed.

**The mass splitting is low by more than either moment is.** Ar − Kr comes out at 1.687
meV from the table above against a harmonic 1.757, −4.0%, because it is a difference of
two numbers whose errors are anti-correlated: a draw that puts extra weight at the top of
the argon band tends to do the same to krypton, so the two move together and the gap
between them moves twice as far. `test_the_species_are_separated_by_the_amount_euphonic_predicts`
compares it against the ensemble scatter of the splitting itself rather than against a
fixed percentage, which is the only fair comparison for a difference of correlated
quantities.

**These are small, cold, ordered crystals.** Nothing here tests a liquid, a molecular
system, or anything anharmonic — by construction, since the reference is harmonic. The
molecular-dynamics tests in `tests/test_md.py` cover the exactly-solvable Einstein
crystal; agreement with a real experiment is milestone M5.

**Euphonic's C extension was unavailable in the environment that produced these
figures**, so the reference fell back to the pure-Python path. That affects speed only,
not the result.

**Ten seeds is not many.** The error bars above are themselves estimated from ten
numbers, so they are good to about 24% each, and a 2σ result at this ensemble size should
not be read as a firm detection of anything. The systematic offset is more convincing
than any individual row precisely because it appears in all six.

**Running for longer will not fix the offset.** Because the mode occupations are frozen
by the initial draw, a longer trajectory converges each run more tightly onto its own
answer without moving the ensemble mean. Extending the runs narrows the inter-segment
spread and leaves the run-to-run scatter essentially unchanged — measured, 4.23% at 41 ps
against 4.35% at 163.8 ps. More seeds, not more picoseconds, is what buys precision here.

`docs/interactive-validation.ipynb` is the interactive counterpart to this page: the same
comparison with the trajectory length, Welch segment length, window, overlap, temperature
and NVE seed under sliders, so the sensitivity of each number here can be read off
directly. It is worth opening for the seed control alone — ⟨E⟩ spans 5.31 to 6.24 meV
across the ten runs, and one of them sits 13σ from the harmonic reference by its own
quoted Welch error bar.

## Reproducing

```bash
uv sync --extra euphonic --group docs
uv run pytest tests/test_euphonic.py -m slow      # 22 tests, about 1 minute
uv run python docs/make_validation_figures.py --refresh
uv run python docs/make_interactive_data.py --refresh   # ~6 min, for the notebook
```

The notebook is committed with its outputs, so it reads without a kernel; the cache is
only needed to re-execute it.
