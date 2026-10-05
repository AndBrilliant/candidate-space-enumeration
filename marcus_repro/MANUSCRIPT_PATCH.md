# Suggested manuscript patch for Sec. VII B

This file is wording guidance, not an automated edit of the manuscript.

## Replace "six-slot event" around the historical five-factor calculation

Preferred structure:

> The conditional replay freezes five primary interval constraints and one
> derived product constraint. We report the five-primary joint and, separately,
> the six-slot joint obtained by additionally imposing the observed
> (m_c m_b/G^2) product window. The latter is evaluated by integrating the
> allowed (c,b) rectangle intersected with the product strip; it is not treated
> as an automatic consequence of the two one-mass constraints.

## Clarify the conditioning measure

Preferred structure:

> The baseline conditional ensemble is the original three-lepton log-uniform
> generator conditioned on the finite observed Koide window
> (|Q_ell-2/3|<2.2	imes10^{-6}). The conditional average is evaluated without
> replacing this measure by a uniform parameterization of the exact Koide
> surface. Exact-surface samplers are retained only as null-shape sensitivity
> variants.

## Clarify P_proc

Add one explicit sentence:

> (P_{m proc}) is neither a Standard-Model p-value nor a global
> look-elsewhere significance for the historical discovery process; it is only
> the repeated-sample frequency of the specified event for the frozen
> algorithm-generator pair ((A,G)).

## Grammar-size robustness

Replace categorical wording such as

> "no plausible discrete enlargement rescues or destroys any conclusion"

with

> "Within the nested anchor enlargements tested here, the accidental-cell rate
> grows approximately linearly with the number of declared comparisons."

## Continuous-parameter rhetoric

Do not identify (ln(10)/0.001) with an "effective" trials factor. Adjacent
grid points are correlated. If retained, describe 2300 only as the number of
0.1%-spaced grid locations in one logarithmic decade, not as 2300 independent
trials.

## Repository statement

Only say that "all seeds, code, outputs and validation paths are archived" once
the cited release/commit actually contains every driver used for the values in
the manuscript.
