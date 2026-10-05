# Statistical derivation notes

## 1. Conditional lepton law

Baseline generator:

[
(ell_1,ell_2,ell_3)stackrel{iid}{sim}mathrm{LogUniform}(L,H).
]

Condition on

[
E_K:quad |Q(ell_1,ell_2,ell_3)-2/3|<epsilon.
]

For fixed (ell_1,ell_2), write (y=sqrt{ell_3}),
(A=sqrt{ell_1}+sqrt{ell_2}), and (B=ell_1+ell_2). The boundary
(Q=q) satisfies

[
(1-q)y^2-2qAy+(B-qA^2)=0.
]

Thus the finite Koide window can be represented exactly as a union of
intervals in (y). Since (ell_3=y^2) and the original law is uniform in
(logell_3), an allowed interval ([y_a,y_b]) has probability

[
rac{2log(y_b/y_a)}{log(H/L)}.
]

The code partitions the physical (y)-range at all four roots from
(q=2/3pmepsilon), tests interval midpoints, and sums the allowed log-width.
This gives (w(ell_1,ell_2)=P(E_K|ell_1,ell_2)).

Sampling (ell_3) uniformly in log mass from the allowed interval union and
weighting by (w) yields

[
E[f|E_K]=rac{E_{\ell_1,\ell_2}[w,E(f|\ell_1,\ell_2,E_K)]}
{E_{\ell_1,\ell_2}[w]}.
]

No exact-surface parameterization or surface-measure assumption is required.

## 2. Cell-oriented interval constraints

Every observed threshold is defined from the manuscript convention

[
r_{m cell}=V/T-1.
]

Examples:

[
|F^2/m_s-1|le t_1
Rightarrow
rac{F^2}{1+t_1}le m_slerac{F^2}{1-t_1},
]

[
|mu_*m_d/m_s^2-1|le t_2
Rightarrow
sqrt{rac{mu_*m_d}{1+t_2}}le m_sle
sqrt{rac{mu_*m_d}{1-t_2}}.
]

Their intersection is the exact shared-(m_s) window.

## 3. True sixth-slot constraint

Slots 4 and 5 define a rectangle in

[
x=log(m_c/c_0),qquad y=log(m_b/b_0),
]

where (c_0=3alpha_Kmu_*), (b_0=mu_*/(2alpha_K)), and therefore
(c_0b_0=G^2=3mu_*^2/2).

The sixth condition is

[
log(1-t_6)le x+ylelog(1+t_6).
]

The exact allowed probability is the area of the slot-4/5 rectangle
intersected by this diagonal strip, divided by the square of the quark
log-span. For a rectangle ([a,b]	imes[c,d]), the area below (x+y=z) is

[
A(z)=rac12[(z-a-c)_+^2-(z-b-c)_+^2-(z-a-d)_+^2+(z-b-d)_+^2].
]

Hence the sixth-slot allowed area is (A(z_{m hi})-A(z_{m lo})).

This is why the product slot costs a finite factor (~0.68 in the current
corrected replay) rather than being an automatic consequence of the two
one-mass slots.
