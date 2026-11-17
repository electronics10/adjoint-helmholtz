# Modes

**Setup.** A straight guide along $x$: $\varepsilon_{r,ij} = \varepsilon_j$ for every $i$. We look for a source-free solution of $Au = 0$ of the form
$$u_{ij} = e_j\,q^i ,$$
with $e$ an $n_y$-vector (the profile) and $q$ a complex number (the step factor along $x$).

**Eq. 1.** Away from the $x$-PML, $D_{xx}$ is the plain stencil $(1, -2, 1)/h^2$. Write $M = D_{yy} + k_0^2\,\mathrm{diag}(\varepsilon)$, of shape $n_y\times n_y$. Row $(i,j)$ of $Au = 0$ is
$$e_j\,\frac{q^{i-1}-2q^i+q^{i+1}}{h^2} + q^i\,(Me)_j = 0 .$$
Divide by $q^i$:
$$(Me)_j = -\frac{q-2+1/q}{h^2}\,e_j .$$

**Eq. 2, the transverse eigenproblem.** Define $\beta^2 = -(q-2+1/q)/h^2$. Then
$$Me = \beta^2 e .$$
This is a 1D eigenproblem in $y$ alone. $M$ is built from the *same* $D_{yy}$ as $A$, PML rows included, so the profile fits the 2D operator exactly.
- **Guided modes** have $k_0^2\varepsilon_\text{clad} < \beta^2 < k_0^2\varepsilon_\text{core}$ and an almost real $\beta^2$, because they barely reach the PML.
- **Spurious modes.** The PML adds eigenvectors that live inside it and have complex $\beta^2$. We keep only $|\mathrm{Im}\,\beta^2| < 10^{-6}\,|\mathrm{Re}\,\beta^2|$. (The prototype found one ranked above a real guided mode.)
- **Normalization.** An eigenvector has an arbitrary scale and phase. We fix them with $\sum_j|e_j|^2h = 1$ and a real positive peak.

**Eq. 3, the discrete wavenumber.** With $q = e^{i\beta_dh}$, we have $q-2+1/q = 2\cos\beta_dh - 2$, so
$$\frac{2-2\cos\beta_dh}{h^2} = \beta^2\quad\Longrightarrow\quad \beta_d = \frac{1}{h}\arccos\Big(1-\frac{\beta^2h^2}{2}\Big).$$
This is the same grid dispersion relation as fdfd1d.md, Eq. 6, now for the guide. $\beta_d\approx\beta(1+\beta^2h^2/24)$. Use $\beta_d$, not $\beta$: only then is $u = e\,e^{i\beta_dx}$ an *exact* solution of the discrete equations, which Eq. 5 needs.

**Eq. 4, the incident field.** The mode crosses grid line $i_0$ and travels in direction $d = \pm1$:
$$u^\text{inc}_{ij} = e_j\,\exp\big(i\,d\,\beta_d\,(i-i_0)h\big).$$

**Eq. 5, the one-way source (total field / scattered field).** We want the field to be $u^\text{inc}$ on the side it travels into and zero behind the source. Let $Q$ be the diagonal 0/1 mask of that side (line $i_0$ included). Then
$$A(Qu^\text{inc}) = (AQ-QA)\,u^\text{inc} + Q\,Au^\text{inc}.$$
Outside the $x$-PML, $Au^\text{inc} = 0$ (Eqs. 1 to 3). So with
$$b = (AQ-QA)\,u^\text{inc},$$
the solution of $Au = b$ is $u = Qu^\text{inc}$: the mode on one side and nothing on the other.
- $AQ - QA$ is zero wherever the stencil does not straddle the edge of $Q$. So $b$ lives on just two lines, $i_0-1$ and $i_0$.
- **Physical picture:** the two lines form a source sheet that cancels the mode's field behind it and creates it ahead of it.
- **Approximation:** inside the $x$-PML, $u^\text{inc}$ is not a solution, so the true field there differs from $Qu^\text{inc}$. What comes back from that difference is the PML's reflection (measured at $4\times10^{-7}$).

**Eq. 6, the monitor.** To read the amplitude $a$ of the profile $e$ on a line, put
$$c = \frac{e}{e^He}$$
on that line and zero elsewhere. Then for $u = a\,e$ on the line, $c^Hu = a\,\frac{e^He}{e^He} = a$. The source launches $a = 1$, so the **transmission** is $T = |c^Hu|^2$.
- **Approximation:** other modes on the same line add $e^He_n/(e^He)$. That is zero for modes of a symmetric guide with opposite parity, and close to zero for any two guided modes, whose profiles are real apart from their tails in the PML.

**Staircasing (approximation).** $\varepsilon_r$ is sampled at the nodes, so the width of the core on the grid is a whole number of cells and changes with $h$. The measured $n_\text{eff}$ of a 0.4 µm Si guide is 3.198 at $h = 0.04$, 3.195 at $0.025$ and 3.214 at $0.02$ µm. This does not affect the tests, which compare the solver against itself, but it sets the accuracy of any absolute number.