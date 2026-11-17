# Topology Optimization

**Setup.** The design region has $n_d\times n_d$ pixels; vectors over it are flattened in C order (conventions.md, item 3). Each pixel has a design variable $\rho_k\in[0,1]$: 0 means SiO₂, 1 means Si. We let $\rho$ be continuous so that $T$ has a gradient, but the final device must have only two materials. The chain is
$$\rho\ \xrightarrow{\ W\ }\ \tilde\rho\ \xrightarrow{\ P_\beta\ }\ \hat\rho\ \longrightarrow\ \varepsilon_r\ \longrightarrow\ T .$$

**Eq. 1, the filter.** $\tilde\rho = W\rho$, where
$$W_{kl} = \frac{w_{kl}}{\sum_m w_{km}},\qquad w_{kl} = \max(0,\ R - d_{kl}),$$
$d_{kl}$ is the distance between pixels $k$ and $l$ (in pixels), and $R$ is the filter radius.
- **Why filter.** Each $\tilde\rho_k$ is a weighted average over a disk of radius $R$, so no feature can be much smaller than the disk. Without the filter, the optimizer makes pixel-sized features and checkerboards: things the grid cannot resolve and a fab cannot make.
- **The edges.** There are no pixels outside the design region, so edge rows average over fewer pixels (each row is renormalized). This makes $W$ **non-symmetric**, so the chain rule needs $W^T$, not $W$.
- Density filter: Bruns & Tortorelli 2001; Bourdin 2001.

**Eq. 2, the projection.** The filter makes everything gray, so we push it back towards 0 and 1 with a smoothed step at threshold $\eta = 0.5$:
$$\hat\rho = P_\beta(\tilde\rho) = \frac{\tanh(\beta\eta)+\tanh\big(\beta(\tilde\rho-\eta)\big)}{\tanh(\beta\eta)+\tanh\big(\beta(1-\eta)\big)} .$$
Check: $P_\beta(0) = 0$ and $P_\beta(1) = 1$ for every $\beta$, and $P_\beta(0.5) = 0.5$. As $\beta\to\infty$, $P_\beta$ becomes a step at 0.5. (Wang, Lazarov & Sigmund 2011, Eq. 3.)

**Eq. 3, its derivative.** Since $\frac{d}{dx}\tanh x = 1-\tanh^2x$,
$$P_\beta'(\tilde\rho) = \frac{\beta\,\big(1-\tanh^2(\beta(\tilde\rho-\eta))\big)}{\tanh(\beta\eta)+\tanh\big(\beta(1-\eta)\big)} .$$

**Eq. 4, the material.** In the design region,
$$\varepsilon_{r,k} = \varepsilon_{min} + \hat\rho_k\,(\varepsilon_{max}-\varepsilon_{min});$$
outside it, $\varepsilon_r$ is the fixed background. **Approximation:** a gray $\hat\rho$ is not a real material. It only exists to make the problem smooth, which is why $T$ is reported both continuous and after thresholding at 0.5.

**Eq. 5, the chain rule.** Our adjoint gives $\partial T/\partial\varepsilon_{r,k}$ (adjoint.md, Eq. 5); call its design-region part $g_\varepsilon$. Each link is a small change:
$$dT = \sum_k\frac{\partial T}{\partial\varepsilon_{r,k}}\,d\varepsilon_{r,k},\qquad d\varepsilon_{r,k} = (\varepsilon_{max}-\varepsilon_{min})\,P_\beta'(\tilde\rho_k)\,d\tilde\rho_k,\qquad d\tilde\rho = W\,d\rho .$$
Put them together. With $D = \mathrm{diag}\big((\varepsilon_{max}-\varepsilon_{min})\,P_\beta'(\tilde\rho)\big)$,
$$dT = g_\varepsilon^T D\,W\,d\rho\quad\Longrightarrow\quad\frac{dT}{d\rho} = W^T D\,g_\varepsilon = W^T\Big[(\varepsilon_{max}-\varepsilon_{min})\,P_\beta'(\tilde\rho)\odot g_\varepsilon\Big].$$
This is reverse mode again (electrodynamics.md, "Cost and limits"): pull the gradient back through each link, last link first. It costs one sparse product with $W^T$, nothing more.

**β-continuation.** Start at $\beta = 1$, where $P_\beta$ is almost linear and the problem is smooth. Optimize, then double $\beta$ and continue from where you stopped, up to $\beta = 128$.
- **Why not start at β = 128.** There, $P'_\beta\approx0$ wherever $|\tilde\rho-0.5|\gg1/\beta$, so most pixels get no gradient (Eq. 5), and the optimizer freezes at its first guess.
- **Why restart the optimizer at each β.** The objective changes with β, so the curvature L-BFGS-B has stored no longer applies.

**The optimizer: L-BFGS-B** (a quasi-Newton method with box bounds). It enforces $0\le\rho\le1$ natively, it builds curvature from gradients alone, and it has no step size to tune. Its limit is that it handles only box bounds. Constraints such as a volume limit, or the worst case over several wavelengths, need **MMA** (Svanberg 1987), the field's standard.

**Measuring "binary".** Measure it on $\hat\rho$, which is what becomes $\varepsilon_r$, not on $\rho$. The fraction of pixels with $\hat\rho\le0.05$ or $\hat\rho\ge0.95$ is the binary fraction.