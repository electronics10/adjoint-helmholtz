# Adjoint gradient

**Setup.** $A(\varepsilon_r)\,u = b$ with $A$ of shape $N\times N$. The objective is $J = |a|^2$, where $a = c^Hu$ is a complex scalar and $c$ is an $N$-vector of probe weights (conventions.md, item 4). For example, $c = e_m$ reads $u_m$.

**Eq. 1.** Only the diagonal of $A$ depends on $\varepsilon_{r,i}$, through $k_0^2\varepsilon_{r,i}$ (fdfd1d.md, Eq. 2; the PML terms do not depend on $\varepsilon_r$):
$$\frac{\partial A}{\partial\varepsilon_{r,i}} = k_0^2\,e_ie_i^T .$$

**Eq. 2.** $J$ is real, so (conventions.md, item 6)
$$dJ = 2\,\mathrm{Re}\Big(\sum_j g_j\,du_j\Big) = 2\,\mathrm{Re}(g^T du),\qquad g_j = \frac{\partial J}{\partial u_j}\ (\bar u \text{ held fixed}).$$

**Eq. 3.** $J = a\bar a$ with $a = \sum_j\bar c_ju_j$ and $\bar a = \sum_j c_j\bar u_j$. Holding $\bar u$ fixed, only $a$ changes, so $\partial J/\partial u_j = \bar a\,\bar c_j$:
$$g = \bar a\,\bar c .$$

**Change of $u$.** Differentiate $Au = b$; $b$ does not depend on $\varepsilon_r$:
$$dA\,u + A\,du = 0\quad\Longrightarrow\quad du = -A^{-1}dA\,u .$$

**Eq. 4, the adjoint field.** Put $du$ into Eq. 2:
$$g^Tdu = -g^TA^{-1}dA\,u = -(A^{-T}g)^T dA\,u .$$
Define $\lambda$ by
$$A^T\lambda = g .$$
Then $g^Tdu = -\lambda^T dA\,u$. One solve gives $\lambda$, whatever $dA$ is.

**Eq. 5, the gradient.** With $dA = k_0^2e_ie_i^T\,d\varepsilon_{r,i}$ (Eq. 1), $\lambda^TdA\,u = k_0^2\lambda_iu_i\,d\varepsilon_{r,i}$, so
$$\frac{dJ}{d\varepsilon_{r,i}} = -2k_0^2\,\mathrm{Re}(\lambda_i u_i).$$
Cost: one LU factorization, then one forward and one transposed solve, for all $N$ gradient entries.

**Physical picture.** $\lambda$ is the field of a source $g$ placed at the probe, sent back into the device. $\lambda_iu_i$ is large where the forward wave and the wave from the probe overlap.

**What goes wrong with the wrong transpose (measured).**
- **$A$ instead of $A^T$.** The gradient is wrong *only in the PML*. Write $A^{-1} = (SA)^{-1}S$ with $SA$ symmetric. Then $A^{-T} = S(SA)^{-1}$, so $\lambda_{A^T} = S\,\lambda_A$ whenever $Sg = g$ (the probe is outside the PML). $S = 1$ outside the PML, so the two agree there. Lesson: the finite-difference check must sample some PML nodes, or this bug passes. In 2D, a probe line that reaches into the PML makes the bug show up in the physical region too.
- **$A^H$ instead of $A^T$.** The gradient is wrong everywhere (110% error in our test). Conjugating flips the PML from absorbing to amplifying, and the solve is global.

**Finite differences.** Use central differences:
$$\frac{J(\varepsilon_{r,i}+\Delta)-J(\varepsilon_{r,i}-\Delta)}{2\Delta}.$$
The truncation error is $O(\Delta^2)$, against $O(\Delta)$ for forward differences. The V plot ![Finite difference vs adjoint](figures/fd_check_1d.png) shows central differences reaching $10^{-10}$ at $\Delta\approx10^{-4}$, where forward differences only reach $10^{-7}$. The test uses $\Delta = 10^{-4}$ and a tolerance of $10^{-6}$ relative to $\max_i|dJ/d\varepsilon_{r,i}|$.