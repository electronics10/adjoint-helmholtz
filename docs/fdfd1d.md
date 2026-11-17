# 1D FDFD with PML

**Setup.** Nodes $x_i = ih$, $i = 0,\dots,N-1$; $u_i$ and $\varepsilon_{r,i}$ live on the nodes (conventions.md, item 3). $h$ is in any length unit and $k_0$ in its inverse (µm in photonics; the tests use wavelength 1). $s_i$ is the stretch factor at $x_i$, and $s_{i\pm1/2}$ is the stretch factor at $x_{i\pm1/2}$.

**Eq. 1, continuous equation with the stretch:**
$$\frac{1}{s}\frac{d}{dx}\Big(\frac{1}{s}\frac{du}{dx}\Big) + k_0^2\varepsilon_r u = b,\qquad s = 1 + i\alpha(x),\quad \alpha = \sigma/\omega \text{ (dimensionless)}.$$

**Eq. 2, row $i$.** Discretize each $\frac{1}{s}\frac{d}{dx}$ as a difference divided by its complex length (electrodynamics.md, "What a change of variables means on a grid"). The first differences sit on half nodes and use $s_{i\pm1/2}h$; the outer difference sits on node $i$ and uses $s_ih$:
$$\frac{1}{s_ih}\Big[\frac{u_{i+1}-u_i}{s_{i+1/2}h}-\frac{u_i-u_{i-1}}{s_{i-1/2}h}\Big]+k_0^2\varepsilon_{r,i}u_i=b_i .$$
Collect the coefficients of $u_{i-1}, u_i, u_{i+1}$:
$$\ell_i = \frac{1}{s_i\,s_{i-1/2}\,h^2},\qquad r_i = \frac{1}{s_i\,s_{i+1/2}\,h^2},\qquad d_i = k_0^2\varepsilon_{r,i}-\ell_i-r_i ,$$
$$\ell_i u_{i-1} + d_i u_i + r_i u_{i+1} = b_i,\qquad u_{-1}=u_N=0 .$$
Rows are **not** multiplied by $s_i$ (conventions.md, item 5), so $A$ is not symmetric, but $SA$ with $S=\mathrm{diag}(s_i)$ is.

**Eq. 3, profile.** The PML is the outermost $n_\text{pml}$ cells on each side, with thickness $L = n_\text{pml}h$. Let $\delta(x)\ge 0$ be the depth into the PML, capped at $L$. Then
$$\alpha(x) = \alpha_{max}\Big(\frac{\delta}{L}\Big)^m,\qquad m=3 .$$
The polynomial ramp makes $\alpha$ start at zero at the entrance, which the grid can resolve (electrodynamics.md, "What it does, in numbers": a too-steep ramp reflects).

**Eq. 4, choosing $\alpha_{max}$.** Inside the layer, $\tilde x = \int s\,dx = x + i\int\alpha\,dx$, so the outgoing wave is
$$e^{ik\tilde x} = e^{ikx}\,e^{-k\int\alpha\,dx}.$$
It crosses the layer, reflects off the mirror $u_N = 0$, and crosses back, so its amplitude returns multiplied by
$$R = \exp\Big(-2k\int_0^L\alpha\,d\delta\Big) = \exp\Big(-2k\,\alpha_{max}\frac{L}{m+1}\Big)
\quad\Longrightarrow\quad
\alpha_{max} = \frac{(m+1)\ln(1/R)}{2kL}.$$
We take $R = 10^{-8}$ and $k = k_0$. **Approximations:** this $R$ is for the continuum at normal incidence; the grid adds its own reflection, which the tests measure. If the medium next to the PML has $n>1$, the wave decays faster than designed, which is harmless.

**Eq. 5, point source:** $b_j = -1/h$ (conventions.md, item 2).

**Eq. 6, exact solution on the infinite grid.** In vacuum without a PML, multiply row $i$ by $h^2$:
$$u_{i-1} + a\,u_i + u_{i+1} = h^2 b_i,\qquad a = k_0^2h^2 - 2 .$$
Try $u_i = C\,q^{|i-j|}$.
- For $i\neq j$ the right side is 0, so $q + 1/q = -a = 2 - k_0^2h^2$. With $q = e^{iKh}$ this is $\cos Kh = 1 - \tfrac12 k_0^2h^2$, the grid's dispersion relation. Take $0<Kh<\pi$; then $q^{|i-j|}$ moves away from the source on both sides, so it is outgoing.
- For $i=j$, both neighbours equal $Cq$: $C(2q + a) = h^2\cdot(-1/h) = -h$. Since $a = -(q+1/q)$, we have $2q + a = q - 1/q = 2i\sin Kh$.

So
$$u_i = \frac{-h\,e^{iKh|i-j|}}{2i\sin Kh}.$$
**Check:** as $h\to0$, $\sin Kh\approx k_0h$, so $u_i\to\frac{-1}{2ik_0}e^{ik_0|x_i-x_j|} = \frac{i}{2k_0}e^{ik_0|x_i-x_j|}$, the continuous $G$. The tests use Eq. 6 to measure the PML alone, with no dispersion error mixed in.