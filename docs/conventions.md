# Conventions

This file is the single reference for the code. `electrodynamics.md` and `convolution.md` are background notes written for learning; where their signs or scalings differ (e.g. they use a $+\delta$ source), this file wins.


1. Time dependence: $e^{-i\omega t}$, so $e^{ikx}$ moves right, and $\text{Im}\varepsilon_r>0$ means loss.
2. Operator and source: $A=\nabla^2 + k_0 ^2\varepsilon_r$. A physical current enters as $b=-i\omega \mu_0 J$. A unit point source is $b=-\delta$, so the Green's functions are $\frac{i}{2k}e^{ik|x|}$ in 1D and $\frac{i}{4}H_0^{(1)}(kr)$ in 2D ($k=k_0\sqrt{\varepsilon_r}$). On the grid, the point source is a single entry $b_j=-1/h$ in 1D and $-1/h^2$ in 2D, so that $\sum_i b_i h=-1$.
3. Grid, boundary, and units
   1. Nodes: $x_i = ih$ for $i=0,\cdots, N-1$, with $u_i = u(x_i)$. $\varepsilon_{r,i}$ lives on the same nodes as $u$.
   2. Boundary: $u=0$ just outside the grid $(u_{-1}=u_N = 0)$. This is a mirror, but it sits behind the PML, so waves are absorbed before they reach it.
   3. PML: the outermost $n_{pml}$ cells on each side, inside the grid. The stretch $s=1+i\sigma/\omega$ is needed at nodes ($s_i$) and at half nodes ($s_{i\pm 1/2}$).
   4. Units: lengths in $\mu m$, and $k_0=2\pi/\lambda_0$ in $\mu m^{-1}$. We drop the physical factor $-i\omega \mu_0$ and work with unit sources. The problem is linear, so that factor only rescales every field.
   5. 2D layout: arrays have shape `(nx, ny)`, and `u[i, j]` $= u(x_i, y_j)$ with $x_i = ih$, $y_j = jh$ (the same $h$ on both axes). Flattened vectors use C order: $k = i\,n_y + j$ (NumPy's default `ravel` and `reshape`).
4. Inner product: $\langle a, b \rangle=a^H b=\sum_i \bar{a}_ib_i$, with the conjugate on the first vector (NumPy: `np.vdot(a,b)`). A measurement is $c^H u$. The discrete product leaves out the factor $h$ that $\int\bar{a}b\,dx \approx ha^H b$ would carry.
5. Adjoint solve: we don't scale the PML rows, so $A$ is not symmetric. $A^T$ is the correct adjoint whether or not $A$ is symmetric, so the code never depends on the PML being written in symmetric form. The adjoint solve always uses $A^T$, with no conjugate. Symmetry of the scaled matrix (reciprocity) becomes a test, not something the code relies on.
6. Differentiating a real $J$ through a complex field: we use the Wirtinger derivative. $\partial J/\partial x$ is taken with $\bar{x}$ held fixed, and for a real parameter $dJ = 2\,\mathrm{Re}\big(\tfrac{\partial J}{\partial x}\,dx\big)$. For example, $J=|x_m|^2$ gives $\partial J/\partial x_m=\bar{x}_m$, see [electrodynamics.md](electrodynamics.md).