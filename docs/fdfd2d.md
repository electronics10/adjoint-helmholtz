# 2D FDFD with PML

**Setup.** $u_{ij}$ and $\varepsilon_{r,ij}$ live on the nodes $(ih, jh)$, $i < n_x$, $j < n_y$. A field is stored as an array of shape $(n_x, n_y)$ and flattened as $k = i\,n_y + j$ (conventions.md, item 3). $D_{xx}$, of shape $n_x\times n_x$, is the 1D stretched second difference along $x$ (fdfd1d.md, Eq. 2 without $k_0^2\varepsilon_r$). $D_{yy}$, of shape $n_y\times n_y$, is the same along $y$. $I_x$ and $I_y$ are identities.

**Eq. 1, the Laplacian.** In 2D the stretched equation is
$$\frac{1}{s_x}\partial_x\Big(\frac{1}{s_x}\partial_xu\Big) + \frac{1}{s_y}\partial_y\Big(\frac{1}{s_y}\partial_yu\Big) + k_0^2\varepsilon_ru = b .$$
$s_x$ depends only on $x$ and $s_y$ only on $y$ (electrodynamics.md, "In 2D and 3D"). So the $x$ term is the 1D operator applied along every column $j$, and the $y$ term is the 1D operator along every row $i$:
$$L = D_{xx}\otimes I_y + I_x\otimes D_{yy}.$$
Why the Kronecker product does this: $(P\otimes Q)_{(i,j),(i',j')} = P_{ii'}Q_{jj'}$ when the pair $(i,j)$ is flattened to $i\,n_y+j$. So
$$[(D_{xx}\otimes I_y)u]_{ij} = \sum_{i',j'}(D_{xx})_{ii'}\,\delta_{jj'}\,u_{i'j'} = \sum_{i'}(D_{xx})_{ii'}\,u_{i'j},$$
which is $D_{xx}$ acting along $x$ at fixed $j$. $L$ is the 5-point stencil, with PML rows near every edge. Both stretches act in the corners.

**Eq. 2, the operator:** $A = L + \mathrm{diag}(k_0^2\varepsilon_r)$, with $\varepsilon_r$ flattened the same way as $u$.

**Eq. 3, symmetry (reciprocity).** In 1D, $S_xD_{xx}$ is symmetric (fdfd1d.md, Eq. 2). Let $S = S_x\otimes S_y$, the diagonal matrix with entries $s_x(i)\,s_y(j)$. By the mixed-product rule $(P\otimes Q)(R\otimes T) = PR\otimes QT$:
$$S(D_{xx}\otimes I_y) = (S_xD_{xx})\otimes S_y,\qquad S(I_x\otimes D_{yy}) = S_x\otimes(S_yD_{yy}).$$
Each is a Kronecker product of a symmetric matrix and a diagonal one, so each is symmetric. $S\,\mathrm{diag}(k_0^2\varepsilon_r)$ is diagonal. So $SA$ is complex symmetric, and $A$ itself is not (conventions.md, item 5).

**Eq. 4, point source:** $b = -1/h^2$ at one node, so that $\sum b\,h^2 = -1$. The response is $G = \frac{i}{4}H_0^{(1)}(kr)$ (conventions.md, item 2).

**Approximations.**
- $H_0^{(1)}$ has a logarithmic singularity at $r=0$ that no grid resolves, so we compare only at $r\ge\lambda/2$.
- The PML is tuned for normal incidence, so waves that hit it at a grazing angle are absorbed less (electrodynamics.md, "In 2D and 3D"). The tests measure the result.