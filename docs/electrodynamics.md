**Note:** this article stems from my conversations with Claude Opus 5.5

## What are physical laws?

Take a field $u(\mathbf{r},t)$ as an example. Most of the $u$ that you can write down never appear in nature. For instance, when you switch on a lamp here, the field far away cannot change at the same instant; light needs time to get there. To discriminate which $u$ are allowed and which aren't, we need a split, which is why most physical laws are in the form of constraints (equations, sometimes inequalities) $$F[u]=0.$$

## Linearity

Now suppose superposition holds: if $u_1$ and $u_2$ are both allowed, so is $au_1+bu_2$. (Two light beams cross and go on as if the other were never there.) So the allowed fields form a vector space. In linear algebra, such a space is always the null space of some linear map. For example, the plane $x+y+z=0$ is a vector space, and it is the null space of the matrix $\begin{pmatrix}1&1&1\end{pmatrix}$. So we can write the law as $Lu=0$ with $L$ linear.

Real experiments, however, often contain what we control (e.g. antenna current, speaker, pulse). Call it the source $s$. Superposition still holds with sources: if $s_1$ gives $u_1$ and $s_2$ gives $u_2$, then $s_1+s_2$ gives $u_1+u_2$. The simplest law with this property is $$Lu=s,$$ which is just $A\mathbf{x}=\mathbf{b}$ for matrices.

What do the solutions look like? Say $u_p$ is one solution, and $u_0$ is anything in the null space ($Lu_0=0$). Then $L(u_p+u_0)=s+0=s$, so $u_p+u_0$ is also a solution. Conversely, if $u_1$ and $u_2$ both solve $Lu=s$, then $L(u_1-u_2)=0$, so they differ by something in the null space. Therefore every solution looks like $$u=u_p+u_0.$$ In other words, the field you measure is the field caused by your source, plus free waves that would be there anyway. Those free waves come from elsewhere or from earlier times. This is why "no source" and "with a source" are the two cases worth separating: they are the two pieces of every solution.

Note that this also means $Lu=s$ alone does not give a unique answer, since we can always add a free wave. To pin down one answer, we need extra conditions: initial conditions, boundary conditions, or "no waves come in from infinity" (radiation condition). Keep this in mind; it comes back when we try to invert $L$.

> [!NOTE] What counts as a source? 
> How do we tell the source $s$ apart from the field $u$? The rule is simple: any term that depends on $u$ belongs to the operator side, and only terms we prescribe independently stay in $s$.
> 
> Glass is a good example. Glass contains electrons. When the field pushes them, they move and form a current. That current is also a source of field, but it is proportional to the field itself, so by the rule it goes into $L$. This is all $\varepsilon_r$ is: a bookkeeping move, the currents of bound electrons moved from the right side to the left side. The antenna current stays on the right because we drive it. So $u$ versus $s$ is really a modeling choice. If the antenna is affected by the field it radiates, its current also depends on $u$, and part of $s$ must move into $L$.
> 
> Two consequences are worth noting now.
> 
> - The electrons do not respond instantly; they keep moving a little after the push. So in time, the bound current is a convolution of the field with some response (see [[convolution]]). At a single frequency, that convolution becomes a plain number, $\varepsilon_r(\omega)$. This is one big reason to work one frequency at a time.
> - When we later design a device by changing $\varepsilon_r(\mathbf{r})$, we are changing $L$, not $s$. So optimizing a device means asking how the solution changes when $L$ changes.

## Why not write the law as $u=Ks$?

How do we get $u$ from $s$? The tempting answer is: if $L$ can be inverted, then $u=L^{-1}s$. Relating back to my previous [[convolution]] note, $L^{-1}$ plays the role of $K$ in $g=Kf$. Its kernel is the [[Green's function]]: the response to a unit point source.

But wait. We just saw that $L$ has a null space, and a matrix with a null space has no inverse. So $L^{-1}$ does not exist yet. It only exists after we add the extra conditions above, because they remove the freedom to add free waves and pick one answer. How the radiation condition does this is a story for later.

Then why don't we write physical laws in the "inverse" form $u=Ks$ directly? In matrix terms, $L$ is sparse while $K$ is dense. A derivative at a point only uses values at neighbouring points (see the next section on how derivatives become matrices), so each row of $L$ has only a few nonzero entries. On the other hand, a source at one point eventually affects every other point, since the effect passes from neighbour to neighbour until it fills the space. So $K$ has almost no zeros. (This holds at one frequency. In space _and_ time, $K$ is zero outside the light cone, because the effect needs time to travel: the lamp example.) The rule is short; its consequence is long, and it depends on the whole geometry and on the extra conditions. In short, nature's laws are short in the $L$ form, and their consequences are written in the $K$ form.

## Locality

So far $L$ is just "some linear operator". Why are physical laws almost always differential equations?

The most general linear rule we can write is an integral, just like in my [[convolution]] note: $$(Lu)(x)=\int W(x,x')\,u(x')\,dx'.$$ (I use $W$ here, since $K$ is already taken by the inverse.) Here $W(x,x')$ says how much the value at $x'$ matters for the rule at $x$.

Now add an assumption: the field here only "listens" to its immediate neighbourhood. Then $W(x,x')$ is nonzero only when $x'$ is very close to $x$. (This is not the same as "nothing happens instantly far away" from the lamp example. The heat equation is local, yet heat spreads to every distance instantly. Finite speed is a separate property of $L$; the wave equation will have it.) Since only $x'$ near $x$ matters, we can Taylor-expand $u(x')$ around $x$: $$u(x')=u(x)+(x'-x)\,u'(x)+\tfrac12(x'-x)^2u''(x)+\dots$$ Put this into the integral, and each term pulls out one derivative of $u$ at $x$: $$(Lu)(x)=m_0(x)\,u(x)+m_1(x)\,u'(x)+m_2(x)\,u''(x)+\dots,\quad m_n(x)=\int W(x,x')\,\frac{(x'-x)^n}{n!}\,dx'.$$ If $W$ only reaches a tiny distance $\epsilon$, then $m_n$ is roughly of size $\epsilon^n$, so the higher terms quickly become negligible. Therefore, a local linear rule is a differential operator with only a few low derivatives. In other words, differential equations are integral equations with very narrow kernels. Here "negligible" means $\epsilon$ is much smaller than the length over which $u$ changes, e.g. the wavelength. When that fails, the higher terms matter; for materials this is spatial dispersion (see "In a medium"). The rigorous version of "local means differential" is **Peetre's theorem**.

Note that if the system is also homogeneous (the same here or there), the $m_n$ cannot depend on $x$, so the coefficients are constant. Then $L$ is a [[convolution]] with a very narrow kernel: $u'=\delta'*u$, $u''=\delta''*u$. Everything in the convolution note applies.

> [!NOTE] **Derivatives as matrices.**
> 
> To see what "sparse" means, put the function on a grid. Pick points $x_1,\dots,x_N$ with spacing $h$, and store $u_i=u(x_i)$: $$\mathbf{u}=(u_1,u_2,\dots,u_N).$$ The function is now a vector, and a linear operator becomes an $N\times N$ matrix.
> 
> What does the second derivative look like on a grid? Taylor-expand one step forward and one step back: $$u_{i\pm1}=u_i\pm h\,u_i'+\tfrac12h^2u_i''\pm\tfrac16h^3u_i'''+\dots$$ Add the two, and the odd terms cancel: $$u_{i+1}+u_{i-1}-2u_i=h^2u_i''+O(h^4)\quad\Longrightarrow\quad u_i''\approx\frac{u_{i-1}-2u_i+u_{i+1}}{h^2}.$$ This is the finite-difference formula, and the pattern $(1,-2,1)$ is called a stencil.
> 
> Row $i$ of the matrix computes $u_i''$, so it only needs $u_{i-1},u_i,u_{i+1}$. Everything else in the row is zero. For $N=5$: $$D_2=\frac{1}{h^2}\begin{pmatrix}-2&1&0&0&0\\1&-2&1&0&0\\0&1&-2&1&0\\0&0&1&-2&1\\0&0&0&1&-2\end{pmatrix}$$ Each row has 3 nonzeros, however large $N$ is. That is exactly what sparse means: entry $(i,j)$ is zero because point $j$ does not directly affect point $i$. Note also that every row is the same stencil shifted by one. This is the grid version of shift invariance (a Toeplitz matrix), and the stencil $(1,-2,1)/h^2$ is the grid version of $\delta''$.
> 
> The 1D Helmholtz operator $u''+k_0^2\varepsilon_r(x)\,u$ is then easy. The second term just multiplies each $u_i$ by its own number, which is a diagonal matrix: $$L=D_2+\mathrm{diag}\big(k_0^2\varepsilon_{r,1},\dots,k_0^2\varepsilon_{r,N}\big).$$ Still sparse. Solving $L\mathbf{u}=\mathbf{s}$ is now literally $A\mathbf{x}=\mathbf{b}$, and this is all a 1D FDFD solver is.
> 
> Now, why is the inverse dense? Take the $3\times3$ version (with $h=1$): $$A=\begin{pmatrix}-2&1&0\\1&-2&1\\0&1&-2\end{pmatrix},\qquad A^{-1}=-\frac14\begin{pmatrix}3&2&1\\2&4&2\\1&2&3\end{pmatrix}.$$ (Multiply them to check.) $A$ has zeros in the corners, but $A^{-1}$ has none. Column $j$ of $A^{-1}$ is the response to a unit source at point $j$, i.e. the Green's function on the grid. The rule only connects neighbours, but the response reaches everywhere.
> 
> Note that this $A$ _is_ invertible, even though I said $L$ has a null space. The trick is hidden in the first and last rows. Writing the matrix this way silently assumes $u=0$ just outside the grid. That is a boundary condition, i.e. one of the extra conditions, and it is what kills the null space. (For this matrix, at least. With the Helmholtz term added, a closed box can still resonate; see the tied rope in "Which solution does nature pick?".) In fact, every finite matrix has made such a choice in its boundary rows. Choosing it well is the whole problem of the radiation condition, which comes later.
> 
> Finally, as $h\to0$ the vector becomes the function again, and a matrix–vector product $\sum_jA_{ij}u_j$ becomes an integral $\int W(x,x')\,u(x')\,dx'$. A sparse matrix becomes a kernel concentrated at $x'=x$ (a derivative of $\delta$), and a dense matrix becomes an ordinary spread-out kernel. So the matrix picture and the integral picture are the same thing.

> [!NOTE] **Locality, stated precisely.** 
> A linear law $Lu=s$ is _local_ if the value of $Lu$ at a point $x$ depends only on $u$ in an arbitrarily small neighbourhood of $x$. Equivalently: if $u=v$ near $x$, then $Lu=Lv$ near $x$. One consequence is used again and again below. Whether $u$ solves the law can be checked point by point. So if we cut the domain into pieces, $u$ solves the law on the whole domain exactly when it solves it inside each piece and across each seam. The conditions across a seam are not new assumptions. They come from integrating the law over a thin strip around the seam (for example, “$u$ is continuous and $u'$ jumps by the source strength”). On a grid this is plain bookkeeping: each row involves only neighbouring values, so a vector solves the system exactly when every row holds. Those rows split into rows inside each piece and the rows at the seams. This is why we may solve piece by piece and glue the results, and why everything outside a region can be replaced by a condition on its edge.

## Homogeneity

So far we know the shape of a linear law, $Lu=s$, where $L$ is a differential operator. We start with the source-free case $s=0$ and treat the general case afterward. Our goal is to first find the free waves: the solutions of $Lu=0$, that is, the null space of $L$.

Linear algebra tells us how to find a null space. Suppose a matrix is diagonal, say $A=\mathrm{diag}(2,0,5)$. Then $A\mathbf{x}=0$ is easy. Each component is multiplied by its own number, so every component with a nonzero number must vanish. The null space is just the second axis. So the strategy is to find a basis in which $L$ is diagonal, that is, to find the eigenvectors of $L$. Then every solution is a combination of eigenvectors, and solving becomes bookkeeping.

For a general $L$, we do not know how to find its eigenvectors. So we add an assumption: the system is physically **homogeneous**, meaning it behaves the same here or there (in both time and space). This assumption is enough. For a homogeneous $L$, the eigenvectors are the exponentials $e^{ikx}$. We give two proofs of this claim. The first is a direct computation. The second shows why homogeneity is what makes it work.

> [!NOTE] 
> Without homogeneity, the eigenvectors still exist, but each operator has different ones. Finding them is as hard as the original problem.
> 
> The eigenvector equation is itself a differential equation: $Lf=\lambda f$. When $L$ has constant coefficients, the guess $e^{ikx}$ solves it at once. When $L$ has coefficients like $m(x)$, $Lf=\lambda f$ is a variable-coefficient equation. Usually it has no solution in elementary functions. So we are back where we started: to solve $Lu=0$, we must first solve $Lf=\lambda f$, which is no easier.
> 
> Some examples (established):
> 
> - $L=-\frac{d^2}{dx^2}+x^2$ (the quantum harmonic oscillator). Its eigenvectors are the **Hermite functions**.
> - $L=\frac{d^2}{dx^2}-x$. Its eigenvectors are shifted **Airy functions**.
> - Bessel's operator gives the **Bessel functions**.
> 
> Each case needs its own special functions, found by its own analysis.

### Proof 1: direct computation

For a differential operator, homogeneity (uniform space, not to be conflated with mathematically homogeneous equations) means constant coefficients. (Proof 2 shows why.) So we take $L=p\!\left(\frac{d}{dx}\right)$ for a polynomial $p$.

Take $e^{ikx}$ and differentiate it: $$\frac{d}{dx}e^{ikx}=ik\,e^{ikx},\qquad \frac{d^2}{dx^2}e^{ikx}=(ik)^2e^{ikx}=-k^2e^{ikx}.$$ The result is the same function times a number. So $e^{ikx}$ is an eigenvector of $\frac{d^n}{dx^n}$ for every $n$, and therefore of $L$: $$L\,e^{ikx}=p(ik)\,e^{ikx}.$$

The coefficients must be constant. If a coefficient depends on $x$, say $m(x)\,\frac{d}{dx}$, then we get $ik\,m(x)\,e^{ikx}$. This is not a number times $e^{ikx}$.

### Proof 2: commuting with shifts

This proof rests on a fact from linear algebra: **commuting operators share eigenvectors** (when eigenvalues are distinct). For example, let $S=\mathrm{diag}(1,2)$. If $AS=SA$, then $A$ is diagonal too, so $A$ has the same eigenvectors as $S$. If we know the eigenvectors of $S$, we get those of $A$ for free.

Here the role of $S$ is played by the shifts. Write $T_a$ for the shift, $(T_a u)(x)=u(x+a)$. In operator language, homogeneous means shift-invariant: $L$ commutes with every shift: $LT_a=T_aL$.

For a differential operator, commuting with shifts is the same as having constant coefficients. Indeed, $T_a\big(m(x)\,u'\big)=m(x+a)\,u'(x+a)$, while $m(x)\,(T_a u)'=m(x)\,u'(x+a)$. These agree for every $u$ only if $m(x+a)=m(x)$ for every $a$, that is, only if $m$ is constant. The same computation works with $u$ or any higher derivative in place of $u'$, so it covers every term of $L$.

The eigenvectors of the shifts are the exponentials: $T_a e^{ikx}=e^{ika}e^{ikx}$. Since $L$ commutes with shifts, $Le^{ikx}$ is an eigenvector of every shift with the same eigenvalue: $$T_a\big(Le^{ikx}\big)=L\big(T_a e^{ikx}\big)=e^{ika}\,Le^{ikx}.$$ The only functions $f$ with $T_a f=e^{ika}f$ for every $a$ are multiples of $e^{ikx}$. (Set $x=0$ to get $f(a)=e^{ika}f(0)$.) So $Le^{ikx}=\lambda\,e^{ikx}$ for some number $\lambda$.

This proof never uses the formula for $L$. It only uses the fact that $L$ commutes with shifts. That is why it also covers convolutions, as in my [[convolution]] note.

### Back to $Lu=0$

Since $Le^{ikx}=p(ik)\,e^{ikx}$, the exponential $e^{ikx}$ is a free wave exactly when $p(ik)=0$. This condition on $k$ is called the dispersion relation (see the following sections). Combinations of these exponentials give the free waves. (The shift argument also allows complex $k$. All it needs is $f(a)=c(a)f(0)$ with $c(a+b)=c(a)c(b)$, and $c(a)=e^{za}$ works for any complex $z$. A complex $k$ gives a wave that decays or grows along $x$: an **evanescent wave**. We meet them again in "What does a free wave look like?".)

(Strictly, the $e^{ikx}$ are not a basis in the finite-dimensional sense. The claim that every function is a combination of them is the content of the [[Fourier transform]], where the sum becomes an integral over $k$.)

The payoff: $L$ acts on each component separately, $$Lu=\int_{-\infty}^{\infty}A(k)\,p(ik)\,e^{ikx}\,dk.$$ So in the $k$ picture, $L$ is just multiplication by $p(ik)$, called the _symbol_ of $L$.

## Space, time, and free waves

### Adding time

The field depends on $t$ too, and the system is also the same now or later. So now there are two kinds of shifts: $T_a$ in space, $(T_a u)(x,t)=u(x+a,t)$, and $S_b$ in time, $(S_b u)(x,t)=u(x,t+b)$. Note that they commute with each other (shifting in space then time is the same as time then space), and $L$ commutes with both.

So let us run the same argument as before, with both shifts at once. Which functions are eigenvectors of every $T_a$ and every $S_b$? Suppose $$f(x+a,\,t+b)=e^{ika}\,e^{-i\omega b}\,f(x,t)\quad\text{for every }a,b.$$ Set $x=0$ and $t=0$, and we get $f(a,b)=e^{ika}e^{-i\omega b}f(0,0)$. So the only such functions are multiples of $$e^{ikx}\,e^{-i\omega t}=e^{i(kx-\omega t)}.$$ So the product form is not an assumption; it is forced. (I write $-i\omega t$ instead of $+i\omega t$ on purpose; the reason comes in a moment.) And by the same argument as before, since $L$ commutes with both shifts, $L\,e^{i(kx-\omega t)}$ is again a number times $e^{i(kx-\omega t)}$.

We can also check it directly. $\partial_x$ only sees the $x$ part, and $\partial_t$ only sees the $t$ part: $$\partial_x\,e^{i(kx-\omega t)}=ik\,e^{i(kx-\omega t)},\qquad \partial_t\,e^{i(kx-\omega t)}=-i\omega\,e^{i(kx-\omega t)}.$$

Is this family big enough? Yes: do Fourier twice. First expand in $x$ at each time, $u(x,t)=\int A(k,t)\,e^{ikx}\,dk$. Then expand each coefficient $A(k,t)$ in time. Together, $$u(x,t)=\iint A(k,\omega)\,e^{i(kx-\omega t)}\,dk\,d\omega.$$ So, at this point, $k$ and $\omega$ are two independent labels. Nothing connects them yet.

What does one basis function look like? Its crests sit where $kx-\omega t$ is constant, i.e. at $x=(\omega/k)\,t+\text{const}$. So it is a **wave**, with crests moving at speed $\omega/k$ (the phase velocity)! The number $k$ counts crests per length ($k=2\pi/\text{wavelength}$), and $\omega$ counts oscillations per time ($\omega=2\pi/\text{period}$). With the minus sign, $k>0$ and $\omega>0$ mean the wave moves to $+x$. This is called a plane wave.

### The symbol and the dispersion relation

So far, linearity, locality and homogeneity in space and time tell us that every free solution is a combination of plane waves. Now we ask which plane waves a given law allows.

$L$ is a polynomial in two derivatives, $L=p(\partial_x,\partial_t)$, and on a plane wave $$L\,e^{i(kx-\omega t)}=p(ik,\,-i\omega)\,e^{i(kx-\omega t)}=D(k,\omega)\,e^{i(kx-\omega t)}.$$ $D(k,\omega)$ is the **symbol** again, now with two labels. For example, $$L=a\,\partial_t^2+b\,\partial_x^2+c\,\partial_x+d\,\partial_t+e \quad\Longrightarrow\quad D(k,\omega)=-a\omega^2-bk^2+ic\,k-id\,\omega+e.$$

What does each term do? Here is the simplest case: drop space, and take $L=\partial_t+\gamma$ with $\gamma>0$. Then $D=-i\omega+\gamma$, and $D=0$ gives $\omega=-i\gamma$. A complex frequency! The wave is $e^{-i\omega t}=e^{-\gamma t}$, which decays. So a first-order time derivative means loss: the field dies out, like a swing with friction.

This also tells us when such terms are absent. For a single real scalar field, a law without loss looks the same when the movie is played backward ($t\to-t$). Under $t\to-t$, $\partial_t$ flips sign but $\partial_t^2$ does not. So such a law can only have even time derivatives. In the same way, a law that looks the same in a mirror ($x\to-x$) can only have even space derivatives, so $c=0$. (Both rules need the restriction to a single real scalar field. For systems of fields, such as Maxwell's equations, they fail; see the note under "Check 3: energy".) (For example, light in vacuum has both properties, so for it, only $\partial_t^2$, $\partial_x^2$ and a constant can appear. We will come back to loss with $\varepsilon_r$: it is exactly what a complex $\varepsilon_r$ describes.)

Before any law, a function of $x$ and $t$ can wiggle at any rate in space and at any rate in time, independently. For example, $\cos(3x)\cos(5t)$ is a perfectly good function, built from $(k,\omega)=(\pm3,\pm5)$. So in the expansion, every pair $(k,\omega)$ is allowed.

Apply $L$ to a general field: $$Lu=\iint A(k,\omega)\,D(k,\omega)\,e^{i(kx-\omega t)}\,dk\,d\omega.$$ For a free wave, this must be zero. The plane waves are independent (like the axes of a basis), so each component must vanish separately: $$A(k,\omega)\,D(k,\omega)=0\quad\text{for every }(k,\omega).$$ So $A(k,\omega)$ can only be nonzero where $$D(k,\omega)=0.$$ Picture the $(k,\omega)$ plane. A general function spreads its $A$ over the whole plane. The set $D=0$ is a few curves in that plane, and a free wave puts its $A$ only on those curves. This condition is called the **dispersion relation**.

> [!NOTE] 
> There is a problem here. If $A$ is an ordinary function that is nonzero only on a curve, the double integral gives zero, since a curve has no area. So $A$ must be concentrated on the curve, like a $\delta$. For each $k$, the equation $D(k,\omega)=0$ has only a few roots in $\omega$. Call them $\omega_1(k),\omega_2(k),\dots$ (if $D$ is second order in $\omega$, there are two). Then $$A(k,\omega)=A_1(k)\,\delta\big(\omega-\omega_1(k)\big)+A_2(k)\,\delta\big(\omega-\omega_2(k)\big)+\dots$$ Put this into the expansion. The $\omega$ integral just picks out each root, since $\int\delta(\omega-\omega_1)\,e^{-i\omega t}\,d\omega=e^{-i\omega_1t}$. So $$u(x,t)=\int\Big[A_1(k)\,e^{i(kx-\omega_1(k)t)}+A_2(k)\,e^{i(kx-\omega_2(k)t)}+\dots\Big]\,dk.$$ This is the general free wave: for each $k$, a few waves, each with its own amplitude.
> 
> Note that the amplitudes $A_1(k),A_2(k),\dots$ are still free. This is the null space from before: there are many free waves, and a law alone does not choose among them. The extra conditions choose. For example, with two roots, there are two amplitudes for each $k$, so we need two pieces of information, such as the field and its time derivative at $t=0$. We will see this concretely for light in the next section.

## Light in free space

Now we use one experimental fact: in vacuum, light of every colour travels at the same speed $c$.

### Light in 1D: the wave equation

In 1D, it can go right or left. So the allowed plane waves are exactly those with $$\omega=ck \quad \text{or} \quad \omega=-ck.$$ We want one polynomial $D$ that vanishes on both. The product does it: $$D(k,\omega)=(\omega-ck)(\omega+ck)=\omega^2-c^2k^2.$$ Now swap the numbers back into derivatives, $\omega^2\leftrightarrow-\partial_t^2$ and $k^2\leftrightarrow-\partial_x^2$: $$\partial_x^2u-\frac{1}{c^2}\,\partial_t^2u=0.$$ This is the wave equation. We did not guess it; it is what "every colour at the same speed, in both directions" looks like as a law. One honest caveat: any multiple $D\,q(k,\omega)$ also vanishes on both lines. The wave equation is the _lowest-order_ law with exactly these free waves, and choosing the lowest order is an assumption.

Now we can write down every free wave of light in 1D. The two roots are $\omega_1(k)=ck$ and $\omega_2(k)=-ck$, so $$u(x,t)=\int\Big[A_1(k)\,e^{ik(x-ct)}+A_2(k)\,e^{ik(x+ct)}\Big]\,dk=F(x-ct)+G(x+ct),$$ where $F(\xi)=\int A_1(k)\,e^{ik\xi}\,dk$ and $G(\xi)=\int A_2(k)\,e^{ik\xi}\,dk$. So a free wave is any shape moving right at speed $c$, plus any shape moving left at speed $c$. (This is called d'Alembert's solution.)

How do the initial conditions fix the amplitudes? Say we know $u(x,0)=f(x)$ and $\partial_tu(x,0)=g(x)$, and expand both: $f(x)=\int\hat f(k)\,e^{ikx}\,dk$, and the same for $g$. Setting $t=0$ in $u$ and in $\partial_tu$ gives, for each $k$, $$\hat f(k)=A_1(k)+A_2(k),\qquad \hat g(k)=-ick\,A_1(k)+ick\,A_2(k).$$ Two equations, two unknowns. Solving, $$A_1(k)=\tfrac12\Big(\hat f(k)+\frac{i\,\hat g(k)}{ck}\Big),\qquad A_2(k)=\tfrac12\Big(\hat f(k)-\frac{i\,\hat g(k)}{ck}\Big).$$ So the two initial conditions pick exactly one free wave out of the null space, as promised. (At $k=0$ the two roots meet and the formula divides by zero. That mode is the uniform motion $u=a+bt$, and $\hat g(0)$ sets $b$.)

(Incidentally, the factorization also works with derivatives: $(\partial_t-c\,\partial_x)(\partial_t+c\,\partial_x)u=0$, one factor for each direction.)

### Going to 3D: richer geometry

In 3D, there are three directions to shift in, so there are three labels, $k_x,k_y,k_z$, and the plane waves are $$e^{i(k_xx+k_yy+k_zz-\omega t)}=e^{i(\mathbf{k}\cdot\mathbf{r}-\omega t)}.$$ The vector $\mathbf{k}=(k_x,k_y,k_z)$ points in the direction the wave travels, and its length $|\mathbf{k}|$ plays the role of $k$ in 1D. On such a wave, $\partial_x\to ik_x$, $\partial_y\to ik_y$, $\partial_z\to ik_z$.

The experimental fact is now: same speed $c$ in every direction. So the allowed waves satisfy $\omega=\pm c|\mathbf{k}|$ for every direction of $\mathbf{k}$, and the same product trick gives $$D(\mathbf{k},\omega)=\omega^2-c^2|\mathbf{k}|^2=\omega^2-c^2\big(k_x^2+k_y^2+k_z^2\big).$$ Swap back, with $k_x^2\leftrightarrow-\partial_x^2$ and so on: $$\big(\partial_x^2+\partial_y^2+\partial_z^2\big)u-\frac{1}{c^2}\,\partial_t^2u=0.$$ The combination $\partial_x^2+\partial_y^2+\partial_z^2$ is called the Laplacian, $\nabla^2$. So this is all the Laplacian is: the derivative version of squared length, $|\mathbf{k}|^2$. That is where the geometry comes in. Rotating your axes does not change a length, so $\nabla^2$ does not depend on how you orient $x,y,z$. That is exactly what "same speed in every direction" needs.

Why second order, and not first? A first-order term like $\partial_x$ would put $k_x$ alone into $D$, which singles out the $x$ direction. The simplest combination that treats all directions alike is $|\mathbf{k}|^2$, which is second order. (This argument assumes $u$ is a single number at each point. Light is not; see the Maxwell section.)

## Light in mediums

### Nonlocal in time

In glass, light is slower by a factor $n$, so $\omega=\pm c|\mathbf{k}|/n$, and the wave equation becomes $$\nabla^2u-\frac{n^2}{c^2}\,\partial_t^2u=0.$$ In real glass, however, $n$ depends on colour (that is how a prism works). Then the allowed waves satisfy $\omega=\pm c|\mathbf{k}|/n(\omega)$, and $D$ is no longer a polynomial in $\omega$. So $L$ is no longer a differential operator in time.

What happened to locality? Recall the bound electrons from earlier: they keep moving a little after the push. The glass has **memory**, so its response at time $t$ depends on the field at earlier times. This is a convolution in time with a kernel of real width. But a convolution in time becomes plain multiplication at a single frequency. So at one $\omega$, the whole memory of the glass is just one number, $$n(\omega)^2=\varepsilon_r(\omega).$$ This is the deep reason textbooks work **one frequency at a time**. In time, the material is a convolution; at one frequency, it is a number.

### Nonhomogeneous in space

Finally, let $n$ change with position, $n(\mathbf{r})$: a lens, a waveguide, a device we want to design. Now the system is not the same here or there, so the exponentials in space are no longer eigenvectors (recall the $m(x)\,\frac{d}{dx}$ example).

The system is still the same now or later, so the exponentials in time still work. So we keep only those. Pick one frequency, and write the field as $u(\mathbf{r})\,e^{-i\omega t}$. Then $\partial_t^2\to-\omega^2$, but $\nabla^2$ has to stay as it is. With $k_0=\omega/c$: $$\nabla^2u+k_0^2\,\varepsilon_r(\mathbf{r})\,u=s.$$ This is the Helmholtz equation. It is algebra in time and a differential equation in space, and which half is which is decided by the symmetry that survives. On the grid, it is the sparse matrix $D_2+\mathrm{diag}(k_0^2\varepsilon_{r,i})$ from before, which is exactly what an FDFD solver builds.

Putting $\varepsilon_r(\mathbf{r})$ inside the scalar equation was a guess. It works when $\varepsilon_r$ changes slowly compared with the wavelength. At a sharp interface it can fail, because light is not a scalar: it has two polarizations, so it is a vector. Maxwell's equations carry this extra structure, and they will tell us exactly when the scalar equation is right.

## Light is a vector: Maxwell's equations

The section "Nonhomogeneous in space" ended with a problem. Light has two polarizations, so the field is a vector, and the scalar equation was only a guess. We need the law for a vector field.

Here I stop deriving and take the law from experiment, just as we took "same speed $c$" from experiment before. A deeper argument for its form exists (see the note at the end of this section), but nothing later depends on it. What matters is that every consequence we use below can be checked.

The law uses two vector fields, the electric field $\mathbf E$ and the magnetic field $\mathbf B$. The sources are the charge density $\rho$ and the current density $\mathbf J$. In SI units, with $\mathbf D=\varepsilon_0\mathbf E$ and $\mathbf H=\mathbf B/\mu_0$ in vacuum, there are four equations, each backed by an experiment:

|Equation|Experiment|
|---|---|
|$\nabla\cdot\mathbf D=\rho$|Coulomb: a charge makes an $\mathbf E$ field that spreads out from it|
|$\nabla\cdot\mathbf B=0$|no magnetic charge has ever been found|
|$\nabla\times\mathbf E=-\partial_t\mathbf B$|Faraday: a changing magnetic field induces a voltage around a loop|
|$\nabla\times\mathbf H=\mathbf J+\partial_t\mathbf D$|Ampère: a current makes a magnetic field; plus Maxwell's term $\partial_t\mathbf D$|

Maxwell's term $\partial_t\mathbf D$ is the only piece that no experiment of his time showed directly. It is forced by charge conservation. Take the divergence of the last equation. Since $\nabla\cdot\nabla\times=0$, $$0=\nabla\cdot\mathbf J+\partial_t(\nabla\cdot\mathbf D)=\nabla\cdot\mathbf J+\partial_t\rho .$$ Without $\partial_t\mathbf D$, the law would demand $\nabla\cdot\mathbf J=0$, so charge could never pile up anywhere. That is false: a capacitor charges. With $\partial_t\mathbf D$, the law demands exactly $\partial_t\rho+\nabla\cdot\mathbf J=0$, which experiment confirms. In the language of the Linearity section, charge conservation is the condition for the source $s$ to be compatible with $L$. If $\rho$ and $\mathbf J$ are given separately and violate it, the four equations have no solution at all. (This is the Fredholm alternative: $Lu=s$ is solvable only if $s$ is orthogonal to the null space of the adjoint $L^*$; see the adjoint note in [[convolution]].) At one frequency the picture changes: the curl-curl equation has a solution for any $\mathbf J$, and Gauss's law then defines $\rho$ from $\mathbf J$ (see "Gauss's law comes for free").

### Check 1: light is a wave with speed $c$

In vacuum, $\rho=0$ and $\mathbf J=0$. Take the curl of Faraday's law and insert Ampère's law: $$\nabla\times\nabla\times\mathbf E=-\partial_t\nabla\times\mathbf B=-\mu_0\varepsilon_0\,\partial_t^2\mathbf E .$$ Use the identity $\nabla\times\nabla\times\mathbf E=\nabla(\nabla\cdot\mathbf E)-\nabla^2\mathbf E$ and $\nabla\cdot\mathbf E=0$: $$\nabla^2\mathbf E-\mu_0\varepsilon_0\,\partial_t^2\mathbf E=0 .$$ Each component obeys the scalar wave equation of the 3D section, with $c=1/\sqrt{\mu_0\varepsilon_0}$. Maxwell computed this number from static measurements of $\varepsilon_0$ and $\mu_0$, and it matched the measured speed of light (established). So the scalar wave equation was a shadow of this one, one component at a time.

### Check 2: two polarizations

Take a plane wave $\mathbf E=\mathbf E_0\,e^{i(\mathbf k\cdot\mathbf r-\omega t)}$. In vacuum $\nabla\cdot\mathbf E=0$, which becomes $\mathbf k\cdot\mathbf E_0=0$: the field is perpendicular to the direction of travel. That leaves two directions, the two polarizations, as observed. Faraday's law gives $\mathbf B_0=\mathbf k\times\mathbf E_0/\omega$, also perpendicular. The scalar equation knew about neither constraint. This is what it was missing.

### Check 3: energy

Dot Faraday's law with $\mathbf H$ and Ampère's law with $\mathbf E$, subtract, and use $\nabla\cdot(\mathbf E\times\mathbf H)=\mathbf H\cdot\nabla\times\mathbf E-\mathbf E\cdot\nabla\times\mathbf H$: $$\partial_t\Big(\tfrac12\varepsilon_0|\mathbf E|^2+\tfrac12\mu_0|\mathbf H|^2\Big)+\nabla\cdot(\mathbf E\times\mathbf H)=-\mathbf J\cdot\mathbf E .$$ The energy at a point changes in only two ways: it flows out with the **Poynting vector** $\mathbf S=\mathbf E\times\mathbf H$, or it does work on the current, $\mathbf J\cdot\mathbf E$. We will use the one-frequency version of this for objective functions.

> [!NOTE] **Why "a single real scalar field"?** 
> In "The symbol and the dispersion relation", the even-derivative rule was limited to a single real scalar field. Here is why. Maxwell's equations are first order in time, yet Check 3 shows they lose nothing: in vacuum, energy only flows. They are unchanged under $t\to-t$ if also $\mathbf B\to-\mathbf B$. The even-derivative rule holds only for a single real scalar field. The general test for "lossless" is: every root $\omega(\mathbf k)$ is real for real $\mathbf k$. Also, lossless and time-reversal symmetric are different things. A magnetized ferrite is lossless but not time-reversal symmetric (it is non-reciprocal).

### The divergence equations are constraints

Count first. There are 6 unknowns (the components of $\mathbf E$ and $\mathbf B$) but 8 equations, so the law looks overdetermined. The two curl equations give $\partial_t\mathbf E$ and $\partial_t\mathbf B$: 6 equations, already enough to step the fields forward in time. The two divergence equations contain no $\partial_t$, so they predict nothing. Could they contradict the curl equations? Take the divergence of the two curl equations. Faraday gives $\partial_t(\nabla\cdot\mathbf B)=0$. Ampère gives $\partial_t(\nabla\cdot\mathbf D)=-\nabla\cdot\mathbf J=\partial_t\rho$. So if the two divergence equations hold at one moment, the curl equations keep them true forever. They constrain the starting state; they do not drive anything.

A small example has the same logic. The motion $\dot x=y$, $\dot y=-x$ keeps $x^2+y^2$ constant. Demanding $x^2+y^2=1$ is not a second law of motion; it only picks a starting point on the circle, and the motion keeps you there.

Then why keep the divergence equations at all? Because the curl equations only _keep_ the divergence; they cannot _choose_ it. Add any static gradient field $\nabla\phi$ to a solution: it has no curl, so the curl equations do not notice it, and the sum is still a solution. Without Gauss's law, such a field could sit in empty space with no charge to produce it. Gauss's law forbids this, by tying $\nabla\cdot\mathbf D$ to real charge. In the circle example, it is the statement that the radius is 1 and not 2. Two more reasons. In statics ($\partial_t=0$), the curl equations say almost nothing, and the divergence equations are most of the physics (electrostatics is Gauss's law). And $\nabla\cdot\mathbf B=0$ is itself an experimental fact (no magnetic charge).

In practice, the two curl equations are the law of motion, and they are what a time-domain solver steps forward. The divergence equations only need to hold in the initial fields. This comes back at one frequency.

> [!NOTE] **The deeper reason (optional).** 
> Why curls, and why two fields? A short answer exists, though I don't use it later. A law that is linear, local, the same in every direction and first order in time can only use $\nabla\times$, $\nabla\cdot$ and $\nabla$ on vector fields. With one real vector field, $\partial_t\mathbf u\propto\nabla\times\mathbf u$ gives imaginary $\omega$: growth and decay, no waves. Two fields that drive each other with opposite signs give exactly the light cone, and the minus sign turns out to be energy conservation. So Maxwell's equations are, in a precise sense, the simplest such law. The same idea gives Dirac's equation for the electron. Keywords: symbol of a first-order system, isotropic tensors, Dirac's trick.

## Maxwell's equations: a revisit

### In a medium

The section "Light in mediums" said the material's memory is a convolution in time. Now we can make that precise. Matter contains bound charges, and $\mathbf E$ pushes them into a current $\mathbf J_b$. Split $\mathbf J=\mathbf J_s+\mathbf J_b$: the part we drive and the part that responds. By the rule from the "What counts as a source?" note, $\mathbf J_b$ goes into $L$.

What can $\mathbf J_b$ look like? We apply the same assumptions as before, now to the material. It is linear (weak fields), time-invariant (the material does not change), causal, and local in space (an atom responds to the field where it sits). Write $\mathbf J_b=\partial_t\mathbf P$, where the polarization $\mathbf P$ is the displacement of bound charge per volume. Then $$\mathbf P(\mathbf r,t)=\varepsilon_0\int_0^\infty\chi(\mathbf r,\tau)\,\mathbf E(\mathbf r,t-\tau)\,d\tau .$$ This is a causal convolution in time (see [[convolution]]), local in space. $\chi$ depends on $\mathbf r$ when the material is not uniform. Ampère's law becomes $\nabla\times\mathbf H=\mathbf J_s+\partial_t\mathbf D$ with $\mathbf D=\varepsilon_0\mathbf E+\mathbf P$.

A toy model shows where $\chi$ comes from. Take an electron on a spring with friction: $m(\ddot x+\gamma\dot x+\omega_0^2x)=-eE$. At one frequency, $x=-(e/m)E/(\omega_0^2-\omega^2-i\gamma\omega)$. With $N$ electrons per volume, $P=-Nex$, so $$\chi(\omega)=\frac{\omega_p^2}{\omega_0^2-\omega^2-i\gamma\omega},\qquad \omega_p^2=\frac{Ne^2}{\varepsilon_0m}.$$ This is the Lorentz model. Its imaginary part comes from the friction $\gamma$, the first-order time derivative from the dispersion-relation section. At one frequency, $\mathbf D=\varepsilon_0\varepsilon_r(\mathbf r,\omega)\mathbf E$ with $\varepsilon_r=1+\chi$.

Three consequences:

- **Loss.** $\varepsilon_r$ is complex. With the $e^{-i\omega t}$ convention, $\operatorname{Im}\varepsilon_r>0$ means loss (proved with the Poynting theorem below, under "At one frequency"). The Lorentz model agrees: $\operatorname{Im}\chi>0$ when $\gamma>0$.
- **Dispersion is unavoidable.** $\chi(\tau)$ is real, so $\varepsilon_r(-\omega)=\overline{\varepsilon_r(\omega)}$. $\chi$ is causal, so $\operatorname{Re}\varepsilon_r$ and $\operatorname{Im}\varepsilon_r$ determine each other ([[Kramers-Kronig relations|Kramers–Kronig relations]]). One consequence: a material with loss at any frequency must have an $\varepsilon_r$ that depends on $\omega$. A constant $\varepsilon_r\neq1$ is always a narrow-band approximation.
- **Spatial locality is an approximation.** It holds when the atomic (or meta-atom) scale is much smaller than the wavelength. It fails in nanoscale plasmonics and near metamaterial resonances, where $\varepsilon_r$ also depends on $\mathbf k$ (spatial dispersion).

(Bound charges can also circulate, giving a magnetization current $\nabla\times\mathbf M$ and $\mu_r$. In optics $\mu_r\approx1$. In RF, ferrites need $\mu_r$, and a biased ferrite needs a tensor $\mu$ that is non-reciprocal. I set $\mu_r=1$ from here on.)

### At one frequency

As in "Light in mediums", we keep the time exponentials and work at one $\omega$. From here on, $\mathbf J$ means the driven current $\mathbf J_s$. Real fields have $\hat u(-\omega)=\overline{\hat u(\omega)}$, so at a single frequency the physical field is $\operatorname{Re}\big(\mathbf E(\mathbf r)\,e^{-i\omega t}\big)$.

With $\partial_t\to-i\omega$: $$\nabla\times\mathbf E=i\omega\mu_0\mathbf H,\qquad \nabla\times\mathbf H=-i\omega\varepsilon_0\varepsilon_r(\mathbf r)\mathbf E+\mathbf J .$$ Eliminate $\mathbf H$: $$\nabla\times\nabla\times\mathbf E-k_0^2\,\varepsilon_r(\mathbf r)\,\mathbf E=i\omega\mu_0\mathbf J,\qquad k_0=\omega/c .$$ This is $Lu=s$ from the Linearity section, with $L=\nabla\times\nabla\times-k_0^2\varepsilon_r$, $u=\mathbf E$ and $s=i\omega\mu_0\mathbf J$. It is the vector version of the Helmholtz equation. (We can also keep both $\mathbf E$ and $\mathbf H$ as unknowns, as the Yee grid and FDTD do. The physics is the same; the matrices differ.)

**Gauss's law comes for free.** Take the divergence; $\nabla\cdot\nabla\times=0$ gives $-k_0^2\nabla\cdot(\varepsilon_r\mathbf E)=i\omega\mu_0\nabla\cdot\mathbf J$. With charge conservation, $\nabla\cdot\mathbf J=i\omega\rho$, this becomes $\nabla\cdot(\varepsilon_0\varepsilon_r\mathbf E)=\rho$. So at $\omega\neq0$ the divergence constraint is already inside the curl-curl equation.

> [!NOTE] 
> What happens as $\omega\to0$? The term $k_0^2\varepsilon_r$ vanishes, and every gradient $\nabla\phi$ is in the null space of $\nabla\times\nabla\times$ (a gradient has no curl). Only the divergence equation could control these fields, and at $\omega=0$ it is no longer implied. So at low frequency $L$ becomes nearly singular. Solvers call this **low-frequency breakdown**, and FEM and MoM codes need special fixes for it. A related point: a discretization must keep "curl of a gradient is zero" exactly, or false (spurious) modes appear. The Yee grid and edge elements do keep it. This comes back in "The Yee grid".

### The symbol is now a matrix

Run the homogeneity argument once more, on the curl-curl equation in vacuum ($\varepsilon_r=1$, no source). On a plane wave $\mathbf E_0e^{i\mathbf k\cdot\mathbf r}$, $\nabla\to i\mathbf k$, so $$\nabla\times\nabla\times\mathbf E\ \to\ -\mathbf k\times(\mathbf k\times\mathbf E_0)=|\mathbf k|^2\mathbf E_0-\mathbf k\,(\mathbf k\cdot\mathbf E_0).$$ So on a plane wave, $L$ is not multiplication by a number but by a $3\times3$ matrix: $$D(\mathbf k)=\big(|\mathbf k|^2-k_0^2\big)I-\mathbf k\mathbf k^{T}.$$ A free wave needs $D\mathbf E_0=0$ with $\mathbf E_0\neq0$. Check the two kinds of direction.

- $\mathbf E_0\perp\mathbf k$: $D\mathbf E_0=(|\mathbf k|^2-k_0^2)\mathbf E_0$. This is zero exactly when $|\mathbf k|=k_0$, the same dispersion relation as the scalar case.
- $\mathbf E_0\parallel\mathbf k$: $D\mathbf E_0=-k_0^2\mathbf E_0$. This is never zero at $\omega\neq0$.

So free waves need $|\mathbf k|=k_0$ and $\mathbf E_0\perp\mathbf k$. The scalar story gave only the first condition. The matrix also fixes the direction of the amplitude: this is polarization (Check 2), now coming straight from the dispersion relation. Note also the entry $-k_0^2$ along $\mathbf k$. It goes to zero as $\omega\to0$, which is the low-frequency breakdown from the note above.

### When is the scalar Helmholtz equation right?

Now we can answer the question from "Nonhomogeneous in space". Take a structure that does not change along $z$, so $\varepsilon_r=\varepsilon_r(x,y)$, as in any 2D simulation (where sources do not change along $z$ either). Fields then do not depend on $z$ either. With $\partial_z=0$, Maxwell's equations split into two independent groups. $E_z$ only produces $H_x,H_y$ (since $(\nabla\times\mathbf E)_x=\partial_yE_z$ and $(\nabla\times\mathbf E)_y=-\partial_xE_z$), and $H_x,H_y$ only produce $E_z$ (since $(\nabla\times\mathbf H)_z=\partial_xH_y-\partial_yH_x$). So $\{E_z,H_x,H_y\}$ is closed, and likewise $\{H_z,E_x,E_y\}$. A scalar $\varepsilon_r$ does not turn vectors, so it does not mix the groups. Every 2D solution is a sum of one from each group, so we can study them separately. (The split fails for a tensor $\varepsilon_r$, a source that varies along $z$, or oblique incidence, $k_z\neq0$.) The two cases are:

- **$\mathbf E=E_z(x,y)\,\hat{\mathbf z}$.** Then $\nabla\cdot(\varepsilon_r\mathbf E)=\partial_z(\varepsilon_rE_z)=0$ automatically, and $\nabla\times\nabla\times(E_z\hat{\mathbf z})=-\nabla^2E_z\,\hat{\mathbf z}$. So $$\nabla^2E_z+k_0^2\varepsilon_r(x,y)\,E_z=-i\omega\mu_0J_z .$$ This is exactly the scalar Helmholtz equation. There is no approximation, even at sharp interfaces.
- **$\mathbf H=H_z(x,y)\,\hat{\mathbf z}$.** Eliminating $\mathbf E$ instead of $\mathbf H$ gives (without sources) $$\nabla\cdot\big(\varepsilon_r^{-1}\nabla H_z\big)+k_0^2H_z=0 .$$ Where $\varepsilon_r$ is constant, this is the Helmholtz equation again. At an interface, however, the $1/\varepsilon_r$ sits inside the derivative. It makes $\varepsilon_r^{-1}\partial_nH_z$ continuous across the interface (this is the tangential $\mathbf E$), and the scalar guess gets exactly this wrong.

So in 2D, one polarization is exactly scalar, and the other is scalar with a different operator. In 3D with a general $\varepsilon_r(\mathbf r)$, the components are coupled, and no scalar equation is exact. (These two cases are called TE and TM, but the labels are swapped between communities. Check which field is along $z$ before trusting a label.)

### Power

Products of real fields need care. The time average of $\operatorname{Re}(ae^{-i\omega t})\operatorname{Re}(be^{-i\omega t})$ is $\tfrac12\operatorname{Re}(a\bar b)$. So the time-averaged Poynting vector is $\langle\mathbf S\rangle=\tfrac12\operatorname{Re}(\mathbf E\times\overline{\mathbf H})$. Use $\nabla\cdot(\mathbf E\times\overline{\mathbf H})=\overline{\mathbf H}\cdot\nabla\times\mathbf E-\mathbf E\cdot\nabla\times\overline{\mathbf H}$ and the two equations above, and take the real part: $$\nabla\cdot\langle\mathbf S\rangle=-\tfrac12\operatorname{Re}\big(\mathbf E\cdot\overline{\mathbf J}\big)-\tfrac{\omega\varepsilon_0}{2}\operatorname{Im}(\varepsilon_r)\,|\mathbf E|^2 .$$ Read it as: power leaving a point = power supplied by the source $-$ power absorbed there. So $\operatorname{Im}\varepsilon_r>0$ is absorption, as claimed. (The imaginary part of the same identity describes stored, reactive energy: the near field of an antenna and its Q.)

This is what objective functions are made of. Transmitted power is $\oint\langle\mathbf S\rangle\cdot d\mathbf A$ over a surface. Absorbed power is $\int\tfrac{\omega\varepsilon_0}{2}\operatorname{Im}(\varepsilon_r)|\mathbf E|^2\,d^3r$.

## Where we are, and what is missing

The scalar part of the note is one chain. Each step adds one assumption and gets one result.

1. **Linearity** (superposition holds). The law is $Lu=s$. Every solution is $u=u_p+u_0$: the field caused by the source plus free waves. So the law alone never gives a unique answer; extra conditions must choose.
2. **Why not $u=Ks$?** $L$ is the short rule, and $K=L^{-1}$ is its long consequence. $K$ only exists after the extra conditions, and its kernel is the Green's function.
3. **Locality** (a point only feels its neighbours). $L$ is a differential operator. On a grid, it is a sparse matrix, the same stencil in every row.
4. **Homogeneity** (the same here and there, now and later). $L$ commutes with shifts, so the plane waves $e^{i(\mathbf k\cdot\mathbf r-\omega t)}$ are its eigenvectors. On a plane wave, $L$ is just multiplication by its symbol $D(\mathbf k,\omega)$.
5. **Dispersion relation.** Free waves live only where $D(\mathbf k,\omega)=0$.
6. **Light in free space** (same speed $c$ for every colour and direction). $D=\omega^2-c^2|\mathbf k|^2$, which is the wave equation.
7. **Light in mediums.** The material has memory, so in time it is a convolution, and at one frequency it is one number $\varepsilon_r(\omega)$. The material also varies in space, so plane waves in space are lost. What is left is the Helmholtz equation, which we solve numerically.

Maxwell's equations do not start a new story. They run the same chain on a vector field. Here is each step side by side:

|Step|Scalar (1–7)|Maxwell|
|---|---|---|
|1. Linearity|$Lu=s$|the same; $u=\mathbf E$, $s=i\omega\mu_0\mathbf J$|
|3. Locality|derivatives only|the same: curls and divergences|
|4. Homogeneity|plane waves are eigenvectors; the symbol is a number|the same plane waves, but each carries a vector amplitude, so the symbol is a $3\times3$ matrix|
|5. Dispersion|$D=0$|$D\mathbf E_0=0$: again $\lvert\mathbf k\rvert=k_0$, plus $\mathbf E_0\perp\mathbf k$|
|6. Free space|the wave equation|each component obeys it (Check 1)|
|7. Mediums|$\varepsilon_r(\omega)$ from memory; Helmholtz in space|$\varepsilon_r(\omega)$ derived from bound currents; Helmholtz becomes curl-curl; the 2D $E_z$ case is exactly step 7|

What is genuinely new, with no scalar counterpart:

- **Polarization.** The amplitude has a direction, and the matrix symbol fixes it ($\mathbf E_0\perp\mathbf k$).
- **Constraints.** Equations without $\partial_t$ (the divergences) restrict the starting state; the curl equations keep them true.
- **A constraint on the source.** $\rho$ and $\mathbf J$ must together conserve charge, or the four equations have no solution at all.
- **Physical power.** $\tfrac12\operatorname{Re}(\mathbf E\times\overline{\mathbf H})$ and absorption $\propto\operatorname{Im}\varepsilon_r$. Objective functions are built from these.

Two ideas run through everything.

- **Symmetry decides the basis.** Whatever shifts $L$ commutes with, those exponentials diagonalize it. A material keeps the time symmetry but breaks the space symmetry. So we use exponentials in time (one frequency at a time) and numerics in space.
- **The null space does not go away by itself.** Free waves can always be added. The finite matrix in "Derivatives as matrices" hid this in its boundary rows.

### What is missing

Look at what the chain has actually produced so far.

For the free waves, it produced a _form_. In vacuum at one frequency, every free wave is $$u(\mathbf r)=\int A(\mathbf k)\,e^{i\mathbf k\cdot\mathbf r}\,d\mathbf k,\qquad A\text{ nonzero only on the sphere }|\mathbf k|=k_0 .$$ But $A$ is unknown. We never said which data fixes it, and we never computed a single free wave from data. That is why the free waves can still feel like playing with symbols.

For a field with a source, the chain produced even less. We know every solution is $u=u_p+u_0$, but we never found any $u_p$. And the law does not say which of the many solutions is the real one.

So two gaps are open:

1. What does a free wave actually look like, and which data picks one?
2. What field does a source make, and which of the many solutions does nature pick?

There is also a second aim from here on. A computer cannot hold a function, only a list of numbers. So every idea below appears twice: once in symbols (functions, derivatives, integrals) and once in numbers (vectors, stencils, sums). The two pictures often look unrelated, though they are the same idea. Saying exactly how they match is half the work of this part.

The rest of the note answers one question at a time, and each answer leaves the next question open:

1. **What does a free wave look like?** Data on one plane fixes it, and an FFT computes it.
2. **What field does a source make?** The Green's function. Computing it means dividing by the symbol, which fails exactly at the free waves.
3. **Which solution does nature pick?** The outgoing one. This also answers what a boundary condition really is.
4. **What does the Green's function say physically?** A small antenna, its radiated power, and reciprocity.
5. **What if a device sits in the way?** There is no formula anymore, so we must compute. There are two routes.
6. **A grid has an edge. How does it say "nothing comes in"?** The PML, which is a change of variables.
7. **Can we trust the grid, and what changes for a real 3D device?** The grid's own dispersion, the Yee grid, and time versus frequency.
8. **How do we improve a design?** The adjoint method.

One small example runs through steps 2 to 8: the 1D equation $u''+k_0^2u=s$ on a grid with $k_0h=0.5$, which is about 12.6 points per wavelength.

## What does a free wave look like?

We stay at one frequency, in vacuum, with the scalar equation $(\nabla^2+k_0^2)\,u=0$. To keep pictures simple, we work in 2D with $u(x,z)$. This is exactly the $E_z$ case from "When is the scalar Helmholtz equation right?".

### The data that picks one free wave

In time, we already know the answer: the initial conditions. A free wave of the 1D wave equation is fixed by $u$ and $\partial_tu$ at $t=0$ ("Light in 1D: the wave equation"). There were two roots, $\omega=\pm ck$, so we needed two pieces of data for each $k$.

Do the same in space. Let $z$ play the role of time, and suppose we know the field on the line $z=0$. Expand it in $x$: $$u(x,0)=\int a(k_x)\,e^{ik_xx}\,dk_x .$$ Each component $e^{ik_xx}$ has to grow into a free wave $e^{i(k_xx+k_zz)}$. The dispersion relation $k_x^2+k_z^2=k_0^2$ allows two roots: $$k_z=\pm\sqrt{k_0^2-k_x^2}.$$ Two roots again, so the line alone is not enough. Just as in time, we would also need $\partial_zu$ on the line.

Usually, however, we know something else instead. Suppose all sources lie below the line. Then above the line, waves can only travel _away_ from the sources, that is, upward ($e^{ik_{z}z-\omega t}$). This keeps one root for each $k_x$, and the field above is $$u(x,z)=\int a(k_x)\,e^{ik_xx}\,e^{ik_zz}\,dk_x,\qquad k_z=+\sqrt{k_0^2-k_x^2}\qquad(z\ge0).$$ This is called the **angular spectrum** of the field. So the field on one line, plus the statement "waves go away from the sources", fixes the free wave everywhere above. The unknown $A(\mathbf k)$ from "What is missing" is now known: for $|k_x|<k_0$, it is $a(k_x)$, placed on the upper half of the circle $|\mathbf k|=k_0$. (The components with $|k_x|>k_0$ have a complex $k_z$; they come next.)

The second statement is not part of the law. It is information about where the sources are. This is our first meeting with the radiation condition.

### Two kinds of components

- **$|k_x|<k_0$.** Here $k_z$ is real. The component is a plane wave, tilted from the $z$ axis by the angle $\theta$ with $\sin\theta=k_x/k_0$. Moving up by $z$ only multiplies it by $e^{ik_zz}$, which has modulus 1. It travels forever.
- **$|k_x|>k_0$.** Here $k_z$ is imaginary: $k_z=i\kappa$ with $\kappa=\sqrt{k_x^2-k_0^2}$. So $e^{ik_zz}=e^{-\kappa z}$: the component wiggles along $x$ but dies along $z$. This is the **evanescent wave** from "Back to $Lu=0$". (We take the root $+i\kappa$, which decays upward. The other root grows upward, and with no source above, nothing could feed it.)

The second case has a striking consequence. A component with $|k_x|>k_0$ has a period $2\pi/k_x$, which is shorter than the wavelength $\lambda=2\pi/k_0$. So any detail of $u(x,0)$ finer than about one wavelength is carried only by evanescent waves, and it dies within about a wavelength of the line. A lens or a camera far away never receives it. This is the **diffraction limit**, and it is nothing more than the dispersion relation (established).

A slit of width $w$ shows both kinds at once. Its spectrum $a(k_x)$ has a width of about $2\pi/w$. If $w\gg\lambda$, almost all of it is propagating and close to $k_x=0$, so the light goes nearly straight on and spreads by an angle of about $\lambda/w$. If $w<\lambda$, most of it is evanescent, and the little that survives spreads in all directions.

### The same thing in numbers

Now compute it. The computer stores the field on the line at $N$ points $x_j=jh$, $j=0,\dots,N-1$, as a vector $(u_0,\dots,u_{N-1})$. The Fourier integral becomes a finite sum, the discrete Fourier transform, which the FFT computes quickly. The whole solver is four steps:

1. Sample $u(x,0)$. This gives a vector of $N$ numbers.
2. Apply the FFT. This gives $N$ amplitudes $a_m$, one for each $k_m=2\pi m/(Nh)$, with $m=-N/2,\dots,N/2-1$.
3. Multiply each $a_m$ by $e^{ik_{z,m}z}$, where $k_{z,m}=\sqrt{k_0^2-k_m^2}$ (the $+i$ root when it is imaginary).
4. Apply the inverse FFT. This gives the field on the line at height $z$.

Notice what is missing: no derivative is ever taken. In the FFT basis, $\partial_x^2$ is just multiplication by $-k_m^2$, so the law splits into $N$ separate one-line problems. This is the Homogeneity section turned into an algorithm: find the basis in which $L$ is diagonal, and solving becomes bookkeeping. It is also exact for every $k_m$ the grid can store.

The grid adds two limits, and both matter later.

- **The shortest wiggle.** The largest $k_m$ is $\pi/h$. A vector of samples cannot hold a period shorter than $2h$. Propagating waves need $|k_x|\le k_0$, so $h\le\lambda/2$ keeps all of them. Finer detail is lost, but it was evanescent anyway.
- **A hidden boundary condition.** Every $e^{ik_mx_j}$ repeats after $N$ points. So the FFT treats the line as a ring: point $N-1$ is the neighbour of point $0$. The computed field is really the field of an endless row of copies of your data. As $z$ grows, the copies spread sideways and run into each other, and a wave leaving at the right edge comes back in at the left. Nobody asked for this; it came with the basis. The usual fix is to pad the data with zeros, so that the copies stay far apart over the range of $z$ you care about.

The second limit is the first example of a theme. A finite computation must say something about what lies outside its window. If we do not choose, the method chooses for us. Here the choice was "outside is a copy of inside". We meet this again in "What is a boundary condition, really?".

So the free waves are no longer just symbols. One line of data plus one direction fixes them, and four lines of code compute them. But where does the data on the line come from? From a source somewhere below it. That is the next question.

## Turning on a source: the Green's function

Now add the source back: $(\nabla^2+k_0^2)\,u=s$, still at one frequency and in vacuum. We want $u_p$, the field the source makes.

### Why one point source is enough

We want to solve $Lu=s$, not for one source, but for any source we might use: a dipole, a line current, an incoming beam. Superposition gives a shortcut.

On the grid, a source is a vector $\mathbf b$, and every vector is a sum of single-point sources: $\mathbf b=\sum_jb_j\,\mathbf e_j$, where $\mathbf e_j$ is 1 at point $j$ and 0 elsewhere. Suppose we know the response to each single point, i.e. the solution $\mathbf g_j$ of $A\mathbf g_j=\mathbf e_j$. Then by superposition, the response to $\mathbf b$ is $$\mathbf u=\sum_jb_j\,\mathbf g_j .$$ The $\mathbf g_j$ are exactly the columns of $A^{-1}$, as we saw in "Derivatives as matrices". So we solve for point sources once, and every other source is just a weighted sum.

In the continuous version, $\mathbf e_j$ becomes $\delta(x-x')$, and $s(x)=\int s(x')\,\delta(x-x')\,dx'$. The response to $\delta(x-x')$ is the **Green's function** $G(x,x')$: the field at $x$ due to a unit point source at $x'$. Then $$u(x)=\int G(x,x')\,s(x')\,dx' .$$ This is the $K$ from "Why not write the law as $u=Ks$?". On the grid, $G(x_i,x_j)$ is the entry $(A^{-1})_{ij}$.

Homogeneity saves even more work. In vacuum, moving the source only moves the response, so $G(x,x')=G(x-x')$. One function $G(x)$, the response to a point source at the origin, then answers every source at once: $$u(x)=\int G(x-x')\,s(x')\,dx' .$$ This is a convolution, and $G$ is the unit point source response $k$ of the [[convolution]] note. That is why we want $G$: it is the whole solution operator in one function.

### Dividing by the symbol

How do we find $G$? In vacuum, the trick from "The same thing in numbers" still works. In the $k$ picture, $L$ is multiplication by its symbol $k_0^2-|\mathbf k|^2$. So solving $LG=\delta$ should be division.

Take 1D. Write the point source as a sum of plane waves, all with the same amplitude: $$\delta(x)=\frac{1}{2\pi}\int e^{ikx}\,dk .$$ $L$ multiplies each component by $k_0^2-k^2$. To get $\delta$ out of $L$, $G$ must therefore divide each component by the same number: $$G(x)=\frac{1}{2\pi}\int\frac{e^{ikx}}{k_0^2-k^2}\,dk .$$ But the denominator is zero at $k=\pm k_0$, and the integral runs straight through these two points. As written, it has no value. In 3D, the bad set is the whole sphere $|\mathbf k|=k_0$.

The bad points are not an accident. They are exactly the free waves of the previous section: the components that $L$ multiplies by zero. If $u_p$ is a solution, so is $u_p+(\text{any free wave})$, because $L$ cannot see the free wave. In the $k$ picture, "$L$ cannot see it" means "$L$ multiplies it by zero", so dividing by $L$ is undefined exactly there. That is to say, three questions are really one:

- Which particular solution is the right one?
- Which free wave do we add?
- What do we do at the zeros of the symbol?

The same thing happens on a computer. Take the 5-point matrix from "Derivatives as matrices", now with the Helmholtz term and $k_0h=0.5$. After multiplying by $h^2$, each row is $(1,\,a,\,1)$ with $a=k_0^2h^2-2=-1.75$. Its eigenvalues are $$-3.48,\quad-2.75,\quad-1.75,\quad-0.75,\quad-0.018 .$$ The last one is almost zero: one standing wave (see the following sections) very nearly fits the box. Solving divides each eigen-component by its eigenvalue, so that component is amplified about $1/0.018\approx56$ times. A unit point source in the middle gives $|u|$ up to $18.9$. This is the tied rope of the next section, in numbers. At a slightly different $k_0$, this eigenvalue is exactly zero, and $A$ has no inverse at all.

To be honest, the grid made a choice here that the continuous problem never made: the rows at the two ends silently say $u=0$ just outside the box. A different choice gives different eigenvalues. So the real question is which choice is _right_, and that is a question about the world, not about the law.

## Which solution does nature pick?

### One dimension: shaking a rope

Picture a very long rope. Shake one point, $x=0$, up and down at frequency $\omega$, and keep shaking. After a while, the rope settles into a steady pattern at that frequency. That pattern is the 1D Green's function. (Small rope motions obey the 1D wave equation, so at one frequency they obey $u''+k_0^2u=\delta(x)$, with $k_0=\omega/v$ and $v$ the wave speed on the rope.)

Away from your hand, nothing pushes the rope, so the equation is source-free there, and on each side $$u=a\,e^{ik_0x}+b\,e^{-ik_0x}.$$ That is 4 unknowns: $a,b$ on the left and $a,b$ on the right. Your hand gives 2 conditions at $x=0$. The rope does not break, so $u$ is continuous. And your hand pushes there, so the rope has a kink: the slope $u'$ jumps (by 1 here, found by integrating the equation across $0$). So 2 unknowns are still free. These are the two free waves, one moving right and one moving left, as the null space predicts.

The equation cannot fix them, so ask the rope. With the time factor $e^{-i\omega t}$, $e^{ik_0x}$ moves right and $e^{-ik_0x}$ moves left. Waves run away from your hand in both directions, and nothing comes back from the far ends. So on the right we keep only $e^{ik_0x}$, and on the left only $e^{-ik_0x}$. These are the 2 missing conditions, and they give $$G(x)=\frac{e^{ik_0|x|}}{2ik_0}.$$ (Check: $G'$ jumps from $-\tfrac12$ to $+\tfrac12$ at $x=0$.)

Compare a rope tied at both ends. Now waves do come back, they pile up into standing waves, and at certain frequencies the rope resonates without limit. This is the nearly zero eigenvalue of the previous section. The equation is the same; only the far ends differ. So what fixes the free waves is not the law but what happens far away.

> [!NOTE] 
> Integrate $u''+k_0^2u=\delta(x)$ from $-\epsilon$ to $\epsilon$: $$u'(\epsilon)-u'(-\epsilon)+k_0^2\int_{-\epsilon}^{\epsilon}u\,dx=1 .$$ $u$ is continuous, so its integral goes to 0 as $\epsilon\to0$. So $u'$ jumps by exactly 1. Why is $u$ continuous? If $u$ jumped, $u'$ would contain a $\delta$, and $u''$ would contain $\delta'$. The right side has no $\delta'$. Picture: the hand makes a kink in the rope, not a break.

### The radiation condition

"Far from all sources, there are only outgoing waves" is the **radiation condition**. It is the extra condition for open space, just as initial conditions are the extra condition in time. In 3D, its precise form is the Sommerfeld condition, $r\,(\partial_ru-ik_0u)\to0$ as $r\to\infty$.

Mindful readers might ask why "outgoing" is the right choice, and not just a taste. There are two honest answers, and they agree.

- **Causality.** In time, the response must come after the source (see the causality note in [[convolution]]). An incoming wave would have to start at infinity before the source was switched on. At one frequency, the retarded response is exactly the outgoing one.
- **A little loss.** Every real medium absorbs a little: $\operatorname{Im}\varepsilon_r>0$. So replace $k_0$ by $\kappa=k_0\sqrt{1+i\eta}$ with a tiny $\eta>0$; then $\operatorname{Im}\kappa>0$. Now $e^{i\kappa x}$ decays to the right and blows up to the left, and $e^{-i\kappa x}$ does the opposite. A bounded solution must keep only the decaying wave on each side, which is exactly the outgoing choice. Then let $\eta\to0$. This is called the **limiting absorption principle** (established).

The second answer also repairs "Dividing by the symbol". With $\kappa$ in place of $k_0$, the denominator $\kappa^2-k^2$ is never zero for real $k$, so the integral for $G$ has a value. Evaluating it is a contour integral (keyword: residue theorem), and the result is $e^{i\kappa|x|}/(2i\kappa)$. As $\eta\to0$, this becomes the rope's answer. So a trace of loss tells us how to step around the zeros of the symbol.

So the null space is not removed by a mathematical trick. It is removed by physics we had already: causality, or equivalently a trace of loss.

### What is a boundary condition, really?

Now we can answer this in general.

The law is local: at each point, it only relates the field to its neighbours ("Locality"). Inside a region, the law ties every point to the next. At the edge of the region, however, the law needs neighbours that lie outside, and it has nothing to say about them. A boundary condition is this missing information. It is a statement about the world outside the region, written as a condition on the region's edge.

The grid makes this literal. Take the 1D Helmholtz rows on $N$ points, $$u_{i-1}+a\,u_i+u_{i+1}=h^2b_i,\qquad i=1,\dots,N .$$ The rows for $i=2,\dots,N-1$ use only points inside: these rows are the law. That is $N-2$ equations for $N$ unknowns, so 2 unknowns are left free. We can even see them. Without a source, the interior rows are solved by $u_i=q^{\,i}$ and by $u_i=q^{-i}$, where $q+1/q=-a$, i.e. $q=e^{iKh}$. These are the grid's right-going and left-going free waves. (On the grid, the wave number is $K$, slightly different from $k_0$; see "The grid has its own dispersion relation".) The rows for $i=1$ and $i=N$ need $u_0$ and $u_{N+1}$, which lie outside. Whatever we put there is the boundary condition, and it picks one point in the 2-dimensional null space. So $$\text{number of boundary conditions}=\text{dimension of the null space of the interior law}.$$ Here both are 2, the same count as on the rope. In 2D or 3D, the edge is a whole curve or surface of points, and we need one condition for each edge point.

Every extra condition in this note is a statement of this kind:

|Condition|Region|What it says about the outside|
|---|---|---|
|Initial conditions|$t>0$|the past: how the field was set up|
|Mirror, $u=0$ just outside|a box|a perfect metal wall: everything reflects|
|Periodic (the FFT ring)|one period|outside is a copy of inside|
|Radiation condition|all of space|empty space forever, with no sources in it: nothing comes in|

One warning follows. A wrong boundary condition does not produce an error message. It produces the correct solution of a different problem. The mirror box gives the field inside a metal box. The FFT ring gives the field of an endless array of copies. Only the radiation condition describes an antenna in open space.

### Three dimensions

For a point source in 3D, the field depends only on the distance $r$. Away from $r=0$, a short computation gives $\nabla^2u=\tfrac1r(ru)''$. So the equation becomes $(ru)''+k_0^2(ru)=0$, and $ru=e^{\pm ik_0r}$. The radiation condition keeps $e^{+ik_0r}$. Near $r=0$, $e^{ik_0r}\approx1$, and $\nabla^2\tfrac{1}{4\pi r}=-\delta$ (Coulomb's law) fixes the constant: $$G(\mathbf r)=-\frac{e^{ik_0r}}{4\pi r}.$$ This is a spherical wave going out. Its amplitude falls like $1/r$, so its power falls like $1/r^2$, and the total power through any large sphere stays the same, as energy conservation requires. (In 2D, $G$ is a Hankel function; see [[Green's function]].)

For the vector curl-curl equation, $G$ becomes a $3\times3$ matrix (the dyadic Green's function). It is built from the same $e^{ik_0r}/r$ and its derivatives. The next section shows a simpler way to handle vector sources.

## What the Green's function says: a small antenna

The Green's function was built from symbols. Here we check that it says things an engineer already knows. On the way, we find one new fact, reciprocity, which the adjoint method will need.

### The field of a short antenna

Take a short wire along $z$, of length $\ell\ll\lambda$, carrying a uniform current $I$ at frequency $\omega$. Charge piles up at the two ends, as charge conservation demands. This is a **Hertzian dipole**, the simplest antenna.

The source is now a vector, so the scalar $G$ does not apply directly. But a standard trick reduces the vector problem to the scalar one: the vector potential.

> [!NOTE] **Potentials at one frequency.** 
> Since $\nabla\cdot\mathbf B=0$, we can write $\mathbf B=\nabla\times\mathbf A$ (a field with no divergence is a curl; established). Faraday's law, $\nabla\times\mathbf E=i\omega\mathbf B$, then says $\nabla\times(\mathbf E-i\omega\mathbf A)=0$. A field with no curl is a gradient, so $\mathbf E=i\omega\mathbf A-\nabla\phi$ for some $\phi$. Put both into Ampère's law, $\nabla\times\mathbf B=\mu_0\mathbf J-i\omega\mu_0\varepsilon_0\mathbf E$, and use $\nabla\times\nabla\times\mathbf A=\nabla(\nabla\cdot\mathbf A)-\nabla^2\mathbf A$: $$\nabla(\nabla\cdot\mathbf A)-\nabla^2\mathbf A=\mu_0\mathbf J+k_0^2\mathbf A+i\omega\mu_0\varepsilon_0\nabla\phi .$$ The pair $(\mathbf A,\phi)$ is not unique. Replace $\mathbf A$ by $\mathbf A+\nabla\chi$ and $\phi$ by $\phi+i\omega\chi$, and neither $\mathbf B$ nor $\mathbf E$ changes. We spend this freedom on the choice $\nabla\cdot\mathbf A=i\omega\mu_0\varepsilon_0\,\phi$ (the **Lorenz gauge**). Then the two gradient terms cancel, and $$(\nabla^2+k_0^2)\,\mathbf A=-\mu_0\mathbf J .$$ Each component of $\mathbf A$ obeys the scalar Helmholtz equation. So the scalar $G$ solves the vector problem, one component at a time.

The current points only along $z$, so only $A_z$ is nonzero. Seen from far away ($r\gg\ell$), the wire is a point source of strength $I\ell$: $J_z=I\ell\,\delta(\mathbf r)$. With $G$ from "Three dimensions", $$A_z=-\mu_0I\ell\,G(\mathbf r)=\frac{\mu_0I\ell}{4\pi}\,\frac{e^{ik_0r}}{r}.$$

Far from the antenna ($k_0r\gg1$), the field is locally a plane wave moving outward, so $\mathbf E$ must be perpendicular to the direction $\hat{\mathbf r}$ (Check 2). The $-\nabla\phi$ term removes exactly the radial part of $i\omega\mathbf A$ there (established). So only the part of $\mathbf A$ perpendicular to $\hat{\mathbf r}$ radiates. In spherical coordinates, $\hat{\mathbf z}=\cos\theta\,\hat{\mathbf r}-\sin\theta\,\hat{\boldsymbol\theta}$, so that part is $A_\theta=-A_z\sin\theta$, and $$E_\theta=-i\omega\mu_0I\ell\,\sin\theta\,\frac{e^{ik_0r}}{4\pi r},\qquad H_\phi=\frac{E_\theta}{\eta_0},\qquad \eta_0=\sqrt{\mu_0/\varepsilon_0}\approx377\ \Omega .$$

Three things can be read off.

- **The pattern.** The factor $\sin\theta$ is largest broadside and zero along the wire. The scalar $G$ only gave $e^{ik_0r}/r$; the direction dependence comes from the vector nature. Only the part of the current perpendicular to your line of sight radiates toward you.
- **The power.** The time-averaged Poynting vector is $\langle S_r\rangle=|E_\theta|^2/(2\eta_0)$, which falls like $1/r^2$. Integrate it over a sphere, using $\int\sin^2\theta\,d\Omega=8\pi/3$ and $\omega\mu_0=k_0\eta_0$: $$P=\frac{\eta_0k_0^2I^2\ell^2}{12\pi}.$$ Define the radiation resistance by $P=\tfrac12I^2R_{\text{rad}}$. Then $$R_{\text{rad}}=\frac{2\pi}{3}\,\eta_0\Big(\frac{\ell}{\lambda}\Big)^2\approx790\,\Big(\frac{\ell}{\lambda}\Big)^2\ \Omega .$$ A wire of length $\lambda/10$ has about $8\ \Omega$. This is the textbook reason why short antennas are hard to match (established).
- **The near field.** The formula above keeps only the $1/r$ part. The exact fields also have $1/r^2$ and $1/r^3$ terms, which come from derivatives of $e^{ik_0r}/r$. Close to the antenna ($k_0r\ll1$) these dominate, and the field looks like that of a static dipole. These terms store energy that sloshes back and forth each cycle instead of leaving. This is the reactive power from the imaginary part of the identity in "Power".

So $G$ is not only a symbol. It is the $e^{ik_0r}/r$ that every antenna formula carries.

### Reciprocity

Look back at $A^{-1}$ in "Derivatives as matrices": $$A^{-1}=-\frac14\begin{pmatrix}3&2&1\\2&4&2\\1&2&3\end{pmatrix}.$$ It is symmetric. Entry $(1,3)$ is the field at point 1 due to a source at point 3. Entry $(3,1)$ is the field at point 3 due to a source at point 1. They are equal. If we swap the source and the probe, the reading does not change. This is **reciprocity**.

It is not a coincidence of this small matrix.

- **On the grid.** $D_2$ is symmetric, and the Helmholtz term is diagonal. So $A^T=A$, even when $\varepsilon_r$ is complex (lossy) and changes from cell to cell. The inverse of a symmetric matrix is symmetric, so $G(x_i,x_j)=G(x_j,x_i)$ for any structure.
- **In symbols.** The continuous version is integration by parts. Let $L=\nabla^2+k_0^2\varepsilon_r(\mathbf r)$, and take any two fields $u$ and $v$. The $\varepsilon_r$ terms cancel in $vLu-uLv$, and the divergence theorem gives (Green's second identity) $$\int\big(v\,Lu-u\,Lv\big)\,dV=\oint\big(v\,\partial_nu-u\,\partial_nv\big)\,dA .$$ Now let both $u$ and $v$ be outgoing. On a large sphere, both look like $e^{ik_0r}/r$, so $\partial_ru\approx ik_0u$ and $\partial_rv\approx ik_0v$. The two surface terms then cancel, and $\int v\,Lu=\int u\,Lv$. Finally, take $u=G(\cdot,\mathbf r_1)$ and $v=G(\cdot,\mathbf r_2)$, so that $Lu=\delta(\mathbf r-\mathbf r_1)$ and $Lv=\delta(\mathbf r-\mathbf r_2)$. The identity becomes $$G(\mathbf r_1,\mathbf r_2)=G(\mathbf r_2,\mathbf r_1).$$

Two details are worth noting.

- No complex conjugate appears anywhere. That is why reciprocity survives loss. $A$ is complex symmetric, $A^T=A$, but it is not Hermitian, $A^H\neq A$.
- Reciprocity fails when $\varepsilon_r$ (or $\mu_r$) is a matrix that is not symmetric, such as the magnetized ferrite from "In a medium". Isolators and circulators exist exactly to break it.

Engineers know reciprocity in other forms: an antenna has the same pattern when transmitting and receiving, and $S_{21}=S_{12}$ for a passive network without ferrites. Here it is simply a property of $G$. The adjoint section rests on it.

## Adding a structure: two routes

Everything so far used homogeneity. Now put a device in the way, so that $\varepsilon_r(\mathbf r)$ varies. Two things break at once.

- Plane waves are no longer eigenvectors of $L$ (recall the $m(x)\,\frac{d}{dx}$ example). So the trick of "The same thing in numbers", diagonalize with an FFT and divide, no longer works.
- $G(\mathbf r,\mathbf r')$ no longer depends only on $\mathbf r-\mathbf r'$, and there is no formula for it. A few shapes have series solutions, such as the sphere (Mie theory). A general device has none.

This is where computing stops being optional. There are two ways to go, and they are the $K$ and $L$ forms from "Why not write the law as $u=Ks$?".

### Route K: keep the vacuum Green's function

Use the rule from "What counts as a source?" in reverse: move the structure to the source side, $$\nabla^2u+k_0^2u=s-k_0^2(\varepsilon_r-1)\,u .$$ The left side is the vacuum operator, whose outgoing $G$ we know. So $$u=G*\big(s-k_0^2(\varepsilon_r-1)\,u\big).$$ This is an integral equation for $u$, called the Lippmann–Schwinger equation. The field appears on both sides, so it is not yet a solution, but it is a law that only involves the structure.

To compute with it, chop the structure into $M$ small cells of volume $V$, and let the unknowns be the field values $u_j$ in the cells. The equation at cell $i$ reads $$u_i+\sum_jk_0^2(\varepsilon_{r,j}-1)\,G(\mathbf r_i-\mathbf r_j)\,V\,u_j=(G*s)_i .$$ Every cell talks to every other cell through $G$, so the $M\times M$ matrix is dense. In exchange, the unknowns live only where $\varepsilon_r\neq1$, and the outgoing condition is automatic, because $G$ is outgoing. (The diagonal entry needs care, because $G$ blows up at $\mathbf r=\mathbf r'$.) For metal objects, the unknowns can even live only on the surface, as surface currents. This is the **method of moments** (MoM).

### Route L: keep the differential law on a grid

Put the whole space on a grid and keep $L$ as a stencil, as in "Derivatives as matrices". This gives a sparse system $A\mathbf x=\mathbf b$.

This route never computes $G$. In principle, $G$ is all the columns of $A^{-1}$. But we never need all of them; we need $A^{-1}\mathbf b$ for one $\mathbf b$. So we _solve_; we do not _invert_. The difference is huge. Take a 3D grid of $100^3=10^6$ points. $A$ has about $7\times10^6$ nonzero entries (7 per row for the 3D stencil). $A^{-1}$ has $10^{12}$ entries, about 16 terabytes as complex numbers. A solver gets $A^{-1}\mathbf b$ without ever forming $A^{-1}$. Either it factorizes $A=LU$ with Gaussian elimination that avoids filling in the zeros (keyword: sparse direct solver), or it iterates, using only products $A\mathbf v$ (keywords: Krylov methods, GMRES).

### So how do we solve with the Green's function?

Now the question has a precise answer.

- **Route K** uses $G$ as a tool. It is known by a formula and evaluated between every pair of cells.
- **Route L** never builds $G$. $G$ is implicitly $A^{-1}$, which is far too big to store, so we compute one response $A^{-1}\mathbf b$ at a time.
- **In both routes**, $G$ is the concept that says what a solution _is_: the outgoing response to point sources, added up.

Which route is better depends on the problem. A large metal antenna in empty space favours route K: the unknowns sit only on the metal surface, and no empty space needs a grid. A complicated dielectric volume favours route L: the volume has to be meshed anyway, and sparsity wins.

Route L, however, has an open problem that route K does not have. The grid ends somewhere, and the outgoing condition that $G$ carried is gone. Since most volume solvers (FDFD, FDTD, FEM) and the adjoint method below take route L, the note follows it from here.

## A grid has an edge: the PML

Recall derivatives as matrices in the Locality section. We knew that the second derivative on the grid is $$u''\approx\frac{u_{i-1}-2u_i+u_{i+1}}{h^2}.$$ We want the Helmholtz equation $u''+k_0^2\varepsilon_{r}u=b$, which results in
$$\frac{u_{i-1}+(k_0^2\varepsilon_{r}h^2-2)\,u_i+u_{i+1}}{h^2}=b_i .$$

So $-2$ is the second derivative, and $k_0^2\varepsilon_{r}h^2$ is the frequency. Look at the last row of the 1D Helmholtz matrix. Multiply it by $h^2$ and write $a=k_0^2h^2-2$ (vacuum, $\varepsilon_r=1$): $$u_{N-1}+a\,u_N+u_{N+1}=h^2b_N .$$ The row needs $u_{N+1}$, but $u_{N+1}$ is not an unknown; it lies outside the grid. This is exactly the situation of "What is a boundary condition, really?". Whatever we put in its place is a statement about the world beyond the edge. The real world has no edge, only the grid does. So "which solution?" has become "what do we pretend lies outside?". Here are three answers, each a small change to the same tridiagonal matrix.

### Attempt 1: a mirror

Set $u_0=u_{N+1}=0$, as the matrix in "Derivatives as matrices" silently did. With 5 unknowns and a point source in the middle: $$\begin{pmatrix}a&1&&&\\1&a&1&&\\&1&a&1&\\&&1&a&1\\&&&1&a\end{pmatrix}\begin{pmatrix}u_1\\u_2\\u_3\\u_4\\u_5\end{pmatrix}=\begin{pmatrix}0\\0\\1\\0\\0\end{pmatrix}.$$ This pretends the outside is a perfect metal wall (for $E_z$, tangential $\mathbf E=0$). It is the tied rope: waves bounce back, pile up into standing waves, and $A$ is singular at the box's resonant frequencies. With $k_0h=0.5$, the solution has $|u|=9.1,\ 16,\ 18.9,\ 16,\ 9.1$. The true outgoing wave has constant $|u|$, because in 1D a wave does not weaken as it travels. We get a large standing wave instead, because this box is close to a resonance (the eigenvalue $-0.018$ from "Dividing by the symbol").

### Attempt 2: a one-way edge

Outside the grid, where we pretend there is vacuum and no source. The device and the source stay inside the grid, so the inside can be as complicated as you like while the outside is simple. Beyond the edge, $\varepsilon_r=1$ and $s=0$. There the equation is $u''+k_0^2u=0$, the same everywhere. So the plane-wave method works again, and every solution outside has the form  
$$u=A\,e^{ik_0x}+B\,e^{-ik_0x}.$$  
"Nothing comes in from outside" sets $B=0$. So outside, $u=A\,e^{ik_0x}$.

The amplitude $A$ depends on everything inside: the device, the source, all of it. But the ratio of two neighbouring values does not depend on $A$:  
$$\frac{u_{N+1}}{u_N}=\frac{A,e^{ik_0(N+1)h}}{A,e^{ik_0Nh}}=e^{ik_0h}=p .$$  
The last row only needs this ratio. So the edge knows nothing about the device, yet it is correct for any device.

All in all, on the right, only the outgoing wave $e^{ik_0x}$ should exist. Then each step multiplies $u$ by the same number, so $u_{N+1}=p\,u_N$ with $p=e^{ik_0h}$. On the left, only $e^{-ik_0x}$ exists, and the same reasoning gives $u_0=p\,u_1$. Put these into the first and last rows: $$\begin{pmatrix}a{+}p&1&&&\\1&a&1&&\\&1&a&1&\\&&1&a&1\\&&&1&a{+}p\end{pmatrix}\begin{pmatrix}u_1\\u_2\\u_3\\u_4\\u_5\end{pmatrix}=\begin{pmatrix}0\\0\\1\\0\\0\end{pmatrix}.$$ Only the two corners changed. Now $|u|\approx1.03$ in every cell. This pretends the outside is empty space that goes on forever. (The exact answer _on the grid_ is $1.033$, not $1$. The difference is the grid's own dispersion; see "The grid has its own dispersion relation". If we use the grid's own wave number in $p$, the edge becomes exact.)

So in 1D the problem is solved. In 2D it is not. A wave $e^{i(k_xx+k_yy)}$ arriving at the right edge needs $u_{N+1}=e^{ik_xh}u_N$, and $k_x=k_0\cos\theta$ depends on the angle of arrival. A real field arrives at many angles at once, and one number in the corner cannot match all of them. Such a row is exact only at one angle. (Keywords: absorbing boundary conditions, Engquist–Majda, Mur.)

### The idea: absorb the wave before it reaches the edge

Keep the mirror, but make sure nothing arrives at it. So we want a layer around the grid that does two things:

1. It absorbs, so the wave is nearly zero when it reaches the mirror.
2. It does not reflect at its entrance.

An ordinary lossy material, $\operatorname{Im}\varepsilon_r>0$, does the first but not the second. It changes the wave's impedance, so the entrance itself reflects.

### The trick: a complex change of variables

Recall what loss _does_ to a wave. In "The radiation condition", loss made $k$ complex, and $e^{i\kappa x}$ decayed. We can get the same decay by making $x$ complex instead of $k$. Inside the layer, replace the real coordinate $x$ by $$\tilde x=x+\frac{i}{\omega}\int_0^x\sigma(x')\,dx',$$ where $\sigma(x)$ is a function we choose: zero outside the layer, positive inside. It is not a property of any real material. Its units are 1/s, so $\sigma/\omega$ is a pure number. (The name comes from conductivity: Bérenger's original PML was a made-up material with electric and magnetic conductivity.) Then $d\tilde x/dx=s(x)$ with $$s=1+\frac{i\sigma}{\omega}.$$

Now demand that the _same_ law hold in $\tilde x$: $\frac{d^2u}{d\tilde x^2}+k_0^2u=0$. Its outgoing solution is $e^{ik_0\tilde x}$, and $$e^{ik_0\tilde x}=e^{ik_0x}\cdot e^{-\frac{k_0}{\omega}\int_0^x\sigma\,dx'} .$$ The first factor has modulus 1. In the second, $k_0$, $\omega$ and $\sigma$ are positive, so the wave decays inside the layer. That is requirement 1.

Requirement 2 comes for free. At the entrance, $\sigma=0$, so $\tilde x=x$. The single function $e^{ik_0\tilde x}$ solves the law on both sides of the entrance and joins smoothly there. No reflected wave is needed, so none appears. This is the continuous statement, though. On a grid there is a small reflection, which we will measure. (Established: Bérenger 1994. Chew and Weedon showed in the same year that it is a coordinate stretch.)

### What a change of variables means on a grid

Inside the layer, the law is written in $\tilde x$. But the grid lives on the real $x$: its points $x_i=ih$ are equally spaced. How do we write a law in $\tilde x$ on a grid in $x$?

In symbols, the answer is the chain rule. Since $d/d\tilde x=(1/s)\,d/dx$, the law becomes $$\frac{1}{s}\frac{d}{dx}\Big(\frac{1}{s}\frac{du}{dx}\Big)+k_0^2\varepsilon_r\,u=0 .$$ (Here $s$ is the stretch factor, not the source. The source sits outside the layer, so I drop it.)

In numbers, the answer is simpler, and it is worth seeing directly. A finite difference is a difference of values divided by a distance. The values do not care which coordinate we use: $u_i$ is the field at point $i$, whatever we call that point. Only the distances change. Measured in $\tilde x$, the grid points are no longer equally spaced. Between points $i$ and $i+1$, the distance is $$\tilde x_{i+1}-\tilde x_i=\int_{x_i}^{x_{i+1}}s\,dx\approx s_{i+1/2}\,h,$$ a complex length. So the first difference that sits between $i$ and $i+1$ is $$\frac{u_{i+1}-u_i}{s_{i+1/2}\,h}.$$ The second difference at $i$ is the difference of two neighbouring first differences, divided by the $\tilde x$-length of cell $i$, which is $s_i\,h$: $$\frac{1}{s_i h}\Big[\frac{u_{i+1}-u_i}{s_{i+1/2}\,h}-\frac{u_i-u_{i-1}}{s_{i-1/2}\,h}\Big]+k_0^2\varepsilon_{r,i}\,u_i=b_i .$$ Multiply by $h^2$ and collect the terms. Row $i$ has three entries: $$\ell_i=\frac{1}{s_i\,s_{i-1/2}},\qquad r_i=\frac{1}{s_i\,s_{i+1/2}},\qquad d_i=k_0^2h^2\varepsilon_{r,i}-\ell_i-r_i .$$ With $s=1$, this is $(1,\,a,\,1)$ again.

In short, in symbols a change of variables rewrites derivatives by the chain rule. On a grid, it keeps the values and changes the lengths the stencil divides by. **The PML is a grid whose cells have the complex length $s\,h$.**

> [!NOTE] **Change of variables, in numbers.** 
> The same picture works for integrals. In symbols, $\int f\,d\tilde x=\int f\,\frac{d\tilde x}{dx}\,dx$. On a grid, an integral is a sum of value times cell length, $\sum_if_i\,\Delta\tilde x_i$, and $\Delta\tilde x_i=s_i\,h$. So the Jacobian $d\tilde x/dx$ is just the ratio "length of cell $i$ in the new coordinate, over its length in the old one". In a derivative the same ratio sits in the denominator, which is the $1/s$ of the chain rule. Nothing else changes: the samples stay where they are.
>
> A complex length is not a longer cell. Write $\Delta\tilde x_i=s_ih=h+i\,\frac{\sigma_i}{\omega}h$. The real part is exactly $h$, so $\operatorname{Re}\tilde x_i=x_i$: every sample keeps its real position, and $u_i$ is the field at $x_i$. The two parts do different jobs on a wave $e^{ik_0\tilde x}$. The real part turns the phase by $k_0h$ per cell, as in vacuum, so the wavelength does not change. The imaginary part multiplies the amplitude by $e^{-k_0h\sigma_i/\omega}$ per cell. So $|s_i|>1$ does not mean the cell is longer in space. It means the wave also loses amplitude there. Example: with $\sigma/\omega=1$ and $k_0h=0.5$, $s=1+i$. Each cell turns the phase by 0.5 rad and multiplies the amplitude by $e^{-0.5}\approx0.61$.
> 
> A complex length is not a longer cell, though. Write $\Delta\tilde x_i=s_ih=h+i\,\frac{\sigma_i}{\omega}h$. The real part is exactly $h$, so $\operatorname{Re}\tilde x_i=x_i$: every sample keeps its real position, and $u_i$ is the field at $x_i$. The two parts do different jobs on a wave $e^{ik_0\tilde x}$. The real part turns the phase by $k_0h$ per cell, as in vacuum, so the wavelength does not change. The imaginary part multiplies the amplitude by $e^{-k_0h\sigma_i/\omega}$ per cell. So $|s_i|>1$ does not mean the cell is longer in space. It means the wave also loses amplitude there. Example: with $\sigma/\omega=1$ and $k_0h=0.5$, $s=1+i$. Each cell turns the phase by 0.5 rad and multiplies the amplitude by $e^{-0.5}\approx0.61$.


> [!NOTE] **A trap.** 
> Do not simplify the stretched law to $\frac{1}{s^2}u''$. Because $s$ varies with $x$, $$\frac1s\Big(\frac1s u'\Big)'=\frac{1}{s^2}u''-\frac{s'}{s^3}u' .$$ Dropping the second term changes the law, and the layer then reflects wherever $s$ changes. The staggered form above keeps that term automatically, because each difference divides by its own local length. Staggering is also the idea of the Yee grid (see "The Yee grid"): $\frac1s\frac{du}{dx}$ at half points plays the role of $H$, and $u$ at whole points plays the role of $E$.

### The matrix, and the three attempts compared

Here is the form with 2 PML cells on each side of 3 physical cells: $$\begin{pmatrix} d_1&r_1&&&&&\\ \ell_2&d_2&r_2&&&&\\ &1&a&1&&&\\ &&1&a&1&&\\ &&&1&a&1&\\ &&&&\ell_6&d_6&r_6\\ &&&&&\ell_7&d_7 \end{pmatrix}$$ The corners keep the mirror ($u_0=u_8=0$). Compare the three attempts:

- **Mirror:** every $s=1$. No extra rows; plain $a$ in the corners.
- **One-way:** no extra rows. Only the corners change, to $a+p$.
- **PML:** extra rows at each end, where $(1,\,a,\,1)$ becomes $(\ell_i,\,d_i,\,r_i)$. The corners stay a mirror.

So all three are the same tridiagonal matrix, and they differ only in the rows near the ends. The PML does not fix the boundary row. It adds a buffer of rows so that the boundary row stops mattering.

Two quick checks for code. Every row except the two corner rows sums to $k_0^2h^2\varepsilon_{r,i}$, since $d_i=k_0^2h^2\varepsilon_{r,i}-\ell_i-r_i$. And deep in the layer the off-diagonal entries shrink, about like $1/s^2$: a stretched cell acts like a longer cell, and longer cells couple more weakly.

### What it does, in numbers

The “no reflection” argument holds for the exact continuous law. A real PML loses that guarantee in two ways: the layer is finite, and the grid only approximates the law. The following measures both effects.

Take the 5 physical cells from before, add 6 PML cells on each side (17 unknowns), and keep $k_0h=0.5$. Ramp $\sigma$ up smoothly, $\sigma/\omega=4\,(\text{depth}/\text{thickness})^3$. Then $|u|$ across the 17 cells is $$0.11,\ 0.38,\ 0.73,\ 0.95,\ 1.02,\ 1.03\ \big|\ \underbrace{1.033\ \ 1.033\ \ 1.033\ \ 1.033\ \ 1.033}_{\text{physical}}\ \big|\ 1.03,\ 1.02,\ 0.95,\ 0.73,\ 0.38,\ 0.11 .$$ The wave fades in the layer and hits the mirror at about 11% amplitude. In the physical cells, the answer matches the exact grid solution to $7\times10^{-4}$. (For a 1D test: with the right-hand side used above, the grid's exact outgoing solution is $u_j=q^{|j|}/(a+2q)$, with $q=e^{iKh}$ and $\cos Kh=1-\tfrac12k_0^2h^2$.)

The error in the physical cells, for different thicknesses and peak values of $\sigma/\omega$:

|PML cells|peak 1|peak 2|peak 4|peak 8|peak 30|
|---|---|---|---|---|---|
|3|0.41|0.13|0.019|0.016|0.10|
|6|0.28|0.033|7e-4|5e-4|6e-3|
|20|7e-3|3e-5|5e-6|1e-5|4e-5|

There are two ways to fail. If $\sigma$ is too weak, the wave survives the trip to the mirror and back. If $\sigma$ is too strong, the ramp is too steep for the grid to resolve, and the entrance reflects after all. A thicker layer fixes both, at the price of more unknowns. In practice people use about 10 to 20 cells.

### In 2D and 3D

Stretch only the coordinate normal to each face: $x$ in the layers on the left and right, $y$ at the top and bottom, and both in the corners. A wave $e^{i(k_xx+k_yy)}$ entering the right layer becomes $$e^{i(k_x\tilde x+k_yy)}=e^{i(k_xx+k_yy)}\,e^{-\frac{k_x}{\omega}\int\sigma\,dx}.$$ Every wave with $k_x>0$ decays, at every angle. This is what the one-way edge could not do.

To be honest about its limits:

- **Grazing waves.** The decay rate is set by $k_x=k_0\cos\theta$, not by $k_0$. A wave travelling almost parallel to the layer ($k_x\to0$) decays slowly.
- **Evanescent waves.** If $k_x=i\kappa$, then $e^{ik_x\tilde x}=e^{-\kappa x}\,e^{-i(\kappa/\omega)\int\sigma}$. The stretch only adds a phase, not extra decay. Variants such as the CFS-PML fix this.
- **The mirror is still there.** The reflection after a round trip is small, not zero.
- **Time domain.** $s=1+i\sigma/\omega$ depends on $\omega$, so in time, dividing by $s$ is a convolution, exactly like the memory of glass in "Nonlocal in time". FDTD codes handle it with extra auxiliary fields (keyword: convolutional PML, CPML).

### Reciprocity survives the PML

The matrix is still sparse, now complex, and invertible: the layer absorbs, so the null space is gone. Where $s=1$, it is symmetric as before. In the layer it is not, as written, because each row divides by its own cell length $s_ih$. Here is the right end of the 7-cell example (rows and columns 5–7; cell 5 is physical, and $s_{11/2}=1$ at the entrance): $$\begin{pmatrix} a & 1 & 0 \\ \dfrac{1}{s_6} & d_{6} & \dfrac{1}{s_6,s_{13/2}}\\ 0 & \dfrac{1}{s_7,s_{13/2}} & d_7 \end{pmatrix} \;\xrightarrow{\ \text{row }i\,\times\,s_i}\; \begin{pmatrix} a & 1 & 0\\ 1 & s_6d_6 & \dfrac{1}{s_{13/2}}\\ 0 & \dfrac{1}{s_{13/2}} & s_7d_7 \end{pmatrix}$$ On the left, the entries across the diagonal differ: $1$ against $1/s_6$, and $\frac{1}{s_6s_{13/2}}$ against $\frac{1}{s_7s_{13/2}}$. Multiplying row $i$ by $s_i$, and $b_i$ too, removes each row’s own cell length. Scaling an equation by a nonzero number does not change its solution, and $s_i\neq0$ because $\operatorname{Re}s_i=1$. The scaling also leaves the physical rows ($s_i=1$) and the right-hand side ($b_i=0$ in the layer) unchanged. On the right, each pair matches, and $A^T=A$. So the grid with a PML is still reciprocal, as “Reciprocity” requires of any honest model of open space. It is complex symmetric but not Hermitian, and the adjoint method uses exactly this symmetry.

So both routes now have a complete answer to "which solution?". Route K builds the outgoing condition into $G$. Route L builds it into extra rows at the ends of the matrix: a PML.

## From the 1D toy to a real solver

The 1D toy now gives numbers we can trust. A real device, however, is 3D, the field is a vector, and the answer is often wanted over a band of frequencies. So two questions remain before the toy becomes a solver. How do we know the grid's numbers are right? And what changes for a real device? Each answer below reuses an idea from earlier in the note.

### The grid has its own dispersion relation

The grid is shift-invariant too, but only by whole cells. So the Homogeneity section applies to the grid itself: $e^{ikx_i}$ is an eigenvector of the stencil. Its eigenvalue is $$\frac{e^{ikh}-2+e^{-ikh}}{h^2}=-\frac{4}{h^2}\sin^2\frac{kh}{2}\;\approx\;-k^2+\frac{k^4h^2}{12}.$$ This is the symbol of the grid. Compare it with the true symbol $-k^2$. For small $kh$ they agree, but not exactly. The extra term $k^4h^2/12$ is the first Taylor term that the 3-point stencil threw away (see "Locality").

So the grid's dispersion relation, $\tfrac{4}{h^2}\sin^2\tfrac{Kh}{2}=k_0^2$, gives a wave number $K$ slightly larger than $k_0$. So grid waves are slightly too slow. In the running example, $k_0h=0.5$ gives $Kh=0.5054$, about 1.1% too large. This is the $1.033$ instead of $1$ in "Attempt 2": on the grid, the amplitude of the outgoing wave is $1/(2\sin Kh)$, while the continuous one is $1/(2k_0h)$.

The relative error in $K$ is about $(k_0h)^2/24$, and it adds up over distance as a phase error. With 20 points per wavelength, $k_0h\approx0.31$, and the phase error is about 0.26 rad after 10 wavelengths. With 10 points per wavelength, it is about 1.0 rad. This is **numerical dispersion**. It is why people quote "10 to 20 points per wavelength", and why large domains need more. When $k_0h>2$, the equation has no real $K$ at all, and the grid cannot carry the wave.

Why did the FFT method of "The same thing in numbers" need only 2 points per wavelength? The FFT uses the exact symbol $-k^2$. Written back as a matrix on the grid, that exact symbol couples every point to every other point: the matrix is dense. The stencil is local, so it is sparse, but a 3-point Taylor expansion is only accurate for small $kh$. So there is a trade-off. Local means sparse with an approximate symbol, and global means dense with an exact symbol. (Keyword: spectral methods. Higher-order stencils sit in between.)

### The Yee grid

The field is a vector, and grid solvers for Maxwell's equations almost all use the Yee grid. Its idea is to place each field component at a different point of the cell. Take the 2D $E_z$ case: put $E_z$ at the grid nodes $(i,j)$, $H_x$ at the edge midpoints $(i,j+\tfrac12)$, and $H_y$ at $(i+\tfrac12,j)$. Then every derivative in the curl equations is a centred difference over one cell. For example, Faraday's law gives $$H_x(i,j+\tfrac12)=\frac{E_z(i,j+1)-E_z(i,j)}{i\omega\mu_0h}.$$ Put the $H$ values into Ampère's law, and the 2D $E_z$ equation becomes the 5-point stencil $$\frac{E_{i+1,j}+E_{i-1,j}+E_{i,j+1}+E_{i,j-1}-4E_{i,j}}{h^2}+k_0^2\varepsilon_{r,ij}E_{i,j}=-i\omega\mu_0J_{i,j}.$$ This is the 1D matrix from "Derivatives as matrices", applied once in $x$ and once in $y$. With PML rows on all four sides, it is the whole 2D $E_z$ FDFD solver.

Why the staggering? It keeps "the curl of a gradient is zero" _exactly_ on the grid. On the grid, a gradient lives on edges: it is the difference of the two end values. A curl lives on faces: it is the sum of the edge values around the face, with signs (Stokes' theorem on one cell). Go around a face, and every node value is added once and subtracted once, so the result is exactly zero. This is what "The divergence equations are constraints" needs on the grid. The curl equations then keep the discrete divergence fixed, and the false (spurious) modes of the low-frequency note cannot appear. Edge elements in FEM are built on the same idea.

> [!NOTE] **Half indices in code.** 
> The $\tfrac12$ lives in the convention, not in the array index. Each component gets its own integer array with a fixed offset: `Ez[i][j]` stores $E_z(i,j)$, `Hx[i][j]` stores $H_x(i,j+\tfrac12)$, and `Hy[i][j]` stores $H_y(i+\tfrac12,j)$. Faraday's law then reads `Hx[i][j] = (Ez[i][j+1] - Ez[i][j]) / ($1j*omega*mu0*h$)`. (`Hx` has one fewer entry than `Ez` along $j$, since there is one fewer midpoint than nodes.) The position still matters wherever a value is evaluated at a point. Material and PML profiles must be sampled where the quantity lives: the $s_{i+1/2}$ in $\frac{u_{i+1}-u_i}{s_{i+1/2}h}$ is $\sigma$ at $x=(i+\tfrac12)h$, and sampling it at $x_i$ makes the layer lopsided and more reflective. A probe on $H_x$ reads at $(i,j+\tfrac12)$, not at node $(i,j)$. A plot of `Hx[i][j]` belongs at $y=(j+\tfrac12)h$, or it looks shifted by half a cell. Off-by-one errors in these offsets are the most common Yee bug.

### Frequency or time?

So far we solved at one frequency: one large sparse system $A\mathbf x=\mathbf b$. The other option is to step the two curl equations forward in time on the Yee grid. This is FDTD.

FDTD works like stepping Newton's law in a physics class. $\mathbf H$ changes at a rate set by $\nabla\times\mathbf E$, and $\mathbf E$ changes at a rate set by $\nabla\times\mathbf H$, just as $x$ changes by $v$ and $v$ changes by $a$. Each new value is a formula of values already known. In 1D, with $c=1$ (the sign depends on which pair of components we use): $$E_i^{\,n+1}=E_i^{\,n}+\frac{\Delta t}{h}\Big(H_{i+1/2}^{\,n+1/2}-H_{i-1/2}^{\,n+1/2}\Big).$$ Everything on the right is known, so there is nothing to solve. Such a scheme is called **explicit**. Note the half steps in time: $\mathbf E$ lives at whole steps and $\mathbf H$ at half steps. This is the staggering of the Yee grid, now in time, and it is called **leapfrog** (or the Verlet method). Updating both at the same instant (forward Euler) adds a little energy to every oscillation at each step, and the run slowly blows up. Leapfrog keeps the energy bounded (established).

- **At one frequency**, a dispersive material is just a number $\varepsilon_r(\omega)$ ("Nonlocal in time"). The price is one large sparse solve per frequency. But once $A$ is factorized, every new source costs almost nothing, which the adjoint will exploit.
- **In time**, there is no linear solve at all, only cheap local updates, so memory grows only like the number of cells. One short pulse covers a whole band of frequencies. The price is that the material's memory must be computed in time. The usual way is to step the Lorentz oscillator equation alongside the fields (the auxiliary differential equation method). The PML also becomes a convolution in time (CPML, from "In 2D and 3D").

Explicit stepping has one more cost, and it follows from locality. Each update uses only neighbouring cells, so a disturbance moves at most one cell per step. If the real wave moves farther in one step ($c\,\Delta t>h$ in 1D), the scheme cannot follow it, and errors grow exponentially. Put a plane wave $e^{i(\mathbf k\cdot\mathbf x-\omega t)}$ into the update, i.e. take the symbol again. This gives the grid's dispersion relation in space and time: $$\sin\frac{\omega\Delta t}{2}=\frac{c\,\Delta t}{h}\sqrt{\textstyle\sum_{j=1}^{d}\sin^2\frac{k_jh}{2}} .$$ The left side is at most 1 for real $\omega$, and the root reaches $\sqrt d$ for the checkerboard wave. So $$c\,\Delta t\le\frac{h}{\sqrt d}$$ in $d$ dimensions. Otherwise some $\omega$ is complex, and that mode grows without bound. This is the **CFL condition** (Courant, Friedrichs and Lewy, 1928). A single spring shows the same limit: leapfrog needs $\omega\,\Delta t\le2$. An **implicit** scheme (Crank–Nicolson, ADI-FDTD) also puts new-time values on the right side. It removes the CFL limit, but each step then couples all the new unknowns and needs a solve. Standard FDTD stays explicit, because avoiding solves is its whole point.

Why, then, is FDTD so common? Mostly because of large 3D problems. A direct solve of $A\mathbf x=\mathbf b$ is Gaussian elimination. Eliminating one unknown couples all its neighbours to each other, so new nonzeros appear (**fill-in**). In 1D, the matrix stays tridiagonal. In 3D, cut an $n\times n\times n$ grid in half. The cut is a plane of $n^2$ cells, and eliminating one side couples all of them: a dense $n^2\times n^2$ block. Even with the best ordering (nested dissection), memory grows like $N^{4/3}$ and work like $N^2$, with $N=n^3$ (established). For $N=10^6$, that is about $10^8$ stored values and $10^{12}$ operations, while one FDTD step touches each of the $10^6$ cells once. Iterative solvers avoid fill-in, but they converge slowly here. The Helmholtz matrix has eigenvalues of both signs, some near zero (the $-0.018$ of "Dividing by the symbol"), and this gets worse at higher frequency (established). In 2D, the cut is only a line, so a direct solve is fast. (My judgment: FDFD with a direct solver for 2D and small 3D devices; FDTD for large 3D devices and broadband runs.)

The two domains are, in some sense, closer than they look. Run FDTD with a single-frequency source until the field settles. The final field is the FDFD solution. So time stepping acts as an iterative solver for $A\mathbf x=\mathbf b$, with time as the iteration count (established). This is also why a high-Q resonator is slow in FDTD: the "iteration" needs about $Q$ periods to settle.

### The solvers, as four choices

Every solver turns the law into a finite problem. Strip away the brand names, and each one makes four choices, all of which we have now met.

1. **Form.** Use $L$ (sparse, needs a PML) or $K$ (dense, outgoing waves built in, unknowns only on the structure). See "Adding a structure: two routes".
2. **Domain.** Work at one frequency and solve $A\mathbf x=\mathbf b$, or work in time and step the curl equations forward.
3. **Basis.** Decide how a function becomes a vector: values on a grid, small local patches (elements), or global functions (plane waves, modes). A grid is simple, but a curved surface becomes a staircase. Elements follow curved surfaces and can be small where detail is needed, at the cost of a more complex code. Global functions are exact for their own symbol but give dense matrices, as the FFT showed.
4. **Solve.** At one frequency, solve $A\mathbf x=\mathbf b$ directly (LU factorization, reusable for every new right-hand side) or iteratively (GMRES, BiCGSTAB; needs a good preconditioner, and every new right-hand side, the adjoint included, is a new run). In time, step explicitly (no solve, but the CFL limit) or implicitly (no CFL limit, but a solve per step).

The common solvers are just different picks:

|Solver|Form|Domain|Basis|Solve|Typical use (my judgment; they all overlap)|
|---|---|---|---|---|---|
|FDFD|$L$|frequency|grid (Yee)|direct (2D), iterative (large 3D)|small devices, inverse design|
|FDTD|$L$|time|grid (Yee)|explicit|broadband runs, nanophotonics (Meep, Lumerical)|
|FIT|$L$|time (or frequency)|grid, integral form|explicit (in time)|CST|
|FEM|$L$|frequency, mostly|elements (edge elements)|direct or iterative|curved, multiscale geometry (HFSS, COMSOL)|
|MoM / BEM|$K$|frequency|patches on surfaces|direct (small), iterative with fast multipole (large)|metal antennas, scattering in open space (FEKO)|
|Modal / RCWA|mixed|frequency|modes, plane waves|small dense solves|periodic and layered structures|

### How to trust a result

Each check below is a principle from earlier in this note.

- **Convergence.** Halve $h$. A second-order scheme should cut the error by about 4 ("The grid has its own dispersion relation"). The phase error grows with distance, about $k_0L\,(k_0h)^2/24$ after a length $L$. So check convergence at the far side of the domain, not only near the source. Large domains need more points per wavelength.
- **Energy balance.** Power from the source = power out through a surface + power absorbed ("Power"). In a lossless structure, all absorption should happen in the PML.
- **Reciprocity.** Swap the source and the probe. The answer should not change, because $A^T=A$ ("Reciprocity").
- **The PML.** Make it thicker. The answer should not change ("What it does, in numbers").
- **Settling (time domain).** Run longer. A single-frequency answer should not change once the field has settled.
- **A known case.** Compare with an analytic Green's function, such as $e^{ik_0|x|}/(2ik_0)$ in 1D, or with the dipole's radiation resistance.

## The adjoint: every gradient for one extra solve

We can now compute the field for one design. The "What counts as a source?" note promised the next step: designing a device means changing $L$, not $s$, and asking how the solution responds. So the last question is: how should we change the design to make it better?

### The problem

Divide the design region into $N$ cells, and let the design parameters be their permittivities $p_i=\varepsilon_{r,i}$. $N$ is easily $10^4$ to $10^6$. The forward problem is the FDFD system $$A(\mathbf p)\,\mathbf x=\mathbf b,\qquad A=D_2+\mathrm{diag}\big(k_0^2\varepsilon_{r,1},\dots\big)\quad(\text{plus PML rows}).$$ We judge a design by one real number $J(\mathbf x)$, the **objective**: for example, the power reaching a probe. To improve the design, we want the gradient $dJ/dp_i$ for every $i$. The naive way changes one $p_i$, solves again, and sees how $J$ changes. That costs $N+1$ solves per gradient, which is hopeless for $N=10^5$. The adjoint method needs two. It rests on one identity, and it needs no symmetry, only a linear system we can solve.

### One measurement, many sources

Let $A\mathbf x=\mathbf b$ with any invertible $A$. We measure the response with a fixed set of weights $\mathbf g$: the reading is $M=\mathbf g^T\mathbf x$. For example, $\mathbf g=\mathbf e_m$ reads the field at cell $m$. Then $$M=\mathbf g^T\big(A^{-1}\mathbf b\big)=\big(A^{-T}\mathbf g\big)^T\mathbf b .$$ The left side says: push the source through the system, then measure. The right side says: pull the measurement back through the system, then overlap it with the source. Both give the same number. Define the **adjoint field** $\boldsymbol\lambda$ by $$A^T\boldsymbol\lambda=\mathbf g,\qquad\text{so that}\qquad M=\boldsymbol\lambda^T\mathbf b\ \text{ for every }\mathbf b .$$ One solve now answers every source at once: $\lambda_j$ is the reading per unit source at cell $j$. So $\boldsymbol\lambda$ is the measurement pulled back into the device, a map of how much each point matters to the reading.

A $2\times2$ check, deliberately not symmetric: $$A=\begin{pmatrix}2&1\\0&3\end{pmatrix},\qquad \mathbf g=\begin{pmatrix}0\\1\end{pmatrix}.$$ Forward, one source at a time: $A^{-1}\mathbf e_1=(\tfrac12,0)$ reads $0$, and $A^{-1}\mathbf e_2=(-\tfrac16,\tfrac13)$ reads $\tfrac13$. That is two solves. Backward: $A^T\boldsymbol\lambda=\mathbf g$ gives $\boldsymbol\lambda=(0,\tfrac13)$, which contains both readings. That is one solve.

### A design change is a source

Now change a design parameter $p_i$ by a small $\delta$. Then $A$ changes by $\delta A=\frac{\partial A}{\partial p_i}\delta$. Use the rule from "What counts as a source?" once more and move the change to the right side. To first order, $$A\,(\mathbf x+\delta\mathbf x)=\mathbf b-\delta A\,\mathbf x .$$ So the change acts as a new source $-\delta A\,\mathbf x$. It sits where $A$ changed, and the forward field there sets its strength. Physically, the extra permittivity is extra bound charge, and the forward field shakes it. By the identity above, the reading changes by $$\delta M=\boldsymbol\lambda^T\big(-\delta A\,\mathbf x\big)=-\delta\,\boldsymbol\lambda^T\frac{\partial A}{\partial p_i}\mathbf x .$$ The same $\boldsymbol\lambda$ serves every $i$. So one forward solve gives $\mathbf x$, one adjoint solve gives $\boldsymbol\lambda$, and each gradient entry is a cheap local product. In FDFD, $\partial A/\partial\varepsilon_{r,i}=k_0^2\,\mathbf e_i\mathbf e_i^T$, so $$\frac{dM}{d\varepsilon_{r,i}}=-k_0^2\,\lambda_i\,x_i .$$

With reciprocity, this has a physical picture. Then $A^T=A$ ("Reciprocity"), so the adjoint system is the same device, and $\boldsymbol\lambda$ is a real wave: the field of a source placed at the probe. The gradient at cell $i$ is the product of two waves meeting there. $x_i$ says how strongly the forward field shakes cell $i$. $\lambda_i$ says how well cell $i$ is heard at the probe. A cell matters when light from the source reaches it and light from it reaches the probe.

### Objectives that are not linear

A real objective, such as the power $J=|x_m|^2$, is not a linear reading. But a small change of it is: $\delta J=2\operatorname{Re}\big(\mathbf g^T\delta\mathbf x\big)$ with $\mathbf g=\partial J/\partial\mathbf x$ (see the note below). So near the current design, $J$ acts like a measurement with weights $\mathbf g$, and nothing else changes: $$A^T\boldsymbol\lambda=\mathbf g,\qquad\frac{dJ}{dp_i}=-2\operatorname{Re}\Big(\boldsymbol\lambda^T\frac{\partial A}{\partial p_i}\mathbf x\Big).$$ For $J=|x_m|^2$, $\mathbf g=\overline{x_m}\,\mathbf e_m$, and in FDFD $dJ/d\varepsilon_{r,i}=-2k_0^2\operatorname{Re}(\lambda_ix_i)$. The factor $\overline{x_m}$ is a phase reference. $J$ is the squared length of the arrow $x_m$ in the complex plane. A small change $\delta x_m$ makes it longer only if it points along $x_m$. With $x_m=1$, a change of $0.1$ gives $J=1.21$, a change of $0.1i$ gives $1.01$ (no first-order change), and a change of $-0.1$ gives $0.81$.

> [!NOTE] **Real objectives of complex fields.** 
> $J$ is real, but $\mathbf x$ is complex, so $J$ depends on both $\mathbf x$ and $\overline{\mathbf x}$. Treat them as independent variables. Then $$dJ=\sum_j\Big(\frac{\partial J}{\partial x_j}dx_j+\frac{\partial J}{\partial\overline{x_j}}d\overline{x_j}\Big).$$ Since $J$ is real, the second term is the complex conjugate of the first, so $dJ=2\operatorname{Re}\big(\sum_jg_j\,dx_j\big)$ with $g_j=\partial J/\partial x_j$. Example: $J=|x_m|^2=x_m\overline{x_m}$ gives $g_m=\overline{x_m}$, and $dJ=2\operatorname{Re}(\overline{x_m}\,dx_m)$. (Keyword: Wirtinger derivatives. Conventions differ: some texts solve $A^H\boldsymbol\mu=\overline{\mathbf g}$ and use $\boldsymbol\mu^H$, which gives $\boldsymbol\mu=\overline{\boldsymbol\lambda}$. The gradient is the same if you stay consistent. Swapping $A^T$ for $A^H$ alone is wrong: with a PML, $A^H=\overline A$ turns the absorbing layer into an emitting one.)

### Without symmetry

Nothing above required $A^T=A$. Symmetry is a convenience. With it, the same solver, even a black-box one, computes the adjoint, and $\boldsymbol\lambda$ is a wave we can picture. Without it, we solve with $A^T$ instead. With a direct solver this costs almost nothing extra: from $A=LU$, $A^T=U^TL^T$ reuses the same factors. An iterative solver needs products $A^T\mathbf v$, which are easy when we have the matrix. Physically, $A^T$ is often still a real device. For a magnetized ferrite, it is the same device with the bias reversed (established).

### The continuous version

The identity carries over, with integrals in place of dot products. A measurement is $M=\int g\,u\,dV$. The transpose becomes the operator $L^T$ defined by integration by parts: $$\int\lambda\,(Lu)\,dV=\int(L^T\lambda)\,u\,dV+\text{boundary terms}.$$ The adjoint field solves $L^T\lambda=g$, and its boundary condition is whatever makes the boundary terms vanish. Then $M=\int\lambda\,s\,dV$, as on the grid.

- **One frequency.** For $L=\nabla^2+k_0^2\varepsilon_r$, the volume terms match, so $L^T=L$. The boundary term is the one from "Reciprocity", and it vanishes when $\lambda$ is outgoing, like $u$. So the adjoint problem is the same device with the source $g$ at the probe. On the grid, the edge rows play the role of the boundary term, which is why the PML rows must be scaled to keep $A^T=A$.
- **In time.** For the wave equation, integrate $\lambda\,\partial_t^2u$ by parts in time: $$\int_0^T\lambda\,\partial_t^2u\,dt=\big[\lambda\,\partial_tu-(\partial_t\lambda)\,u\big]_0^T+\int_0^T(\partial_t^2\lambda)\,u\,dt .$$ The forward field starts at rest, $u=\partial_tu=0$ at $t=0$, so that end vanishes. At $t=T$, $u$ is whatever it is, so the term vanishes only if $\lambda=\partial_t\lambda=0$ at $t=T$. So the adjoint field has a *final* condition, and we step it backward from $t=T$. This is why adjoint FDTD runs backward in time (the adjoint-state method; established).

### Cost and limits

- **Two solves, one factorization.** With a direct solver, the factorization of $A$ is the expensive part, and the adjoint solve reuses it, with or without symmetry. So the gradient costs barely more than the forward solve ("Frequency or time?").
- **Time domain.** The adjoint simulation runs backward in time, and the gradient needs the product of both fields at every time. So we must store (or recompute) the forward fields. The same method is used in seismic imaging.
- **Reverse and forward mode.** Pulling the measurement back is the chain rule evaluated from the output side, called reverse mode. It wins for one output and many inputs. Many outputs and one input favour the other order, forward mode. With $k$ objectives, we need $k$ adjoint solves. Backpropagation in neural networks is reverse mode (established).
- **A gradient is only a direction.** Everything here is first order. The design problem is not convex, so gradient steps find a good local design, not the best one. Real devices also need constraints: $\varepsilon_r$ must end up as one of the two real materials, and features must be large enough to fabricate. The usual tools are filtering and projection (keyword: **topology optimization**).

### Check it against finite differences

The adjoint gradient must agree with the naive one. For a few cells, change $p_i$ by a small step $\Delta$, solve again, and compare $\big(J(p_i+\Delta)-J(p_i)\big)/\Delta$ with the adjoint value.

Here is the running example: the 1D grid with PML from "What it does, in numbers", with 13 physical cells, a slab of 5 cells with $\varepsilon_r=2$, the source on one side, and the probe on the other. For the first cell of the slab:

| |$dJ/d\varepsilon_r$|
|---|---|
|adjoint|$-0.093618$|
|finite difference, $\Delta=10^{-1}$|$-0.0998$|
|finite difference, $\Delta=10^{-4}$|$-0.093624$|
|finite difference, $\Delta=10^{-7}$|$-0.093618$|
|finite difference, $\Delta=10^{-12}$|$-0.0928$|

If $\Delta$ is too large, the finite-difference formula is wrong (truncation error). If $\Delta$ is too small, round-off takes over. In between there is a valley of agreement. If the two never agree, the usual culprits are a missing $2\operatorname{Re}$, $A^T$ versus $A^H$, a missing conjugate in $\mathbf g$, or PML rows that were not scaled to make $A$ symmetric.

> [!NOTE] **The same idea elsewhere, and why it is not a standard EE course.** 
> The identity is old, and several fields found it on their own (established). Circuit simulators compute the sensitivity of one output to every component value with one extra solve on the **adjoint network** (Director and Rohrer, 1969, using Tellegen's theorem). In optimal control it became the adjoint-state method in the 1970s, and in machine learning it became backpropagation in the 1980s. In reactor physics, $\boldsymbol\lambda$ is the *importance function*: how much a neutron at each point contributes to a detector. In quantum mechanics, the first-order energy shift is $\delta E=\langle\psi|\delta H|\psi\rangle$ (Hellmann–Feynman). There the operator is self-adjoint and the problem is an eigenproblem, so the "adjoint field" is the state itself, and no extra solve is needed. My guess (speculation) is that courses teach analysis, "given a design, find the field", and leave design to experience. The adjoint only pays off when a design has thousands of parameters, which became common in EM only with fast solvers and inverse design in the 2010s.

In one sentence: to learn how every input affects one output, ask the output once, backward.

## What this note leaves out

This note followed one thread, from "what is a linear law?" to a solver and its gradient. Textbooks cover much more. Below are the topics I would add next, each with the place it attaches to this note.

- **Modes: $Lu=\lambda u$.** Every problem here was driven: given a source, find the field. The other half of linear algebra is the eigenproblem. Fix $\omega$ and look for waves $e^{i\beta z}$ along a waveguide: $\beta$ is an eigenvalue and the mode profile is an eigenvector (a mode solver). Or fix the geometry and look for the $\omega$ at which a source-free solution exists: cavity resonances, like the tied rope and the nearly zero eigenvalue in "Dividing by the symbol". With loss or a PML, these $\omega$ become complex (quasi-normal modes), the same complex frequency as the friction example in "The symbol and the dispersion relation". In photonic design, the objective is often the overlap of the field with a mode, so the adjoint source is the mode itself. Keywords: mode solver, quasi-normal modes, mode overlap.
- **Interfaces.** Integrate Maxwell's equations across a thin layer at an interface. Tangential $\mathbf E$ and $\mathbf H$ are continuous; normal $\mathbf D$ and $\mathbf B$ jump by any surface charge or current. The $H_z$ case was one example. For a plane wave hitting a flat interface, the system is still shift-invariant _along_ the interface, so $k_x$ cannot change (Proof 2). That is Snell's law. Matching the fields gives the Fresnel coefficients. When $k_x>n_2k_0$, the transmitted wave is evanescent, as in "Two kinds of components": total internal reflection, which is what guides light in fibres and waveguides. Keywords: Fresnel equations, Brewster angle, total internal reflection.
- **The equivalence principle and the far field of real antennas.** The field outside a closed surface can be reproduced by currents on that surface alone: $\mathbf J_s=\hat{\mathbf n}\times\mathbf H$ and a (fictitious) magnetic current $\mathbf M_s=-\hat{\mathbf n}\times\mathbf E$. This is "What counts as a source?" once more: fields on a surface can be traded for sources. MoM uses it. Simulations use it to inject a clean plane wave (total-field/scattered-field sources) and to predict an antenna's far-field pattern from a small near-field simulation. Directivity, gain and input impedance follow, as $R_{\text{rad}}$ did for the short dipole. Keywords: Love's equivalence principle, TF/SF source, near-to-far-field transform.
- **Fourier optics.** "What does a free wave look like?" is its first chapter. A lens performs a Fourier transform, an imaging system is a filter on $k_x$, and the image of a point is the point spread function ($K^*K$ in the convolution note). Keywords: numerical aperture, point spread function, near-field imaging.

More briefly:

- **The finite element method.** The "elements" basis in depth: the weak form, meshing, and edge elements.
- **Large solves.** Preconditioners for iterative solvers, and fast algorithms that make route K's dense matrices cheap (keyword: fast multipole method).
- **Statics and potentials.** $\omega=0$, where low-frequency breakdown lives. The potentials from the short-antenna note in full, gauge choice, and the Laplace and Poisson equations. Capacitance and inductance are the circuit limit. MoM codes are usually written with potentials (the mixed-potential integral equation).
- **Transmission lines and ports.** 1D waves of voltage and current. S-parameters are the matrix from incoming to outgoing waves at the ports: the radiation condition, written for a circuit. This is the daily language of RF.
- **Scattering.** Cross sections; Mie theory (the exact solution for a sphere, and a standard test case for solvers); the Born approximation (one step of the Lippmann–Schwinger equation).
- **Further out.** Momentum and radiation pressure (the Maxwell stress tensor); relativity (the speed of light is the same for every observer, and $\mathbf E$ and $\mathbf B$ become one object); anisotropic and gyrotropic media (biased ferrites); nonlinear optics (where Linearity fails); quantum optics (photons).