---
tags:
  - "#physics"
  - "#electrodynamics"
  - numerical
---
**Note:** this article stems from my conversation with Claude Opus 5.5

## What are physical laws?

Take a field $u(\mathbf{r},t)$ as an example. Most of the $u$ that you can write down never appear in nature. For instance, when you switch on a lamp here, the field far away cannot change at the same instant; light needs time to get there. To discriminate which $u$ are allowed and which aren't, we need a split, which is why most physical laws are in the form of constraints (equations, sometimes inequalities)
$$F[u]=0.$$

## Linearity

Now suppose superposition holds: if $u_1$ and $u_2$ are both allowed, so is $au_1+bu_2$. (Two light beams cross and go on as if the other were never there.) That is to say, the allowed fields form a vector space. In linear algebra, such a space is always the null space of some linear map. For example, the plane $x+y+z=0$ is a vector space, and it is the null space of the matrix $(1\;1\;1)$. So we can write the law as $Lu=0$ with $L$ linear.

Real experiments, however, often contain what we control (e.g. antenna current, speaker, pulse). Call it the source $s$. Superposition still holds with sources: if $s_1$ gives $u_1$ and $s_2$ gives $u_2$, then $s_1+s_2$ gives $u_1+u_2$. The simplest law with this property is
$$Lu=s,$$
which is just $A\mathbf{x}=\mathbf{b}$ for matrices.

What do the solutions look like? Say $u_p$ is one solution, and $u_0$ is anything in the null space ($Lu_0=0$). Then $L(u_p+u_0)=s+0=s$, so $u_p+u_0$ is also a solution. Conversely, if $u_1$ and $u_2$ both solve $Lu=s$, then $L(u_1-u_2)=0$, so they differ by something in the null space. Therefore every solution looks like
$$u=u_p+u_0.$$
That is to say, the field you measure is the field caused by your source, plus free waves that would be there anyway. Those free waves come from elsewhere or from earlier times. This is why "no source" and "with a source" are the two cases worth separating: they are the two pieces of every solution.

Note that this also means $Lu=s$ alone does not give a unique answer, since we can always add a free wave. To pin down one answer, we need extra conditions: initial conditions, boundary conditions, or "no waves come in from infinity" (radiation condition). Keep this in mind; it comes back when we try to invert $L$.

> [!NOTE]
> What counts as a source? Mindful readers might be curious how we discriminate the source $s$ from the field $u$. The rule is simple: any term that depends on $u$ belongs to the operator side, and only terms we prescribe independently stay in $s$.
> 
> Glass is a good example. Glass contains electrons. When the field pushes them, they move and form a current. That current is also a source of field, but it is proportional to the field itself, so by the rule it goes into $L$. This is all $\varepsilon_r$ is: a bookkeeping move, the currents of bound electrons moved from the right side to the left side. The antenna current stays on the right because we drive it. To be honest, then, $u$ versus $s$ is a modeling choice. If the antenna is affected by the field it radiates, its current also depends on $u$, and part of $s$ must move into $L$.
> 
> Two consequences are worth noting now.
> - The electrons do not respond instantly; they keep moving a little after the push. So in time, the bound current is a convolution of the field with some response (see [convolution](convolution.md)). At a single frequency, that convolution becomes a plain number, $\varepsilon_r(\omega)$. This is one big reason to work one frequency at a time.
> - When we later design a device by changing $\varepsilon_r(\mathbf{r})$, we are changing $L$, not $s$. So optimizing a device means asking how the solution changes when $L$ changes.

## Why not write the law as $u=Ks$?

How do we get $u$ from $s$? The tempting answer is: if $L$ can be inverted, then $u=L^{-1}s$. Relating back to my previous [convolution](convolution.md) note, $L^{-1}$ plays the role of $K$ in $g=Kf$. Its kernel is the [Green's function](Green's%20function.md): the response to a unit point source.

But wait. We just saw that $L$ has a null space, and a matrix with a null space has no inverse. So $L^{-1}$ does not exist yet. It only exists after we add the extra conditions above, because they remove the freedom to add free waves and pick one answer. How the radiation condition does this is a story for later.

Then why don't we write physical laws in the "inverse" form $u=Ks$ directly? In matrix terms, $L$ is sparse while $K$ is dense. A derivative at a point only uses values at neighbouring points (see the next section on how derivatives become matrices), so each row of $L$ has only a few nonzero entries. On the other hand, a source at one point eventually affects every other point, since the effect passes from neighbour to neighbour until it fills the space. So $K$ has almost no zeros. (This holds at one frequency. In space *and* time, $K$ is zero outside the light cone, because the effect needs time to travel: the lamp example.) The rule is short; its consequence is long, and it depends on the whole geometry and on the extra conditions. That is to say, nature's laws are short in the $L$ form, and their consequences are written in the $K$ form.

## Locality

So far $L$ is just "some linear operator". Why are physical laws almost always differential equations?

The most general linear rule we can write is an integral, just like in my [convolution](convolution.md) note:
$$(Lu)(x)=\int W(x,x')\,u(x')\,dx'.$$
(I use $W$ here, since $K$ is already taken by the inverse.) Here $W(x,x')$ says how much the value at $x'$ matters for the rule at $x$.

Now add an assumption: the field here only "listens" to its immediate neighbourhood. That is to say, $W(x,x')$ is nonzero only when $x'$ is very close to $x$. (This is not the same as "nothing happens instantly far away" from the lamp example. The heat equation is local, yet heat spreads to every distance instantly. Finite speed is a separate property of $L$; the wave equation will have it.) Since only $x'$ near $x$ matters, we can Taylor-expand $u(x')$ around $x$:
$$u(x')=u(x)+(x'-x)\,u'(x)+\tfrac12(x'-x)^2u''(x)+\dots$$
Put this into the integral, and each term pulls out one derivative of $u$ at $x$:
$$(Lu)(x)=m_0(x)\,u(x)+m_1(x)\,u'(x)+m_2(x)\,u''+\dots,\quad m_n(x)=\int W(x,x')\,\frac{(x'-x)^n}{n!}\,dx'.$$
If $W$ only reaches a tiny distance $\epsilon$, then $m_n$ is roughly of size $\epsilon^n$, so the higher terms quickly become negligible. Therefore, a local linear rule is a differential operator with only a few low derivatives. In other words, differential equations are integral equations with very narrow kernels. To be honest, "negligible" means $\epsilon$ is much smaller than the length over which $u$ changes, e.g. the wavelength. When that fails, the higher terms matter; for materials this is spatial dispersion (see "In a medium"). The rigorous version of "local means differential" is **Peetre's theorem**.

Note that if the system is also homogeneous (the same here or there), the $m_n$ cannot depend on $x$, so the coefficients are constant. Then $L$ is a [convolution](convolution.md) with a very narrow kernel: $u'=\delta'*u$, $u''=\delta''*u$. Everything in the convolution note applies.

> [!NOTE]
> **Derivatives as matrices:**
> 
> To see what "sparse" means, put the function on a grid. Pick points $x_1,\dots,x_N$ with spacing $h$, and store $u_i=u(x_i)$:
> $$\mathbf{u}=(u_1,u_2,\dots,u_N).$$
> The function is now a vector, and a linear operator becomes an $N\times N$ matrix.
> 
> What does the second derivative look like on a grid? Taylor-expand one step forward and one step back:
> $$u_{i\pm1}=u_i\pm h\,u_i'+\tfrac12h^2u_i''\pm\tfrac16h^3u_i'''+\dots$$
> Add the two, and the odd terms cancel:
> $$u_{i+1}+u_{i-1}-2u_i=h^2u_i''+O(h^4)\quad\Longrightarrow\quad u_i''\approx\frac{u_{i-1}-2u_i+u_{i+1}}{h^2}.$$
> This is the finite-difference formula, and the pattern $(1,-2,1)$ is called a stencil.
> 
> Row $i$ of the matrix computes $u_i''$, so it only needs $u_{i-1},u_i,u_{i+1}$. Everything else in the row is zero. For $N=5$:
> $$D_2=\frac{1}{h^2}\begin{pmatrix}-2&1&0&0&0\\1&-2&1&0&0\\0&1&-2&1&0\\0&0&1&-2&1\\0&0&0&1&-2\end{pmatrix}$$
> Each row has 3 nonzeros, however large $N$ is. That is exactly what sparse means: entry $(i,j)$ is zero because point $j$ does not directly affect point $i$. Note also that every row is the same stencil shifted by one. This is the grid version of shift invariance (a Toeplitz matrix), and the stencil $(1,-2,1)/h^2$ is the grid version of $\delta''$.
> 
> The 1D Helmholtz operator $u''+k_0^2\varepsilon_r(x)\,u$ is then easy. The second term just multiplies each $u_i$ by its own number, which is a diagonal matrix:
> $$L=D_2+\mathrm{diag}\big(k_0^2\varepsilon_{r,1},\dots,k_0^2\varepsilon_{r,N}\big).$$
> Still sparse. Solving $L\mathbf{u}=\mathbf{s}$ is now literally $A\mathbf{x}=\mathbf{b}$, and this is all a 1D FDFD solver is.
> 
> Now, why is the inverse dense? Take the $3\times3$ version (with $h=1$):
> $$A=\begin{pmatrix}-2&1&0\\1&-2&1\\0&1&-2\end{pmatrix},\qquad A^{-1}=-\frac14\begin{pmatrix}3&2&1\\2&4&2\\1&2&3\end{pmatrix}.$$
> (Multiply them to check.) $A$ has zeros in the corners, but $A^{-1}$ has none. Column $j$ of $A^{-1}$ is the response to a unit source at point $j$, i.e. the Green's function on the grid. The rule only connects neighbours, but the response reaches everywhere.
> 
> Mindful readers might notice that this $A$ *is* invertible, even though I said $L$ has a null space. The trick is hidden in the first and last rows. Writing the matrix this way silently assumes $u=0$ just outside the grid. That is a boundary condition, i.e. one of the extra conditions, and it is what kills the null space. (For this matrix, at least. With the Helmholtz term added, a closed box can still resonate; see the tied rope in "Choosing the outgoing solution".) To be honest, every finite matrix has made such a choice in its boundary rows. Choosing it well is the whole problem of the radiation condition, which comes later.
> 
> Finally, as $h\to0$ the vector becomes the function again, and a matrix–vector product $\sum_jA_{ij}u_j$ becomes an integral $\int W(x,x')\,u(x')\,dx'$. A sparse matrix becomes a kernel concentrated at $x'=x$ (a derivative of $\delta$), and a dense matrix becomes an ordinary spread-out kernel. So the matrix picture and the integral picture are the same thing.

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

Take $e^{ikx}$ and differentiate it:
$$\frac{d}{dx}e^{ikx}=ik\,e^{ikx},\qquad \frac{d^2}{dx^2}e^{ikx}=(ik)^2e^{ikx}=-k^2e^{ikx}.$$
The result is the same function times a number. So $e^{ikx}$ is an eigenvector of $\frac{d^n}{dx^n}$ for every $n$, and therefore of $L$:
$$L\,e^{ikx}=p(ik)\,e^{ikx}.$$

The coefficients must be constant. If a coefficient depends on $x$, say $m(x)\,\frac{d}{dx}$, then we get $ik\,m(x)\,e^{ikx}$. This is not a number times $e^{ikx}$.

### Proof 2: commuting with shifts

This proof rests on a fact from linear algebra: **commuting operators share eigenvectors** (when eigenvalues are distinct). For example, let $S=\mathrm{diag}(1,2)$. If $AS=SA$, then $A$ is diagonal too, so $A$ has the same eigenvectors as $S$. If we know the eigenvectors of $S$, we get those of $A$ for free.

Here the role of $S$ is played by the shifts. Write $T_a$ for the shift, $(T_a u)(x)=u(x+a)$. In operator language, homogeneous means shift-invariant: $L$ commutes with every shift: $LT_a=T_aL$.

For a differential operator, commuting with shifts is the same as having constant coefficients. Indeed, $T_a\big(m(x)\,u'\big)=m(x+a)\,u'(x+a)$, while $m(x)\,(T_a u)'=m(x)\,u'(x+a)$. These agree for every $u$ only if $m(x+a)=m(x)$ for every $a$, that is, only if $m$ is constant. The same computation works with $u$ or any higher derivative in place of $u'$, so it covers every term of $L$.

The eigenvectors of the shifts are the exponentials: $T_a e^{ikx}=e^{ika}e^{ikx}$. Since $L$ commutes with shifts, $Le^{ikx}$ is an eigenvector of every shift with the same eigenvalue:
$$T_a\big(Le^{ikx}\big)=L\big(T_a e^{ikx}\big)=e^{ika}\,Le^{ikx}.$$
The only functions $f$ with $T_a f=e^{ika}f$ for every $a$ are multiples of $e^{ikx}$. (Set $x=0$ to get $f(a)=e^{ika}f(0)$.) So $Le^{ikx}=\lambda\,e^{ikx}$ for some number $\lambda$.

This proof never uses the formula for $L$. It only uses the fact that $L$ commutes with shifts. That is why it also covers convolutions, as in my [convolution](convolution.md) note.

### Back to $Lu=0$

Since $Le^{ikx}=p(ik)\,e^{ikx}$, the exponential $e^{ikx}$ is a free wave exactly when $p(ik)=0$. This condition on $k$ is called the dispersion relation (see the following sections). Combinations of these exponentials give the free waves. (The shift argument also allows complex $k$. All it needs is $f(a)=c(a)f(0)$ with $c(a+b)=c(a)c(b)$, and $c(a)=e^{za}$ works for any complex $z$. A complex $k$ gives a wave that decays or grows along $x$: an **evanescent wave**. We meet them again in the PML section.)

(Strictly, the $e^{ikx}$ are not a basis in the finite-dimensional sense. The claim that every function is a combination of them is the content of the [Fourier transform](Fourier%20transform.md) , where the sum becomes an integral over $k$.)

The payoff: $L$ acts on each component separately,
$$Lu=\int_{-\infty}^{\infty}A(k)\,p(ik)\,e^{ikx}\,dk.$$
So in the $k$ picture, $L$ is just multiplication by $p(ik)$, called the *symbol* of $L$.

## Space, time, and free waves

### Adding time

The field depends on $t$ too, and the system is also the same now or later. So now there are two kinds of shifts: $T_a$ in space, $(T_a u)(x,t)=u(x+a,t)$, and $S_b$ in time, $(S_b u)(x,t)=u(x,t+b)$. Note that they commute with each other (shifting in space then time is the same as time then space), and $L$ commutes with both.

So let us run the same argument as before, with both shifts at once. Which functions are eigenvectors of every $T_a$ and every $S_b$? Suppose
$$f(x+a,\,t+b)=e^{ika}\,e^{-i\omega b}\,f(x,t)\quad\text{for every }a,b.$$
Set $x=0$ and $t=0$, and we get $f(a,b)=e^{ika}e^{-i\omega b}f(0,0)$. That is to say, the only such functions are multiples of
$$e^{ikx}\,e^{-i\omega t}=e^{i(kx-\omega t)}.$$
So the product form is not an assumption; it is forced. (I write $-i\omega t$ instead of $+i\omega t$ on purpose; the reason comes in a moment.) And by the same argument as before, since $L$ commutes with both shifts, $L\,e^{i(kx-\omega t)}$ is again a number times $e^{i(kx-\omega t)}$.

We can also check it directly. $\partial_x$ only sees the $x$ part, and $\partial_t$ only sees the $t$ part:
$$\partial_x\,e^{i(kx-\omega t)}=ik\,e^{i(kx-\omega t)},\qquad \partial_t\,e^{i(kx-\omega t)}=-i\omega\,e^{i(kx-\omega t)}.$$

Is this family big enough? Yes: do Fourier twice. First expand in $x$ at each time, $u(x,t)=\int A(k,t)\,e^{ikx}\,dk$. Then expand each coefficient $A(k,t)$ in time. Together,
$$u(x,t)=\iint A(k,\omega)\,e^{i(kx-\omega t)}\,dk\,d\omega.$$
So, at this point, $k$ and $\omega$ are two independent labels. Nothing connects them yet.

What does one basis function look like? Its crests sit where $kx-\omega t$ is constant, i.e. at $x=(\omega/k)\,t+\text{const}$. So it is a **wave**, with crests moving at speed $\omega/k$ (the phase velocity)! The number $k$ counts crests per length ($k=2\pi/\text{wavelength}$), and $\omega$ counts oscillations per time ($\omega=2\pi/\text{period}$). With the minus sign, $k>0$ and $\omega>0$ mean the wave moves to $+x$. This is called a plane wave.

### The symbol and the dispersion relation

So far, linearity, locality and homogeneity in space and time tell us that every free solution is a combination of plane waves. Now we ask which plane waves a given law allows.

$L$ is a polynomial in two derivatives, $L=p(\partial_x,\partial_t)$, and on a plane wave
$$L\,e^{i(kx-\omega t)}=p(ik,\,-i\omega)\,e^{i(kx-\omega t)}=D(k,\omega)\,e^{i(kx-\omega t)}.$$
$D(k,\omega)$ is the **symbol** again, now with two labels. For example,
$$L=a\,\partial_t^2+b\,\partial_x^2+c\,\partial_x+d\,\partial_t+e
\quad\Longrightarrow\quad
D(k,\omega)=-a\omega^2-bk^2+ic\,k-id\,\omega+e.$$

Mindful readers might wonder what each term does. Here is the simplest case: drop space, and take $L=\partial_t+\gamma$ with $\gamma>0$. Then $D=-i\omega+\gamma$, and $D=0$ gives $\omega=-i\gamma$. A complex frequency! The wave is $e^{-i\omega t}=e^{-\gamma t}$, which decays. That is to say, a first-order time derivative means loss: the field dies out, like a swing with friction.

This also tells us when such terms are absent. Laws without loss look the same when the movie is played backward ($t\to-t$). Under $t\to-t$, $\partial_t$ flips sign but $\partial_t^2$ does not. So a lossless law can only have even time derivatives. In the same way, a law that looks the same in a mirror ($x\to-x$) can only have even space derivatives, so $c=0$. (To be honest, both rules hold only for a single real scalar field; see the correction under "Check 3: energy".) (For example, light in vacuum has both properties, so for it, only $\partial_t^2$, $\partial_x^2$ and a constant can appear. We will come back to loss with $\varepsilon_r$: it is exactly what a complex $\varepsilon_r$ describes.)

Before any law, a function of $x$ and $t$ can wiggle at any rate in space and at any rate in time, independently. For example, $\cos(3x)\cos(5t)$ is a perfectly good function, built from $(k,\omega)=(\pm3,\pm5)$. So in the expansion, every pair $(k,\omega)$ is allowed.

Apply $L$ to a general field:
$$Lu=\iint A(k,\omega)\,D(k,\omega)\,e^{i(kx-\omega t)}\,dk\,d\omega.$$
For a free wave, this must be zero. The plane waves are independent (like the axes of a basis), so each component must vanish separately:
$$A(k,\omega)\,D(k,\omega)=0\quad\text{for every }(k,\omega).$$
That is to say, $A(k,\omega)$ can only be nonzero where
$$D(k,\omega)=0.$$
Picture the $(k,\omega)$ plane. A general function spreads its $A$ over the whole plane. The set $D=0$ is a few curves in that plane, and a free wave puts its $A$ only on those curves. This condition is called the **dispersion relation**. 

> [!NOTE]
> Mindful readers might notice a problem. If $A$ is an ordinary function that is nonzero only on a curve, the double integral gives zero, since a curve has no area. So $A$ must be concentrated on the curve, like a $\delta$. For each $k$, the equation $D(k,\omega)=0$ has only a few roots in $\omega$. Call them $\omega_1(k),\omega_2(k),\dots$ (if $D$ is second order in $\omega$, there are two). Then
> $$A(k,\omega)=A_1(k)\,\delta\big(\omega-\omega_1(k)\big)+A_2(k)\,\delta\big(\omega-\omega_2(k)\big)+\dots$$
> Put this into the expansion. The $\omega$ integral just picks out each root, since $\int\delta(\omega-\omega_1)\,e^{-i\omega t}\,d\omega=e^{-i\omega_1t}$. So
> $$u(x,t)=\int\Big[A_1(k)\,e^{i(kx-\omega_1(k)t)}+A_2(k)\,e^{i(kx-\omega_2(k)t)}+\dots\Big]\,dk.$$
> This is the general free wave: for each $k$, a few waves, each with its own amplitude.
> 
> Note that the amplitudes $A_1(k),A_2(k),\dots$ are still free. This is the null space from before: there are many free waves, and a law alone does not choose among them. The extra conditions choose. For example, with two roots, there are two amplitudes for each $k$, so we need two pieces of information, such as the field and its time derivative at $t=0$. We will see this concretely for light in the next section.

## Light in free space

Now we use one experimental fact: in vacuum, light of every colour travels at the same speed $c$. 

### Light in 1D: the wave equation

In 1D, it can go right or left. So the allowed plane waves are exactly those with
$$\omega=ck \quad \text{or} \quad \omega=-ck.$$
We want one polynomial $D$ that vanishes on both. The product does it:
$$D(k,\omega)=(\omega-ck)(\omega+ck)=\omega^2-c^2k^2.$$
Now swap the numbers back into derivatives, $\omega^2\leftrightarrow-\partial_t^2$ and $k^2\leftrightarrow-\partial_x^2$:
$$\partial_x^2u-\frac{1}{c^2}\,\partial_t^2u=0.$$
This is the wave equation. We did not guess it; it is what "every colour at the same speed, in both directions" looks like as a law. One honest caveat: any multiple $D\,q(k,\omega)$ also vanishes on both lines. The wave equation is the *lowest-order* law with exactly these free waves, and choosing the lowest order is an assumption. 

Now we can write down every free wave of light in 1D. The two roots are $\omega_1(k)=ck$ and $\omega_2(k)=-ck$, so
$$u(x,t)=\int\Big[A_1(k)\,e^{ik(x-ct)}+A_2(k)\,e^{ik(x+ct)}\Big]\,dk=F(x-ct)+G(x+ct),$$
where $F(\xi)=\int A_1(k)\,e^{ik\xi}\,dk$ and $G(\xi)=\int A_2(k)\,e^{ik\xi}\,dk$. That is to say, a free wave is any shape moving right at speed $c$, plus any shape moving left at speed $c$. (This is called d'Alembert's solution.)

How do the initial conditions fix the amplitudes? Say we know $u(x,0)=f(x)$ and $\partial_tu(x,0)=g(x)$, and expand both: $f(x)=\int\hat f(k)\,e^{ikx}\,dk$, and the same for $g$. Setting $t=0$ in $u$ and in $\partial_tu$ gives, for each $k$,
$$\hat f(k)=A_1(k)+A_2(k),\qquad \hat g(k)=-ick\,A_1(k)+ick\,A_2(k).$$
Two equations, two unknowns. Solving,
$$A_1(k)=\tfrac12\Big(\hat f(k)+\frac{i\,\hat g(k)}{ck}\Big),\qquad A_2(k)=\tfrac12\Big(\hat f(k)-\frac{i\,\hat g(k)}{ck}\Big).$$
So the two initial conditions pick exactly one free wave out of the null space, as promised. (At $k=0$ the two roots meet and the formula divides by zero. That mode is the uniform motion $u=a+bt$, and $\hat g(0)$ sets $b$.)

(Incidentally, the factorization also works with derivatives: $(\partial_t-c\,\partial_x)(\partial_t+c\,\partial_x)u=0$, one factor for each direction.)

### Going to 3D: richer geometry

In 3D, there are three directions to shift in, so there are three labels, $k_x,k_y,k_z$, and the plane waves are
$$e^{i(k_xx+k_yy+k_zz-\omega t)}=e^{i(\mathbf{k}\cdot\mathbf{r}-\omega t)}.$$
The vector $\mathbf{k}=(k_x,k_y,k_z)$ points in the direction the wave travels, and its length $|\mathbf{k}|$ plays the role of $k$ in 1D. On such a wave, $\partial_x\to ik_x$, $\partial_y\to ik_y$, $\partial_z\to ik_z$.

The experimental fact is now: same speed $c$ in every direction. So the allowed waves satisfy $\omega=\pm c|\mathbf{k}|$ for every direction of $\mathbf{k}$, and the same product trick gives
$$D(\mathbf{k},\omega)=\omega^2-c^2|\mathbf{k}|^2=\omega^2-c^2\big(k_x^2+k_y^2+k_z^2\big).$$
Swap back, with $k_x^2\leftrightarrow-\partial_x^2$ and so on:
$$\big(\partial_x^2+\partial_y^2+\partial_z^2\big)u-\frac{1}{c^2}\,\partial_t^2u=0.$$
The combination $\partial_x^2+\partial_y^2+\partial_z^2$ is called the Laplacian, $\nabla^2$. So this is all the Laplacian is: the derivative version of squared length, $|\mathbf{k}|^2$. That is where the geometry comes in. Rotating your axes does not change a length, so $\nabla^2$ does not depend on how you orient $x,y,z$. That is exactly what "same speed in every direction" needs.

Why second order, and not first? A first-order term like $\partial_x$ would put $k_x$ alone into $D$, which singles out the $x$ direction. The simplest combination that treats all directions alike is $|\mathbf{k}|^2$, which is second order. (This argument assumes $u$ is a single number at each point. Light is not; see the Maxwell section.)

## Light in mediums

### Nonlocal in time

In glass, light is slower by a factor $n$, so $\omega=\pm c|\mathbf{k}|/n$, and the wave equation becomes
$$\nabla^2u-\frac{n^2}{c^2}\,\partial_t^2u=0.$$
In real glass, however, $n$ depends on colour (that is how a prism works). Then the allowed waves satisfy $\omega=\pm c|\mathbf{k}|/n(\omega)$, and $D$ is no longer a polynomial in $\omega$. So $L$ is no longer a differential operator in time.

What happened to locality? Recall the bound electrons from earlier: they keep moving a little after the push. The glass has **memory**, so its response at time $t$ depends on the field at earlier times. This is a convolution in time with a kernel of real width. But a convolution in time becomes plain multiplication at a single frequency. So at one $\omega$, the whole memory of the glass is just one number,
$$n(\omega)^2=\varepsilon_r(\omega).$$
This is the deep reason textbooks work **one frequency at a time**. In time, the material is a convolution; at one frequency, it is a number.

### Nonhomogeneous in space

Finally, let $n$ change with position, $n(\mathbf{r})$: a lens, a waveguide, a device we want to design. Now the system is not the same here or there, so the exponentials in space are no longer eigenvectors (recall the $m(x)\,\frac{d}{dx}$ example).

The system is still the same now or later, so the exponentials in time still work. So we keep only those. Pick one frequency, and write the field as $u(\mathbf{r})\,e^{-i\omega t}$. Then $\partial_t^2\to-\omega^2$, but $\nabla^2$ has to stay as it is. With $k_0=\omega/c$:
$$\nabla^2u+k_0^2\,\varepsilon_r(\mathbf{r})\,u=s.$$
This is the Helmholtz equation. It is algebra in time and a differential equation in space, and which half is which is decided by the symmetry that survives. On the grid, it is the sparse matrix $D_2+\mathrm{diag}(k_0^2\varepsilon_{r,i})$ from before, which is exactly what an FDFD solver builds.

To be honest, putting $\varepsilon_r(\mathbf{r})$ inside the scalar equation was a guess. It works when $\varepsilon_r$ changes slowly compared with the wavelength. At a sharp interface it can fail, because light is not a scalar: it has two polarizations, so it is a vector. Maxwell's equations carry this extra structure, and they will tell us exactly when the scalar equation is right.

## Light is a vector: Maxwell's equations

The section "Nonhomogeneous in space" ended with a problem. Light has two polarizations, so the field is a vector, and the scalar equation was only a guess. We need the law for a vector field.

Here I stop deriving and take the law from experiment, just as we took "same speed $c$" from experiment before. To be honest, a deeper argument for its form exists (see the note at the end of this section), but nothing later depends on it. What matters is that every consequence we use below can be checked.

The law uses two vector fields, the electric field $\mathbf E$ and the magnetic field $\mathbf B$. The sources are the charge density $\rho$ and the current density $\mathbf J$. In SI units, with $\mathbf D=\varepsilon_0\mathbf E$ and $\mathbf H=\mathbf B/\mu_0$ in vacuum, there are four equations, each backed by an experiment:

| Equation | Experiment |
|---|---|
| $\nabla\cdot\mathbf D=\rho$ | Coulomb: a charge makes an $\mathbf E$ field that spreads out from it |
| $\nabla\cdot\mathbf B=0$ | no magnetic charge has ever been found |
| $\nabla\times\mathbf E=-\partial_t\mathbf B$ | Faraday: a changing magnetic field induces a voltage around a loop |
| $\nabla\times\mathbf H=\mathbf J+\partial_t\mathbf D$ | Ampère: a current makes a magnetic field; plus Maxwell's term $\partial_t\mathbf D$ |

Maxwell's term $\partial_t\mathbf D$ is the only piece that no experiment of his time showed directly. It is forced by charge conservation. Take the divergence of the last equation. Since $\nabla\cdot\nabla\times=0$,
$$0=\nabla\cdot\mathbf J+\partial_t(\nabla\cdot\mathbf D)=\nabla\cdot\mathbf J+\partial_t\rho .$$
Without $\partial_t\mathbf D$, the law would demand $\nabla\cdot\mathbf J=0$, so charge could never pile up anywhere. That is false: a capacitor charges. With $\partial_t\mathbf D$, the law demands exactly $\partial_t\rho+\nabla\cdot\mathbf J=0$, which experiment confirms. In the language of the Linearity section, charge conservation is the condition for the source $s$ to be compatible with $L$. A source that violates it has no solution at all. (This is the Fredholm alternative: $Lu=s$ is solvable only if $s$ is orthogonal to the null space of the adjoint $L^*$; see the adjoint note in [convolution](convolution.md).)

### Check 1: light is a wave with speed $c$

In vacuum, $\rho=0$ and $\mathbf J=0$. Take the curl of Faraday's law and insert Ampère's law:
$$\nabla\times\nabla\times\mathbf E=-\partial_t\nabla\times\mathbf B=-\mu_0\varepsilon_0\,\partial_t^2\mathbf E .$$
Use the identity $\nabla\times\nabla\times\mathbf E=\nabla(\nabla\cdot\mathbf E)-\nabla^2\mathbf E$ and $\nabla\cdot\mathbf E=0$:
$$\nabla^2\mathbf E-\mu_0\varepsilon_0\,\partial_t^2\mathbf E=0 .$$
Each component obeys the scalar wave equation of the 3D section, with $c=1/\sqrt{\mu_0\varepsilon_0}$. Maxwell computed this number from static measurements of $\varepsilon_0$ and $\mu_0$, and it matched the measured speed of light (established). So the scalar wave equation was a shadow of this one, one component at a time.

### Check 2: two polarizations

Take a plane wave $\mathbf E=\mathbf E_0\,e^{i(\mathbf k\cdot\mathbf r-\omega t)}$. In vacuum $\nabla\cdot\mathbf E=0$, which becomes $\mathbf k\cdot\mathbf E_0=0$: the field is perpendicular to the direction of travel. That leaves two directions, the two polarizations, as observed. Faraday's law gives $\mathbf B_0=\mathbf k\times\mathbf E_0/\omega$, also perpendicular. The scalar equation knew about neither constraint. This is what it was missing.

### Check 3: energy

Dot Faraday's law with $\mathbf H$ and Ampère's law with $\mathbf E$, subtract, and use $\nabla\cdot(\mathbf E\times\mathbf H)=\mathbf H\cdot\nabla\times\mathbf E-\mathbf E\cdot\nabla\times\mathbf H$:
$$\partial_t\Big(\tfrac12\varepsilon_0|\mathbf E|^2+\tfrac12\mu_0|\mathbf H|^2\Big)+\nabla\cdot(\mathbf E\times\mathbf H)=-\mathbf J\cdot\mathbf E .$$
The energy at a point changes in only two ways: it flows out with the **Poynting vector** $\mathbf S=\mathbf E\times\mathbf H$, or it does work on the current, $\mathbf J\cdot\mathbf E$. We will use the one-frequency version of this for objective functions.

> [!NOTE]
> **Correction.** Earlier I said a lossless law can only have even time derivatives. Maxwell's equations are first order in time, yet Check 3 shows they lose nothing: in vacuum, energy only flows. They are unchanged under $t\to-t$ if also $\mathbf B\to-\mathbf B$. The even-derivative rule holds only for a single real scalar field. The general test for "lossless" is: every root $\omega(\mathbf k)$ is real for real $\mathbf k$. Also, lossless and time-reversal symmetric are different things. A magnetized ferrite is lossless but not time-reversal symmetric (it is non-reciprocal).

### The divergence equations are constraints

Count first. There are 6 unknowns (the components of $\mathbf E$ and $\mathbf B$) but 8 equations, so the law looks overdetermined. The two curl equations give $\partial_t\mathbf E$ and $\partial_t\mathbf B$: 6 equations, already enough to step the fields forward in time. The two divergence equations contain no $\partial_t$, so they predict nothing. Could they contradict the curl equations? Take the divergence of the two curl equations. Faraday gives $\partial_t(\nabla\cdot\mathbf B)=0$. Ampère gives $\partial_t(\nabla\cdot\mathbf D)=-\nabla\cdot\mathbf J=\partial_t\rho$. So if the two divergence equations hold at one moment, the curl equations keep them true forever. They constrain the starting state; they do not drive anything.

A small example has the same logic. The motion $\dot x=y$, $\dot y=-x$ keeps $x^2+y^2$ constant. Demanding $x^2+y^2=1$ is not a second law of motion; it only picks a starting point on the circle, and the motion keeps you there.

Then why keep the divergence equations at all? Because the curl equations only *keep* the divergence; they cannot *choose* it. Add any static gradient field $\nabla\phi$ to a solution: it has no curl, so the curl equations do not notice it, and the sum is still a solution. Without Gauss's law, such a field could sit in empty space with no charge to produce it. Gauss's law forbids this, by tying $\nabla\cdot\mathbf D$ to real charge. In the circle example, it is the statement that the radius is 1 and not 2. Two more reasons. In statics ($\partial_t=0$), the curl equations say almost nothing, and the divergence equations are most of the physics (electrostatics is Gauss's law). And $\nabla\cdot\mathbf B=0$ is itself an experimental fact (no magnetic charge).

In practice, the two curl equations are the law of motion, and they are what a time-domain solver steps forward. The divergence equations only need to hold in the initial fields. This comes back at one frequency.

> [!NOTE]
> **The deeper reason (optional).** Why curls, and why two fields? A short answer exists, though I don't use it later. A law that is linear, local, the same in every direction and first order in time can only use $\nabla\times$, $\nabla\cdot$ and $\nabla$ on vector fields. With one real vector field, $\partial_t\mathbf u\propto\nabla\times\mathbf u$ gives imaginary $\omega$: growth and decay, no waves. Two fields that drive each other with opposite signs give exactly the light cone, and the minus sign turns out to be energy conservation. So Maxwell's equations are, in a precise sense, the simplest such law. The same idea gives Dirac's equation for the electron. Keywords: symbol of a first-order system, isotropic tensors, Dirac's trick.

## Maxwell's equations: a revisit

### In a medium

The section "Light in mediums" said the material's memory is a convolution in time. Now we can make that precise. Matter contains bound charges, and $\mathbf E$ pushes them into a current $\mathbf J_b$. Split $\mathbf J=\mathbf J_s+\mathbf J_b$: the part we drive and the part that responds. By the rule from the "What counts as a source?" note, $\mathbf J_b$ goes into $L$.

What can $\mathbf J_b$ look like? We apply the same assumptions as before, now to the material. It is linear (weak fields), time-invariant (the material does not change), causal, and local in space (an atom responds to the field where it sits). Write $\mathbf J_b=\partial_t\mathbf P$, where the polarization $\mathbf P$ is the displacement of bound charge per volume. Then
$$\mathbf P(\mathbf r,t)=\varepsilon_0\int_0^\infty\chi(\mathbf r,\tau)\,\mathbf E(\mathbf r,t-\tau)\,d\tau .$$
This is a causal convolution in time (see [convolution](convolution.md)), local in space. $\chi$ depends on $\mathbf r$ when the material is not uniform. Ampère's law becomes $\nabla\times\mathbf H=\mathbf J_s+\partial_t\mathbf D$ with $\mathbf D=\varepsilon_0\mathbf E+\mathbf P$.

A toy model shows where $\chi$ comes from. Take an electron on a spring with friction: $m(\ddot x+\gamma\dot x+\omega_0^2x)=-eE$. At one frequency, $x=-(e/m)E/(\omega_0^2-\omega^2-i\gamma\omega)$. With $N$ electrons per volume, $P=-Nex$, so
$$\chi(\omega)=\frac{\omega_p^2}{\omega_0^2-\omega^2-i\gamma\omega},\qquad \omega_p^2=\frac{Ne^2}{\varepsilon_0m}.$$
This is the Lorentz model. Its imaginary part comes from the friction $\gamma$, the first-order time derivative from the dispersion-relation section. At one frequency, $\mathbf D=\varepsilon_0\varepsilon_r(\mathbf r,\omega)\mathbf E$ with $\varepsilon_r=1+\chi$.

Three consequences:

- **Loss.** $\varepsilon_r$ is complex. With the $e^{-i\omega t}$ convention, $\operatorname{Im}\varepsilon_r>0$ means loss (proved with the Poynting theorem below, under "At one frequency"). The Lorentz model agrees: $\operatorname{Im}\chi>0$ when $\gamma>0$.
- **Dispersion is unavoidable.** $\chi(\tau)$ is real, so $\varepsilon_r(-\omega)=\overline{\varepsilon_r(\omega)}$. $\chi$ is causal, so $\operatorname{Re}\varepsilon_r$ and $\operatorname{Im}\varepsilon_r$ determine each other ([Kramers–Kronig relations](Kramers-Kronig%20relations.md)). One consequence: a material with loss at any frequency must have an $\varepsilon_r$ that depends on $\omega$. A constant $\varepsilon_r\neq1$ is always a narrow-band approximation.
- **Spatial locality is an approximation.** It holds when the atomic (or meta-atom) scale is much smaller than the wavelength. It fails in nanoscale plasmonics and near metamaterial resonances, where $\varepsilon_r$ also depends on $\mathbf k$ (spatial dispersion).

(Bound charges can also circulate, giving a magnetization current $\nabla\times\mathbf M$ and $\mu_r$. In optics $\mu_r\approx1$. In RF, ferrites need $\mu_r$, and a biased ferrite needs a tensor $\mu$ that is non-reciprocal. I set $\mu_r=1$ from here on.)

### At one frequency

As in "Light in mediums", we keep the time exponentials and work at one $\omega$. From here on, $\mathbf J$ means the driven current $\mathbf J_s$. Real fields have $\hat u(-\omega)=\overline{\hat u(\omega)}$, so at a single frequency the physical field is $\operatorname{Re}\big(\mathbf E(\mathbf r)\,e^{-i\omega t}\big)$.

With $\partial_t\to-i\omega$:
$$\nabla\times\mathbf E=i\omega\mu_0\mathbf H,\qquad \nabla\times\mathbf H=-i\omega\varepsilon_0\varepsilon_r(\mathbf r)\mathbf E+\mathbf J .$$
Eliminate $\mathbf H$:
$$\nabla\times\nabla\times\mathbf E-k_0^2\,\varepsilon_r(\mathbf r)\,\mathbf E=i\omega\mu_0\mathbf J,\qquad k_0=\omega/c .$$
This is $Lu=s$ from the Linearity section, with $L=\nabla\times\nabla\times-k_0^2\varepsilon_r$, $u=\mathbf E$ and $s=i\omega\mu_0\mathbf J$. It is the vector version of the Helmholtz equation. (We can also keep both $\mathbf E$ and $\mathbf H$ as unknowns, as the Yee grid and FDTD do. The physics is the same; the matrices differ.)

**Gauss's law comes for free.** Take the divergence; $\nabla\cdot\nabla\times=0$ gives $-k_0^2\nabla\cdot(\varepsilon_r\mathbf E)=i\omega\mu_0\nabla\cdot\mathbf J$. With charge conservation, $\nabla\cdot\mathbf J=i\omega\rho$, this becomes $\nabla\cdot(\varepsilon_0\varepsilon_r\mathbf E)=\rho$. So at $\omega\neq0$ the divergence constraint is already inside the curl-curl equation.

> [!NOTE]
> Mindful readers might ask what happens as $\omega\to0$. The term $k_0^2\varepsilon_r$ vanishes, and every gradient $\nabla\phi$ is in the null space of $\nabla\times\nabla\times$ (a gradient has no curl). Only the divergence equation could control these fields, and at $\omega=0$ it is no longer implied. So at low frequency $L$ becomes nearly singular. Solvers call this **low-frequency breakdown**, and FEM and MoM codes need special fixes for it. A related point: a discretization must keep "curl of a gradient is zero" exactly, or false (spurious) modes appear. The Yee grid and edge elements do keep it. This comes back with discretization.

### The symbol is now a matrix

Run the homogeneity argument once more, on the curl-curl equation in vacuum ($\varepsilon_r=1$, no source). On a plane wave $\mathbf E_0e^{i\mathbf k\cdot\mathbf r}$, $\nabla\to i\mathbf k$, so
$$\nabla\times\nabla\times\mathbf E\ \to\ -\mathbf k\times(\mathbf k\times\mathbf E_0)=|\mathbf k|^2\mathbf E_0-\mathbf k\,(\mathbf k\cdot\mathbf E_0).$$
So on a plane wave, $L$ is not multiplication by a number but by a $3\times3$ matrix:
$$D(\mathbf k)=\big(|\mathbf k|^2-k_0^2\big)I-\mathbf k\mathbf k^{T}.$$
A free wave needs $D\mathbf E_0=0$ with $\mathbf E_0\neq0$. Check the two kinds of direction.
- $\mathbf E_0\perp\mathbf k$: $D\mathbf E_0=(|\mathbf k|^2-k_0^2)\mathbf E_0$. This is zero exactly when $|\mathbf k|=k_0$, the same dispersion relation as the scalar case.
- $\mathbf E_0\parallel\mathbf k$: $D\mathbf E_0=-k_0^2\mathbf E_0$. This is never zero at $\omega\neq0$.

So free waves need $|\mathbf k|=k_0$ and $\mathbf E_0\perp\mathbf k$. The scalar story gave only the first condition. The matrix also fixes the direction of the amplitude: this is polarization (Check 2), now coming straight from the dispersion relation. Note also the entry $-k_0^2$ along $\mathbf k$. It goes to zero as $\omega\to0$, which is the low-frequency breakdown from the note above.

### When is the scalar Helmholtz equation right?

Now we can answer the question from "Nonhomogeneous in space". Take a structure that does not change along $z$, so $\varepsilon_r=\varepsilon_r(x,y)$, as in any 2D simulation. Fields then do not depend on $z$ either. With $\partial_z=0$, Maxwell's equations split into two independent groups. $E_z$ only produces $H_x,H_y$ (since $(\nabla\times\mathbf E)_x=\partial_yE_z$ and $(\nabla\times\mathbf E)_y=-\partial_xE_z$), and $H_x,H_y$ only produce $E_z$ (since $(\nabla\times\mathbf H)_z=\partial_xH_y-\partial_yH_x$). So $\{E_z,H_x,H_y\}$ is closed, and likewise $\{H_z,E_x,E_y\}$. A scalar $\varepsilon_r$ does not turn vectors, so it does not mix the groups. Every 2D solution is a sum of one from each group, so we can study them separately. (The split fails for a tensor $\varepsilon_r$, a source that varies along $z$, or oblique incidence, $k_z\neq0$.) The two cases are:

- **$\mathbf E=E_z(x,y)\,\hat{\mathbf z}$.** Then $\nabla\cdot(\varepsilon_r\mathbf E)=\partial_z(\varepsilon_rE_z)=0$ automatically, and $\nabla\times\nabla\times(E_z\hat{\mathbf z})=-\nabla^2E_z\,\hat{\mathbf z}$. So
$$\nabla^2E_z+k_0^2\varepsilon_r(x,y)\,E_z=-i\omega\mu_0J_z .$$
This is exactly the scalar Helmholtz equation. There is no approximation, even at sharp interfaces.
- **$\mathbf H=H_z(x,y)\,\hat{\mathbf z}$.** Eliminating $\mathbf E$ instead of $\mathbf H$ gives (without sources)
$$\nabla\cdot\big(\varepsilon_r^{-1}\nabla H_z\big)+k_0^2H_z=0 .$$
Where $\varepsilon_r$ is constant, this is the Helmholtz equation again. At an interface, however, the $1/\varepsilon_r$ sits inside the derivative. It makes $\varepsilon_r^{-1}\partial_nH_z$ continuous across the interface (this is the tangential $\mathbf E$), and the scalar guess gets exactly this wrong.

So in 2D, one polarization is exactly scalar, and the other is scalar with a different operator. In 3D with a general $\varepsilon_r(\mathbf r)$, the components are coupled, and no scalar equation is exact. (These two cases are called TE and TM, but the labels are swapped between communities. Check which field is along $z$ before trusting a label.)

### Power

Products of real fields need care. The time average of $\operatorname{Re}(ae^{-i\omega t})\operatorname{Re}(be^{-i\omega t})$ is $\tfrac12\operatorname{Re}(a\bar b)$. So the time-averaged Poynting vector is $\langle\mathbf S\rangle=\tfrac12\operatorname{Re}(\mathbf E\times\overline{\mathbf H})$. Use $\nabla\cdot(\mathbf E\times\overline{\mathbf H})=\overline{\mathbf H}\cdot\nabla\times\mathbf E-\mathbf E\cdot\nabla\times\overline{\mathbf H}$ and the two equations above, and take the real part:
$$\nabla\cdot\langle\mathbf S\rangle=-\tfrac12\operatorname{Re}\big(\mathbf E\cdot\overline{\mathbf J}\big)-\tfrac{\omega\varepsilon_0}{2}\operatorname{Im}(\varepsilon_r)\,|\mathbf E|^2 .$$
Read it as: power leaving a point = power supplied by the source $-$ power absorbed there. So $\operatorname{Im}\varepsilon_r>0$ is absorption, as claimed. (The imaginary part of the same identity describes stored, reactive energy: the near field of an antenna and its Q.)

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

| Step           | Scalar (1–7)                                            | Maxwell                                                                                                               |
| -------------- | ------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------- |
| 1. Linearity   | $Lu=s$                                                  | the same; $u=\mathbf E$, $s=i\omega\mu_0\mathbf J$                                                                    |
| 3. Locality    | derivatives only                                        | the same: curls and divergences                                                                                       |
| 4. Homogeneity | plane waves are eigenvectors; the symbol is a number    | the same plane waves, but each carries a vector amplitude, so the symbol is a $3\times3$ matrix                       |
| 5. Dispersion  | $D=0$                                                   | $D\mathbf E_0=0$: again $\lvert\mathbf k\rvert=k_0$, plus $\mathbf E_0\perp\mathbf k$                                 |
| 6. Free space  | the wave equation                                       | each component obeys it (Check 1)                                                                                     |
| 7. Mediums     | $\varepsilon_r(\omega)$ from memory; Helmholtz in space | $\varepsilon_r(\omega)$ derived from bound currents; Helmholtz becomes curl-curl; the 2D $E_z$ case is exactly step 7 |

What is genuinely new, with no scalar counterpart:
- **Polarization.** The amplitude has a direction, and the matrix symbol fixes it ($\mathbf E_0\perp\mathbf k$).
- **Constraints.** Equations without $\partial_t$ (the divergences) restrict the starting state; the curl equations keep them true.
- **A constraint on the source.** $\mathbf J$ must conserve charge, or $Lu=s$ has no solution at all.
- **Physical power.** $\tfrac12\operatorname{Re}(\mathbf E\times\overline{\mathbf H})$ and absorption $\propto\operatorname{Im}\varepsilon_r$. Objective functions are built from these.

Two ideas run through everything.
- **Symmetry decides the basis.** Whatever shifts $L$ commutes with, those exponentials diagonalize it. A material keeps the time symmetry but breaks the space symmetry. So we use exponentials in time (one frequency at a time) and numerics in space.
- **The null space does not go away by itself.** Free waves can always be added. The finite matrix in "Derivatives as matrices" hid this in its boundary rows.

So the problem is now one sparse system $A\mathbf x=\mathbf b$ on a grid, with one open question: what goes into the boundary rows? At real $\omega$, $L$ on all of space still has a null space; in vacuum, it is all plane waves with $|\mathbf k|=k_0$. So $L^{-1}$ does not exist yet. Physically, we know which solution we want: the field the source *sends out*, with nothing coming in from infinity. The radiation condition states this precisely.

## Choosing the outgoing solution

We work at one frequency, in vacuum, with the scalar equation $(\nabla^2+k_0^2)\,u=s$. This is exact for the 2D $E_z$ case, and the vector case changes little (see the end of this section).

### Why one point source is enough

We want to solve $Lu=s$, not for one source, but for any source we might use: a dipole, a line current, an incoming beam. Superposition gives a shortcut.

On the grid, a source is a vector $\mathbf b$, and every vector is a sum of single-point sources: $\mathbf b=\sum_jb_j\,\mathbf e_j$, where $\mathbf e_j$ is 1 at point $j$ and 0 elsewhere. Suppose we know the response to each single point, i.e. the solution $\mathbf g_j$ of $A\mathbf g_j=\mathbf e_j$. Then by superposition, the response to $\mathbf b$ is
$$\mathbf u=\sum_jb_j\,\mathbf g_j .$$
The $\mathbf g_j$ are exactly the columns of $A^{-1}$, as we saw in "Derivatives as matrices". So we solve for point sources once, and every other source is just a weighted sum.

In the continuous version, $\mathbf e_j$ becomes $\delta(x-x')$, and $s(x)=\int s(x')\,\delta(x-x')\,dx'$. The response to $\delta(x-x')$ is the **Green's function** $G(x,x')$, and
$$u(x)=\int G(x,x')\,s(x')\,dx' .$$
This is the $K$ from "Why not write the law as $u=Ks$?".

Homogeneity saves even more work. In vacuum, moving the source only moves the response, so $G(x,x')=G(x-x')$. One function $G(x)$, the response to a point source at the origin, then answers every source at once:
$$u(x)=\int G(x-x')\,s(x')\,dx' .$$
This is a convolution, and $G$ is the unit point source response $k$ of the [convolution](convolution.md) note. That is why we want $G$: it is the whole solution operator in one function.

### One dimension: shaking a rope

Picture a very long rope. Shake one point, $x=0$, up and down at frequency $\omega$, and keep shaking. After a while, the rope settles into a steady pattern at that frequency. That pattern is the 1D Green's function. (Small rope motions obey the 1D wave equation, so at one frequency they obey $u''+k_0^2u=\delta(x)$, with $k_0=\omega/v$ and $v$ the wave speed on the rope.)

Away from your hand, nothing pushes the rope, so the equation is source-free there, and on each side
$$u=a\,e^{ik_0x}+b\,e^{-ik_0x}.$$
That is 4 unknowns: $a,b$ on the left and $a,b$ on the right. Your hand gives 2 conditions at $x=0$. The rope does not break, so $u$ is continuous. And your hand pushes there, so the rope has a kink: the slope $u'$ jumps (by 1 here, found by integrating the equation across $0$). So 2 unknowns are still free. These are the two free waves, one moving right and one moving left, as the null space predicts.

The equation cannot fix them, so ask the rope. With the time factor $e^{-i\omega t}$, $e^{ik_0x}$ moves right and $e^{-ik_0x}$ moves left. Waves run away from your hand in both directions, and nothing comes back from the far ends. So on the right we keep only $e^{ik_0x}$, and on the left only $e^{-ik_0x}$. These are the 2 missing conditions, and they give
$$G(x)=\frac{e^{ik_0|x|}}{2ik_0}.$$
(Check: $G'$ jumps from $-\tfrac12$ to $+\tfrac12$ at $x=0$.)

Compare a rope tied at both ends. Now waves do come back, they pile up into standing waves, and at certain frequencies the rope resonates without limit. The equation is the same; only the far ends differ. So what fixes the free waves is not the law but what happens far away. We meet the tied rope again with the PML.

### The same problem in the $k$ picture

The free waves show up in the $k$ picture too. In vacuum, $L$ is homogeneous, so on a plane wave $e^{i\mathbf k\cdot\mathbf r}$ it is multiplication by its symbol $k_0^2-|\mathbf k|^2$. Inverting $L$ should then be division:
$$\hat u(\mathbf k)=\frac{\hat s(\mathbf k)}{k_0^2-|\mathbf k|^2}.$$
But the denominator is zero on the sphere $|\mathbf k|=k_0$. These are exactly the free waves. So the freedom we found on the rope appears here as a division by zero, and only there. Everywhere else in $\mathbf k$, inverting $L$ is easy.

### The radiation condition

"Far from all sources, there are only outgoing waves" is the **radiation condition**. It is the extra condition for open space, just as initial conditions are the extra condition in time. In 3D, its precise form is the Sommerfeld condition, $r\,(\partial_ru-ik_0u)\to0$ as $r\to\infty$.

Mindful readers might ask why "outgoing" is the right choice, and not just a taste. There are two honest answers, and they agree.
- **Causality.** In time, the response must come after the source (see the causality note in [convolution](convolution.md)). An incoming wave would have to start at infinity before the source was switched on. At one frequency, the retarded response is exactly the outgoing one.
- **A little loss.** Every real medium absorbs a little: $\operatorname{Im}\varepsilon_r>0$. So replace $k_0$ by $\kappa=k_0\sqrt{1+i\eta}$ with a tiny $\eta>0$; then $\operatorname{Im}\kappa>0$. Now $e^{i\kappa x}$ decays to the right and blows up to the left, and $e^{-i\kappa x}$ does the opposite. A bounded solution must keep only the decaying wave on each side, which is exactly the outgoing choice. In the $k$ picture, the denominator $\kappa^2-k^2$ is never zero for real $k$, so the division is well defined. Then let $\eta\to0$. This is called the **limiting absorption principle** (established).

So the null space is not removed by a mathematical trick. It is removed by physics we had already: causality, or equivalently a trace of loss.

### Three dimensions

For a point source in 3D, the field depends only on the distance $r$. Away from $r=0$, a short computation gives $\nabla^2u=\tfrac1r(ru)''$. So the equation becomes $(ru)''+k_0^2(ru)=0$, and $ru=e^{\pm ik_0r}$. The radiation condition keeps $e^{+ik_0r}$. Near $r=0$, $e^{ik_0r}\approx1$, and $\nabla^2\tfrac{1}{4\pi r}=-\delta$ (Coulomb's law) fixes the constant:
$$G(\mathbf r)=-\frac{e^{ik_0r}}{4\pi r}.$$
This is a spherical wave going out. Its amplitude falls like $1/r$, so its power falls like $1/r^2$, and the total power through any large sphere stays the same, as energy conservation requires. (In 2D, $G$ is a Hankel function; see [Green's function](Green's%20function.md).)

For the vector curl-curl equation, $G$ becomes a $3\times3$ matrix (the dyadic Green's function). It is built from the same $e^{ik_0r}/r$ and its derivatives. The extra derivative terms dominate near the source (the near field) and fade far away (the far field).

### A structure, and the $K$ form

In vacuum, the convolution $u=G*s$ from "Why one point source is enough" now uses the outgoing $G$. So this $K$ form has the radiation condition built in.

With a structure $\varepsilon_r(\mathbf r)$, there is no such formula, but the $K$ form survives in another way. Use the rule from "What counts as a source?" in reverse: move the structure to the source side,
$$\nabla^2u+k_0^2u=s-k_0^2(\varepsilon_r-1)\,u ,$$
so that
$$u=G*\big(s-k_0^2(\varepsilon_r-1)\,u\big).$$
This is an integral equation for $u$ (the Lippmann–Schwinger equation). Its unknowns live only where $\varepsilon_r\neq1$, and outgoing waves are automatic. The price is a dense matrix. Solvers of this family (volume integral equations, and MoM for metal surfaces) take this route.

### A finite grid: the PML

The differential route keeps $L$ sparse, but a grid ends somewhere. Now we can say what goes into the boundary rows: rows that let outgoing waves leave without reflection. The obvious choices fail.
- **$u=0$ at the edge** (what the matrix in "Derivatives as matrices" did) is a perfect mirror. This is the tied rope: the grid becomes a closed box, waves bounce back, and $L$ is singular at the box's resonant frequencies.
- **An "outgoing only" boundary row.** In 1D, $u'=ik_0u$ at the right edge is exact, since only one wave goes out. In 2D and 3D, waves arrive at every angle, and such a condition is exact only at normal incidence.

The **perfectly matched layer (PML)** uses the loss idea again, but only in a layer around the grid. Inside the layer, replace the real coordinate $x$ by a complex one:
$$\tilde x=x+\frac{i}{\omega}\int_0^x\sigma(x')\,dx',\qquad \frac{d}{d\tilde x}=\frac1{s(x)}\frac{d}{dx},\quad s=1+\frac{i\sigma}{\omega},$$
with $\sigma=0$ outside the layer and $\sigma>0$ inside. An outgoing wave becomes
$$e^{ik\tilde x}=e^{ikx}\,e^{-(k/\omega)\int\sigma\,dx},$$
which decays inside the layer, for every $k>0$, and hence for every angle. It does not reflect at the entrance, because nothing physical changes there. The stretched solution is the same function, evaluated at complex $x$. (Established: Bérenger 1994; Chew and Weedon showed it is a coordinate stretch.)

To be honest about its limits:
- The layer ends in a mirror, so a tiny reflection comes back after a round trip through the layer. It is small, not zero.
- On a grid, a sudden jump in $\sigma$ reflects. So $\sigma$ is ramped up gradually over 10 or more cells.
- Evanescent waves ($k=i\kappa$) only gain a phase in the layer, not extra decay. Variants such as the CFS-PML fix this.

On the 1D grid, the PML changes only the rows inside the layer. Row $i$ becomes
$$\frac{1}{s_i h}\Big[\frac{u_{i+1}-u_i}{s_{i+1/2}\,h}-\frac{u_i-u_{i-1}}{s_{i-1/2}\,h}\Big]+k_0^2\varepsilon_{r,i}\,u_i=b_i ,$$
where $b_i$ is the source and $s_i,s_{i\pm1/2}$ are the stretch factors. The matrix is still sparse, now complex, and invertible: the layer absorbs, so the null space is gone. One more detail matters later for the adjoint. After multiplying row $i$ by $s_i$, the matrix is symmetric ($A^T=A$) but not Hermitian. This is reciprocity on the grid.

So both routes now have a complete answer to "which solution?". The $K$ route builds the outgoing condition into $G$. The $L$ route builds it into the boundary rows, with a PML.

## From law to solver

Every solver turns the law into a finite problem a computer can do. Strip away the brand names, and each one makes three choices, all of which we have already met.

1. **Form.** Use $L$ (the differential law: sparse, needs boundary rows such as a PML) or $K$ (the Green's function: dense, outgoing waves built in, unknowns only on the structure).
2. **Domain.** Work at one frequency and solve $A\mathbf x=\mathbf b$, or work in time and step the curl equations forward.
3. **Basis.** Decide how a function becomes a vector: values on a grid, small local patches (elements), or global functions (plane waves, modes).

The common solvers are just different picks:

| Solver | Form | Domain | Basis | Typical use (my judgment; they all overlap) |
|---|---|---|---|---|
| FDFD | $L$ | frequency | grid (Yee) | small devices, inverse design |
| FDTD | $L$ | time | grid (Yee) | broadband runs, nanophotonics (Meep, Lumerical) |
| FIT | $L$ | time (or frequency) | grid, integral form | CST |
| FEM | $L$ | frequency, mostly | elements (edge elements) | curved, multiscale geometry (HFSS, COMSOL) |
| MoM / BEM | $K$ | frequency | patches on surfaces | metal antennas, scattering in open space (FEKO) |
| Modal / RCWA | mixed | frequency | modes, plane waves | periodic and layered structures |

### What each choice costs

- **$L$ or $K$.** $L$ gives a sparse matrix over the whole volume, plus a PML. $K$ gives a dense matrix over the structure only, with no PML. So a large metal antenna in empty space favours $K$: the unknowns sit only on the metal surface. A complicated dielectric volume favours $L$: the volume has to be meshed anyway, and sparsity wins.
- **Frequency or time.** At one frequency, a dispersive material is just a number $\varepsilon_r(\omega)$ ("Nonlocal in time"). The price is one large sparse solve per frequency. But once $A$ is factorized, every new source costs almost nothing, which the adjoint will exploit. In time, there is no linear solve at all, only cheap local updates, and one pulse covers a whole band of frequencies. The price is that the material's memory must be computed in time. The usual way is to step the Lorentz oscillator equation alongside the fields (the auxiliary differential equation method).
- **Basis.** A grid is simple, but a curved surface becomes a staircase. Elements follow curved surfaces and can be small where detail is needed, at the cost of a more complex code.

One more cost of the time domain follows from locality. Each update uses only neighbouring cells, so information moves at most about one cell per step. The physical wave must not outrun that, so the time step is limited: $c\,\Delta t\le h/\sqrt d$ in $d$ dimensions (the **CFL condition**). The exact constant comes from putting a plane wave into the update, i.e. the symbol again.

### The Yee grid

Grid solvers for Maxwell's equations almost all use the Yee grid. Its idea is to place each field component at a different point of the cell. Take the 2D $E_z$ case: put $E_z$ at the grid nodes $(i,j)$, $H_x$ at the edge midpoints $(i,j+\tfrac12)$, and $H_y$ at $(i+\tfrac12,j)$. Then every derivative in the curl equations is a centred difference over one cell. For example, Faraday's law gives
$$H_x(i,j+\tfrac12)=\frac{E_z(i,j+1)-E_z(i,j)}{i\omega\mu_0h}.$$
Put the $H$ values into Ampère's law, and the 2D $E_z$ equation becomes the 5-point stencil
$$\frac{E_{i+1,j}+E_{i-1,j}+E_{i,j+1}+E_{i,j-1}-4E_{i,j}}{h^2}+k_0^2\varepsilon_{r,ij}E_{i,j}=-i\omega\mu_0J_{i,j}.$$
This is the 1D matrix from "Derivatives as matrices" applied once in $x$ and once in $y$. With the PML rows, it is the whole 2D $E_z$ FDFD solver.

Why the staggering? It keeps "the curl of a gradient is zero" *exactly* on the grid. On the grid, a gradient lives on edges: it is the difference of the two end values. A curl lives on faces: it is the sum of the edge values around the face, with signs (Stokes' theorem on one cell). Go around a face, and every node value is added once and subtracted once, so the result is exactly zero. This is what the "Divergence equations are constraints" section needs on the grid. The curl equations then keep the discrete divergence fixed, and the false (spurious) modes of the low-frequency note cannot appear. Edge elements in FEM are built on the same idea.

### The grid has its own dispersion relation

The grid is shift-invariant too, but only by whole cells. So the homogeneity argument still works: $e^{ikx_i}$ is an eigenvector of the stencil. Its eigenvalue is
$$\frac{e^{ikh}-2+e^{-ikh}}{h^2}=-\frac{4}{h^2}\sin^2\frac{kh}{2}\;\approx\;-k^2+\frac{k^4h^2}{12}.$$
This is the symbol of the grid. Compare it with the true symbol $-k^2$. For small $kh$ they agree, but not exactly. So the grid's dispersion relation, $\tfrac{4}{h^2}\sin^2\tfrac{kh}{2}=k_0^2$, gives a slightly larger $k$ than $k_0$. Grid waves are slightly too slow. The relative error is about $(k_0h)^2/24$, and it adds up over distance as a phase error. With 20 points per wavelength, $k_0h\approx0.31$, and the phase error is about 0.25 rad after 10 wavelengths. This is **numerical dispersion**. It is why people quote "10 to 20 points per wavelength", and why large domains need more. When $k_0h>2$, the grid cannot carry the wave at all.

### How to trust a result

Each check below is a principle from earlier in this note.
- **Convergence.** Halve $h$. A second-order scheme should cut the error by about 4.
- **Energy balance (Poynting).** Power from the source = power out through a surface + power absorbed. In a lossless structure, all absorption should happen in the PML.
- **Reciprocity.** Swap the source and the probe. The answer should not change, because $A^T=A$.
- **The PML.** Make it thicker. The answer should not change.
- **A known case.** Compare with an analytic Green's function, such as $e^{ik_0|x|}/(2ik_0)$ in 1D.

## The adjoint: every gradient for one extra solve

The "What counts as a source?" note promised this: designing a device means changing $L$, not $s$, and asking how the solution responds. Now we can answer it.

### The problem

Divide the design region into $N$ cells, and let the design parameters be their permittivities $p_i=\varepsilon_{r,i}$. $N$ is easily $10^4$ to $10^6$. The forward problem is the FDFD system
$$A(\mathbf p)\,\mathbf x=\mathbf b .$$
We judge a design by one real number $J(\mathbf x)$, the **objective**: for example, the power reaching a probe, or the Poynting flux through a port. To improve the design, we want the gradient, $dJ/dp_i$ for every $i$.

The naive way is to change one $p_i$ a little, solve again, and see how $J$ changes. That costs $N+1$ solves per gradient. With $N=10^5$, this is hopeless.

### The chain rule, and where to put the brackets

Differentiate $A\mathbf x=\mathbf b$ with respect to $p_i$ ($\mathbf b$ does not depend on $p_i$):
$$\frac{\partial A}{\partial p_i}\mathbf x+A\frac{d\mathbf x}{dp_i}=0\quad\Longrightarrow\quad \frac{d\mathbf x}{dp_i}=-A^{-1}\frac{\partial A}{\partial p_i}\mathbf x .$$
$J$ is real but $\mathbf x$ is complex, so $J$ depends on both $\mathbf x$ and $\overline{\mathbf x}$. For real $p_i$ this gives $dJ=2\operatorname{Re}\big(\mathbf g^T d\mathbf x\big)$, where $\mathbf g^T=\partial J/\partial\mathbf x$ is a row vector. So
$$\frac{dJ}{dp_i}=-2\operatorname{Re}\Big(\mathbf g^T A^{-1}\,\frac{\partial A}{\partial p_i}\,\mathbf x\Big).$$
Everything now depends on where we put the brackets. Computing $A^{-1}\big(\tfrac{\partial A}{\partial p_i}\mathbf x\big)$ means one solve for each $i$: the naive method again. But $\mathbf g^TA^{-1}$ is the same for every $i$. So compute it once: define $\boldsymbol\lambda$ by
$$A^T\boldsymbol\lambda=\mathbf g .$$
Then every component of the gradient is
$$\frac{dJ}{dp_i}=-2\operatorname{Re}\Big(\boldsymbol\lambda^T\frac{\partial A}{\partial p_i}\mathbf x\Big).$$
That is the whole adjoint method: one forward solve for $\mathbf x$, one **adjoint solve** for $\boldsymbol\lambda$, and then all $N$ derivatives come from cheap products. The cost no longer depends on $N$. (In general, one output and many inputs favour this "brackets from the left" order, called reverse mode. Many outputs and one input favour the other order, called forward mode. Backpropagation in neural networks is reverse mode. Established.)

For FDFD, $p_i$ appears only in the diagonal term $k_0^2\varepsilon_{r,i}$ of row $i$, so $\partial A/\partial p_i=k_0^2\,\mathbf e_i\mathbf e_i^T$. The gradient becomes
$$\frac{dJ}{d\varepsilon_{r,i}}=-2k_0^2\operatorname{Re}\big(\lambda_i\,x_i\big).$$
The gradient at each cell is just the forward field times the adjoint field at that cell.

### What the adjoint field is

The adjoint equation is $A^T\boldsymbol\lambda=\mathbf g$. With the PML rows scaled as in "A finite grid: the PML", $A^T=A$ (reciprocity). So $\boldsymbol\lambda$ solves *the same problem* as the forward field, with a different source, $\mathbf g$. And $\mathbf g=\partial J/\partial\mathbf x$ is nonzero only where the objective looks: at the probe or the port. So the adjoint field is the field you get by putting a source at the probe and letting it radiate back into the device.

Here is a concrete case. Let the objective be the power at one probe cell $m$: $J=|x_m|^2=x_m\overline{x_m}$. Then $\partial J/\partial x_m=\overline{x_m}$, so $\mathbf g=\overline{x_m}\,\mathbf e_m$: a single point source at the probe.

This also explains *why* the product $\lambda_ix_i$ is the sensitivity. Change $\varepsilon_r$ at cell $i$ by a small $\delta$. By the rule in "What counts as a source?", run in reverse, this is a new small source at cell $i$, proportional to $\delta\,x_i$: the forward field arriving there. That small source radiates, and some of it reaches the probe. By reciprocity, the amount that travels from cell $i$ to the probe equals the amount that travels from the probe to cell $i$, which is $\lambda_i$. So the sensitivity is "field arriving from the source" times "field arriving from the probe". This is the adjoint from the convolution note at work: integrating over observers, i.e. sending the information back from the observation point to the sources. (Conventions differ: some texts use $A^H$ and get $\overline{\boldsymbol\lambda}$. The gradient is the same if you stay consistent.)

### The cost, honestly

- **Two solves, one factorization.** With a direct solver, the factorization of $A$ is the expensive part. Since $A^T=A$, the adjoint solve reuses it, so the gradient costs barely more than the forward solve ("What each choice costs").
- **Time domain.** The same idea works in FDTD: the adjoint simulation runs backward in time, the time reversal from the convolution note. The price is storing (or recomputing) the forward fields. In continuous form, this is the **adjoint-state method**, also used in seismic imaging.
- **A gradient is only a direction.** The design problem is not convex, so gradient steps find a good local design, not the best one. Real devices also need constraints: $\varepsilon_r$ must end up as one of the two real materials, and features must be large enough to fabricate. The usual tools are filtering and projection (keyword: **topology optimization**).

### Check it against finite differences

The adjoint gradient must agree with the naive one. For a few cells, change $p_i$ by a small $h$, solve again, and compare $\big(J(p_i+h)-J(p_i)\big)/h$ with the adjoint value. They should agree to several digits. If $h$ is too large, the finite-difference formula is wrong (truncation error); if $h$ is too small, round-off takes over. A scan over $h$ shows a valley of agreement in between. If the two disagree, the usual culprits are a missing $2\operatorname{Re}$, $A^T$ versus $A^H$, a missing conjugate in $\mathbf g$, or PML rows that were not scaled to make $A$ symmetric.

## What this note leaves out

This note followed one thread, from "what is a linear law?" to a solver and its gradient. Textbooks cover much more. Below are the topics I would add next, each with the place it attaches to this note.

- **Modes: $Lu=\lambda u$.** Every problem here was driven: given a source, find the field. The other half of linear algebra is the eigenproblem. Fix $\omega$ and look for waves $e^{i\beta z}$ along a waveguide: $\beta$ is an eigenvalue and the mode profile is an eigenvector (a mode solver). Or fix the geometry and look for the $\omega$ at which a source-free solution exists: cavity resonances, like the tied rope. With loss or a PML, these $\omega$ become complex (quasi-normal modes), the same complex frequency as the friction example in "The symbol and the dispersion relation". In photonic design, the objective is often the overlap of the field with a mode, so the adjoint source is the mode itself. Keywords: mode solver, quasi-normal modes, mode overlap.
- **Fourier optics and the diffraction limit.** Write a field in the plane $z=0$ as a sum of $e^{i(k_xx+k_yy)}$. The dispersion relation $|\mathbf k|=k_0$ then fixes $k_z=\sqrt{k_0^2-k_x^2-k_y^2}$ for each one. When $k_x^2+k_y^2>k_0^2$, $k_z$ is imaginary: the wave is evanescent and dies within about a wavelength. So detail with a period shorter than $\lambda$ never reaches a far detector. This is the diffraction limit, and it is the symbol again. Moving a field from one plane to the next is multiplication by $e^{ik_zz}$ (the angular spectrum method), and the image of a point is the point spread function ($K^*K$ in the convolution note). Keywords: angular spectrum, diffraction limit, numerical aperture, near-field imaging.
- **The equivalence principle and the far field.** The field outside a closed surface can be reproduced by currents on that surface alone: $\mathbf J_s=\hat{\mathbf n}\times\mathbf H$ and a (fictitious) magnetic current $\mathbf M_s=-\hat{\mathbf n}\times\mathbf E$. This is "What counts as a source?" once more: fields on a surface can be traded for sources. MoM uses it (the unknowns are surface currents). Simulations use it to inject a clean plane wave (total-field/scattered-field sources). And it lets a small near-field simulation predict an antenna's pattern far away: far from the source, the field $G*s$ of any source becomes $e^{ik_0r}/r$ times a function of direction, which is the radiation pattern. Directivity, gain, radiation resistance and input impedance follow from it and from the Poynting vector. Keywords: Love's equivalence principle, TF/SF source, near-to-far-field transform.
- **Interfaces.** Integrate Maxwell's equations across a thin layer at an interface. Tangential $\mathbf E$ and $\mathbf H$ are continuous; normal $\mathbf D$ and $\mathbf B$ jump by any surface charge or current. The $H_z$ case was one example. For a plane wave hitting a flat interface, the system is still shift-invariant *along* the interface, so $k_x$ cannot change (Proof 2). That is Snell's law. Matching the fields gives the Fresnel coefficients. When $k_x>n_2k_0$, the transmitted wave is evanescent: total internal reflection, which is what guides light in fibres and waveguides. Keywords: Fresnel equations, Brewster angle, total internal reflection.

More briefly:
- **Statics and potentials.** $\omega=0$, where low-frequency breakdown lives. The potentials $\phi$ and $\mathbf A$, gauge choice, and the Laplace and Poisson equations. Capacitance and inductance are the circuit limit. MoM codes are usually written with potentials (the mixed-potential integral equation).
- **Transmission lines and ports.** 1D waves of voltage and current. S-parameters are the matrix from incoming to outgoing waves at the ports: the radiation condition, written for a circuit. This is the daily language of RF.
- **Scattering.** Cross sections; Mie theory (the exact solution for a sphere, and a standard test case for solvers); the Born approximation (one step of the Lippmann–Schwinger equation).
- **Further out.** Momentum and radiation pressure (the Maxwell stress tensor); relativity (the speed of light is the same for every observer, and $\mathbf E$ and $\mathbf B$ become one object); anisotropic and gyrotropic media (biased ferrites); nonlinear optics (where Linearity fails); quantum optics (photons).
