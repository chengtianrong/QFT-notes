# QFT

## Contents

- [0. Mathematical Foundations](#qft-ch0)
  - [0.1 Fourier Transform](#fourier-transform)
  - [0.2 $\gamma$ Matrices: Weyl Representation](#gamma-matrices)
  - [0.3 IO(1,3) Group](#classical-poincare)
- [1. Classical Field Theory](#qft-ch1)
  - [1.1 Definition of Classical Fields](#classical-field)
  - [1.2 Dynamics of Classical Fields](#classical-dynamics)
    - [Noether's Theorem](#classical-noether)
  - [1.3 Dirac Constraint](#dirac-constraint)
  - [1.4 The Dirac Field $\psi$](#dirac-field)
- [2. Operator Quantization](#qft-ch2)
  - [2.1 Canonical Quantization](#canonical-quantization)
- [3. Path Integrals](#qft-ch3)
- [4. The Standard Model](#qft-ch4)
  - [4.1 Yang–Mills Theory](#yang-mills)
  - [4.2 Quantum Chromodynamics](#qcd)
- [5. New Physics Beyond the Standard Model](#qft-ch5)
- [6. Scattering Theory](#qft-ch6)
  - [6.2 Perturbation Theory](#perturbation-theory)

---

<a id="qft-ch0"></a>
## 0. Mathematical Foundations

<a id="fourier-transform"></a>
### 0.1 Fourier Transform

We use the following convention for the four-dimensional Fourier transform, with the factor $(2\pi)^{-4}$ included in the inverse transform:

$$
\begin{aligned}
\widetilde f(p)&=\int d^4x\,e^{+ip\cdot x}f(x),\\
f(x)&=\int\frac{d^4p}{(2\pi)^4}\,e^{-ip\cdot x}\widetilde f(p).
\end{aligned}
\tag{0.1}
$$

Here $x^\mu=(t,\mathbf x)$ and $p^\mu=(p^0,\mathbf p)$, and the phase is given by

$$
p\cdot x=p_\mu x^\mu=p^0t-\mathbf p\cdot\mathbf x.
\tag{0.2}
$$

The Fourier representations of the Dirac $\delta$ function are

$$
\begin{aligned}
\int\frac{d^4p}{(2\pi)^4}\,e^{-ip\cdot(x-y)}&=\delta^{(4)}(x-y),\\[6pt]
\int d^4x\,e^{i(p-q)\cdot x}&=(2\pi)^4\delta^{(4)}(p-q).
\end{aligned}
\tag{0.3}
$$


<br><br>

<a id="gamma-matrices"></a>
### 0.2 $\gamma$ Matrices: Weyl Representation

The Pauli matrices are

$$
\sigma^1=\begin{pmatrix}0&1\\1&0\end{pmatrix},\qquad
\sigma^2=\begin{pmatrix}0&-i\\i&0\end{pmatrix},\qquad
\sigma^3=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
\tag{0.4}
$$

<br><br>

**Properties of the Pauli Matrices**

For each $i=1,2,3$,

$$
(\sigma^i)^\dagger=\sigma^i,\qquad
(\sigma^i)^2=I_2,\qquad
(\sigma^i)^{-1}=\sigma^i,\qquad
\operatorname{tr}\sigma^i=0,\qquad
\det\sigma^i=-1.
\tag{0.5}
$$

Thus, each Pauli matrix is both Hermitian and unitary, with eigenvalues $+1$ and $-1$.

The product identity is

$$
\sigma^i\sigma^j=\delta^{ij}I_2+i\epsilon^{ijk}\sigma^k.
\tag{0.6}
$$

This gives the commutation and anticommutation relations, as well as the trace orthogonality relation:

$$
[\sigma^i,\sigma^j]=2i\epsilon^{ijk}\sigma^k,\qquad
\{\sigma^i,\sigma^j\}=2\delta^{ij}I_2,\qquad
\operatorname{tr}(\sigma^i\sigma^j)=2\delta^{ij}.
\tag{0.7}
$$

For any $\mathbf a,\mathbf b\in\mathbb R^3$, the product identity can also be written as

$$
(\mathbf a\cdot\boldsymbol\sigma)(\mathbf b\cdot\boldsymbol\sigma)
=(\mathbf a\cdot\mathbf b)I_2
+i(\mathbf a\times\mathbf b)\cdot\boldsymbol\sigma,
\qquad \boldsymbol\sigma=(\sigma^1,\sigma^2,\sigma^3).
\tag{0.8}
$$

<br><br>

**$\gamma$ matrices**

Pauli four-vectors

$$
\sigma^\mu:=(I_2,\boldsymbol\sigma),\qquad
\bar\sigma^\mu:=(I_2,-\boldsymbol\sigma),\qquad \mu=0,1,2,3.
\tag{0.9}
$$

Using them, the $\gamma$ matrices in the Weyl representation are defined as:

$$
\gamma^\mu:=
\begin{pmatrix}
0&\sigma^\mu\\
\bar\sigma^\mu&0
\end{pmatrix}_{4\times4}.
\tag{0.10}
$$

These matrices satisfy 

$$
\{\gamma^\mu,\gamma^\nu\}=2\eta^{\mu\nu}I_4.
\tag{0.11}
$$

$$
(\gamma^0)^2=I_4,\qquad
(\gamma^i)^2=-I_4,\\
(\gamma^0)^\dagger=\gamma^0,\qquad
(\gamma^i)^\dagger=-\gamma^i.
\tag{0.12}
$$

Hence each $\gamma^\mu$ is unitary:

$$
(\gamma^\mu)^\dagger\gamma^\mu=I_4.
\tag{0.13}
$$

Define

$$
\gamma^5:=i\gamma^0\gamma^1\gamma^2\gamma^3
=\begin{pmatrix}-I_2&0\\0&I_2\end{pmatrix}
.
\tag{0.14}
$$

It satisfies

$$
(\gamma^5)^2=I_4,\qquad
\{\gamma^5,\gamma^\mu\}=0,\qquad
(\gamma^5)^\dagger = (\gamma^5).

\tag{0.15}
$$

<br><br>

<a id="classical-poincare"></a>
### 0.3 IO(1,3) Group

**The Group O(1,3)** 

To distinguish a matrix $\Lambda$ from the corresponding tensor $\Lambda$ of coordinate transformation, we denote the matrix by $[\Lambda]$ and the tensor simply by $\Lambda$. We use the same notation for $\eta$.

The matrix $[\Lambda]$

$$
O(1,3):=\{[\Lambda]\in GL(4,\mathbb R)\mid [\Lambda]^T[\eta][\Lambda]=[\eta]\}.
\tag{0.16}
$$

where $[\eta]=\operatorname{diag}(1,-1,-1,-1)$. 

The tensor $\Lambda$ 

Given a specific diffeomorphism of Minkowskispacetime $\phi:x\mapsto x'$, we obtain a tensor

$$
\Lambda^\mu{}_\nu:=\frac{\partial x'^\mu}{\partial x^\nu}.
\tag{0.17}
$$

By the way, we can prove that the isometries of Minkowski spacetime must be linear up to a translation.

The entry of $[\Lambda]$ in row $\mu$ and column $\nu$ is the corresponding component of the tensor $\Lambda$:

$$
[\Lambda]_{\mu\nu}=\Lambda^\mu{}_\nu.
\tag{0.18}
$$

Similarly,

$$
[\eta]_{\mu\nu}=\eta_{\mu\nu}.
\tag{0.19}
$$

Therefore, equivalently, we obtain the following relationship:

$$
\eta_{\mu\nu}\Lambda^\mu{}_\rho\Lambda^\nu{}_\sigma=\eta_{\rho\sigma}.
\tag{0.20}
$$

where we have used the relation $([\Lambda]^T)_{\mu\nu}=\Lambda^\nu{}_\mu$. This is exactly the condition for $\phi$ to be an isometry of Minkowski spacetime.


The inverse and transpose of $[\Lambda]$ satisfy

$$
[\Lambda]^{-1}=[\eta][\Lambda]^T[\eta],\qquad
[\Lambda][\eta][\Lambda]^T=[\eta],
\tag{0.21}
$$

so $[\Lambda]^{-1}$ and $[\Lambda]^T$ also belong to $O(1,3)$. 

In components, the first relation reads

$$
(\Lambda^{-1})^\mu{}_\nu=\Lambda_\nu{}^\mu.
\tag{0.22}
$$

<br><br>

**Connectedness**

Taking the determinant of $[\Lambda]^T[\eta][\Lambda]=[\eta]$ gives

$$
(\det[\Lambda])^2=1\quad\Longrightarrow\quad\det[\Lambda]=\pm1.
\tag{0.23}
$$

Setting $\rho=\sigma=0$ in $\eta_{\mu\nu}\Lambda^\mu{}_\rho\Lambda^\nu{}_\sigma=\eta_{\rho\sigma}$ gives

$$
(\Lambda^0{}_0)^2=1+\sum_i(\Lambda^i{}_0)^2\ge1\quad\Longrightarrow\quad\Lambda^0{}_0\ge1\ \text{or}\ \Lambda^0{}_0\le-1.
\tag{0.24}
$$

We can prove that $O(1,3)$ has exactly these four connected components.

Define the time reversal $T$ and the space inversion $P$ by

$$
T:(t,\mathbf x)\mapsto(-t,\mathbf x),\qquad
P:(t,\mathbf x)\mapsto(t,-\mathbf x),
\tag{0.25}
$$

whose matrices are

$$
[T]=\operatorname{diag}(-1,1,1,1),\qquad
[P]=\operatorname{diag}(1,-1,-1,-1).
\tag{0.26}
$$

<div style="display: flex; justify-content: center;">

| Connected component | $\det[\Lambda]$ | $\Lambda^0{}_0$ | Representative |
|:---:|:---:|:---:|:---:|
| $\mathcal L^\uparrow_+$ | $+1$ | $\ge1$ | $I$ |
| $\mathcal L^\downarrow_+$ | $+1$ | $\le-1$ | $[PT]$ |
| $\mathcal L^\uparrow_-$ | $-1$ | $\ge1$ | $[P]$ |
| $\mathcal L^\downarrow_-$ | $-1$ | $\le-1$ | $[T]$ |

</div>

$SO(1,3)$ consists of the two classes with $\det[\Lambda]=+1$:

$$
SO(1,3)=\mathcal L^\uparrow_+\cup\mathcal L^\downarrow_+.
\tag{0.27}
$$

Furthermore, $SO(1,3)^+=\mathcal L^\uparrow_+$ is a subgroup of $O(1,3)$.

<br><br>

**Lie Algebra of O(1,3)**

Consider a one-parameter subgroup $e^{tX}$, $t\in\mathbb R$, of $O(1,3)$. It satisfies

$$
(e^{tX})^T[\eta]\,e^{tX}=[\eta].
\tag{0.28}
$$

Expanding the left-hand side to first order in $t$ gives

$$
(I+tX^T)[\eta](I+tX)+O(t^2)=[\eta]+t\left(X^T[\eta]+[\eta]X\right)+O(t^2),
\tag{0.29}
$$

so the first-order term must vanish:

$$
X^T[\eta]+[\eta]X=0.
\tag{0.30}
$$

Therefore,

$$
\mathfrak o(1,3):=\{X\in M_4(\mathbb R)\mid X^T[\eta]+[\eta]X=0\}.
\tag{0.31}
$$

Furthermore, this relation can be written more explicitly as

$$
[\eta]X=-X^T[\eta]=-([\eta]X)^T,
\tag{0.32}
$$

which means that the matrix $[\eta]X$ is antisymmetric. 

An antisymmetric $4\times4$ matrix has $4\times3/2=6$ independent real parameters. Hence the general form of $X$ is

$$
X=\begin{pmatrix}
0&b^1&b^2&b^3\\
b^1&0&-\theta^3&\theta^2\\
b^2&\theta^3&0&-\theta^1\\
b^3&-\theta^2&\theta^1&0
\end{pmatrix}
=\sum_i\theta^iJ^i+\sum_ib^iK^i
\tag{0.33}
$$

In physics, the following notation is commonly used:

$$
e^{-\frac i2\omega^{\mu\nu}J^{\mu\nu}}
\tag{0.34}
$$

where $\omega_{\mu\nu}$ is antisymmetric in its indices, and

$$
-iJ^{ij}=\epsilon_{ijk}J^k,\qquad -iJ^{0i}=K^i.
\tag{0.35}
$$

Each $J^{\mu\nu}$ is simply a $4\times4$ matrix, with $\alpha,\beta$ labelling its rows and columns. Its entries are

$$
(J^{\mu\nu})_{\alpha\beta}=i\left(\eta^{\mu\alpha}\delta^{\nu\beta}-\eta^{\nu\alpha}\delta^{\mu\beta}\right).
\tag{0.36}
$$

Explicitly (row index $\alpha$, column index $\beta$; the remaining components follow from $J^{\nu\mu}=-J^{\mu\nu}$), the rotation generators $J^{ij}=i\epsilon_{ijk}J^k$ are Hermitian:

$$
J^{23}=\begin{pmatrix}0&0&0&0\\0&0&0&0\\0&0&0&-i\\0&0&i&0\end{pmatrix},\quad
J^{31}=\begin{pmatrix}0&0&0&0\\0&0&0&i\\0&0&0&0\\0&-i&0&0\end{pmatrix},\quad
J^{12}=\begin{pmatrix}0&0&0&0\\0&0&-i&0\\0&i&0&0\\0&0&0&0\end{pmatrix},
\tag{0.37}
$$

and the boost generators $J^{0i}=iK^i$ are anti-Hermitian:

$$
J^{01}=\begin{pmatrix}0&i&0&0\\i&0&0&0\\0&0&0&0\\0&0&0&0\end{pmatrix},\quad
J^{02}=\begin{pmatrix}0&0&i&0\\0&0&0&0\\i&0&0&0\\0&0&0&0\end{pmatrix},\quad
J^{03}=\begin{pmatrix}0&0&0&i\\0&0&0&0\\0&0&0&0\\i&0&0&0\end{pmatrix}.
\tag{0.38}
$$





---

<a id="qft-ch1"></a>
## 1. Classical Field Theory

<a id="classical-field"></a>
### 1.1 Definition of Classical Fields

**Definition**

$$
\text{Field}:\text{Spacetime}\longrightarrow\text{Vector Space},
\tag{1.1}
$$

$$
\boxed{\Psi:\mathbb R^4\longrightarrow\mathbb R^n\ \text{or}\ \mathbb C^n
}
\tag{1.2}
$$

Additional regularity conditions, such as $C^\infty$ smoothness, may be imposed.

The space of fields $\{\Psi\}$ carries a representation of $ISO(1,3)^+$. For $(\Lambda,a)\in ISO(1,3)^+$,

$$
\Psi\xrightarrow{\;(\Lambda,a)\;}U(\Lambda,a)\Psi,
\tag{1.3}
$$

where $U$ is an infinite-dimensional representation that can be constructed from the finite-dimensional representation $D(\Lambda)$:

$$
\boxed{
\begin{aligned}
[U(\Lambda,a)\Psi]_\alpha(x)
&:=D_{\alpha\beta}(\Lambda)
[(\Lambda,a)_*\Psi]_\beta(x)\\
&=D_{\alpha\beta}(\Lambda)
\Psi_\beta(\Lambda^{-1}x-\Lambda^{-1}a).
\end{aligned}
}
\tag{1.4}
$$

**infinite-dimensional representations of $SO(1,3)$**


<br><br>

<a id="classical-dynamics"></a>
### 1.2 Dynamics of Classical Fields

Having defined classical fields, we now turn to their dynamics. In the Lagrangian formalism, this amounts to constructing a Lagrangian density for the fields.
In classical mechanics, the dynamical variables are the generalized coordinates $q_i(t)$, which depend only on time. In field theory, the dynamical variables are the field components $\Psi_\alpha(x)$, which depend on the spacetime point $x=(t,\mathbf x)$.

There are two ways to interpret $x$:

1. At fixed time $t$, the pair $(\alpha,\mathbf x)$ plays the role of the label $i$ in $q_i(t)$, so $\Psi_\alpha(t,\mathbf x)$ corresponds to $q_i(t)$.
2. More physically, the whole spacetime point $x$ plays the role of $t$ in classical mechanics, so $\Psi_\alpha(x)$ corresponds to $q_i(t)$.


**Conditions on the Lagrangian Density**

The following assumptions are commonly used to construct local relativistic field theories in Minkowski spacetime, subject to the qualifications listed below.

| Condition | Reason |
| --- | --- |
| Locality: dependence only on fields and finitely many derivatives at the same spacetime point | Yields local field equations; causal propagation requires additional conditions. |
| At most first derivatives (a standard choice) | Yields equations of motion of at most second order. The Ostrogradsky instability arises from nondegenerate dependence on higher time derivatives; higher derivatives are not excluded in general. |
| No explicit dependence on $x$ | Ensures translation invariance and local energy-momentum conservation; conservation of the total charges also requires suitable boundary conditions. |
| Lorentz scalar, up to a total divergence | Ensures Lorentz covariance of the field equations; boundary terms must respect the chosen boundary conditions. |
| Reality (Hermiticity for operator-valued fields) | Provides the standard starting point for a Hermitian Hamiltonian and unitary time evolution after consistent quantization. |
| Correct signs for the kinetic terms of physical degrees of freedom | Avoids negative kinetic energies. A lower bound on the total energy also constrains the interactions; for ordinary scalar fields, the potential must be bounded below. |
| Internal and gauge symmetries; renormalizability | Symmetry requirements depend on the theory. Power-counting renormalizability is optional in an effective field theory within its regime of validity. |

<br><br>
<a id="classical-noether"></a>
#### Noether's Theorem

Let $\phi_\lambda$ be a one-parameter family of transformations with $\phi_0=\mathrm{id}$, acting on spacetime as $x\mapsto x'=\phi_\lambda(x)$, with $x'^\mu=:x^\mu+\delta x^\mu$, and on the field as $\Psi\mapsto\Psi'$. If the action is invariant under $\phi_\lambda$ for every integration region $\Omega$ in (1.8), then the quantity

$$
\boxed{
\int d^3x\left(
\frac{\partial\mathcal L}{\partial\dot\Psi_\alpha}
\frac{\partial\,\delta x^\nu}{\partial\lambda}
\partial_\nu\Psi_\alpha
-\frac{\partial\mathcal L}{\partial\dot\Psi_\alpha}
\frac{\partial\,\delta_0\Psi_\alpha}{\partial\lambda}
-\mathcal L\frac{\partial\,\delta x^0}{\partial\lambda}
\right),
}
\tag{1.5}
$$

is conserved when the field $\Psi_\alpha$ satisfies the Euler–Lagrange equations and, together with its first derivatives, falls off fast enough at spatial infinity for the integral in (1.5) to converge and for the surface term in (1.12) to vanish. Here $\delta_0\Psi_\alpha$ is defined in (1.6), and the derivatives with respect to $\lambda$ are taken at $\lambda=0$.

**Proof**

**(A) Coordinate Transformation**

Let $x'=\phi_\lambda(x)$, with $\lambda$ infinitesimal. Define two field variations:

$$
\begin{aligned}
\delta_0\Psi_\alpha(x)
&:=\Psi'_\alpha(x')-\Psi_\alpha(x)
=D_{\alpha\beta}(\Lambda)\Psi_\beta(x)-\Psi_\alpha(x),\\
\delta\Psi_\alpha(x)
&:=\Psi'_\alpha(x)-\Psi_\alpha(x)
=[U(\Lambda,a)\Psi]_\alpha(x)-\Psi_\alpha(x).
\end{aligned}
\tag{1.6}
$$

Here $\Psi'$ is the transformed field. Using the pushforward $(\phi_\lambda)_*\Psi:=\Psi\circ\phi_\lambda^{-1}$, we write $\Psi'_\alpha=D_{\alpha\beta}(\Lambda)\,[(\phi_\lambda)_*\Psi]_\beta$. Evaluating at $x'$ gives $\Psi'_\alpha(x')=D_{\alpha\beta}(\Lambda)\Psi_\beta(x)$, which is the last equality in the first line, so the first variation, $\delta_0\Psi_\alpha$, describes the transformation of the field itself. For $\phi_\lambda=(\Lambda,a)$, $\Psi'=U(\Lambda,a)\Psi$ by (1.4), which is the last equality in the second line. If $\phi_\lambda$ is not a Poincaré transformation, $D(\Lambda)$ is replaced by the identity and $\delta_0\Psi_\alpha=0$.

To first order in $\lambda$, the two variations are related by

$$
\begin{aligned}
\delta\Psi_\alpha(x)
&=\bigl[\Psi'_\alpha(x')-\Psi_\alpha(x)\bigr]
-\bigl[\Psi'_\alpha(x')-\Psi'_\alpha(x)\bigr]\\
&=\delta_0\Psi_\alpha(x)-\delta x^\mu\partial_\mu\Psi'_\alpha(x)\\
&=\delta_0\Psi_\alpha(x)-\delta x^\mu\partial_\mu\Psi_\alpha(x),
\end{aligned}
\tag{1.7}
$$

where the last step uses $\Psi'_\alpha=\Psi_\alpha+O(\lambda)$.

Over a spacetime region $\Omega$, the actions before and after the transformation are

$$
\begin{aligned}
I&=\int_\Omega d^4x\,
\mathcal L\bigl(\Psi_\alpha(x),\partial_\mu\Psi_\alpha(x)\bigr),\\
I'&=\int_{\phi_\lambda(\Omega)}d^4x'\,
\mathcal L\bigl(\Psi'_\alpha(x'),\partial'_\mu\Psi'_\alpha(x')\bigr),
\end{aligned}
\tag{1.8}
$$

where $\partial'_\mu:=\partial/\partial x'^\mu$. The same function $\mathcal L$ is used in $I'$: the transformation acts on the field configuration and on the region $\Omega$, not on the Lagrangian. No property of $\mathcal L$ is assumed here.

Change the integration variable in $I'$ from $x'$ back to $x$, so that $\phi_\lambda(\Omega)$ becomes $\Omega$. To first order in $\lambda$,

$$
\begin{aligned}
d^4x'
&=\det\!\left(\frac{\partial x'^\mu}{\partial x^\nu}\right)d^4x
=\bigl(1+\partial_\mu\delta x^\mu\bigr)d^4x,\\
\partial'_\mu
&=\frac{\partial x^\nu}{\partial x'^\mu}\partial_\nu
=\partial_\mu-(\partial_\mu\delta x^\nu)\partial_\nu,\\
\Psi'_\alpha(x')
&=\Psi_\alpha(x)+\delta_0\Psi_\alpha(x),\\
\partial'_\mu\Psi'_\alpha(x')
&=\partial_\mu\Psi_\alpha(x)
+\partial_\mu\delta_0\Psi_\alpha(x)
-(\partial_\mu\delta x^\nu)\partial_\nu\Psi_\alpha(x),
\end{aligned}
\tag{1.9}
$$

where the third line is (1.6). Hence

$$
I'=\int_\Omega d^4x\,
\bigl(1+\partial_\mu\delta x^\mu\bigr)\,
\mathcal L\bigl(
\Psi_\alpha+\delta_0\Psi_\alpha,\ 
\partial_\mu\Psi_\alpha
+\partial_\mu\delta_0\Psi_\alpha
-(\partial_\mu\delta x^\nu)\partial_\nu\Psi_\alpha
\bigr).
\tag{1.10}
$$

Expanding (1.10) to first order in $\lambda$ gives the off-shell identity

$$
\begin{aligned}
I'-I
&=\int_\Omega d^4x\,\biggl[
\mathcal L\,\partial_\mu\delta x^\mu
+\frac{\partial\mathcal L}{\partial\Psi_\alpha}\delta_0\Psi_\alpha\\
&\qquad
+\frac{\partial\mathcal L}{\partial(\partial_\mu\Psi_\alpha)}
\bigl(\partial_\mu\delta_0\Psi_\alpha-(\partial_\mu\delta x^\nu)\partial_\nu\Psi_\alpha\bigr)
\biggr]\\
&=\int_\Omega d^4x\,\biggl\{
\biggl[
\frac{\partial\mathcal L}{\partial\Psi_\alpha}
-\partial_\mu\frac{\partial\mathcal L}{\partial(\partial_\mu\Psi_\alpha)}
\biggr]\delta\Psi_\alpha\\
&\qquad
+\partial_\mu\biggl[
\mathcal L\,\delta x^\mu
+\frac{\partial\mathcal L}{\partial(\partial_\mu\Psi_\alpha)}\delta\Psi_\alpha
\biggr]
\biggr\}.
\end{aligned}
\tag{1.11}
$$

The second equality substitutes $\delta_0\Psi_\alpha=\delta\Psi_\alpha+\delta x^\nu\partial_\nu\Psi_\alpha$ from (1.7) and uses the product rule together with $\partial_\nu\mathcal L=\frac{\partial\mathcal L}{\partial\Psi_\alpha}\partial_\nu\Psi_\alpha+\frac{\partial\mathcal L}{\partial(\partial_\mu\Psi_\alpha)}\partial_\nu\partial_\mu\Psi_\alpha$; the latter holds because $\mathcal L$ has no explicit $x$-dependence.

By hypothesis, the action is invariant: $I'=I$. Take $\Omega=[t_1,t_2]\times V$ and let $\Psi_\alpha$ satisfy the Euler–Lagrange equations, so that the Euler–Lagrange term in (1.11) vanishes. Gauss's theorem then gives

$$
\begin{aligned}
0
&=\int_\Omega d^4x\,\partial_\mu\biggl[
\mathcal L\,\delta x^\mu
+\frac{\partial\mathcal L}{\partial(\partial_\mu\Psi_\alpha)}\delta\Psi_\alpha
\biggr]\\
&=\biggl[\int_V d^3x\,\biggl(
\mathcal L\,\delta x^0
+\frac{\partial\mathcal L}{\partial\dot\Psi_\alpha}\delta\Psi_\alpha
\biggr)\biggr]_{t_1}^{t_2}\\
&\quad
+\int_{t_1}^{t_2}dt\oint_{\partial V}dS_i\,\biggl(
\mathcal L\,\delta x^i
+\frac{\partial\mathcal L}{\partial(\partial_i\Psi_\alpha)}\delta\Psi_\alpha
\biggr).
\end{aligned}
\tag{1.12}
$$

By the fall-off condition, the surface term vanishes as $V\to\mathbb R^3$. Since $t_1$ and $t_2$ are arbitrary,

$$
\int d^3x\,\biggl(
\mathcal L\,\delta x^0
+\frac{\partial\mathcal L}{\partial\dot\Psi_\alpha}\delta\Psi_\alpha
\biggr)=\text{const}.
\tag{1.13}
$$

Substituting $\delta\Psi_\alpha=\delta_0\Psi_\alpha-\delta x^\nu\partial_\nu\Psi_\alpha$ from (1.7), differentiating with respect to $\lambda$, and reversing the overall sign gives (1.5) for coordinate transformations.

**(B) Internal Transformation**

Now $\phi_\lambda$ acts trivially on spacetime, and on the field by a constant matrix $g(\lambda)$ with $g(0)=I_n$:

$$
x'^\mu=x^\mu,\qquad
\Psi'_\alpha(x)=g_{\alpha\beta}(\lambda)\,\Psi_\beta(x).
\tag{1.14}
$$

Hence $\delta x^\mu=0$, and the two variations in (1.6) coincide: $\delta\Psi_\alpha(x)=\delta_0\Psi_\alpha(x)=g_{\alpha\beta}(\lambda)\Psi_\beta(x)-\Psi_\alpha(x)$. Repeating (1.8)–(1.13) with $\delta x^\mu=0$: in (1.9), $d^4x'=d^4x$ and $\partial'_\mu=\partial_\mu$; the terms containing $\delta x^\mu$ drop out of (1.10)–(1.12); and (1.13) becomes

$$
\int d^3x\,\frac{\partial\mathcal L}{\partial\dot\Psi_\alpha}\,\delta_0\Psi_\alpha=\text{const}.
\tag{1.15}
$$

Differentiating with respect to $\lambda$ and reversing the overall sign gives (1.5) with $\delta x^\mu=0$. $\square$











<br><br>

<a id="dirac-constraint"></a>
### 1.3 Dirac Constraint













<br><br>
<a id="dirac-field"></a>
### 1.4 The Dirac Field $\psi$

The Lagrangian of the Dirac Field has the following form:

$$
\boxed{
\mathcal L=\bar\psi\bigl(i\gamma^\mu\partial_\mu-m\bigr)\psi,
}
\tag{1.16}
$$

where $\psi:\mathbb R^4\to\mathbb C^4$ is the Dirac field, $\bar\psi:=\psi^\dagger\gamma^0$ is its Dirac adjoint, the $\gamma^\mu$ are given by (0.10), and $m$ is the mass. 

Treating $\psi$ and $\bar\psi$ as independent, the Euler–Lagrange equation for $\bar\psi$ is the Dirac equation:

$$
\boxed{
\bigl(i\gamma^\mu\partial_\mu-m\bigr)\psi=0.
}
\tag{1.17}
$$

<br><br>

**Weyl Representation**

For the Dirac field, $D(\Lambda)$ in (1.4) is

$$
\Lambda_{\frac12}=\exp\Bigl(-\frac i2\omega_{\mu\nu}S^{\mu\nu}\Bigr),
\qquad
S^{\mu\nu}:=\frac i4[\gamma^\mu,\gamma^\nu],
\tag{1.18}
$$

with the same $\omega_{\mu\nu}$ as in (0.34). In the Weyl representation every $S^{\mu\nu}$ is block diagonal:

$$
S^{0i}=-\frac i2\begin{pmatrix}\sigma^i&0\\0&-\sigma^i\end{pmatrix},
\qquad
S^{ij}=\frac12\epsilon^{ijk}\begin{pmatrix}\sigma^k&0\\0&\sigma^k\end{pmatrix}.
\tag{1.19}
$$

Hence $\Lambda_{\frac12}$ is block diagonal: writing $\psi=\begin{pmatrix}\psi_L\\\psi_R\end{pmatrix}$, the left-handed and right-handed Weyl spinors $\psi_L,\psi_R$ transform independently. For an infinitesimal transformation, with $\omega_{ij}=\epsilon_{ijk}\theta^k$ and $\omega_{0i}=b^i$ as in (0.33),

$$
\begin{aligned}
\Lambda_{\frac12}
&\simeq I_4-\frac i2\omega_{\mu\nu}S^{\mu\nu}\\
&=I_4-\frac i2\begin{pmatrix}\boldsymbol\theta\cdot\boldsymbol\sigma&0\\0&\boldsymbol\theta\cdot\boldsymbol\sigma\end{pmatrix}
-i\Bigl(-\frac i2\Bigr)\begin{pmatrix}\mathbf b\cdot\boldsymbol\sigma&0\\0&-\mathbf b\cdot\boldsymbol\sigma\end{pmatrix}\\
&=\begin{pmatrix}I_2-i\boldsymbol\theta\cdot\frac{\boldsymbol\sigma}2-\mathbf b\cdot\frac{\boldsymbol\sigma}2&0\\0&I_2-i\boldsymbol\theta\cdot\frac{\boldsymbol\sigma}2+\mathbf b\cdot\frac{\boldsymbol\sigma}2\end{pmatrix},
\end{aligned}
\tag{1.20}
$$

so that

$$
\begin{aligned}
\psi_L&\xrightarrow{\ \Lambda\ }\Bigl(I_2-i\boldsymbol\theta\cdot\frac{\boldsymbol\sigma}2-\mathbf b\cdot\frac{\boldsymbol\sigma}2\Bigr)\psi_L,\\
\psi_R&\xrightarrow{\ \Lambda\ }\Bigl(I_2-i\boldsymbol\theta\cdot\frac{\boldsymbol\sigma}2+\mathbf b\cdot\frac{\boldsymbol\sigma}2\Bigr)\psi_R.
\end{aligned}
\tag{1.21}
$$

In terms of $\psi_L,\psi_R$, the Dirac equation (1.17) reads

$$
\begin{pmatrix}-m&i\sigma^\mu\partial_\mu\\i\bar\sigma^\mu\partial_\mu&-m\end{pmatrix}
\begin{pmatrix}\psi_L\\\psi_R\end{pmatrix}=0
\quad\Longleftrightarrow\quad
\begin{cases}
i\sigma^\mu\partial_\mu\psi_R-m\psi_L=0,\\
i\bar\sigma^\mu\partial_\mu\psi_L-m\psi_R=0.
\end{cases}
\tag{1.22}
$$

For $m=0$ it reduces to the Weyl equations

$$
\begin{cases}
i\bigl(\partial_t-\boldsymbol\sigma\cdot\nabla\bigr)\psi_L=0,\\
i\bigl(\partial_t+\boldsymbol\sigma\cdot\nabla\bigr)\psi_R=0.
\end{cases}
\tag{1.23}
$$

<br><br>

**Plane-Wave Solutions of the Dirac Equation**

Apply $-i\gamma^\mu\partial_\mu-m$ to the left-hand side of (1.17):

$$
\begin{aligned}
\bigl(-i\gamma^\mu\partial_\mu-m\bigr)\bigl(i\gamma^\nu\partial_\nu-m\bigr)\psi
&=\bigl(\gamma^\mu\gamma^\nu\partial_\mu\partial_\nu+m^2\bigr)\psi\\
&=\bigl(\eta^{\mu\nu}\partial_\mu\partial_\nu+m^2\bigr)\psi\\
&=\bigl(\Box+m^2\bigr)\psi,
\end{aligned}
\tag{1.24}
$$

where $\Box:=\eta^{\mu\nu}\partial_\mu\partial_\nu=\partial_t^2-\nabla^2$. Since the left-hand side of (1.24) vanishes by (1.17), every component of $\psi$ satisfies the Klein–Gordon equation

$$
\bigl(\Box+m^2\bigr)\psi=0.
\tag{1.25}
$$

A plane wave $e^{-ip\cdot x}$ satisfies (1.25) only if $p^2=m^2$. In the Weyl representation (0.10), the plane-wave solutions of (1.17) therefore take the form

$$
\psi(x)=\begin{pmatrix}u_L\\u_R\end{pmatrix}e^{-ip\cdot x},
\qquad
p^0=\pm E_{\mathbf p},\quad
E_{\mathbf p}:=\sqrt{\mathbf p^2+m^2},
\tag{1.26}
$$

where $u_L,u_R\in\mathbb C^2$ are constant. Substituting (1.26) into (1.17) gives

$$
\bigl(\not p-m\bigr)\begin{pmatrix}u_L\\u_R\end{pmatrix}=0,
\qquad
\not p:=\gamma\cdot p=p^\mu\gamma_\mu.
\tag{1.27}
$$

**(1) Rest Frame**

Take $m>0$ and work in the rest frame, $\mathbf p=\mathbf 0$; there we denote the constants in (1.26) by $u_L(\mathbf 0)$ and $u_R(\mathbf 0)$. First consider $p^0>0$. Then $p^\mu=(m,\mathbf 0)$, $e^{-ip\cdot x}=e^{-imt}$ and $\not p=m\gamma^0$, so (1.27) becomes

$$
\begin{aligned}
\bigl(m\gamma^0-m\bigr)\begin{pmatrix}u_L(\mathbf 0)\\u_R(\mathbf 0)\end{pmatrix}
&=\begin{pmatrix}-mI_2&mI_2\\mI_2&-mI_2\end{pmatrix}
\begin{pmatrix}u_L(\mathbf 0)\\u_R(\mathbf 0)\end{pmatrix}=0\\[6pt]
\Longrightarrow\quad u_L(\mathbf 0)&=u_R(\mathbf 0).
\end{aligned}
\tag{1.28}
$$

Writing $u:=\begin{pmatrix}u_L\\u_R\end{pmatrix}$, there are two independent solutions

$$
u^s(\mathbf 0)=\sqrt m\begin{pmatrix}\xi^s\\\xi^s\end{pmatrix},\qquad s=1,2,
\tag{1.29}
$$

where $\{\xi^s\}$ is a basis of $\mathbb C^2$ with $\xi^{r\dagger}\xi^s=\delta^{rs}$.

**(2) Boost**

Boost along $z$ with rapidity $\eta$, i.e. $\omega_{03}=\eta$ in (1.18), so that $\mathbf p=(0,0,p^3)$ and $\tanh\eta=v=p^3/E_{\mathbf p}$. Then

$$
\begin{aligned}
\Lambda_{\frac12}=e^{-i\eta S^{03}}
&=\exp\Bigl[-i\eta\Bigl(-\frac i2\Bigr)\begin{pmatrix}\sigma^3&0\\0&-\sigma^3\end{pmatrix}\Bigr]
=\begin{pmatrix}e^{-\frac12\eta\sigma^3}&0\\0&e^{\frac12\eta\sigma^3}\end{pmatrix}\\
&=\operatorname{diag}\bigl(e^{-\eta/2},\,e^{\eta/2},\,e^{\eta/2},\,e^{-\eta/2}\bigr),
\end{aligned}
\tag{1.30}
$$

$$
\tanh\eta=v
\ \Longrightarrow\
e^{\pm\eta}=\sqrt{\frac{1\pm v}{1\mp v}}=\frac{E_{\mathbf p}\pm p^3}{m},
\qquad
e^{\pm\eta/2}=\sqrt{\frac{E_{\mathbf p}\pm p^3}{m}}.
\tag{1.31}
$$

Hence

$$
u^s(p)=\Lambda_{\frac12}u^s(\mathbf 0)
=\operatorname{diag}\Bigl(\sqrt{E_{\mathbf p}-p^3},\,\sqrt{E_{\mathbf p}+p^3},\,\sqrt{E_{\mathbf p}+p^3},\,\sqrt{E_{\mathbf p}-p^3}\Bigr)
\begin{pmatrix}\xi^s\\\xi^s\end{pmatrix}.
\tag{1.32}
$$

With $\xi^1=\begin{pmatrix}1\\0\end{pmatrix}$ and $\xi^2=\begin{pmatrix}0\\1\end{pmatrix}$,

$$
u^1(p)=\begin{pmatrix}\sqrt{E_{\mathbf p}-p^3}\\0\\\sqrt{E_{\mathbf p}+p^3}\\0\end{pmatrix}
\xrightarrow{\ \eta\to\infty\ }\sqrt{2E_{\mathbf p}}\begin{pmatrix}0\\0\\1\\0\end{pmatrix},
\qquad
u^2(p)=\begin{pmatrix}0\\\sqrt{E_{\mathbf p}+p^3}\\0\\\sqrt{E_{\mathbf p}-p^3}\end{pmatrix}
\xrightarrow{\ \eta\to\infty\ }\sqrt{2E_{\mathbf p}}\begin{pmatrix}0\\1\\0\\0\end{pmatrix}.
\tag{1.33}
$$

**(3) Arbitrary Momentum**

Since $\operatorname{diag}\bigl(\sqrt{E_{\mathbf p}-p^3},\sqrt{E_{\mathbf p}+p^3}\bigr)=\sqrt{E_{\mathbf p}-p^3\sigma^3}$ and $\operatorname{diag}\bigl(\sqrt{E_{\mathbf p}+p^3},\sqrt{E_{\mathbf p}-p^3}\bigr)=\sqrt{E_{\mathbf p}+p^3\sigma^3}$, for arbitrary $\mathbf p$

$$
\sqrt m\,\Lambda_{\frac12}=\begin{pmatrix}\sqrt{p\cdot\sigma}&0\\0&\sqrt{p\cdot\bar\sigma}\end{pmatrix},
\qquad
p\cdot\sigma:=p_\mu\sigma^\mu=p^0I_2-\mathbf p\cdot\boldsymbol\sigma,
\quad
p\cdot\bar\sigma:=p_\mu\bar\sigma^\mu=p^0I_2+\mathbf p\cdot\boldsymbol\sigma,
\tag{1.34}
$$

with positive square roots. The positive-energy solutions are

$$
\boxed{
\psi=u^s(p)\,e^{-ip\cdot x},\qquad
u^s(p)=\begin{pmatrix}\sqrt{p\cdot\sigma}\,\xi^s\\\sqrt{p\cdot\bar\sigma}\,\xi^s\end{pmatrix},\qquad
p^0=E_{\mathbf p}.
}
\tag{1.35}
$$

**Verification**

$$
\begin{aligned}
(p\cdot\bar\sigma)(p\cdot\sigma)
&=\bigl(p^0I_2+p^i\sigma^i\bigr)\bigl(p^0I_2-p^j\sigma^j\bigr)
=(p^0)^2-\frac12\{\sigma^i,\sigma^j\}p^ip^j
=(p^0)^2-\mathbf p^2=m^2,\\[6pt]
\bigl(\not p-m\bigr)u^s(p)
&=\begin{pmatrix}-m&p\cdot\sigma\\p\cdot\bar\sigma&-m\end{pmatrix}
\begin{pmatrix}\sqrt{p\cdot\sigma}\,\xi^s\\\sqrt{p\cdot\bar\sigma}\,\xi^s\end{pmatrix}
=\begin{pmatrix}\sqrt{p\cdot\sigma}\,\bigl(\sqrt{(p\cdot\sigma)(p\cdot\bar\sigma)}-m\bigr)\xi^s\\\sqrt{p\cdot\bar\sigma}\,\bigl(\sqrt{(p\cdot\bar\sigma)(p\cdot\sigma)}-m\bigr)\xi^s\end{pmatrix}=0.
\end{aligned}
\tag{1.36}
$$

<br><br>

**Negative-Energy Solutions**

For $p^0=-E_{\mathbf p}$, (1.27) becomes

$$
\bigl(-\gamma^0E_{\mathbf p}-\boldsymbol\gamma\cdot\mathbf p-m\bigr)u(p)=0
\quad\Longleftrightarrow\quad
\bigl(\gamma^0E_{\mathbf p}-\boldsymbol\gamma\cdot(-\mathbf p)+m\bigr)u(p)=0,
\tag{1.37}
$$

which is the positive-energy equation $\bigl(\gamma^0E_{\mathbf p}-\boldsymbol\gamma\cdot\mathbf p-m\bigr)u=0$ with $\mathbf p\to-\mathbf p$ and $m\to-m$. Hence

$$
u^s(p)=\begin{pmatrix}\sqrt{E_{\mathbf p}I_2+\mathbf p\cdot\boldsymbol\sigma}\,\xi^s\\-\sqrt{E_{\mathbf p}I_2-\mathbf p\cdot\boldsymbol\sigma}\,\xi^s\end{pmatrix}
=\begin{pmatrix}\sqrt{-p\cdot\sigma}\,\xi^s\\-\sqrt{-p\cdot\bar\sigma}\,\xi^s\end{pmatrix},
\qquad p^0=-E_{\mathbf p}.
\tag{1.38}
$$

**Antiparticle Solutions**

Define

$$
u^c(p):=\begin{pmatrix}i\sigma^2&0\\0&-i\sigma^2\end{pmatrix}\gamma^0\,u(p),
\qquad
i\sigma^2=\begin{pmatrix}0&1\\-1&0\end{pmatrix}.
\tag{1.39}
$$

For $\mathbf p=(0,0,p^3)$,

$$
i\sigma^2\sqrt{p\cdot\bar\sigma}
=\begin{pmatrix}0&1\\-1&0\end{pmatrix}\begin{pmatrix}\sqrt{E_{\mathbf p}+p^3}&0\\0&\sqrt{E_{\mathbf p}-p^3}\end{pmatrix}
=\begin{pmatrix}0&\sqrt{E_{\mathbf p}-p^3}\\-\sqrt{E_{\mathbf p}+p^3}&0\end{pmatrix}
=\sqrt{p\cdot\sigma}\;i\sigma^2,
\tag{1.40}
$$

and likewise $i\sigma^2\sqrt{p\cdot\sigma}=\sqrt{p\cdot\bar\sigma}\;i\sigma^2$. Applying (1.39) to (1.35) and writing $\eta^s:=i\sigma^2\xi^s$,

$$
u^c(p)=\begin{pmatrix}i\sigma^2\sqrt{p\cdot\bar\sigma}\,\xi^s\\-i\sigma^2\sqrt{p\cdot\sigma}\,\xi^s\end{pmatrix}
=\begin{pmatrix}\sqrt{p\cdot\sigma}\,\eta^s\\-\sqrt{p\cdot\bar\sigma}\,\eta^s\end{pmatrix},
\qquad p^0=E_{\mathbf p},
\tag{1.41}
$$

which is (1.38) at $-p$ with $\xi^s\to\eta^s$. For arbitrary $\mathbf p$, the antiparticle solutions are

$$
\boxed{
\psi=v^s(p)\,e^{+ip\cdot x},\qquad
v^s(p):=\begin{pmatrix}\sqrt{p\cdot\sigma}\,\eta^s\\-\sqrt{p\cdot\bar\sigma}\,\eta^s\end{pmatrix},\qquad
\bigl(\not p+m\bigr)v^s(p)=0,\qquad
p^0=E_{\mathbf p}.
}
\tag{1.42}
$$

<br><br>

**Properties of the Solutions**

**(1) Orthonormality**

From (1.35), $u^\dagger u=2E_{\mathbf p}\,\xi^\dagger\xi$ and $\bar uu=2m\,\xi^\dagger\xi$. With $\xi^{r\dagger}\xi^s=\eta^{r\dagger}\eta^s=\delta^{rs}$,

$$
\boxed{
\begin{aligned}
u^{r\dagger}(p)\,u^s(p)&=2E_{\mathbf p}\,\delta^{rs}, &\qquad \bar u^r(p)\,u^s(p)&=2m\,\delta^{rs},\\
v^{r\dagger}(p)\,v^s(p)&=2E_{\mathbf p}\,\delta^{rs}, &\qquad \bar v^r(p)\,v^s(p)&=-2m\,\delta^{rs},\\
u^{r\dagger}(p)\,v^s(p')&=v^{r\dagger}(p')\,u^s(p)=0, &\qquad p'&:=(p^0,-\mathbf p).
\end{aligned}
}
\tag{1.43}
$$

**(2) Completeness (Spin Sums)**

With $\sum_s\xi^s\xi^{s\dagger}=I_2$,

$$
\boxed{
\sum_s u^s(p)\,\bar u^s(p)=\not p+m,\qquad
\sum_s v^s(p)\,\bar v^s(p)=\not p-m.
}
\tag{1.44}
$$

<br><br>

**Dirac Bilinears**

Bilinears have the form $\bar\psi\Gamma\psi$ with $\Gamma$ a $4\times4$ matrix. A basis of the $4\times4$ matrices consists of 16 matrices:

<div style="display: flex; justify-content: center;">

| $\Gamma$ | Number | Replaced by |
|:---:|:---:|:---:|
| $I_4$ | 1 | |
| $\gamma^\mu$ | 4 | |
| $S^{\mu\nu}$ | 6 | |
| $\gamma^{[\mu}\gamma^\nu\gamma^{\rho]}$ | 4 | $\gamma^5\gamma^\mu$ |
| $\gamma^{[\mu}\gamma^\nu\gamma^\rho\gamma^{\sigma]}$ | 1 | $\gamma^5$ |

</div>

where

$$
\gamma^5=\frac{i}{4!}\epsilon_{\mu\nu\rho\sigma}\gamma^\mu\gamma^\nu\gamma^\rho\gamma^\sigma,
\qquad \epsilon_{0123}=+1,
\tag{1.45}
$$

which agrees with (0.14). Under $\psi\to\Lambda_{\frac12}\psi$, $\bar\psi\to\bar\psi\Lambda_{\frac12}^{-1}$,

$$
\begin{aligned}
\bar\psi\psi&\xrightarrow{\ \Lambda\ }\bar\psi\psi,\\
\bar\psi\gamma^\mu\psi&\xrightarrow{\ \Lambda\ }\bar\psi\Lambda_{\frac12}^{-1}\gamma^\mu\Lambda_{\frac12}\psi=\Lambda^\mu{}_\nu\,\bar\psi\gamma^\nu\psi,\\
\bar\psi S^{\mu\nu}\psi&\xrightarrow{\ \Lambda\ }\Lambda^\mu{}_\rho\Lambda^\nu{}_\sigma\,\bar\psi S^{\rho\sigma}\psi,\\
\bar\psi\gamma^5\gamma^\mu\psi&\xrightarrow{\ \Lambda\ }\det[\Lambda]\,\Lambda^\mu{}_\nu\,\bar\psi\gamma^5\gamma^\nu\psi\qquad\text{(pseudovector)},\\
\bar\psi\gamma^5\psi&\xrightarrow{\ \Lambda\ }\det[\Lambda]\,\bar\psi\gamma^5\psi\qquad\text{(pseudoscalar)}.
\end{aligned}
\tag{1.46}
$$

**Proof** (for $\gamma^5$)

$$
\begin{aligned}
\Lambda_{\frac12}^{-1}\gamma^5\Lambda_{\frac12}
&=\frac{i}{4!}\epsilon_{\mu\nu\rho\sigma}\Lambda^\mu{}_\alpha\Lambda^\nu{}_\beta\Lambda^\rho{}_\kappa\Lambda^\sigma{}_\lambda\,\gamma^\alpha\gamma^\beta\gamma^\kappa\gamma^\lambda\\
&=\frac{i}{4!}\det[\Lambda]\,\epsilon_{\alpha\beta\kappa\lambda}\gamma^\alpha\gamma^\beta\gamma^\kappa\gamma^\lambda
=\det[\Lambda]\,\gamma^5.
\end{aligned}
\tag{1.47}
$$

<br><br>

**Chirality**

$$
\gamma^5\begin{pmatrix}\psi_L\\0\end{pmatrix}=-\begin{pmatrix}\psi_L\\0\end{pmatrix},
\qquad
\gamma^5\begin{pmatrix}0\\\psi_R\end{pmatrix}=\begin{pmatrix}0\\\psi_R\end{pmatrix}.
\tag{1.48}
$$

The left- and right-handed projection operators are

$$
P_L:=\frac{I_4-\gamma^5}2,\qquad P_R:=\frac{I_4+\gamma^5}2
\quad\Longrightarrow\quad
\begin{cases}
P_L+P_R=I_4,\\
P_LP_R=P_RP_L=0,\\
P_L^n=P_L,\quad P_R^n=P_R.
\end{cases}
\tag{1.49}
$$

<br><br>

**Vector and Axial Currents**

$$
j^\mu:=\bar\psi\gamma^\mu\psi,\qquad
j^{\mu5}:=\bar\psi\gamma^\mu\gamma^5\psi.
\tag{1.50}
$$

The Euler–Lagrange equation of (1.16) for $\psi$ is $-i\partial_\mu\bar\psi\gamma^\mu-m\bar\psi=0$. Together with (1.17),

$$
\partial_\mu j^\mu=(\partial_\mu\bar\psi)\gamma^\mu\psi+\bar\psi\gamma^\mu\partial_\mu\psi
=im\bar\psi\psi-im\bar\psi\psi=0
\quad\Longrightarrow\quad
\int d^3x\,\bar\psi\gamma^0\psi=\text{const}.
\tag{1.51}
$$

<br><br>

**Helicity**

$$
h:=\frac{\mathbf p}{|\mathbf p|}\cdot\mathbf S,\qquad
\mathbf S:=\frac12\begin{pmatrix}\boldsymbol\sigma&0\\0&\boldsymbol\sigma\end{pmatrix}.
\tag{1.52}
$$

For $\mathbf p=(0,0,p^3)$ with $p^3>0$, the spinors (1.33) satisfy $h\,u^1=\frac12u^1$ and $h\,u^2=-\frac12u^2$.






<a id="qft-ch2"></a>
## 2. Operator Quantization


The first viewpoint in [Section 1.2](#classical-dynamics) naturally leads to Operator quantization: spatial position labels the degrees of freedom, while time remains the evolution parameter. 



<br><br>
<a id="canonical-quantization"></a>
### 2.1 Canonical Quantization

Consider unconstrained  fields $\phi_a(x)$ with a Lagrangian density

$$
\mathcal L=\mathcal L(\phi_a,\partial_\mu\phi_a),
\qquad
L(t)=\int d^3x\,\mathcal L.
\tag{2.1}
$$

At each fixed time, the field values form an infinite set of generalized coordinates:

$$
i\longleftrightarrow(a,\mathbf x),
\qquad
q_i(t)\longleftrightarrow\phi_a(t,\mathbf x).
\tag{2.2}
$$

**Conjugate momentum**

The momentum density conjugate to $\phi_a(t,\mathbf x)$ is

$$
\pi_a(t,\mathbf x)
\equiv
\frac{\delta L(t)}{\delta\dot\phi_a(t,\mathbf x)}
=
\frac{\partial\mathcal L}{\partial\dot\phi_a(t,\mathbf x)}.
\tag{2.3}
$$

**Poisson brackets and commutation relations**

The classical fields obey the equal-time canonical Poisson brackets

$$
\begin{aligned}
\{\phi_a(t,\mathbf x),\pi_b(t,\mathbf y)\}_{\mathrm{PB}}
&=\delta_{ab}\delta^{(3)}(\mathbf x-\mathbf y),\\
\{\phi_a(t,\mathbf x),\phi_b(t,\mathbf y)\}_{\mathrm{PB}}
&=\{\pi_a(t,\mathbf x),\pi_b(t,\mathbf y)\}_{\mathrm{PB}}=0.
\end{aligned}
\tag{2.4}
$$

Canonical quantization promotes the fields and their conjugate momenta to operators and replaces these fundamental Poisson brackets according to

$$
\{\ ,\ \}_{\mathrm{PB}}
\longrightarrow
\frac{1}{i}[\ ,\ ].
\tag{2.5}
$$

Thus,

$$
\boxed{
\begin{aligned}
[\hat\phi_a(t,\mathbf x),\hat\pi_b(t,\mathbf y)]
&=i\,\delta_{ab}\delta^{(3)}(\mathbf x-\mathbf y),\\
[\hat\phi_a(t,\mathbf x),\hat\phi_b(t,\mathbf y)]
&=[\hat\pi_a(t,\mathbf x),\hat\pi_b(t,\mathbf y)]=0.
\end{aligned}}
\tag{2.6}
$$

Here $\delta_{ab}\delta^{(3)}(\mathbf x-\mathbf y)$ generalizes the Kronecker delta $\delta_{ij}$ to the field labels $(a,\mathbf x)$.

<a id="qft-ch3"></a>
## 3. Path Integrals

<a id="qft-ch4"></a>
## 4. The Standard Model

<a id="yang-mills"></a>
### 4.1 Yang–Mills Theory

**Lie Group**

$$
\begin{aligned}
U(n)&=\{A\in M_n(\mathbb C)\mid A^\dagger A=I\},\\[6pt]
SU(n)&=\{A\in M_n(\mathbb C)\mid \det A=1,\ A^\dagger A=I\}.
\end{aligned}
\tag{4.1}
$$



**$U(n)$**

Consider the one-parameter subgroup

$$
\begin{aligned}
g(\theta)&=e^{i\theta T},\qquad \theta\in\mathbb R,\quad T\in M_n(\mathbb C),\\[6pt]
\left.\frac{d}{d\theta}g(\theta)\right|_{\theta=0}&=iT,\\[6pt]
g^\dagger(\theta)g(\theta)=I\quad\text{for all }\theta\in\mathbb R
&\Longrightarrow i(T-T^\dagger)=0
\Longrightarrow T^\dagger=T.
\end{aligned}
\tag{4.6}
$$

 The Lie algebra is

$$
\mathfrak{u}(n)=\{iT\in M_n(\mathbb C)\mid T^\dagger=T\}.
$$




**$SU(n)$**

The one-parameter subgroup and its generator are

$$
\begin{aligned}
g(\theta)&=e^{i\theta T},\qquad \theta\in\mathbb R,\quad T\in M_n(\mathbb C),\\[6pt]
\left.\frac{d}{d\theta}g(\theta)\right|_{\theta=0}&=iT.
\end{aligned}
\tag{4.7}
$$


$$
0=\left.\frac{d}{d\theta}\bigl[g^\dagger(\theta)g(\theta)\bigr]\right|_{\theta=0}
=i(T-T^\dagger)
\quad\Longrightarrow\quad T^\dagger=T.
$$

We use the theorem $\det(e^{i\theta T})=e^{i\theta\operatorname{tr}T}$, valid for any $T\in M_n(\mathbb C)$ and $\theta\in\mathbb C$.

$$
0=\left.\frac{d}{d\theta}\det g(\theta)\right|_{\theta=0}
=\left.\frac{d}{d\theta}e^{i\theta\operatorname{tr}T}\right|_{\theta=0}
=i\operatorname{tr}T
\quad\Longrightarrow\quad \operatorname{tr}T=0.
$$

Thus $T$ is traceless and Hermitian, and the Lie algebra is

$$
\mathfrak{su}(n)=\{iT\in M_n(\mathbb C)\mid
\operatorname{tr}T=0,\ T^\dagger=T\}.
\tag{4.8}
$$

Scalar multiplication is taken over $\mathbb R$ to preserve anti-Hermiticity; hence $\mathfrak{su}(n)$ and  $\mathfrak{u}(n)$ are both real Lie algebras.

<br><br>

<a id="ym-4"></a>
**Structure Constants**

Let $\{iT^a\}$, with $a=1,\ldots,n^2-1$, be a basis of the real Lie algebra $\mathfrak{su}(n)$.
The structure constants $f^{abc}$ are defined by $[T^a,T^b]=if^{abc}T^c$.

**Reality: $f^{abc}\in\mathbb R$.**

Since the matrices $T^a$ are Hermitian,

$$
\begin{aligned}
([T^a,T^b])^\dagger&=-[T^a,T^b],\qquad
\operatorname{tr}[T^a,T^b]=0,\\[6pt]
[T^a,T^b]&\in\mathfrak{su}(n),\\[6pt]
[T^a,T^b]&=f^{abc}(iT^c),\qquad f^{abc}\in\mathbb R.
\end{aligned}
\tag{4.17}
$$

The expansion coefficients are real because $\{iT^c\}$ is a basis over $\mathbb R$.

<br><br>

<a id="ym-orthonormal"></a>
**Orthonormal Basis of $\mathfrak{su}(n)$**

For $iT_1,iT_2\in\mathfrak{su}(n)$, define

$$
\begin{aligned}
(\cdot,\cdot)&:\mathfrak{su}(n)\times\mathfrak{su}(n)\longrightarrow\mathbb R,\\[6pt]
(iT_1,iT_2)&:=-2\operatorname{tr}\bigl[(iT_1)(iT_2)\bigr]\\
&=2\operatorname{tr}(T_1T_2).
\end{aligned}
\tag{4.18}
$$

This form is bilinear over $\mathbb R$ by linearity of the trace. Since $T_1$ and $T_2$ are Hermitian, it is real-valued:

$$
\begin{aligned}
\bigl[\operatorname{tr}(T_1T_2)\bigr]^*
&=\operatorname{tr}\bigl[(T_1T_2)^\dagger\bigr]\\
&=\operatorname{tr}(T_2T_1)
=\operatorname{tr}(T_1T_2).
\end{aligned}
$$

Cyclicity of the trace also gives symmetry:

$$
\begin{aligned}
(iT_1,iT_2)&=2\operatorname{tr}(T_1T_2)\\
&=2\operatorname{tr}(T_2T_1)=(iT_2,iT_1).
\end{aligned}
$$

Finally, the form is positive definite:

$$
\begin{aligned}
(iT,iT)&=2\operatorname{tr}(T^\dagger T)\\
&=2\sum_{i,j}|T_{ij}|^2\ge0,
\end{aligned}
$$

with equality if and only if $T=0$.

Having established that this form defines an inner product on $\mathfrak{su}(n)$, we may choose an orthonormal basis $\{iT^a\}$ satisfying $(iT^a,iT^b)=\delta^{ab}$. Equivalently,

$$
\boxed{\operatorname{tr}(T^aT^b)=\frac12\delta^{ab}.}
\tag{4.19}
$$

**Total antisymmetry.** With this normalization and cyclicity of the trace, we obtain

$$
\begin{aligned}
\operatorname{tr}([T^a,T^b]T^c)&=\frac{i}{2}f^{abc},\\
f^{abc}&=-2i\operatorname{tr}([T^a,T^b]T^c),\\
f^{bac}&=-f^{abc},\\
f^{acb}
&=-2i\operatorname{tr}(T^aT^cT^b-T^cT^aT^b)\\
&=-2i\operatorname{tr}(T^bT^aT^c-T^aT^bT^c)
=-f^{abc}.
\end{aligned}
\tag{4.20}
$$

$$
\boxed{
f^{abc}=f^{[abc]}\in\mathbb R
}
\tag{4.21}
$$

<br><br>

<a id="ym-5"></a>
**Local Symmetry and Covariant Derivative**

Fiber bundles provide a useful geometric picture of fields in Yang–Mills theory.
In this picture, an internal space is attached to each spacetime point, and the field takes its value in that space.
This makes an additional freedom explicit: we may choose a basis in each internal space,
and this choice can vary from point to point.
This local freedom of choice underlies gauge symmetry.

Building on this idea, we now develop the framework of gauge theory.

The field operator is defined as

$$
\Psi(x):=\begin{pmatrix}
\psi_1(x)\\
\vdots\\
\psi_n(x)
\end{pmatrix}.
\tag{4.22}
$$

Define the Dirac adjoint by

$$
\begin{aligned}
\bar\Psi(x)&:=\bigl(\bar\psi_1(x),\ldots,\bar\psi_n(x)\bigr),\\[6pt]
\bar\psi_a(x)&:=\psi_a^\dagger(x)\gamma^0.
\end{aligned}
\tag{4.23}
$$

The matter Lagrangian density is

$$
\begin{aligned}
\mathcal L_{\mathrm{matter}}
&=\bar\Psi\,i\gamma^\mu I_n\partial_\mu\Psi\\[6pt]
&=\sum_{a=1}^{n}\bar\psi_a\,i\gamma^\mu\partial_\mu\psi_a.
\end{aligned}
\tag{4.24}
$$

<br><br>

<a id="qcd"></a>
### 4.2 Quantum Chromodynamics

**Color Group $SU(3)$**

The gauge group of QCD is $SU(3)$. By (4.8) with $n=3$, $\mathfrak{su}(3)$ has a basis $\{iT^a\}$, $a=1,\ldots,8$, which we take orthonormal in the sense of (4.19). The standard choice is $T^a:=\tfrac12\lambda^a$, where the Gell-Mann matrices are

$$
\begin{aligned}
\lambda^1&=\begin{pmatrix}0&1&0\\1&0&0\\0&0&0\end{pmatrix},&
\lambda^2&=\begin{pmatrix}0&-i&0\\i&0&0\\0&0&0\end{pmatrix},&
\lambda^3&=\begin{pmatrix}1&0&0\\0&-1&0\\0&0&0\end{pmatrix},\\[6pt]
\lambda^4&=\begin{pmatrix}0&0&1\\0&0&0\\1&0&0\end{pmatrix},&
\lambda^5&=\begin{pmatrix}0&0&-i\\0&0&0\\i&0&0\end{pmatrix},&
\lambda^6&=\begin{pmatrix}0&0&0\\0&0&1\\0&1&0\end{pmatrix},\\[6pt]
\lambda^7&=\begin{pmatrix}0&0&0\\0&0&-i\\0&i&0\end{pmatrix},&
\lambda^8&=\frac1{\sqrt3}\begin{pmatrix}1&0&0\\0&1&0\\0&0&-2\end{pmatrix}.
\end{aligned}
\tag{4.25}
$$

Each $\lambda^a$ is Hermitian and traceless, and $\operatorname{tr}(\lambda^a\lambda^b)=2\delta^{ab}$, so $\operatorname{tr}(T^aT^b)=\tfrac12\delta^{ab}$ as in (4.19). By (4.20), $f^{abc}=-2i\operatorname{tr}([T^a,T^b]T^c)$; up to the antisymmetry (4.21), the nonzero structure constants are

$$
\begin{aligned}
f^{123}&=1,\\[4pt]
f^{147}=f^{246}=f^{257}=f^{345}&=\tfrac12,\qquad
f^{156}=f^{367}=-\tfrac12,\\[4pt]
f^{458}=f^{678}&=\tfrac{\sqrt3}{2}.
\end{aligned}
\tag{4.26}
$$

The quadratic Casimir $T^aT^a$ commutes with every $T^b$, since $f^{bac}$ is antisymmetric in $a,c$:

$$
[T^b,T^aT^a]=[T^b,T^a]T^a+T^a[T^b,T^a]=if^{bac}\bigl(T^cT^a+T^aT^c\bigr)=0.
$$

The defining representation of $SU(3)$ is irreducible, so by Schur's lemma $T^aT^a=C_FI_3$. Taking the trace with (4.19),

$$
3C_F=\operatorname{tr}(T^aT^a)=\tfrac12\cdot8
\quad\Longrightarrow\quad
\boxed{T^aT^a=C_FI_3,\qquad C_F=\tfrac43.}
\tag{4.27}
$$

<br><br>

**Quark Fields**

For each flavor $f\in\{u,d,s,c,b,t\}$ the quark field is a color triplet, i.e. the multiplet (4.22) with $n=3$, with Dirac adjoint as in (4.23):

$$
q_f(x):=\begin{pmatrix}\psi_{f,1}(x)\\\psi_{f,2}(x)\\\psi_{f,3}(x)\end{pmatrix},
\qquad
\bar q_f(x):=\bigl(\bar\psi_{f,1}(x),\bar\psi_{f,2}(x),\bar\psi_{f,3}(x)\bigr).
\tag{4.28}
$$

Throughout this section $\gamma^\mu$ acts on the spinor index of each $\psi_{f,a}$ and the $3\times3$ matrices $T^a$ act on the color index, so $\gamma^\mu$ commutes with every matrix built from the $T^a$.

<br><br>

**Local $SU(3)$ Transformation**

Let $\theta^a:\mathbb R^4\to\mathbb R$, $a=1,\ldots,8$, and set

$$
g(x):=e^{i\theta^a(x)T^a}\in SU(3).
\tag{4.29}
$$

This is the internal transformation (1.14) with $g$ now depending on $x$. The quark fields transform as

$$
q_f'(x)=g(x)\,q_f(x),\qquad
\bar q_f'(x)=\bar q_f(x)\,g^\dagger(x).
\tag{4.30}
$$

Since $\partial_\mu(gq_f)=g\,\partial_\mu q_f+(\partial_\mu g)\,q_f$, the derivative term of (4.24) is not invariant under (4.30):

$$
\bar q_f'\,i\gamma^\mu\partial_\mu q_f'
=\bar q_f\,i\gamma^\mu\partial_\mu q_f
+\bar q_f\,i\gamma^\mu g^\dagger(\partial_\mu g)\,q_f.
\tag{4.31}
$$

<br><br>

**Covariant Derivative and Gluon Fields**

Introduce eight real vector fields $A^a_\mu$, the gluon fields, and define

$$
A_\mu:=A^a_\mu T^a,\qquad
D_\mu:=\partial_\mu-ig_sA_\mu,
\tag{4.32}
$$

where $g_s\in\mathbb R$ is the strong coupling constant and $D_\mu$ acts on color triplets. Let $A_\mu'$ denote the transformed gluon field and $D_\mu':=\partial_\mu-ig_sA_\mu'$. We require $D_\mu'q_f'=g\,D_\mu q_f$:

$$
\begin{aligned}
D_\mu'(gq_f)&=g\,\partial_\mu q_f+(\partial_\mu g)\,q_f-ig_sA_\mu'\,g\,q_f,\\[6pt]
g\,D_\mu q_f&=g\,\partial_\mu q_f-ig_s\,gA_\mu\,q_f.
\end{aligned}
\tag{4.33}
$$

Equating the two for all $q_f$ gives $A_\mu'g=gA_\mu-\frac{i}{g_s}\partial_\mu g$, i.e.

$$
\boxed{
A_\mu'=gA_\mu g^\dagger-\frac{i}{g_s}(\partial_\mu g)\,g^\dagger
=gA_\mu g^\dagger+\frac{i}{g_s}\,g\,\partial_\mu g^\dagger.
}
\tag{4.34}
$$

Differentiating $gg^\dagger=I_3$ gives $\bigl((\partial_\mu g)g^\dagger\bigr)^\dagger=-(\partial_\mu g)g^\dagger$, and Jacobi's formula $\partial_\mu\det g=\det g\,\operatorname{tr}(g^\dagger\partial_\mu g)$ with $\det g=1$ gives $\operatorname{tr}\bigl((\partial_\mu g)g^\dagger\bigr)=0$. Hence $A_\mu'$ is again Hermitian and traceless, $A_\mu'=A_\mu'^aT^a$ with real $A_\mu'^a$.

Since $D_\mu'(gq_f)=g\,D_\mu q_f$ holds for all $q_f$, as operators on color triplets

$$
D_\mu'=g\,D_\mu\,g^\dagger.
\tag{4.35}
$$

For infinitesimal $\theta^a$, $g=I_3+i\theta^aT^a+O(\theta^2)$ and (4.34) gives

$$
\begin{aligned}
A_\mu'&=A_\mu+i\theta^b[T^b,A_\mu]+\frac1{g_s}(\partial_\mu\theta^a)T^a+O(\theta^2)\\
&=A_\mu-f^{bca}\theta^bA^c_\mu T^a+\frac1{g_s}(\partial_\mu\theta^a)T^a+O(\theta^2),\\[6pt]
\delta A^a_\mu:=A_\mu'^a-A_\mu^a
&=\frac1{g_s}\partial_\mu\theta^a+f^{abc}A^b_\mu\theta^c+O(\theta^2).
\end{aligned}
\tag{4.36}
$$

<br><br>

**Field Strength**

Define $F_{\mu\nu}$ by $[D_\mu,D_\nu]=:-ig_sF_{\mu\nu}$. Since $[\partial_\mu,A_\nu]q_f=(\partial_\mu A_\nu)q_f$,

$$
\begin{aligned}
[D_\mu,D_\nu]
&=[\partial_\mu-ig_sA_\mu,\ \partial_\nu-ig_sA_\nu]\\
&=-ig_s\bigl(\partial_\mu A_\nu-\partial_\nu A_\mu\bigr)-g_s^2[A_\mu,A_\nu],\\[6pt]
F_{\mu\nu}&=\partial_\mu A_\nu-\partial_\nu A_\mu-ig_s[A_\mu,A_\nu].
\end{aligned}
\tag{4.37}
$$

Expanding in the basis, $F_{\mu\nu}=F^a_{\mu\nu}T^a$ with

$$
\boxed{
F^a_{\mu\nu}=\partial_\mu A^a_\nu-\partial_\nu A^a_\mu+g_sf^{abc}A^b_\mu A^c_\nu.
}
\tag{4.38}
$$

By (4.35), $[D_\mu',D_\nu']=g[D_\mu,D_\nu]g^\dagger$, so

$$
F_{\mu\nu}'=gF_{\mu\nu}g^\dagger,
\qquad
\delta F^a_{\mu\nu}=f^{abc}F^b_{\mu\nu}\theta^c+O(\theta^2).
\tag{4.39}
$$

By cyclicity of the trace and (4.19),

$$
\operatorname{tr}\bigl(F_{\mu\nu}'F'^{\mu\nu}\bigr)
=\operatorname{tr}\bigl(F_{\mu\nu}F^{\mu\nu}\bigr)
=\tfrac12F^a_{\mu\nu}F^{a\mu\nu}.
\tag{4.40}
$$

<br><br>

**The QCD Lagrangian**

$$
\boxed{
\mathcal L_{\mathrm{QCD}}
=\sum_{f}\bar q_f\bigl(i\gamma^\mu D_\mu-m_f\bigr)q_f
-\frac12\operatorname{tr}\bigl(F_{\mu\nu}F^{\mu\nu}\bigr)
=\sum_{f}\bar q_f\bigl(i\gamma^\mu D_\mu-m_f\bigr)q_f
-\frac14F^a_{\mu\nu}F^{a\mu\nu},
}
\tag{4.41}
$$

where the sum runs over the six flavors in (4.28), $m_f$ is the mass of flavor $f$, and $D_\mu$ acts on the $q_f$ to its right. Under (4.30) and (4.34), $D_\mu'q_f'=g\,D_\mu q_f$ by construction, so

$$
\bar q_f'\bigl(i\gamma^\mu D_\mu'-m_f\bigr)q_f'
=\bar q_f\,g^\dagger\bigl(i\gamma^\mu g\,D_\mu-m_f\,g\bigr)q_f
=\bar q_f\bigl(i\gamma^\mu D_\mu-m_f\bigr)q_f.
\tag{4.42}
$$

Together with (4.40), $\mathcal L_{\mathrm{QCD}}$ is invariant under the local transformations (4.30), (4.34).

Substituting (4.32) and (4.38) into (4.41) separates the free and interaction terms:

$$
\begin{aligned}
\mathcal L_{\mathrm{QCD}}
={}&\sum_f\bar q_f\bigl(i\gamma^\mu\partial_\mu-m_f\bigr)q_f
-\frac14\bigl(\partial_\mu A^a_\nu-\partial_\nu A^a_\mu\bigr)\bigl(\partial^\mu A^{a\nu}-\partial^\nu A^{a\mu}\bigr)\\[6pt]
&+g_sA^a_\mu\sum_f\bar q_f\gamma^\mu T^aq_f
-\frac{g_s}{2}f^{abc}\bigl(\partial_\mu A^a_\nu-\partial_\nu A^a_\mu\bigr)A^{b\mu}A^{c\nu}
-\frac{g_s^2}{4}f^{abc}f^{ade}A^b_\mu A^c_\nu A^{d\mu}A^{e\nu}.
\end{aligned}
\tag{4.43}
$$

<br><br>

**Equations of Motion**

Treating $q_f$, $\bar q_f$ and $A^a_\mu$ as independent, the Euler–Lagrange equation for $\bar q_f$ is

$$
\bigl(i\gamma^\mu D_\mu-m_f\bigr)q_f=0.
\tag{4.44}
$$

For $A^a_\nu$, define the color current $j^{a\nu}:=\sum_f\bar q_f\gamma^\nu T^aq_f$. From (4.41) with (4.38) and (4.43),

$$
\begin{aligned}
\frac{\partial\mathcal L_{\mathrm{QCD}}}{\partial(\partial_\mu A^a_\nu)}
&=-\frac12F^{b\rho\sigma}\,\delta^{ab}\bigl(\delta^\mu_\rho\delta^\nu_\sigma-\delta^\mu_\sigma\delta^\nu_\rho\bigr)
=-F^{a\mu\nu},\\[6pt]
\frac{\partial\mathcal L_{\mathrm{QCD}}}{\partial A^a_\nu}
&=g_sj^{a\nu}-\frac12F^{b\rho\sigma}\,g_s\bigl(f^{bad}\delta^\nu_\rho A^d_\sigma+f^{bca}A^c_\rho\delta^\nu_\sigma\bigr)\\
&=g_sj^{a\nu}-g_sf^{abc}A^c_\mu F^{b\mu\nu}\\
&=g_sj^{a\nu}+g_sf^{abc}A^b_\mu F^{c\mu\nu}.
\end{aligned}
\tag{4.45}
$$

The Euler–Lagrange equation $\partial_\mu\bigl(\partial\mathcal L_{\mathrm{QCD}}/\partial(\partial_\mu A^a_\nu)\bigr)=\partial\mathcal L_{\mathrm{QCD}}/\partial A^a_\nu$ then reads

$$
\boxed{
\partial_\mu F^{a\mu\nu}+g_sf^{abc}A^b_\mu F^{c\mu\nu}=-g_sj^{a\nu}.
}
\tag{4.46}
$$

The left-hand side is the $a$-component of $\partial_\mu F^{\mu\nu}-ig_s[A_\mu,F^{\mu\nu}]$, the covariant derivative of $F^{\mu\nu}$ in the adjoint representation; with $j^\nu:=j^{a\nu}T^a$, (4.46) is

$$
D_\mu F^{\mu\nu}:=\partial_\mu F^{\mu\nu}-ig_s[A_\mu,F^{\mu\nu}]=-g_s\,j^\nu.
\tag{4.47}
$$

<a id="qft-ch5"></a>
## 5. New Physics Beyond the Standard Model

<a id="qft-ch6"></a>
## 6. Scattering Theory

<a id="perturbation-theory"></a>
### 6.2 Perturbation Theory

**Interaction Picture**

Let the Hamiltonian be

$$
\hat H=\hat H_0+\hat V,
\tag{6.1}
$$

with $\hat H_0$ and $\hat V$ time independent in the Schrödinger picture. In the Schrödinger picture, states evolve as $|\psi_S(t)\rangle=e^{-i\hat Ht}|\psi_S(0)\rangle$ and operators $\hat O_S$ are fixed; in the Heisenberg picture, $\hat O_H(t)=e^{i\hat Ht}\hat O_Se^{-i\hat Ht}$ and $|\psi_H\rangle=|\psi_S(0)\rangle$. The interaction picture evolves operators with $\hat H_0$ instead of $\hat H$:

$$
\boxed{
\hat O_I(t):=e^{i\hat H_0t}\,\hat O_S\,e^{-i\hat H_0t},
\qquad
|\psi_I(t)\rangle:=e^{i\hat H_0t}\,|\psi_S(t)\rangle .
}
\tag{6.2}
$$

The three pictures coincide at $t=0$, and $\langle\psi_I(t)|\hat O_I(t)|\psi_I(t)\rangle=\langle\psi_S(t)|\hat O_S|\psi_S(t)\rangle$.

Differentiating the first of (6.2),

$$
\frac{d\hat O_I}{dt}=i\,[\hat H_0,\hat O_I].
\tag{6.3}
$$

In particular $\hat H_{0,I}=\hat H_0$, while the interaction term becomes time dependent. Differentiating the second of (6.2) and using $i\partial_t|\psi_S\rangle=\hat H|\psi_S\rangle$,

$$
i\partial_t|\psi_I(t)\rangle
=e^{i\hat H_0t}\bigl(\hat H-\hat H_0\bigr)|\psi_S(t)\rangle
=e^{i\hat H_0t}\hat Ve^{-i\hat H_0t}\,|\psi_I(t)\rangle,
$$

i.e.

$$
\boxed{
i\partial_t|\psi_I(t)\rangle=\hat V_I(t)\,|\psi_I(t)\rangle,
\qquad
\hat V_I(t):=e^{i\hat H_0t}\,\hat V\,e^{-i\hat H_0t}.
}
\tag{6.4}
$$

Substituting $\hat O_S=e^{-i\hat H_0t}\hat O_I(t)e^{i\hat H_0t}$ into the definition of $\hat O_H$,

$$
\hat O_H(t)=\hat U_I(t,0)^\dagger\,\hat O_I(t)\,\hat U_I(t,0),
\qquad
\hat U_I(t,0):=e^{i\hat H_0t}e^{-i\hat Ht}.
\tag{6.5}
$$

For a field operator, (6.5) reads $\hat\phi_H(t,\mathbf x)=\hat U_I(t,0)^\dagger\,\hat\phi_I(t,\mathbf x)\,\hat U_I(t,0)$, and by (6.3) $\hat\phi_I(t,\mathbf x)$ obeys the Heisenberg equation of $\hat H_0$ alone.

<br><br>

**Dyson Series**

Define the interaction-picture evolution operator $\hat U_I(t,t_0)$ by $|\psi_I(t)\rangle=\hat U_I(t,t_0)|\psi_I(t_0)\rangle$. By (6.4),

$$
i\partial_t\hat U_I(t,t_0)=\hat V_I(t)\,\hat U_I(t,t_0),
\qquad
\hat U_I(t_0,t_0)=1 .
\tag{6.6}
$$

By (6.2) and $|\psi_S(t)\rangle=e^{-i\hat H(t-t_0)}|\psi_S(t_0)\rangle$,

$$
\hat U_I(t,t_0)=e^{i\hat H_0t}\,e^{-i\hat H(t-t_0)}\,e^{-i\hat H_0t_0},
\tag{6.7}
$$

which is unitary, satisfies $\hat U_I(t_1,t_2)\hat U_I(t_2,t_3)=\hat U_I(t_1,t_3)$ and $\hat U_I(t,t_0)^\dagger=\hat U_I(t_0,t)$, and reduces to the $\hat U_I(t,0)$ of (6.5) at $t_0=0$. Since $[\hat V_I(t_1),\hat V_I(t_2)]\neq0$ in general, (6.6) is not solved by $e^{-i\int\hat V_I}$. Instead, integrate (6.6) and iterate:

$$
\begin{aligned}
\hat U_I(t,t_0)
&=1-i\int_{t_0}^{t}dt_1\,\hat V_I(t_1)\,\hat U_I(t_1,t_0)\
&=1+(-i)\int_{t_0}^{t}dt_1\,\hat V_I(t_1)
+(-i)^2\int_{t_0}^{t}dt_1\int_{t_0}^{t_1}dt_2\,\hat V_I(t_1)\hat V_I(t_2)+\cdots ,
\end{aligned}
\tag{6.8}
$$

where the $n$-th term is integrated over $t\ge t_1\ge\cdots\ge t_n\ge t_0$, with later times to the left.

Define the time-ordered product

$$
\mathcal T[\hat V_I(t_1)\hat V_I(t_2)]
:=\theta(t_1-t_2)\,\hat V_I(t_1)\hat V_I(t_2)+\theta(t_2-t_1)\,\hat V_I(t_2)\hat V_I(t_1),
\tag{6.9}
$$

and likewise for $n$ operators, ordered with time decreasing from left to right. $\mathcal T[\hat V_I(t_1)\cdots\hat V_I(t_n)]$ is symmetric under permutations of $t_1,\ldots,t_n$; the cube $[t_0,t]^n$ splits into $n!$ regions according to the ordering of the $t_i$, each contributing equally, and on $t_1>\cdots>t_n$ the time-ordered product is the ordered product in (6.8). Hence

$$
\int_{t_0}^{t}dt_1\int_{t_0}^{t_1}dt_2\cdots\int_{t_0}^{t_{n-1}}dt_n\,
\hat V_I(t_1)\cdots\hat V_I(t_n)
=\frac1{n!}\int_{t_0}^{t}dt_1\cdots\int_{t_0}^{t}dt_n\,
\mathcal T[\hat V_I(t_1)\cdots\hat V_I(t_n)],
\tag{6.10}
$$

and (6.8) becomes the Dyson series

$$
\boxed{
\hat U_I(t,t_0)
=\sum_{n=0}^{\infty}\frac{(-i)^n}{n!}\int_{t_0}^{t}dt_1\cdots dt_n\,
\mathcal T[\hat V_I(t_1)\cdots\hat V_I(t_n)]
=:\mathcal T\exp\Bigl(-i\int_{t_0}^{t}dt'\,\hat V_I(t')\Bigr).
}
\tag{6.11}
$$

Alternatively, differentiate the $n$-th term of (6.11) with respect to $t$: each of the $n$ integration variables may take the upper limit $t$, which is then the latest time, so $\mathcal T$ places $\hat V_I(t)$ leftmost; the $n$ contributions are equal and give $-i\hat V_I(t)$ times the $(n-1)$-th term. Thus the series satisfies (6.6), whose solution is unique.

In field theory, $\hat V=\int d^3x\,\hat{\mathcal H}_{\mathrm{int}}\bigl(\hat\phi_S(\mathbf x)\bigr)$. By (6.2), $\hat V_I(t)=\int d^3x\,\hat{\mathcal H}_{\mathrm{int}}\bigl(\hat\phi_I(t,\mathbf x)\bigr)$, and (6.11) becomes

$$
\hat U_I(t,t_0)=\mathcal T\exp\Bigl(-i\int_{t_0}^{t}dt'\int d^3x\,\hat{\mathcal H}_{\mathrm{int}}\bigl(\hat\phi_I(x)\bigr)\Bigr).
\tag{6.12}
$$

For interactions without derivative couplings, $\hat{\mathcal H}_{\mathrm{int}}=-\hat{\mathcal L}_{\mathrm{int}}$.
