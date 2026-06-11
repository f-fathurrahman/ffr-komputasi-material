Here’s a detailed explanation of the mathematics behind the Julia program, covering the Hamiltonian construction, discretization, and antisymmetrization procedure.

---

## 1. The Physical System

We consider $N_{\text{el}}$ spin‑1/2 fermions moving in one dimension. The Hamiltonian (in atomic units, $\hbar = m = 1$) is

$$
\hat H = \sum_{i=1}^{N_{\text{el}}} \left[ -\frac12 \frac{d^2}{dx_i^2} + V_{\text{ext}}(x_i) \right] + \sum_{1\le i<j\le N_{\text{el}}} v_{\text{ee}}(x_i-x_j),
$$

where  
- $V_{\text{ext}}(x) = \frac12 \omega^2 x^2$ (harmonic trap),  
- $v_{\text{ee}}(r) = 1/\sqrt{r^2 + a_{\text{soft}}^2}$ (soft‑Coulomb interaction).

The total wavefunction $\Psi(x_1,\sigma_1,\dots,x_{N_{\text{el}}},\sigma_{N_{\text{el}}})$ must be antisymmetric under exchange of any two particles (including their spins).

---

## 2. Discretization on a Spatial Grid

We restrict the spatial domain to $x\in[-L,L]$ and discretize it into $N$ equally spaced points:

$$
x_j = -L + (j-1)\Delta x,\qquad \Delta x = \frac{2L}{N-1},\qquad j=1,\dots,N.
$$

The single‑particle Hamiltonian (kinetic + external) is represented as an $N\times N$ matrix in this grid basis. The kinetic term $-\frac12\frac{d^2}{dx^2}$ is approximated by a finite‑difference stencil of order $m$:

$$
\left(-\frac12\frac{d^2}{dx^2}\psi\right)(x_j) \approx -\frac12\sum_{k=-s}^{s} c_k \frac{\psi(x_{j+k})}{\Delta x^2},
$$

where $s = (m-1)/2$ and the coefficients $c_k$ are standard finite‑difference weights for the second derivative (e.g., for $m=3$: $c_{-1}=1, c_0=-2, c_1=1$). This yields a sparse symmetric matrix $T$ with elements

$$
T_{j,j+k} = -\frac12 \frac{c_k}{\Delta x^2}.
$$

The external potential is diagonal: $V_{j,j} = V_{\text{ext}}(x_j)$.

Thus the single‑particle matrix is

$$
\mathbf{h} = \mathbf{T} + \mathbf{V}_{\text{ext}}.
$$

---

## 3. Many‑Body Hamiltonian in the Product Basis

For distinguishable particles, we take the basis states as simple products of single‑particle grid states:

$$
|i_1,i_2,\dots,i_{N_{\text{el}}}\rangle = |i_1\rangle \otimes |i_2\rangle \otimes \cdots \otimes |i_{N_{\text{el}}}\rangle,
$$

where each $i_p$ runs from $1$ to $N$. The corresponding wavefunction is $\phi_{i_1,\dots,i_{N_{\text{el}}}}(x_1,\dots,x_{N_{\text{el}}}) = \prod_{p=1}^{N_{\text{el}}} \delta_{x_p, x_{i_p}}$. This basis has dimension $N^{N_{\text{el}}}$.

### 3.1 Non‑interacting part

The operator $\sum_{p=1}^{N_{\text{el}}} \left( -\frac12 \frac{d^2}{dx_p^2} + V_{\text{ext}}(x_p) \right)$ acts on each particle independently. In the product basis, it becomes a sum of Kronecker products:

$$
\mathbf{H}_0 = \sum_{p=1}^{N_{\text{el}}} \underbrace{\mathbf{I} \otimes \cdots \otimes \mathbf{h}}_{\text{position }p} \otimes \cdots \otimes \mathbf{I}.
$$

If we order the basis such that the first particle index varies fastest (column‑major order), the linear index corresponding to $(i_1,i_2,\dots,i_{N_{\text{el}}})$ is

$$
\text{idx} = 1 + \sum_{p=1}^{N_{\text{el}}} (i_p-1) N^{p-1}.
$$

Then the term for particle $p$ can be constructed by taking the Kronecker product in reverse order: start with the factor for particle 1, then multiply from the left by the factor for particle 2, etc. The Julia function `build_term(p)` implements this.

### 3.2 Interaction part

The two‑body interaction $\sum_{i<j} v_{\text{ee}}(x_i-x_j)$ is diagonal in the product basis because the basis states are localized at grid points. Its matrix elements are

$$
\langle i_1,\dots,i_{N_{\text{el}}} | \hat V_{\text{int}} | j_1,\dots,j_{N_{\text{el}}}\rangle = \delta_{i_1,j_1}\cdots\delta_{i_{N_{\text{el}}},j_{N_{\text{el}}}} \sum_{p<q} v_{\text{ee}}(x_{i_p}-x_{i_q}).
$$

Thus the interaction contributes a diagonal matrix $\mathbf{U}$ with entries

$$
U_{\text{idx}} = \sum_{1\le p<q\le N_{\text{el}}} v_{\text{ee}}(x_{i_p}-x_{i_q}).
$$

The total many‑body Hamiltonian in the product basis is

$$
\mathbf{H} = \mathbf{H}_0 + \mathbf{U}.
$$

---

## 4. Exact Diagonalization

We use ARPACK to compute the lowest few eigenvalues and eigenvectors of the sparse matrix $\mathbf{H}$ (real, symmetric). For a given number of computed states $M$, we obtain eigenvalues $E^{(n)}$ and corresponding eigenvectors $\mathbf{c}^{(n)}$ in the product basis. Each eigenvector can be reshaped into a tensor $\Psi^{(n)}(i_1,\dots,i_{N_{\text{el}}})$ of size $N\times\cdots\times N$.

---

## 5. Antisymmetrization

The eigenvectors from the product basis correspond to **distinguishable** particles. To obtain physical fermionic states, we must antisymmetrize the wavefunction including spin.

### 5.1 Spin part

For a given spin configuration specified by the string `electrons` (e.g., `"uud"`), we construct the spin tensor as a product of spinors:

$$
\chi(\sigma_1,\dots,\sigma_{N_{\text{el}}}) = \bigotimes_{p=1}^{N_{\text{el}}} \begin{cases} \begin{pmatrix}1\\0\end{pmatrix} & \text{if spin up},\\ \begin{pmatrix}0\\1\end{pmatrix} & \text{if spin down}. \end{cases}
$$

This yields a tensor of shape $(2,2,\dots,2)$ (dimension $2^{N_{\text{el}}}$) with a single non‑zero entry at the index corresponding to the given spin configuration.

### 5.2 Full wavefunction (spatial + spin)

We combine spatial and spin degrees of freedom into a single tensor with alternating dimensions:

$$
\Phi(x_1,\sigma_1,x_2,\sigma_2,\dots,x_{N_{\text{el}}},\sigma_{N_{\text{el}}}) = \Psi(x_1,\dots,x_{N_{\text{el}}}) \cdot \chi(\sigma_1,\dots,\sigma_{N_{\text{el}}}).
$$

In the discretized representation, this is an array of size $(N,2,N,2,\dots,N,2)$. The element at indices $(i_1,s_1,i_2,s_2,\dots,i_{N_{\text{el}}},s_{N_{\text{el}}})$ equals $\Psi(i_1,\dots,i_{N_{\text{el}}}) \cdot \chi(s_1,\dots,s_{N_{\text{el}}})$.

### 5.3 Permutation sum

To enforce antisymmetry, we sum over all permutations $P$ of the $N_{\text{el}}$ particles, each weighted by the sign of the permutation:

$$
\Phi_{\text{anti}} = \sum_{P} \operatorname{sgn}(P) \; \Phi^{(P)},
$$

where $\Phi^{(P)}$ is the wavefunction with the coordinates of particle $p$ replaced by those of particle $P(p)$. In tensor terms, we permute the particle dimensions accordingly. For example, for $N_{\text{el}}=2$,

$$
\Phi_{\text{anti}}(x_1,\sigma_1,x_2,\sigma_2) = \Phi(x_1,\sigma_1,x_2,\sigma_2) - \Phi(x_2,\sigma_2,x_1,\sigma_1).
$$

The sign $\operatorname{sgn}(P)$ is $+1$ for even permutations, $-1$ for odd ones. In the implementation, we loop over all permutations, compute the required dimension permutation using `permutedims`, and accumulate the result.

### 5.4 Normalization

After summation, we normalize:

$$
\Phi_{\text{anti}} \leftarrow \frac{\Phi_{\text{anti}}}{\sqrt{\sum_{\text{all indices}} |\Phi_{\text{anti}}|^2 \; (\Delta x)^{N_{\text{el}}}}},
$$

where the factor $(\Delta x)^{N_{\text{el}}}$ accounts for the volume element in the continuous inner product.

---

## 6. Final States

The antisymmetrized states $\Phi_{\text{anti}}$ now represent proper fermionic wavefunctions. Duplicate states (identical up to a global phase) are removed. The ground state (or any desired state) is extracted, and its energy is the corresponding eigenvalue from the diagonalization.

---

## 7. Remarks on the Implementation

- **Memory and performance**: The dimension $N^{N_{\text{el}}}$ grows quickly. For $N=60$ and $N_{\text{el}}=2$ it is 3600; for $N=40$ and $N_{\text{el}}=3$ it is 64,000. Sparse matrices and ARPACK make such sizes feasible.
- **Basis ordering**: The code carefully uses column‑major indexing (particle 1 varies fastest) to match the `reshape` and `LinearIndices` conventions.
- **Antisymmetrization loop**: The permutation loop uses `permutedims` with a `Tuple` of the new dimension order. This is efficient for small $N_{\text{el}}$ (the number of permutations is $N_{\text{el}}!$, which is manageable for $N_{\text{el}}\le 4$).
- **Spin**: The spin part is constructed as a tensor product of spinors. For configurations with more than one electron of the same spin, the resulting spin tensor may have multiple non‑zero entries, but the antisymmetrization automatically mixes them to yield the correct spin multiplets.

This approach directly implements the many‑body Schrödinger equation on a real‑space grid, without approximating the interaction by a basis of non‑interacting orbitals. It is conceptually straightforward and yields exact (up to discretization) results for the given grid.