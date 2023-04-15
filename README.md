# Discontinuous Galerkin Linear Wave Test
```bash
make
./a.out
```
```bash
python -m http.server
```
Open http://localhost:8000/index.html

## Method
Solves

$$
\partial_t u + \frac{1}{2}\partial_x u^2 = 0
$$

in the spatial interval [0, 1]. Note the linearized wave speed is $u$ and flux is $u^2/2$.

Expand within each cell centered at $x = a$, $x \in [a - \Delta x/2, a + \Delta x / 2]$

$$
u(x, t) = \sum_n c_n(t) w_n(x)
$$

where the basis functions $w_n$ are Legendre polynomials

$$
w_n(x) = P_n(\frac{2 (x - a)}{\Delta x})
$$

Require the discretized weak form

$$
0 = \int w_m (\partial_t u + \partial_x u^2/2)dx
$$

and integrate the 2nd term by parts, while using orthogonality of $w$'s on the
first term. The coefficient update is obtained as

$$
\frac{\Delta x}{2m + 1} (\partial_t c_m) = -[w_m J_{\text{numerical}}]_-^+ + \frac{1}{2} \sum_{nl} c_n c_l \int w_n w_l (\partial_x w_m)dx
$$

where $J_{\text{numerical}}$ is the numerical flux at cell edges to be found by a Riemann solver.

$J_{\text{numerical}}$ is found by evaluating $u(x, t)$ at cell
edges (where $w_m = \pm 1$) and then using HLLE.

The term $\int w_n w_l (\partial_x w_m)dx$ can be precomputed as a single
array indexed by $m,n,l$.

The method is conservative since only the zeroth Legendre polynomial contributes
to $\int udx$ in a cell, and the update for the zeroth coefficient reduces to

$$
\Delta x (\partial_t c_0) = -[J_\text{numerical}]_-^+
$$
