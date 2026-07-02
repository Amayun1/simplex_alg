simplex_alg
======

This project is a simple C implementation of the simplex algorithm for linear programming. The simplex algorithm 
solves linear programs in the standard form

```math
\begin{align*}
	\text{minimize } &z = c^Tx \\
	\text{subject to} & Ax = b
\end{align*}
```

where $c \in \mathbb{R}^n$, $A$ is an $m \times n$ matrix, and $b \in \mathbb{R}^m$ with $b \geq 0$.
The project may also solve linear programs with constraints of the form $Ax \leq b$ using a wrapper which converts the problem to standard form.

The main algorithm functions are in [`simplex.c`](simplex.c). Files [`simplex_test.c`](simplex_test.c), [`ez_simplex_test.c`](ez_simplex_test.c), and [`klee_minty.c`](klee_minty.c) contain tests programs for the algorithm.

For an overview on the simplex algorithm see <https://en.wikipedia.org/wiki/Simplex_algorithm>.
For further details on the simplex algorithm and general optimization, see [*Linear and Nonlinear Optimization*](https://epubs.siam.org/doi/book/10.1137/1.9780898717730) by Griva, Nash, and Sofer or another suitable optimization reference.
