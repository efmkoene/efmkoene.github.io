---
layout: post
title: "Shaping regularization"
subtitle: ""
tags: [geophysics, Python]
---

I'm going through the paper [Shaping regularization in geophysical estimation problems](https://library.seg.org/doi/abs/10.1190/1.2433716), and in particular the conjugate gradient program to estimate the shaping-regularized solutions.
I struggled a little bit with the paper, and found a slight mistake in the explanation of the algorithm. So, perhaps the following notes will help another person out, at some point.


#### Inverting a system of equations
The setting is the solution to linear (algebra) problems. In their simplest form these are problems of the form

$$ \mathbf{A} \vec{x} = \vec{b}, $$

for which the solution is $\vec{x} = \mathbf{A}^{-1} \vec{b}$. However, that assumes that $\mathbf{A}$ is square and invertible, which isn't always the case.
In case that requirement does not hold, we can solve a closely linked problem instead, which is the *generalized least squares* solution:

#### Least-squares solution
$$ \mathbf{A}^T\mathbf{A} \vec{x} = \mathbf{A}^T\vec{b} \Longrightarrow \hat{\vec{x}} = \left( \mathbf{A}^T\mathbf{A} \right)^{-1} \mathbf{A}^T\vec{b}. $$

The solution $\hat{\vec{x}}$ does not satisfy the first set of equations, but it does satisfy the second set of equations.
The product $\mathbf{A}^T\mathbf{A}$ is hopefully invertible.

#### Regularized (Tikhonov least-squares solution)
However, even the least-squres solution is not always nicely invertible, or may otherwise be too sensitive to small errors in the inputs. To counteract this problem, one can use regularization.
The problem that one solves then is

$$ \left( \lambda^2\mathbf{I} + \mathbf{A}^T\mathbf{A}\right) \vec{x} = \mathbf{A}^T\vec{b} \Longrightarrow \hat{\vec{x}} = \left( \lambda^2\mathbf{I} + \mathbf{A}^T\mathbf{A}\right)^{-1} \mathbf{A}^T\vec{b}. $$

The addition of values along the diagonal ($\lambda^2$ is a scalar) has the effect of dampening the final solution, and making the inverse a lot more stable.

#### Shaping regularization solution
In the linked paper, eq. (12) corresponds to wanting to solve the following problem, and corresponding solution:

$$ \left( \lambda^2\mathbf{S}^{-1} + \mathbf{A}^T\mathbf{A}\right) \vec{x} = \mathbf{A}^T\vec{b} \Longrightarrow \hat{\vec{x}} = \left( \lambda^2\mathbf{S}^{-1} + \mathbf{A}^T\mathbf{A}\right)^{-1} \mathbf{A}^T\vec{b}. $$

In the paper, they consider the case where the equation is pre-multiplied with $\mathbf{S}$ on both sides, which yields the same solution! This formulation shows that the regularization may vary along the entire matrix, as $\mathbf{S}$ or $\mathbf{S}^{-1}$ is not forced to be diagonal or unitary.

