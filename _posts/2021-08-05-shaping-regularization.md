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

$$ \mathbf{A} \mathbf{x} = \mathbf{b}, $$

for which the solution is $\mathbf{x} = \mathbf{A}^{-1} \mathbf{b}$. However, that assumes that $\mathbf{A}$ is square and invertible, which isn't always the case.
In case that requirement does not hold, we can solve a closely linked problem instead, which is the *generalized least squares* solution:

#### Least-squares solution
$$ \mathbf{A}^T\mathbf{A} \mathbf{x} = \mathbf{A}^T\mathbf{b} \Longrightarrow \hat{\mathbf{x}} = \left( \mathbf{A}^T\mathbf{A} \right)^{-1} \mathbf{A}^T\mathbf{b}. $$

The solution $\hat{\mathbf{x}}$ does not satisfy the first set of equations, but it does satisfy the second set of equations.
The product $\mathbf{A}^T\mathbf{A}$ is hopefully invertible.

#### Regularized (Tikhonov least-squares solution)
However, even the least-squres solution is not always nicely invertible, or may otherwise be too sensitive to small errors in the inputs. To counteract this problem, one can use regularization.
The problem that one solves then is

$$ \left( \lambda^2\mathbf{I} + \mathbf{A}^T\mathbf{A}\right) \mathbf{x} = \mathbf{A}^T\mathbf{b} \Longrightarrow \hat{\mathbf{x}} = \left( \lambda^2\mathbf{I} + \mathbf{A}^T\mathbf{A}\right)^{-1} \mathbf{A}^T\mathbf{b}. $$

The addition of values along the diagonal ($\lambda^2$ is a scalar) has the effect of dampening the final solution, and making the inverse a lot more stable.

#### Shaping regularization solution
In the linked paper, eq. (12) corresponds to wanting to solve the following problem, and corresponding solution:

$$ \left( \lambda^2\mathbf{S}^{-1} + \mathbf{A}^T\mathbf{A} - \lambda^2\mathbf{I}\right) \mathbf{x} = \mathbf{A}^T\mathbf{b} \Longrightarrow \hat{\mathbf{x}} = \left( \lambda^2\mathbf{S}^{-1} + \mathbf{A}^T\mathbf{A} - \lambda^2\mathbf{I}\right)^{-1} \mathbf{A}^T\mathbf{b}. $$

In the paper, they consider the case where the equation is pre-multiplied with $\mathbf{S}$ on both sides, which yields the same solution! This formulation shows that the regularization may vary along the entire matrix, as $\mathbf{S}$ or $\mathbf{S}^{-1}$ is not forced to be diagonal or unitary. Note that for $\mathbf{S}^{-1}=2\mathbf{I}$ we recover the Tikhonov solution.

Then, the assumption is made that $\mathbf{S}=\mathbf{S}^{T} = \mathbf{H}\mathbf{H}^T$, i.e., that $\mathbf{S}$ is symmetric. That means that we can rewrite the system as

$$ \left( \lambda^2\mathbf{H}^{-T}\mathbf{H}^{-1} + \mathbf{A}^T\mathbf{A} - \lambda^2\mathbf{I}\right) \mathbf{x} = \mathbf{A}^T\mathbf{b},  $$

and pre-multiply with $\mathbf{H}^T$ to find

$$ \left( \lambda^2\mathbf{H}^{-1} + \mathbf{H}^T\left( \mathbf{A}^T\mathbf{A} - \lambda^2\mathbf{I}\right)\right) \mathbf{x} = \mathbf{H}^T\mathbf{A}^T\mathbf{b},  $$

and we may re-write the terms in the brackets as

$$ \left( \lambda^2\mathbf{I} + \mathbf{H}^T\left(\mathbf{A}^T\mathbf{A}- \lambda^2\mathbf{I}\right)\mathbf{H}\right) \mathbf{H}^{-1} \mathbf{x} = \mathbf{H}^T\mathbf{A}^T\mathbf{b}.  $$

That gives the final equation that is inverted to find a solution, where we use $(M H^{-1})^{-1}=H M^{-1}$,

$$ \hat{\mathbf{x}} = \mathbf{H}\left( \lambda^2\mathbf{I} + \mathbf{H}^T\left(\mathbf{A}^T\mathbf{A}- \lambda^2\mathbf{I}\right)\mathbf{H}\right)^{-1}\mathbf{H}^T\mathbf{A}^T\mathbf{b}. $$

So, that shows how the equation in the paper may be derived. Note, again, that the Tikhonov solution is recovered for $\mathbf{H}=\mathbf{I}/\sqrt{2}$.

#### Implementation
The paper, additionally, provides a conjugate gradient algorithm to compute the shaping regularized solution. Unfortunately, the algorithm is erroneous. The good news is that the fix is simple: every occurrence of $\lambda$ should be replaced with $\lambda^2$, to be in line with the mathematical notation in the paper. A Python implementation is

````python
import numpy
 
def norm(a):
    return a.T @ a
 
def conj_grad_shaping(L,H,d,lambd,tol=1e-20,N=500):
    p = np.zeros(L.shape[1])
    m = np.zeros(L.shape[1])
    r = -d
    for n in range(N):
        gm = L.T @ r  - lambd**2 * m
        gp = H.T @ gm + lambd**2 * p
        gm = H @ gp
        gr = L @ gm
        rho = norm(gp)
        if n==0:
            beta = 0
            rho0 = rho
            sp = gp
            sm = gm
            sr = gr
        else:
            beta = rho/rhocur
            if beta < tol or rho/rho0 < tol:
                return m
            sp = gp + beta*sp
            sm = gm + beta*sm
            sr = gr + beta*sr
        alpha = rho/( norm(sr) + lambd**2*( norm(sp)  - norm(sm) ))
        p = p - alpha*sp
        m = m - alpha*sm
        r = r - alpha*sr
        rhocur = rho
    return m
 
L = np.array([[1, 3 ],
              [2, 4 ],
              [1, 6 ]])
H = np.array([[1  ,0.2],
              [0.2,1  ]])
d = np.array([4, 1, 3])
lambd = 1.9
 
S = H @ H.T
   
m_est_cg = conj_grad_shaping(L,H,d,lambd)
print("conjug grad result=", m_est_cg )
print("error...", norm(  (lambd**2 * np.eye(L.shape[1]) + S@(norm(L) - lambd**2*np.eye(L.shape[1])))@m_est_cg  - S@L.T@d ) )

# The following is only to compare it against a built-in linear algebra solver
from numpy.linalg import inv, solve
m_est_las = solve( lambd**2 * inv(S) + norm(L) - lambd**2 * np.eye(L.shape[1]), L.T @ d )
print("full solver result=", m_est_las )
print("error...", norm(  (lambd**2 * np.eye(L.shape[1]) + S@(norm(L) - lambd**2*np.eye(L.shape[1])))@m_est_las - S@L.T@d ) )
```

