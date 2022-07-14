---
layout: post
title: "The Meijering method in 3D"
subtitle: ""
tags: [python]
---

I recently contributed to my [first Github pull request](https://github.com/scikit-image/scikit-image/pull/6149) for the `scikit-image` project.

The contribution deals with an implementation of a method proposed by Meijering et al. (2004), https://onlinelibrary.wiley.com/doi/10.1002/cyto.a.20022 which was not quite perfect in `scikit-image`.

One of the issues was that the tuneable factor $\alpha$ was not chosen perfectly in the previous implementation, because Meijering et al. (2004) show that $\alpha=-1/3$ is the optimal value.

I'd like to clarify where this value comes from.

#### Preliminaries
To start with, define the image as $f(\mathbf{x})$ with $\mathbf{x}\in\mathbb{R}^n$ to define an $n$-D image, and define any derivative operations [as a convolution with a (similarly) differentiated Gaussian](https://www.crisluengo.net/archives/22/),

$$f_j(\mathbf{x}) = f(\mathbf{x}) * \frac{\partial}{\partial x_j}\frac{e^{-\mathbf{x}^2/(2\sigma^2)}}{(\sqrt{2\pi}\sigma)^n}.$$

and

$$f_{ij}(\mathbf{x}) = f(\mathbf{x}) * \frac{\partial}{\partial x_i}\frac{\partial}{\partial x_j}\frac{e^{-\mathbf{x}^2/(2\sigma^2)}}{(\sqrt{2\pi}\sigma)^n}.$$


Furthermore, we define the Hessian of $f$ as

$$H_f = \left[ \begin{array}{ccc}
f_{xx} & f_{xy} \\
f_{xy} & f_{yy} \end{array} \right].$$

Pre- and post-multiplying this Hessian matrix with a vector $d\in\mathbb{R}^n$ of length 1 gives us the second derivative in the direction $d$ (see, e.g., [here](https://math.stackexchange.com/questions/2573376/second-directional-derivative-and-hessian-matrix)),

$$
\mathbf{d}^T H_f \mathbf{d} = f_{\mathbf{d}\mathbf{d}} = (\mathbf{d}\cdot\nabla)^2 f(\mathbf{x}),
$$

where $f_{\mathbf{d}\mathbf{d}}$ is abuse of notation to indicate the second derivative in direction $\mathbf{d}$. If, for example, $\mathbf{d}^T=[1\quad 0]$, we obtain $\mathbf{d}^TH_f\mathbf{d}=f_{xx}$, so the 2nd derivative in the $x$ (or '1') direction.
A special occasian presents itself when the Hessian $H_f$ has an eigenvalue decomposition for $n$ orthonormal eigenvectors $\mathbf{v}\_i$ (that is, $\mathbf{v}\_i^T\mathbf{v}\_j=0$ for $i\neq j$ and each eigenvectors has length 1), with corresponding eigenvalues $\lambda_i$ (which satisfy $H_f \mathbf{v}\_i = \lambda_i \mathbf{v}\_i$), as in that case we may rewrite the above equation as

$$
\mathbf{v}\_i^T H_f \mathbf{v}\_i = \lambda_i = f_{\mathbf{v}\_i\mathbf{v}\_i} = (\mathbf{v}_i\cdot\nabla)^2 f(\mathbf{x}),
$$

where we [made use of the fact that](https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix#Real_symmetric_matrices) the Hessian $H_f$ is symmetric and real (and we assume that the eigendecomposition exists). For example, for a 2D image if $\mathbf{v}\_i^T=(v_1\quad v_2)$, then $(\mathbf{v}\_i\cdot\nabla)^2=v_1^2\partial^2/\partial x_1^2 + 2v_1v_2 \partial^2/(\partial x_1\partial x_2) + v_2^2\partial^2/\partial x_2^2$.

#### An altered Hessian

