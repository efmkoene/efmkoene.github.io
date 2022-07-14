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

To start with, define the image as $f(\mathbf{x})$ with $\mathbf{x}\in\mathbb{R}^n$ to define an $n$-D image, and define any derivative operations [as a convolution with a (similarly) differentiated Gaussian](https://www.crisluengo.net/archives/22/),
\begin{equation}
  f_j(\mathbf{x}) = f(\mathbf{x}) * \frac{\frac{\partial}{\partial x_j}e^{-\mathbf{x}^2/(2\sigma^2)}}{(\sqrt{2\pi}\sigma)^n}.
\end{equation}
and
\begin{equation}
  f_{ij}(\mathbf{x}) = f(\mathbf{x}) * \frac{\frac{\partial}{\partial x_i}\frac{\partial}{\partial x_j}e^{-\mathbf{x}^2/(2\sigma^2)}}{(\sqrt{2\pi}\sigma)^n}.
\end{equation}

Furthermore, we define the Hessian of $f$ as
\begin{equation}
H_f = \left[ \begin{array}{cc} f_{xx} & f_{xy} \\ f_{xy} & f_{yy} \end{array} \right].
\end{equation}
Pre- and post-multiplying this Hessian matrix with a vector $d\in\mathbb{R}^n$ of length 1 gives us the second derivative in the direction $d$ (see, e.g., [here](https://math.stackexchange.com/questions/2573376/second-directional-derivative-and-hessian-matrix)),
\begin{equation}
\mathbf{d}^T H_f \mathbf{d} = f_{\mathbf{d}\mathbf{d}},
\end{equation}
where $f_{\mathbf{d}\mathbf{d}}$ is abuse of notation to indicate the second derivative in direction $\mathbf{d}$. If, for example, $\mathbf{d}^T=[1\quad 0]$, we obtain $\mathbf{d}^TH_f\mathbf{d}=f_xx$, so the 2nd derivative in the $x$ (or '1') direction.
A special occasian presents itself when the Hessian $H_f$ has an eigenvalue decomposition for $n$ eigenvectors $\mathbf{v}\_i$ and corresponding eigenvalues $\lambda_i$ (which satisfy $H_f \mathbf{v}_i = \lambda_i \mathbf{v}_i$), as in that case we may rewrite the above equation as
\begin{equation}
  \mathbf{v}_i^T H_f \mathbf{v} = \lambda_i = f_{\mathbf{v}_i\mathbf{v}_i},
\end{equation}
where we [made use of the fact that](https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix#Real_symmetric_matrices) the Hessian $H_f$ is symmetric and real.

