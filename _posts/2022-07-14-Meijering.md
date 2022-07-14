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
To start with, define the image as $f(\mathbf{x})$ with $\mathbf{x}\in\mathbb{R}^n$ to define an $n$-D image, and define a derivative operation on images [as a convolution with a (similarly) differentiated Gaussian](https://www.crisluengo.net/archives/22/),

$$f_j(\mathbf{x}) \equiv f(\mathbf{x}) * \frac{\partial}{\partial x_j}\frac{e^{-\mathbf{x}^2/(2\sigma^2)}}{(\sqrt{2\pi}\sigma)^n}.$$

and

$$f_{ij}(\mathbf{x}) = f(\mathbf{x}) * \frac{\partial}{\partial x_i}\frac{\partial}{\partial x_j}\frac{e^{-\mathbf{x}^2/(2\sigma^2)}}{(\sqrt{2\pi}\sigma)^n}.$$


Furthermore, we define the Hessian of $f$ as

$$H_f(\mathbf{x}) = \left[ \begin{array}{ccc}
f_{xx}(\mathbf{x}) & f_{xy}(\mathbf{x}) \\
f_{xy}(\mathbf{x}) & f_{yy}(\mathbf{x}) \end{array} \right].$$

Because the Hessian is [symmetric and real](https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix#Real_symmetric_matrices), it allows an eigendecomposition with $n$ orthonormal eigenvectors $\mathbf{v}\_i$ (that is, $\mathbf{v}\_i^T\mathbf{v}\_j=0$ for $i\neq j$ and each eigenvectors has length 1), with corresponding eigenvalues $\lambda_i$ (which satisfy $H_f \mathbf{v}\_i = \lambda_i \mathbf{v}\_i$). This may be written as

$$
\mathbf{v}\_i^T H_f(\mathbf{x}) \mathbf{v}\_i = \lambda_i.
$$

Interestingly, it also holds that pre- and post-multiplying a Hessian matrix with a vector $d\in\mathbb{R}^n$ of length 1 gives us the second derivative in the direction $d$ (see, e.g., [here](https://math.stackexchange.com/questions/2573376/second-directional-derivative-and-hessian-matrix)),

$$
\mathbf{d}^T H_f(\mathbf{x}) \mathbf{d} = f_{\mathbf{d}\mathbf{d}}(\mathbf{x}) = (\mathbf{d}\cdot\nabla)^2 f(\mathbf{x}),
$$

where $f_{\mathbf{d}\mathbf{d}}(\mathbf{x})$ is abuse of notation to indicate the second derivative in direction $\mathbf{d}$ at $(\mathbf{x})$. We realize that this pre- and post-multiplication pattern could also be seen in the case of the eigendecomposition, and we thus have that:

$$
\mathbf{v}_i^T H_f(\mathbf{x}) \mathbf{v}_i = f_{\mathbf{v}_i\mathbf{v}\_i}(\mathbf{x}) = (\mathbf{v}_i\cdot\nabla)^2 f(\mathbf{x}) = \lambda_i.
$$

For example, for a 2D image if $\mathbf{v}\_i^T=(v_1\quad v_2)$, then $(\mathbf{v}\_i\cdot\nabla)^2=v_1^2\partial^2/\partial x_1^2 + 2v_1v_2 \partial^2/(\partial x_1\partial x_2) + v_2^2\partial^2/\partial x_2^2$ and $(\mathbf{v}\_i\cdot\nabla)^2 f(\mathbf{x})=v_1^2 f_{xx} + 2v_1v_y f_{xy} + v_2^2f_{yy}$. And, remarkably, that derivative will be exactly equal to the eigenvalue $\lambda_i$ at $(\mathbf{x})$!

#### An altered Hessian
Meijering et al. (2004) define an altered Hessian (denoted by a prime $'$), with a tunable parameter $\alpha$ (omitting all $(\mathbf{x})$ dependencies momentarily)

$$H_f'(\mathbf{x}) = \left[ \begin{array}{ccc}
f_{xx}+\alpha f_{yy} & (1-\alpha)f_{xy} \\
(1-\alpha)f_{xy} & f_{yy} + \alpha f_{xx} \end{array} \right].$$

This Hessian also allows an eigendecomposition with normalized eigenvectors $\mathbf{v}\_i'$ and corresponding eigenvalues $\lambda_i'$. It turns out that, neatly, the eigenvalues of this system are not altered at all ($\mathbf{v}\_i'=\mathbf{v}\_i$), and that the eigenvalues are simply related through a simple relation, $\lambda_1'=\lambda_1+\alpha\lambda_2$ and $\lambda_2'=\lambda_2+\alpha\lambda_1$. Thus, we know that the following eigendecomposition holds,

$$ \mathbf{v}_i^T H_f'(\mathbf{x}) \mathbf{v}_i = \lambda_i' = \lambda_i + \alpha \lambda_j,\quad (j\neq i). $$

From the previous section, we know that $\lambda_i$ is simply the second directional derivative in direction $\mathbf{v}\_i$; and similarly $\alpha\lambda_j$ will be the second directional derivative in direction $\mathbf{v}\_j$ scaled with a parameter $\alpha$. Thus,

$$ \mathbf{v}\_i^T H_f'(\mathbf{x}) \mathbf{v}\_i = \lambda_i' = \lambda_i + \alpha \lambda_j = \left((\mathbf{v}_i\cdot\nabla)^2 + \alpha(\mathbf{v}_j\cdot\nabla)^2\right) f(\mathbf{x}),\quad (j\neq i). $$

Now comes the step in which we define $\alpha$: we want that the filter is as 'straight' as possible in the orthogonal direction as possible (such that our filter operation resembles a box car) -- which corresponds to setting its second derivative in this orthogonal direction to zero. 

$$\lim_{\mathbf{x}\to 0} \left[ (\mathbf{v}_j \cdot \nabla)^2 \left(\mathbf{v}_i^T H_f'(\mathbf{x}) \mathbf{v}_i\right)\right] = 0 \quad (j\neq i).$$

Then it becomes a matter of linear algebra and differentiating Gaussians, to work out what that left-hand side corresponds to. Using the $\lambda_i$ and directional derivative relations found above, we find

$$\lim_{\mathbf{x}\to 0} \left[ (\mathbf{v}_j \cdot \nabla)^2 \left((\mathbf{v}_i\cdot\nabla)^2 + \alpha(\mathbf{v}_j\cdot\nabla)^2\right) f(\mathbf{x})\right] = \lim_{\mathbf{x}\to 0}\left[ (\mathbf{v}_i\cdot\nabla)^2(\mathbf{v}_j\cdot\nabla)^2 +\alpha(\mathbf{v}_j\cdot\nabla)^4\right]f(\mathbf{x}) \quad (j\neq i).$$

Before we work out all the terms in brackets, we establish that every term will correspond to some fourth derivative $f_{ijkl}(\mathbf{x})$, so it is helpful to see what limits we will find. Note that the derivative of a Gaussian is typically also a Gaussian, multiplied with some terms, and the Gaussian evaluated at $\mathbf{x}=0$ equals $G(0)=1$. We know that the order of differentiation doesn't matter, thus $f_{xxxy}=f_{xxyx}=f_{xyxx}=f_{yxxx}$ etc. We thus need to establish these five possible variations (use, e.g., [WolframAlpha](https://www.wolframalpha.com/input?i=Lim%5B+D%5BD%5Be%5E%28-%28x%5E2%2By%5E2%29%2F%282*s%5E2%29%29%2C+%7Bx%2C2%7D%5D%2C+%7By%2C2%7D%5D+%2C+%7Bx-%3E0%2C+y-%3E0%7D%5D))

$$
\lim_{\mathbf{x}\to 0}  \left[ \begin{array}{ccc}
f_{xxxx} \\
f_{xxxy} \\
f_{xxyy} \\
f_{xyyy} \\
f_{yyyy}\end{array} \right] = \left[ \begin{array}{ccc}
3/\sigma^2 \\
0 \\
1/\sigma^2 \\
0 \\
3/\sigma^2\end{array} \right] .
$$
