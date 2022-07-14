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
\mathbf{v}_i^T H_f(\mathbf{x}) \mathbf{v}_i = \lambda_i.
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

This Hessian also allows an eigendecomposition with normalized eigenvectors $\mathbf{v}\_i'$ and corresponding eigenvalues $\lambda_i'$. It turns out that, neatly, the eigenvalues of this system are not altered at all ($\mathbf{v}\_i'=\mathbf{v}\_i$), and that the eigenvalues are related through a simple formula, $\lambda_1'=\lambda_1+\alpha\lambda_2$ and $\lambda_2'=\lambda_2+\alpha\lambda_1$. Thus, we know that the following eigendecomposition holds,

$$ \mathbf{v}_i^T H_f'(\mathbf{x}) \mathbf{v}_i = \lambda_i' = \lambda_i + \alpha \lambda_j,\quad (j\neq i). $$

From the previous section, we know that $\lambda_i$ is simply the second directional derivative in direction $\mathbf{v}\_i$; and similarly $\alpha\lambda_j$ will be the second directional derivative in direction $\mathbf{v}\_j$ scaled with a parameter $\alpha$. Thus,

$$ \mathbf{v}_i^T H_f'(\mathbf{x}) \mathbf{v}_i = \lambda_i' = \lambda_i + \alpha \lambda_j = \left((\mathbf{v}_i\cdot\nabla)^2 + \alpha(\mathbf{v}_j\cdot\nabla)^2\right) f(\mathbf{x}),\quad (j\neq i). $$

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
f_{yyyy}\end{array} \right] = f(0)\left[ \begin{array}{ccc}
3/\sigma^2 \\
0 \\
1/\sigma^2 \\
0 \\
3/\sigma^2\end{array} \right] .
$$

Hence, we find that e.g. the first term equals (setting $\mathbf{v}\_i^T=(v_1\quad v_2)$ and $\mathbf{v}\_j=(r_1\quad r_2)$ and simply computing each term and evaluating the derivatives at zero,

$$ \lim_{\mathbf{x}\to 0}\left[ (\mathbf{v}_i\cdot\nabla)^2(\mathbf{v}_j\cdot\nabla)^2 f(\mathbf{x}) \right] = f(0) \frac{3v_1^2 r_1^2 + v_1^2r_2^2 + v_2^2r_1^2 + 3v_2r^2 + 4r_1r_2v_1v_2}{\sigma^4}. $$

We now use two results from linear algebra for orthonormal systems of eigenvectors $\mathbf{v}\_i$:

$$ \mathbf{v}_i \left[ \begin{array}{ccc}
1 & 0 \\
0 & 1 \end{array} \right] \mathbf{v}_j = \begin{cases} 0 \quad (j\neq i) \\ 1 \quad (j=i) \end{cases} $$

and

$$ \mathbf{v}_i \left[ \begin{array}{ccc}
0 & 1 \\
1 & 0 \end{array} \right] \mathbf{v}_j = \pm 1 \quad (j\neq i) $$

In our example that corresponds to the terms $v_1r_1+v_2r_2=0$, or $(v_1r_1+v_2r_2)^2=v_1^2r_1^2+2v_1v_2r_1r_2+v_2^2r_2^2=0$ and $(v_1r_2 - v_2r_1)=\pm 1$, thus $(v_1r_2 - v_2r_1)^2=v_1^2r_2^2 + v_2^2r_1^2 - 2r_1r_2v_1v_2=1$. The sum of these two squared terms corresponds exactly to the terms collected above, thus we find that 

$$ \lim_{\mathbf{x}\to 0}\left[ (\mathbf{v}_i\cdot\nabla)^2(\mathbf{v}_j\cdot\nabla)^2 f(\mathbf{x}) \right] = f(0) \frac{1}{\sigma^4}\quad (j\neq i). $$

Correspondingly, for the other term in the limit, containing the $\alpha$, we obtain for $\mathbf{v}\_j=(r_1\quad r_2)$

$$ \lim_{\mathbf{x}\to 0}\left[ (\alpha(\mathbf{v}_j\cdot\nabla)^4\right]f(\mathbf{x}) = f(0)\frac{3r_1^4 + 6r_1^2r_2^2 + 3r_2^2}{\sigma^2} = f(0)\frac{3}{\sigma^4}, $$

where we used the orthogonality relation as described above. Combining the findings for the two limits, we obtain

$$\lim_{\mathbf{x}\to 0} \left[ (\mathbf{v}_j \cdot \nabla)^2 \left((\mathbf{v}_i\cdot\nabla)^2 + \alpha(\mathbf{v}_j\cdot\nabla)^2\right) f(\mathbf{x})\right] = f(\mathbf{x})\frac{1-3\alpha}{\sigma^4}\equiv 0\quad (j\neq i).$$

It is then easy to see that for $\alpha=-1/3$ we achieve our objective of setting the term to 0. This is what Meijering et al. (2004) also found (using a similar notation).


#### Going to n-D
We establish an augmented Hessian matrix in any dimension as follows. Define $H_f(\mathbf{x})$ as the standard Hessian matrix, then the augmented form is created by computing $H_f-\alpha H_f + \alpha \mathrm{Tr}(H_f)I$ where $\mathrm{Tr}$ is the trace of the matrix, and $I$ is the identity matrix. In 3D that gives the following matrix,

$$H_f'(\mathbf{x}) = H_f-\alpha H_f + \alpha \mathrm{Tr}(H_f)I = \left[ \begin{array}{ccc}
f_{xx}+\alpha (f_{yy}+ f_{zz}) & (1-\alpha)f_{xy} & (1-\alpha)f_{xz} \\
(1-\alpha)f_{xy} & f_{yy}+\alpha(f_{xx}+f_{zz}) & (1-\alpha)f_{yz} \\
(1-\alpha)f_{yz} & (1-\alpha)f_{xz} & f_{zz}+\alpha(f_{xx}+f_{yy}) \end{array} \right].$$

Assume that the original Hessian $H_f(\mathbf{x})$ allowed the eigendecomposition into orthonormal eigenvectors $\mathbf{v}\_i$ with corresponding eigenvalues $\lambda_i$ such that $H_f\mathbf{v}\_i=\lambda_i\mathbf{v}\_i$. Using this identity when multiplying the augmented Hessian with the same eigenvector $\mathbf{v}\_i$ we then find that this eigenvector is also a solution to the augmented system:

$$ H'_f\mathrm{v}\_i = \left(H_f-\alpha H_f + \alpha \mathrm{Tr}(H_f)I\right)\mathbf{v}\_i = \underbrace{\left( \lambda_i - \alpha\lambda_i + \alpha\mathrm{Tr}(H_f)\lambda_i \right)}_{\mathrm{new eigenvalues}}\mathrm{v}\_i $$

From linear algebra, we know that [the trace of a matrix equals the sum of its eigenvalues](https://en.wikipedia.org/wiki/Trace_(linear_algebra)#Trace_as_the_sum_of_eigenvalues), $\mathrm{Tr}(H_f)=\sum_i \lambda_i$, thus we can find that the new eigenvalues are described by a similar rule as in the 2-D case:

$$ \lambda_i' = \lambda_i + \sum_{j\neq i} \alpha \lambda_j. $$

**Just to re-cap what we found so far**: the eigenvalues that satisfy $H_f(\mathbf{x})$ also satisfied the augmented Hessian $H_f'$, and the eigenvalues are described as $\lambda_i=\lambda_i+\sum_{j\neq i}\alpha \lambda_j$. This generalizes the 2-D case to $n$ dimensions.

As previously done, we can ascribe a meaning to each un-primed eigenvalue, $\lambda_i=f_{\mathbf{v}\_i\mathbf{v}\_i}$, that is, the eigenvalue corresponds to the second derivative in the corresponding eigenvector direction. That means that we have

$$ \mathbf{v}_i^T H_f'(\mathbf{x}) \mathbf{v}_i = \lambda_i' = \lambda_i + \sum_{j\neq i}\alpha \lambda_j = \left[(\mathbf{v}_i\cdot\nabla)^2 + \sum_{j\neq i} \alpha (\mathbf{v}_j \cdot \nabla)^2 \right]f(\mathbf{x}). $$

In words, thus, we have that the *primed* $\lambda_i'$ corresponds to the directional derivative in direction $\mathbf{v}\_i$ as well as the sum of directional derivatives in the other orthonormal directions, scaled by $\alpha$. If we want to minimize this term, we want that the filter is maximally 'straight' in the direction corresponding to the smallest eigenvalue. Say, the eigenvalues are ordered from largest to smallest, then $\mathbf{v}\_n$ is the smallest eigenvector (note that the next part of the derivation does not in any way depend on this ordering; I suspect that we just want maximum 'straightness' in all directions...), and we want

$$ \lim_{\mathbf{x}\to 0} \left[ (\mathbf{v}_n\cdot\nabla)^2(\mathbf{v}_i^T H_f'(\mathbf{x}) \mathbf{v}_i) \right]f(\mathbf{x}) = 0$$

We end up doing the same kind of calculations as we did in the 2-D case. We expand the expressions, and [use the following derivative relations](https://www.wolframalpha.com/input?i=Lim%5B+D%5BD%5BD%5BD%5Be%5E%28-%28x%5E2%2By%5E2%2Bz%5E2%2Ba%5E2%29%2F%282*s%5E2%29%29%2C+%7Bx%2C1%7D%5D%2C+%7By%2C1%7D%5D%2C+%7Bz%2C1%7D%5D%2C%7Ba%2C1%7D%5D+%2C+%7Bx-%3E0%2C+y-%3E0%2C+z-%3E0%7D%5D) (realize that at the limit of $\mathbf{x}\to 0$ we have $f_{xxyy}=f_{xxzz}=f_{yyxx}$ etc., so the actual indices are not what's relevant, only the relative occurence of any given derivative relation):

$$
\lim_{\mathbf{x}\to 0}  \left[ \begin{array}{ccc}
f_{xxxx} \\
f_{xxxy} \\
f_{xxyy} \\
f_{xyyy} \\
f_{yyyy} \\
f_{xxyz} \\
f_{xyza}\end{array} \right] = f(0)\left[ \begin{array}{ccc}
3/\sigma^2 \\
0 \\
1/\sigma^2 \\
0 \\
3/\sigma^2 \\ 0 \\ 0\end{array} \right] .
$$

We find that the limit may be expanded into $n$ terms of multiplications of terms $(\mathbf{v}\_j\cdot\nabla)^2(\mathbf{v}\_i\cdot\nabla)^2f(\mathbf{x})$, which follow the same relation as found in the 2-D case:

$$ \lim_{\mathbf{x}\to 0} (\mathbf{v}\_j\cdot\nabla)^2(\mathbf{v}\_i\cdot\nabla)^2 f(\mathbf{x}) = \begin{cases} \frac{f(0)}{\sigma} &\mathrm{if}\ i\neq j, \\ \frac{3f(0)}{\sigma}&\mathrm{if}\ i=j. \end{cases}$$

For example, in 3-D for vectors $\mathbf{v}\_i^T=(v_1\quad v_2\quad v_3)$ and $\mathbf{v}\_j^T=(r_1\quad r_2\quad r_3)$ we expend the terms to find
$$ \lim_{\mathbf{x}\to 0} (\left[ \begin{array}{ccc}
r_1 \\
r_2 \\
r_3 \end{array} \right]\cdot\nabla)^2(\left[ \begin{array}{ccc}
v_1 \\
v_2 \\
v_3 \end{array} \right]\cdot\nabla)^2 f(\mathbf{x}) = f(0)\frac{4 r_2 r_3 v_2 v_3 + 4 r_1 v_1 (r_2 v_2 + r_3 v_3) + r_1^2 (3 v_1^2 + v_2^2 + v_3^2) + r_2^2 (v_1^2 + 3 v_2^2 + v_3^2) +  r_3^2 (v_1^2 + v_2^2 + 3 v_3^2)}{\sigma^2} $$

We now again look for instances of the positive inner product ($r_1v_1+r_2v_2+r_3v_3=0$) and 

$$ \lim_{\mathbf{x}\to 0} \left[ (\mathbf{v}_n\cdot\nabla)^2(\mathbf{v}_i^T H_f'(\mathbf{x}) \mathbf{v}_i) \right]f(\mathbf{x}) = f(0)\frac{1+\alpha-n\alpha}{\sigma^4}$$
