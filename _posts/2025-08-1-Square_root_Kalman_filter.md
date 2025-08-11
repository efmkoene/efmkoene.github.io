---
layout: post
title: "The Square Root Kalman Filter"
subtitle: ""
tags: [atmospheric modeling, kalman filter]
---

## Introduction
Consider we have access to observations $\mathbf{y}\_\mathrm{obs} \in \mathbb{R}^{m}$ and they are due to some linear operator acting on a 'hidden state' $\mathbf{x}\_\mathrm{true}\in\mathbb{R}^{n}$ with added noise $\mathbf{e}$ like
\begin{align}
    \mathbf{y}\_\mathrm{obs} = \mathbf{H}\mathbf{x}\_\mathrm{true} + \mathbf{e},
\end{align}
where $\mathbf{H} \in \mathbb{R}^{n\times m}$, and the error has zero mean $\mathbb{E}[\mathbf{e}]=0$ but a covariance structure like $\mathbb{E}[\mathbf{e}\mathbf{e}^{\mathsf{T}}]=\mathbf{R}$ which we call the *observation error covariance*. For example, if $\mathbf{R}$ is diagonal, its entries simply correspond to the variance (typically denoted as $\sigma^2$) of each measurement.

Typically, we don't know $\mathbf{x}\_\mathrm{true}$ and want to estimate it based on the data. Assuming Gaussian distributed errors we find that we need to minimize a cost function of the form 

\begin{align}
J(\mathbf{x}) = (\mathbf{x} - \mathbf{x}_b) \mathbf{P}_b^{-1}(\mathbf{x} - \mathbf{x}_b)^{\mathsf{T}} + (\mathbf{y}\_\mathrm{obs}-\mathbf{H}\mathbf{x})\mathbf{R}^{-1}(\mathbf{y}\_\mathrm{obs}-\mathbf{H}\mathbf{x})^{\mathsf{T}},\tag{1}\label{eq:costfunction}
\end{align}

where $\mathbf{x}_b$ is a prior model (''background''), $\mathbf{P}_b$ is the *error covariance model* or *background covariance model*, i.e., the uncertainty about the prior model. The minimizer of $J(x)$ in eq. \eqref{eq:costfunction} gives the optimal estimate, also known as the 'posterior' or 'analysis' mean (see, e.g., Fichtner 2021, *Lecture Notes on Inverse Theory*),

\begin{equation}
\mathbf{x}_a = \mathbf{x}_b + \mathbf{K} (\mathbf{y}\_{\mathrm{obs}} - \mathbf{H} \mathbf{x}_b),\tag{2}\label{eq:mean_update}
\end{equation}

where the optimal Kalman gain $\mathbf{K}$ is:

\begin{equation}
\mathbf{K} = \mathbf{P}_b \mathbf{H}^{\mathsf{T}} (\mathbf{H} \mathbf{P}_b \mathbf{H}^{\mathsf{T}} + \mathbf{R})^{-1},\label{eq:kalman_gain}
\end{equation}

and the analysis covariance is

\begin{equation}
\mathbf{P}_a = (\mathbf{I} - \mathbf{K} \mathbf{H}) \mathbf{P}_b.\tag{3}\label{eq:cov_update}
\end{equation}

Equations (\ref{eq:mean_update})--(\ref{eq:cov_update}) form the foundation of the stationary Kalman filter (i.e., those without state-space propagation), but also of analytical inversions or variational approaches.

## Square root approach
The original *square root* formulation of the Kalman filter, attributed to J. Potter, was developed in the context of the Apollo mission to minimize the 'word length' (i.e., the number of bits that a computer could process simultaneously). For instance, if the matrix $\mathbf{R}$ is diagonal with entries $\sigma_i^2$ along the diagonal, and certain measurements possess significantly higher precision than others, then the variations between $\sigma_i^2$ values necessitate greater word lengths compared to the variations between $\sigma_i$ values. However, it was later discovered that the square root formulation is typically also more stable and often the appropriate way to implement a Kalman filter in practice. In particular, we are talking about the square root of the $\mathbf{P}_b$ and $\mathbf{P}_a$ matrices here. 

In short, the square root of a matrix is any matrix $\mathbf{Z}$ for which $\mathbf{Z}\mathbf{Z}^{\mathsf{T}}$ forms that matrix. So, if we can write $\mathbf{P}_b=\mathbf{Z}\mathbf{Z}^{\mathsf{T}}$ then $\sqrt{\mathbf{P}_b}=\mathbf{Z}$. (Square root matrices are not necessarily unique). We can thus re-write the Kalman equations using this substitution,

\begin{align}
    \mathbf{x}_a &= \mathbf{x}_b + \mathbf{K}(\mathbf{y}\_\mathrm{obs} - \mathbf{H}\mathbf{x}_b), \\\\
    \mathbf{K} &= \mathbf{Z}_b\mathbf{Z}_b^{\mathsf{T}} \mathbf{H}^{\mathsf{T}} (\mathbf{H} \mathbf{Z}_b\mathbf{Z}_b^{\mathsf{T}} \mathbf{H}^{\mathsf{T}} + \mathbf{R})^{-1}, \\\\
    \mathbf{Z}_a\mathbf{Z}_a^{\mathsf{T}} &= (\mathbf{I} - \mathbf{K} \mathbf{H}) \mathbf{Z}_b\mathbf{Z}_b^{\mathsf{T}}.\tag{4}\label{eq:kalgain_notyet}
\end{align}

The advantage is not immediately clear, but we can make a substitution $\mathbf{Y}_b=\mathbf{H}\mathbf{Z}_b$ to find

\begin{align}
    \mathbf{K} &= \mathbf{Z}_b\mathbf{Y}_b^{\mathsf{T}} (\mathbf{Y}_b\mathbf{Y}_b^{\mathsf{T}} + \mathbf{R})^{-1},
\end{align}

Although this form is mathematically valid, it still requires us to explicitly compute terms like $\mathbf{Z}_b\mathbf{Z}_b^{\mathsf{T}}$ which means it suffers from the same computational limitations. Hence, we can go one step further and compute one possible square root of eq. \eqref{eq:kalgain_notyet} which is

\begin{align}
    \mathbf{Z}_a = \mathbf{Z}_b\left(\mathbf{I} - \mathbf{Y}_b^{\mathsf{T}}\sqrt{\mathbf{Y}_b\mathbf{Y}_b^{\mathsf{T}}+\mathbf{R}}^{-\top}\left(\sqrt{\mathbf{Y}_b\mathbf{Y}_b^{\mathsf{T}}+\mathbf{R}}+\sqrt{\mathbf{R}}\right)^{-1}\mathbf{Y}_b\right).
\end{align}

This square root formulation was given by Andrews (1968, *'A square root formulation of the Kalman covariance equations'*), and a derivation is given in the appendix here. Thus, the square root formulation of the Kalman filter is

$$
\begin{align}
    \mathbf{x}_a &= \mathbf{x}_b + \mathbf{K}(\mathbf{y}\_\mathrm{obs} - \mathbf{H}\mathbf{x}_b),\tag{3}\label{eq:squaremean} \\
    \mathbf{K} &= \mathbf{Z}_b\mathbf{Y}_b^{\mathsf{T}} (\mathbf{Y}_b\mathbf{Y}_b^{\mathsf{T}} + \mathbf{R})^{-1}, \\
    \mathbf{Z}_a &= \mathbf{Z}_b\left(\mathbf{I} - \mathbf{Y}_b^{\mathsf{T}}\sqrt{\mathbf{Y}_b\mathbf{Y}_b^{\mathsf{T}}+\mathbf{R}}^{-\top}\left(\sqrt{\mathbf{Y}_b\mathbf{Y}_b^{\mathsf{T}}+\mathbf{R}}+\sqrt{\mathbf{R}}\right)^{-1}\mathbf{Y}_b\right).\tag{6}\label{eq:squareZa}
\end{align}
$$

The square root Kalman filter in \eqref{eq:squaremean}--\eqref{eq:squareZa} works exclusively with the square roots of (combinations of) the covariance matrices, which means this system overcomes the main limitations of the 'word length' indicated above.

## Ensemble square root approach
Using the square root form of the Kalman filter improves numerical stability and reduces issues with limited precision. However, when the state vector is very large, computations like $\mathbf{H}\mathbf{Z}_b$ become expensive. The Ensemble Kalman Filter (EnKF), particularly in its square root form, addresses this by approximating $\mathbf{Z}_b$ using a reduced-rank representation.

We define $\mathbf{x}_b = \mathbf{Z}_b \mathbf{G}$, where $\mathbf{G} \in \mathbb{R}^{n \times N}$ is a random Gaussian matrix with zero mean and unit variance (e.g., as produced by \texttt{randn} in MATLAB or \texttt{numpy.random.randn} in Python), and $N<n$, often $N \ll n$. Since $\mathbb{E}[\mathbf{G}\mathbf{G}^{\mathsf{T}}] = \mathbf{I}$, we have:

\begin{align}
    \mathbb{E}[\mathbf{Z}_b\mathbf{G}\mathbf{G}^{\mathsf{T}}\mathbf{Z}_b^{\mathsf{T}}] &= \mathbf{Z}_b \mathbb{E}[\mathbf{G}\mathbf{G}^{\mathsf{T}}] \mathbf{Z}_b^{\mathsf{T}} \\\\
    &= \mathbf{Z}_b\mathbf{Z}_b^{\mathsf{T}},\\\\
    & = \mathbf{P}_b.
\end{align}

Hence, $\mathbf{x}_b$ is a rank-reduced square root approximation of $\mathbf{P}_b$, with shape $\mathbb{R}^{n \times N}$. We then approximate:

\begin{equation}
    \mathbf{P}_b = \mathbf{Z}_b\mathbf{Z}_b^{\mathsf{T}} \approx \frac{1}{N}\mathbf{x}_b\mathbf{x}_b^{\mathsf{T}},
\end{equation}

and analogously for $\mathbf{P}_a$. Thus, we can formulate the square root EnKF by replacing all occurences of $\mathbf{Z}$ in the square root formulation of the Kalman filter with $\mathbf{x}/\sqrt{N}$, knowing that we are making an approximation.\footnote{In the literature, the factor $\mathbf{G}\mathbf{G}^{\mathsf{T}}/(N-1)$ is often used as the \textit{unbiased} estimator of the covariance of $\mathbf{G}$, which tends to the identity matrix for large $N$ (often only when $N\gg n$). However, since $\mathbf{G}$ has zero mean, no degree of freedom is used up, and in my view, from the derivation presented here, the correct factor is $1/N$; the literature is simply mistaken. Rather, the other literature starts off by drawing $N$ members first, subtracting their mean from each member to get to 'state vector deviations'. We skipped any such a step here, although our equations are now identical save for the factor $N-1$. We are not doing a statistical trick and merely approximate one quantity in the Kalman filter. Of course, the practical difference is negligible for $N > 100$.} Thus the bulk implementation of the square root EnKF becomes (using an apostrophe [$'$] to indicate we are approximating the quantities):

\begin{equation}
\begin{aligned}
    \mathbf{x}_b' &= \mathbf{Z}_b\mathbf{G},\\\\
    \mathbf{Y}' &= \mathbf{H} \mathbf{x}_b', \label{eq:bulkensupdate}\\\\
    \mathbf{K}' &= \frac{1}{N} \mathbf{x}_b' \mathbf{Y}'{}^{\mathsf{T}} \left( \frac{1}{N} \mathbf{Y}' \mathbf{Y}'{}^{\mathsf{T}} + \mathbf{R} \right)^{-1}, \\\\
    \mathbf{x}'_a &= \mathbf{x}_b + \mathbf{K}'(\mathbf{y}\_\mathrm{obs} - \mathbf{H}\mathbf{x}_b), \\\\
    \mathbf{x}'_a &= \mathbf{x}_b'\left(\mathbf{I} - \frac{1}{N}\mathbf{Y}'{}^{\mathsf{T}}\sqrt{\frac{1}{N}\mathbf{Y}'\mathbf{Y}'{}^{\mathsf{T}}+\mathbf{R}}^{-\top}\left(\sqrt{\frac{1}{N}\mathbf{Y}'\mathbf{Y}'{}^{\mathsf{T}}+\mathbf{R}}+\sqrt{\mathbf{R}}\right)^{-1}\mathbf{Y}'\right).
\end{aligned}
\end{equation}
