---
layout: post
title: "The Stationary Square Root Kalman Filter"
subtitle: ""
tags: [atmospheric modeling, kalman filter]
---

## Introduction -- the stationary Kalman Filter
Consider we have access to observations $\mathbf{y}\_\mathrm{obs} \in \mathbb{R}^{m}$ and they are due to some linear operator acting on a 'hidden state' $\mathbf{x}\_\mathrm{true}\in\mathbb{R}^{n}$ with added noise $\mathbf{e}$ like
\begin{align}
    \mathbf{y}\_\mathrm{obs} = \mathbf{H}\mathbf{x}\_\mathrm{true} + \mathbf{e},\tag{1}\label{eq:obsmodel}
\end{align}
where $\mathbf{H} \in \mathbb{R}^{n\times m}$ is the so-called 'observation operator', while the error has zero mean $\mathbb{E}[\mathbf{e}]=0$ and a covariance structure like $\mathbb{E}[\mathbf{e}\mathbf{e}^{\mathsf{T}}]=\mathbf{R}$ which we call the *observation error covariance*. For example, if $\mathbf{R}$ is diagonal, its entries simply correspond to the variance (typically denoted as $\sigma^2$) of each measurement.

Typically, we don't know $\mathbf{x}\_\mathrm{true}$ and want to estimate it based on the data. Assuming Gaussian distributed errors we find that we need to minimize a cost function of the form 

\begin{align}
J(\mathbf{x}) = (\mathbf{x} - \mathbf{x}_b) \mathbf{P}_b^{-1}(\mathbf{x} - \mathbf{x}_b)^{\mathsf{T}} + (\mathbf{y}\_\mathrm{obs}-\mathbf{H}\mathbf{x})\mathbf{R}^{-1}(\mathbf{y}\_\mathrm{obs}-\mathbf{H}\mathbf{x})^{\mathsf{T}},\tag{2}\label{eq:costfunction}
\end{align}

where $\mathbf{x}_b$ is a prior model (''background''), $\mathbf{P}_b$ is the *error covariance model* or *background covariance model*, i.e., the uncertainty about the prior model. The minimizer of $J(x)$ in eq. \eqref{eq:costfunction} gives the optimal estimate, also known as the 'posterior' or 'analysis' mean (see, e.g., Fichtner 2021, *Lecture Notes on Inverse Theory*),

\begin{equation}
\mathbf{x}_a = \mathbf{x}_b + \mathbf{K} (\mathbf{y}\_{\mathrm{obs}} - \mathbf{H} \mathbf{x}_b),\tag{3}\label{eq:mean_update}
\end{equation}

where the optimal Kalman gain $\mathbf{K}$ is:

\begin{equation}
\mathbf{K} = \mathbf{P}_b \mathbf{H}^{\mathsf{T}} (\mathbf{H} \mathbf{P}_b \mathbf{H}^{\mathsf{T}} + \mathbf{R})^{-1},\label{eq:kalman_gain}
\end{equation}

and the analysis covariance is

\begin{equation}
\mathbf{P}_a = (\mathbf{I} - \mathbf{K} \mathbf{H}) \mathbf{P}_b.\tag{4}\label{eq:cov_update}
\end{equation}

Equations (\ref{eq:mean_update})--(\ref{eq:cov_update}) form the foundation of the stationary Kalman filter (i.e., those without state-space propagation), but also of analytical inversions or variational approaches.

## Square root approach
The original *square-root* formulation of the Kalman filter -- attributed to J. Potter -- was [developed during the Apollo program](https://ntrs.nasa.gov/api/citations/19860003843/downloads/19860003843.pdf) to address hardware limitations, specifically the *word length* of onboard computers (i.e., the number of bits processed in parallel). The motivation was that numerical ranges arising from covariance values could require more bits for accurate representation than their square roots. For example, if the measurement noise covariance matrix $\mathbf{R}$ is diagonal with entries $\sigma_i^2$ on the diagonal, and certain measurements are orders of magnitude more precise than others, then the spread in $\sigma_i^2$ values can be far larger than the spread in $\sigma_i$ values. Working with the latter therefore reduces the required word length. Although originally a hardware-oriented optimization, the square-root formulation was later found to have a significant numerical advantage: it is inherently more stable and better conditioned (i.e., tends to be more stable in finite-precision arithmetic, and has a better  less prone to round-off errors). As a result, it is now often the preferred way to implement the Kalman filter in practice.

In short, the square root of a matrix is any matrix $\mathbf{Z}$ for which $\mathbf{Z}\mathbf{Z}^{\mathsf{T}}$ forms that matrix. So, if we can write $\mathbf{P}_b=\mathbf{Z}_b\mathbf{Z}_b^{\mathsf{T}}$ then $\sqrt{\mathbf{P}_b}=\mathbf{Z}_b$. (Square root matrices are not necessarily unique). We can thus re-write the Kalman equations using this substitution,

$$
\begin{align}
    \mathbf{x}_a &= \mathbf{x}_b + \mathbf{K}(\mathbf{y}_\mathrm{obs} - \mathbf{H}\mathbf{x}_b), \\
    \mathbf{K} &= \mathbf{Z}_b\mathbf{Z}_b^{\mathsf{T}} \mathbf{H}^{\mathsf{T}} (\mathbf{H} \mathbf{Z}_b\mathbf{Z}_b^{\mathsf{T}} \mathbf{H}^{\mathsf{T}} + \mathbf{R})^{-1}, \\
    \mathbf{Z}_a\mathbf{Z}_a^{\mathsf{T}} &= (\mathbf{I} - \mathbf{K} \mathbf{H}) \mathbf{Z}_b\mathbf{Z}_b^{\mathsf{T}}.\tag{5}\label{eq:kalgain_notyet}
\end{align}
$$
The advantage is not immediately clear, but we can make a substitution $\mathbf{Y}_b=\mathbf{H}\mathbf{Z}_b$ (compare to eq. \eqref{eq:obsmodel}) to find

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
    \mathbf{x}_a &= \mathbf{x}_b + \mathbf{K}(\mathbf{y}_\mathrm{obs} - \mathbf{H}\mathbf{x}_b),\tag{6}\label{eq:squaremean} \\
    \mathbf{K} &= \mathbf{Z}_b\mathbf{Y}_b^{\mathsf{T}} (\mathbf{Y}_b\mathbf{Y}_b^{\mathsf{T}} + \mathbf{R})^{-1}, \\
    \mathbf{Z}_a &= \mathbf{Z}_b\left(\mathbf{I} - \mathbf{Y}_b^{\mathsf{T}}\sqrt{\mathbf{Y}_b\mathbf{Y}_b^{\mathsf{T}}+\mathbf{R}}^{-\top}\left(\sqrt{\mathbf{Y}_b\mathbf{Y}_b^{\mathsf{T}}+\mathbf{R}}+\sqrt{\mathbf{R}}\right)^{-1}\mathbf{Y}_b\right).\tag{7}\label{eq:squareZa}
\end{align}
$$

The square root Kalman filter in \eqref{eq:squaremean}--\eqref{eq:squareZa} works exclusively with the square roots of (combinations of) the covariance matrices, which means this system overcomes the main limitations of the 'word length' indicated above.

## Ensemble square root approach
Using the square root form of the Kalman filter improves numerical stability and reduces issues with limited precision. However, when the state vector is very large, computations like $\mathbf{H}\mathbf{Z}_b$ become expensive. The Ensemble Kalman Filter (EnKF), particularly in its square root form, addresses this by approximating $\mathbf{Z}_b$ using a reduced-rank representation.

We define $\mathbf{X}_b = \mathbf{Z}_b \mathbf{G}$, where $\mathbf{G} \in \mathbb{R}^{n \times N}$ is a random Gaussian matrix with zero mean and unit variance (e.g., as produced by `randn` in MATLAB or `numpy.random.randn` in Python), and $N<n$, often $N \ll n$. Since $\mathbb{E}[\mathbf{G}\mathbf{G}^{\mathsf{T}}] = \mathbf{I}$, we have:

$$
\begin{align}
    \mathbb{E}[\mathbf{Z}_b\mathbf{G}\mathbf{G}^{\mathsf{T}}\mathbf{Z}_b^{\mathsf{T}}] &= \mathbf{Z}_b \mathbb{E}[\mathbf{G}\mathbf{G}^{\mathsf{T}}] \mathbf{Z}_b^{\mathsf{T}} \\
    &= \mathbf{Z}_b\mathbf{Z}_b^{\mathsf{T}},\\
    & = \mathbf{P}_b.
\end{align}
$$

Hence, $\mathbf{X}_b$ with shape $\mathbb{R}^{n \times N}$ allows for a rank-reduce approximation $\mathbf{P}_b \approx \mathbb{E}[\mathbf{X}_b\mathbf{X}_b^{\mathsf{T}}]=\mathbf{X}_b\mathbf{X}_b^{\mathsf{T}}/N$. We then approximate:

\begin{equation}
    \mathbf{P}_b = \mathbf{Z}_b\mathbf{Z}_b^{\mathsf{T}} \approx \frac{1}{N}\mathbf{X}_b\mathbf{X}_b^{\mathsf{T}},
\end{equation}

and analogously for $\mathbf{P}_a$. Thus, we can formulate the square root EnKF by replacing all occurences of $\mathbf{Z}$ in the square root formulation of the Kalman filter with $\mathbf{X}/\sqrt{N}$, knowing that we are making an approximation.<sup id="fnref-1">[1](#fn-1)</sup> Thus the bulk implementation of the square root EnKF becomes (using an apostrophe [$'$] to indicate we are approximating the quantities):

$$
\begin{aligned}
    \mathbf{X}_b' &= \mathbf{Z}_b\mathbf{G},\\
    \mathbf{Y}' &= \mathbf{H} \mathbf{X}_b', \\
    \mathbf{K}' &= \frac{1}{N} \mathbf{x}_b' \mathbf{Y}'{}^{\mathsf{T}} \left( \frac{1}{N} \mathbf{Y}' \mathbf{Y}'{}^{\mathsf{T}} + \mathbf{R} \right)^{-1}, \\
    \mathbf{x}'_a &= \mathbf{x}_b + \mathbf{K}'(\mathbf{y}_\mathrm{obs} - \mathbf{H}\mathbf{x}_b), \\
    \mathbf{X}'_a &= \mathbf{X}_b'\left(\mathbf{I} - \frac{1}{N}\mathbf{Y}'{}^{\mathsf{T}}\sqrt{\frac{1}{N}\mathbf{Y}'\mathbf{Y}'{}^{\mathsf{T}}+\mathbf{R}}^{-\top}\left(\sqrt{\frac{1}{N}\mathbf{Y}'\mathbf{Y}'{}^{\mathsf{T}}+\mathbf{R}}+\sqrt{\mathbf{R}}\right)^{-1}\mathbf{Y}'\right).
\end{aligned}\tag{8}\label{eq:bulkensupdate}
$$

It is worth considering for a moment what is required to implement the ensemble square root Kalman filter as described here (whether through a bulk or sequential implementation). We first build the full error covariance matrix $\mathbf{P}_b$. We compute its square root $\mathbf{Z}_b$ and multiply it with a Gaussian random matrix of reduced shape, $\mathbf{X}_b=\mathbf{Z}_b\mathbf{G}$. We then compute $\mathbf{H}\mathbf{X}_b$ which yields $N$ (rather than $n$) 'ensemble' realizations, i.e., all equally plausible model realizations. Essentially, our computations are cheaper by a factor $N/n$. If $N$ is chosen too small, however, $\mathbf{X}_b\mathbf{X}_b^{\mathsf{T}}/N$ does not tend to $\mathbf{P}_b$ very well, so $N$ must be small-but-not-too-small. If $\mathbf{P}_b$ is purely diagonal, this strategy will not work, because we cannot rank-reduce it. Therefore, we assume (or must make sure) that $\mathbf{P}_b$ exhibits considerable covariance between the various state vector elements to make a rank reduction possible. Finally, once $\mathbf{H}\mathbf{X}_b$ and $\mathbf{H}\mathbf{x}_b$ are computed, we have all the ingredients to compute the (square root form of the) Kalman filter which we compute by computing eq. \eqref{eq:bulkensupdate}.


## Implementations
### Square root filter
#### Bulk implementation square root filter
This algorithm can be implemented in 'bulk' mode (i.e., assimilating all $\mathbf{y}\_\mathrm{obs}$ at once) through a direct implementation of eqs. \eqref{eq:squaremean}--\eqref{eq:squareZa}. That is, one formulates $\mathbf{P}_b$ and using its Cholesky decomposition $\mathbf{P}_b=\mathbf{Z}_b\mathbf{Z}_b^{\mathsf{T}}$ and then computing $\mathbf{Y}_b=\mathbf{H}\mathbf{Z}_b$ one can compute all the relevant parameters. We can construct the posterior covariance (if needed) by computing $\mathbf{Z}_a\mathbf{Z}_a^{\mathsf{T}}$.

#### Sequential implementation of square root filter I -- basic implementation
If $\mathbf{R}$ is a diagonal matrix (i.e., all measurements are independent) we can assimilate the observations one-by-one (i.e., sequentially). 
We then take 
$\left[\mathbf{y}\_\mathrm{obs}\right]\_i$
for observation $i$, 
$[\mathbf{H}]\_i$ 
to select row $i$, 
$\[\mathbf{R}\]\_{ii}$ as the diagonal observation error component, we get the following algorithm. 
Starting with $\mathbf{Z}^{(0)} = \mathbf{Z}_b$ and $\mathbf{x}^{(0)}=\mathbf{x}_b$ we compute

$$
\begin{aligned}
\mathbf{a} &= [\mathbf{H}]_{i} \mathbf{Z}^{(i-1)}&&\in\mathbb{R}^{1\times n}, \\
    b & = \mathbf{a}\mathbf{a}^{\mathsf{T}}+ [\mathbf{R}]_{ii} && \in\mathbb{R},\\
    \alpha & = \left(1+\sqrt{\frac{[\mathbf{R}]_{ii}}{b}}\right)^{-1} &&\in\mathbb{R}, \\
    \mathbf{K} & = \frac{\mathbf{Z}_b^{(i-1)}\mathbf{a}^{\mathsf{T}}}{b} &&\in \mathbb{R}^{n\times 1},\\
    \mathbf{x}^{(i)} & = \mathbf{x}^{(i-1)} +\mathbf{K}([\mathbf{y}_\mathrm{obs}]_i - [\mathbf{H}]_i\mathbf{x}^{(i-1)}) && \in\mathbb{R}^{n\times 1}, \\
    \mathbf{Z}^{(i)} & = \mathbf{Z}^{(i-1)} - \alpha \mathbf{K} \mathbf{a} && \in\mathbb{R}^{n\times n}.
\end{aligned}
\tag{9}\label{eq:Zupdate}
$$
After all $m$ observations are assimilated, $\mathbf{x}^{(m)}=\mathbf{x}_a$ and $\mathbf{Z}^{(m)}=\mathbf{Z}_a$. Note how we substituted the Kalman gain into the expression for the square root of the error covariance matrix. This was possible because terms like $\mathbf{Y}_b\mathbf{Y}_b^{\mathsf{T}}+\mathbf{R}$ are simple scalars when considering single updates.


#### Sequential implementation of square root filter II -- pre-whitening the observations
The sequential update scheme is efficient when many observations are assimilated at once, but it requires the observations to be independent, which is a stringent requirement. We can get around this by `whitening' the observations such that they are guaranteed to be independent. We do this by pre-multiplying our measurement model by $\sqrt{\mathbf{R}}^{-1}$ such that
\begin{equation}
    \sqrt{\mathbf{R}}^{-1}\mathbf{y}\_\mathrm{obs} = \sqrt{\mathbf{R}}^{-1}\mathbf{H}\mathbf{x}\_\mathrm{true} +\sqrt{\mathbf{R}}^{-1}\mathbf{e}.
\end{equation}
The expectation of this noise model remains $\mathbb{E}[\sqrt{\mathbf{R}}^{-1}\mathbf{e}]=0$ while 

$$
\begin{align}
\mathbb{E}[\sqrt{\mathbf{R}}^{-1}\mathbf{e}\mathbf{e}^{\mathsf{T}}\sqrt{\mathbf{R}}^{-1}]&=\sqrt{\mathbf{R}}^{-1}\mathbb{E}[\mathbf{e}\mathbf{e}^T]\sqrt{\mathbf{R}}^{-T}, \\
&=\sqrt{\mathbf{R}}^{-1}\mathbf{R}\sqrt{\mathbf{R}}^{-T}, \\
& =\sqrt{\mathbf{R}}^{-1}\sqrt{\mathbf{R}}\sqrt{\mathbf{R}}^{\mathsf{T}}\sqrt{\mathbf{R}}^{-T}, \\
&=\mathbf{I}.
\end{align}
$$ 

Hence, this pre-multiplication of the data makes the noise uncorrelated and of uniform variance. We can, therefore, simply change our algorithm into the following form

$$
\begin{aligned}
    \mathbf{a} &= [\sqrt{\mathbf{R}}^{-1}\mathbf{H}]_{i} \mathbf{Z}^{(i-1)}&&\in\mathbb{R}^{1\times n}, \\
    b & = \mathbf{a}\mathbf{a}^{\mathsf{T}}+ 1 && \in\mathbb{R},\\
    \alpha & = \left(1+\sqrt{\frac{1}{b}}\right)^{-1} &&\in\mathbb{R}, \\
    \mathbf{K} & = \frac{\mathbf{Z}_b^{(i-1)}\mathbf{a}^{\mathsf{T}}}{b} &&\in \mathbb{R}^{n\times 1},\\
    \mathbf{x}^{(i)} & = \mathbf{x}^{(i-1)} +\mathbf{K}([\sqrt{\mathbf{R}}^{-1}\mathbf{y}_\mathrm{obs}]_i - [\sqrt{\mathbf{R}}^{-1}\mathbf{H}]_i\mathbf{x}^{(i-1)}) && \in\mathbb{R}^{n\times 1}, \\
    \mathbf{Z}^{(i)} & = \mathbf{Z}^{(i-1)} - \alpha \mathbf{K} \mathbf{a} && \in\mathbb{R}^{n\times n}.
\end{aligned}
$$

### Sequential implementation of square root filter III -- independent of H
In the two sequential schemes, we have to apply $\mathbf{H}$ to the updated mean and covariance matrix as each measurement comes in, so $m$ times in serial. We can consider an alternative scheme that is free of such re-application of $\mathbf{H}$ which is, for example, of interest when $\mathbf{H}$ is expensive to compute in serial but cheap(er) to compute in parallel. The scheme rests on the fact that 
$[\mathbf{H}]_i\mathbf{Z}^{(i-1)}=[\mathbf{H}\mathbf{Z}^{(i-1)}]_i$, 
and that by applying $\mathbf{H}$ to eq. \eqref{eq:Zupdate} we have an expression depending only on factors including $\mathbf{H}\mathbf{Z}^{(i-1)}$. The same holds for the mean update. We again start with $\mathbf{Z}^{(0)}=\mathbf{Z}_b$, $\mathbf{x}^{(0)}=\mathbf{x}_b$ but now also write $\mathbf{Y}^{(0)}=\mathbf{H}\mathbf{Z}_b$ and $\mathbf{y}^{(0)}=\mathbf{H}\mathbf{x}_b$. Then,

$$
\begin{aligned}
    y & = [\mathbf{y}^{(i-1)}]_{i} &&\in\mathbb{R},\\
    \mathbf{a} &= [\mathbf{Y}^{(i-1)}]_{i} &&\in\mathbb{R}^{1\times n}, \\
    b & = \mathbf{a}\mathbf{a}^{\mathsf{T}}+ [\mathbf{R}]_{ii} && \in\mathbb{R},\\
    \alpha & = \left(1+\sqrt{\frac{[\mathbf{R}]_{ii}}{b}}\right)^{-1} &&\in\mathbb{R}, \\
    \mathbf{K} & = \frac{\mathbf{Z}_b^{(i-1)}\mathbf{a}^{\mathsf{T}}}{b} &&\in \mathbb{R}^{n\times 1},\\
    \mathbf{V} & = \frac{\mathbf{Y}_b^{(i-1)}\mathbf{a}^{\mathsf{T}}}{b} &&\in \mathbb{R}^{m\times 1}\\
    \mathbf{x}^{(i)} & = \mathbf{x}^{(i-1)} +\mathbf{K}([\mathbf{y}_\mathrm{obs}]_i - y) && \in\mathbb{R}^{n\times 1}, \\
    \mathbf{Z}^{(i)} & = \mathbf{Z}^{(i-1)} - \alpha \mathbf{K} \mathbf{a} && \in\mathbb{R}^{n\times n}, \\
    \mathbf{y}^{(i)} & = \mathbf{y}^{(i-1)} + \mathbf{V}([\mathbf{y}_\mathrm{obs}]_i - y)&&\in\mathbb{R}^{m\times 1},\\
    \mathbf{Y}^{(i)} & = \mathbf{Y}^{(i-1)} - \alpha \mathbf{V} \mathbf{a}&& \in\mathbb{R}^{m\times n}.
\end{aligned}
$$

We note here that $\mathbf{V}=\mathbf{H}\mathbf{K}$, $\mathbf{y}^{(i)}=\mathbf{H}\mathbf{x}^{(i)}$ and $\mathbf{Y}^{(i)}=\mathbf{H}\mathbf{x}^{(i)}$. The pre-whitening trick (sequential scheme II) applies just as well in this scheme for the case of correlated errors. In that case, we simply change $\mathbf{Y}^{(0)}\to \sqrt{\mathbf{R}}^{-1}\mathbf{Y}^{(0)}$ and $\mathbf{y}^{(0)}\to\sqrt{\mathbf{R}}^{-1}\mathbf{y}^{(0)}$, along with 
$[\mathbf{R}]_{ii} \to 1$ and 
$\mathbf{y}\_\mathrm{obs}\to \sqrt{\mathbf{R}}^{-1}\mathbf{y}\_\mathrm{obs}$.

### Ensemble square root filter

#### Serial implementation of ensemble square root filter
It is clear how the various serial schemes generalize to the ensemble case, we will just cover serial implementation III from above, which was the variation where we do not re-apply $\mathbf{H}$ within the scheme but compute all observations ahead of computing the updates. We start with $\mathbf{X}^{(0)}{}'=\mathbf{X}_b=\mathbf{Z}_b\mathbf{G}$, and $\mathbf{x}^{(0)}{}'=\mathbf{x}_b$, along with $\mathbf{Y}^{(0)}{}'=\mathbf{H}\mathbf{X}_b$ and $\mathbf{y}^{(0)}{}'=\mathbf{H}\mathbf{x}_b$. Then, the serial scheme requires one to iterate over the scheme below for each observation:

$$
\begin{aligned}
    y & = [\mathbf{y}^{(i-1)}{}']_{i} &&\in\mathbb{R},\\
    \mathbf{a} &= [\mathbf{Y}^{(i-1)}{}']_{i} &&\in\mathbb{R}^{1\times N}, \\
    b & = \frac{1}{N}\mathbf{a}\mathbf{a}^{\mathsf{T}}+ [\mathbf{R}]_{ii} && \in\mathbb{R},\\
    \alpha & = \left(1+\sqrt{\frac{[\mathbf{R}]_{ii}}{b}}\right)^{-1} &&\in\mathbb{R}, \\
    \mathbf{K}' & = \frac{1}{N}\frac{\mathbf{X}_b^{(i-1)}{}'\mathbf{a}^{\mathsf{T}}}{b} &&\in \mathbb{R}^{n\times 1},\\
    \mathbf{V}' & = \frac{1}{N}\frac{\mathbf{Y}_b^{(i-1)}{}'\mathbf{a}^{\mathsf{T}}}{b} &&\in \mathbb{R}^{m\times 1}\\
    \mathbf{x}^{(i)}{}' & = \mathbf{x}^{(i-1)}{}' +\mathbf{K}'([\mathbf{y}_\mathrm{obs}]_i - y) && \in\mathbb{R}^{n\times 1}, \\
    \mathbf{X}^{(i)}{}' & = \mathbf{X}^{(i-1)}{}' - \alpha \mathbf{K}' \mathbf{a} && \in\mathbb{R}^{n\times M}, \\
    \mathbf{y}^{(i)}{}' & = \mathbf{y}^{(i-1)}{}' + \mathbf{V}'([\mathbf{y}_\mathrm{obs}]_i - y)&&\in\mathbb{R}^{m\times 1},\\
    \mathbf{Y}^{(i)}{}' & = \mathbf{Y}^{(i-1)}{}' - \alpha \mathbf{V}' \mathbf{a}&& \in\mathbb{R}^{m\times N}.
\end{aligned}
$$

At the final iteration, then, $\mathbf{x}^{(m)}{}'\approx \mathbf{x}_b$ and $\mathbf{X}^{(i)}{}'\mathbf{X}^{(i)}{}'{}^{\mathsf{T}}\approx \mathbf{P}_a$. The extension to the pre-whitened case is identical to what was described in the square root Kalman filter section. Note how some of the shapes of variables are now $N$ instead of $n$.

## Appendix: Square root formulation of the covariance update
The definition of the Kalman filter gives the update of the prior covariance matrix 

$$
\begin{equation}
    \mathbf{P}_a = (\mathbf{I} - \overbrace{\mathbf{P}_b\mathbf{H}^T\underbrace{(\mathbf{H}\mathbf{P}_b\mathbf{H}^T+\mathbf{R})^{-1}}_{\mathbf{D}^{-1}}}^{\mathbf{K}}\mathbf{H})\mathbf{P}_b,\tag{A1}\label{eq:firstkalman}
\end{equation}
$$

where $\mathbf{K}$ is the Kalman gain and $\mathbf{D}$ is the innovation covariance matrix. We want to find a square root decomposition of the above expression, making an ansatz

\begin{equation}
    \mathbf{P}_a = (\mathbf{I} - \mathbf{P}_b\mathbf{H}^T\mathbf{W}\mathbf{H})\mathbf{P}_b(\mathbf{I} - \mathbf{P}_b\mathbf{H}^T\mathbf{W}\mathbf{H})^T,\tag{A2}\label{eq:ansatz}
\end{equation}

for an as-of-yet unspecified $\mathbf{W}$ value. Expanding the previous expression yields

$$
\begin{align}
    \mathbf{P}_a = \left(\mathbf{I} - \mathbf{P}_b\mathbf{H}^T\left(\mathbf{W} + \mathbf{W}^T - \mathbf{W}\mathbf{H}\mathbf{P}_b\mathbf{H}^T\mathbf{W}^T\right)\mathbf{H}\right)\mathbf{P}_b.\tag{A3}\label{eq:secondkalman}
\end{align}
$$

Comparing eqs. \eqref{eq:firstkalman} and \eqref{eq:secondkalman} shows we need to solve for a $\mathbf{W}$ that satisfies

$$
\begin{align}
    \mathbf{D}^{-1} & = \mathbf{W} + \mathbf{W}^T -\mathbf{W}\mathbf{H}\mathbf{P}_b\mathbf{H}^T\mathbf{W}^T, \\
    & = \mathbf{W}\left(\mathbf{W}^{-T} + \mathbf{W}^{-1} -\left(\mathbf{D}-\mathbf{R}\right)\right)\mathbf{W}^T,\tag{A4}\label{eq:target}
\end{align}
$$

(recalling innovation covariance matrix $\mathbf{D}=\mathbf{H}\mathbf{P}_b\mathbf{H}^T+\mathbf{R}$, thus $\mathbf{D}-\mathbf{R}=\mathbf{H}\mathbf{P}_b\mathbf{H}^T$). We derive a matrix inverse identity for $\mathbf{D}^{-1}$ using its matrix square root $\mathbf{A}=\sqrt{\mathbf{D}}$,

$$
\begin{align}
    \mathbf{D}^{-1} &= (\mathbf{A}\mathbf{A}^T)^{-1}, \\
    & = \mathbf{A}^{-T} \mathbf{A}^{-1},\\
    &= \mathbf{A}^{-T}(\mathbf{A}+\mathbf{B})^{-1}\left[(\mathbf{A}+\mathbf{B})(\mathbf{A}+\mathbf{B})^T \right](\mathbf{A}+\mathbf{B})^{-T}\mathbf{A}^{-1}, \\
    & = \mathbf{A}^{-T}(\mathbf{A}+\mathbf{B})^{-1}\left[\mathbf{A}(\mathbf{A}+\mathbf{B})^T + (\mathbf{A}+\mathbf{B})\mathbf{A}^T - (\mathbf{A}\mathbf{A}^T-\mathbf{B}\mathbf{B}^T) \right](\mathbf{A}+\mathbf{B})^{-T}\mathbf{A}^{-1}, \\
    & = \mathbf{C}\left[\mathbf{C}^{-T} +\mathbf{C}^{-1} - (\mathbf{A}\mathbf{A}^T - \mathbf{B}\mathbf{B}^T)\right]\mathbf{C}^T,\tag{A5}\label{eq:targetfound}
\end{align}
$$

for $\mathbf{C}=\mathbf{A}^{-T}(\mathbf{A}+\mathbf{B})^{-1}$, using an arbitrary extra matrix $\mathbf{B}$ (where we assume the inverse $\mathbf{A}+\mathbf{B}$ exists), where it still holds that $\mathbf{A}\mathbf{A}^T=\mathbf{D}$. Comparing eqs. \eqref{eq:target} and \eqref{eq:targetfound}, we can see that for $\mathbf{B}\mathbf{B}^T=\mathbf{R}$ we obtain $\mathbf{W}$ by inspection as

$$
\begin{equation}
    \mathbf{W}=\mathbf{A}^{-T}(\mathbf{A}+\mathbf{B}) = \sqrt{\mathbf{D}}^{-T}(\sqrt{\mathbf{D}}+\sqrt{\mathbf{R}})^{-1}.
\end{equation}
$$

Substituting this form for $\mathbf{W}$ into eq. \eqref{eq:ansatz} gives us

$$
\begin{equation}
    \mathbf{P}_a = \left(\mathbf{I} - \mathbf{P}_b\mathbf{H}^T\sqrt{\mathbf{D}}^{-T}(\sqrt{\mathbf{D}}+\sqrt{\mathbf{R}})^{-1}\mathbf{H}\right)\mathbf{P}_b\left(\mathbf{I} - \mathbf{P}_b\mathbf{H}^T\sqrt{\mathbf{D}}^{-T}(\sqrt{\mathbf{D}}+\sqrt{\mathbf{R}})^{-1}\mathbf{H}\right)^T,
\end{equation}
$$

thus, the square root of the covariance update may be written as

$$
\begin{equation}
    \sqrt{\mathbf{P}_a} = \left(\mathbf{I} - \mathbf{P}_b\mathbf{H}^T\sqrt{\mathbf{D}}^{-T}(\sqrt{\mathbf{D}}+\sqrt{\mathbf{R}})^{-1}\mathbf{H}\right)\sqrt{\mathbf{P}_b}\label{eq:endofderivation}
\end{equation}
$$

(as, then $\sqrt{\mathbf{P}_a}\sqrt{\mathbf{P}_a}^T=\mathbf{P}_a$).


<a id="fn-1">[1]</a>: [↩](#fnref-1) In the literature, the factor $\mathbf{G}\mathbf{G}^{\mathsf{T}}/(N-1)$ is often used as the *unbiased* estimator of the covariance of $\mathbf{G}$, which tends to the identity matrix for large $N$ (often only when $N\gg n$). However, since $\mathbf{G}$ has zero mean (we know this *a priori*), no degree of freedom is used up to subtract a data-driven mean from our ensemble (you can see it -- at no point in our algorithm do we add or subtract anything from our ensemble), thus the correct factor is $1/N$; the literature is simply mistaken. Other literature might start off by drawing $N$ members first, subtracting their *mean* from each member to get to 'state vector deviations'; in such a case we should use the unbiased estimator of the covariance. 
Their derivation will then be identical to what we arrived at, save for the factor $(N-1)^{-1}$ rather than $N^{-1}$. 
However, as we are not doing such a statistical trick, and merely approximate $\mathbf{P}_b$ in the Kalman filter, I don't see how a factor $(N-1)^{-1}$ should show its head. See also [the Wikipedia entry on computing the sample covariance in case the population mean is known (rather than estimated), in our case known to be 0](https://en.wikipedia.org/wiki/Covariance#Calculating_the_sample_covariance). Of course, the practical difference is negligible for $N > 100$, but it seems worth pointing out.