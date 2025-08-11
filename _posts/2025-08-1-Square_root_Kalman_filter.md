---
layout: post
title: "The Square Root Kalman Filter"
subtitle: ""
tags: [atmospheric modeling, kalman filter]
---
Consider we have access to observations $\mathbf{y}\_\mathrm{obs} \in \mathbb{R}^{m}$ and they are due to some linear operator acting on a `hidden state' $\mathbf{x}\_\mathrm{true}\in\mathbb{R}^{n}$ with added noise $\mathbf{e}$ like
\begin{align}
    \mathbf{y}\_\mathrm{obs} = \mathbf{H}\mathbf{x}\_\mathrm{true} + \mathbf{e},
\end{align}
where $\mathbf{H} \in \mathbb{R}^{n\times m}$, and the error has zero mean $\mathbb{E}[\mathbf{e}]=0$ but a covariance structure like $\mathbb{E}[\mathbf{e}\mathbf{e}^{\mathsf{T}}]=\mathbf{R}$ which we call the \textit{observation error covariance}. For example, if $\mathbf{R}$ is diagonal, its entries simply correspond to the variance (typically denoted as $\sigma^2$) of each measurement.
