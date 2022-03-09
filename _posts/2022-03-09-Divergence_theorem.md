---
layout: post
title: "The divergence theorem (a numerical explanation)"
subtitle: ""
tags: [atmospheric modeling, python]
---

The first step you do in numerical modelling is forgetting everything you've ever learned about calculus. Well, kind of, anyhow. You replace derivatives with finite-differences, integrations with summations, and so on.
But the interesting thing is that, in some way, you are actually carrying out what these procedures are all about. The derivative is defined as $[f(x)+f(x+\Delta x)]/\Delta x$, just for $\Delta x\to 0$. In numerical modeling you cannot go to 0 (or the computer will complain!), and just keep the number small. But it is how the inventors of calculus (Newton and Leibniz) originally invented these ideas, without the concepts of 'limits' etc. available to them.

On that note, I think that the divergence theorem is a neat theorem that has a numerical proof that is much easier to understand than the mathematical derivation underlying it.

#### The fundamental theorem of calculus, from a numerical algorithm
Consider the following definition of a finite-difference derivative,
\begin{equation}
  \frac{\mathrm{D}f(x)}{\mathrm{D}x} \equiv \frac{f(x+\Delta x)-f(x-\Delta x)}{2\Delta x}.
\end{equation}
And now consider adding this derivative to another derivative right next to it,
\begin{align}
  \frac{\mathrm{D}f(x)}{\mathrm{D}x} + \frac{\mathrm{D}f(x+\Delta x)}{\mathrm{D}x} & = \frac{f(x+\Delta x)-f(x-\Delta x)}{2\Delta x} + \frac{f(x+2\Delta x)-f(x)}{2\Delta x}, \\
  & = \frac{f(x+2\Delta x)+f(x+\Delta x)-f(x)-f(x-\Delta x)}{2\Delta x}.
\end{align}
And now consider adding the derivative right next to it, too:
\begin{align}
  \frac{\mathrm{D}f(x)}{\mathrm{D}x} + \frac{\mathrm{D}f(x+\Delta x)}{\mathrm{D}x} + \frac{\mathrm{D}f(x+2\Delta x)}{\mathrm{D}x} & = \frac{f(x+\Delta x)-f(x-\Delta x)}{2\Delta x} + \frac{f(x+2\Delta x)-f(x)}{2\Delta x} + \frac{f(x+3\Delta x)-f(x+\Delta x)}{2\Delta x}, \\
  & = \frac{f(x+3\Delta x) + f(x+2\Delta x) - f(x) - f(x-\Delta x)}{2\Delta x}.
\end{align}
And the one next to it...
\begin{align}
  \frac{\mathrm{D}f(x)}{\mathrm{D}x} + \frac{\mathrm{D}f(x+\Delta x)}{\mathrm{D}x} + \frac{\mathrm{D}f(x+2\Delta x)}{\mathrm{D}x} + \frac{\mathrm{D}f(x+3\Delta x)}{\mathrm{D}x}  & = \frac{f(x+4\Delta x) + f(x+3\Delta x) - f(x) - f(x-\Delta x)}{2\Delta x}.
\end{align}
The pattern that emerges is
\begin{align}
  \sum_{i=0}^N \frac{\mathrm{D}f(x+i\Delta x)}{\mathrm{D}x} \Delta x = \frac{f(x+(N+1)\Delta x)+f(x+N\Delta x) - f(x)-f(x-\Delta x)}{2}
\end{align}
We can see that the right-hand side simply corresponds to the arithmic mean $(f(x+(N+1)\Delta x)+f(x+N\Delta x))/2)=f(x+(N+1/2)\Delta x)$. So, the pattern is equally
\begin{align}
  \sum_{i=0}^N \frac{\mathrm{D}f(x+i\Delta x)}{\mathrm{D}x} \Delta x = f\left( x+\left(N+\frac{1}{2}\right)\Delta x \right) - f\left( x-\left(\frac{1}{2}\right)\Delta x \right).
\end{align}
This result actually corresponds to one of the two fundamental theorems of calculus,
\begin{align}
  \int_{a}^{b} \frac{\mathrm{d}f(x)}{\mathrm{d}x} \mathrm{d}x = f(b)-f(a).
\end{align}
It says that the area of the graph plotted by the $f'(x)$ corresponds to the difference between the function values at the ends of the original graph of $f(x)$. In words, it makes little sense that such a connection should exist, but numerically it is fairly obvious that the [telescoping sum](https://en.wikipedia.org/wiki/Telescoping_series) (i.e., the series of terms that cancel each other due to opposing signs) gives rise to such a result!

#### The divergence theorem
Once you understand the above example for the fundamental theorem of calculus, the divergence theorem actually becomes entirely trivial to see too. There is a vector field called $\mathbf{F}=F_x\mathbf{i}+F_y\mathbf{j}$ for basis vectors $\mathbf{i}$ and $\mathbf{j}$.
The divergence is defined as
\begin{equation}
  \frac{\mathrm{Di}f(x,y)}{\mathrm{Di}x} \equiv \frac{f(x+\Delta x,y)-f(x-\Delta x,y)}{2\Delta x} + \frac{f(x,y+\Delta y)-f(x,y-\Delta y)}{2\Delta y}.
\end{equation}
If we sum this expression in the $x$ direction, we find what we already obtained for the fundamental theorem of algebra,
\begin{equation}
  \sum_{i=0}^N \frac{\mathrm{Di}f(x+i\Delta x,y)}{\mathrm{Di}x} = \frac{f\left( x+\left(N+\frac{1}{2}\right)\Delta x,y \right) - f\left( x-\left(\frac{1}{2}\right)\Delta x,y \right)}{\Delta x} - \frac{f(xx,y+\Delta y)-f(x,y-\Delta y)}{2\Delta y}.
\end{equation}
And if we do the same in the $y$ direction we obtain
\begin{equation}
  \sum_{j=0}^M\sum_{i=0}^N \frac{\mathrm{Di}f(x+i\Delta x,y+j\Delta y)}{\mathrm{Di}x} = \sum_{j=0}^M\left( \frac{f\left( x+\left(N+\frac{1}{2}\right)\Delta x,y \right) - f\left( x-\left(\frac{1}{2}\right)\Delta x,y+j\Delta y \right)}{\Delta x}\right) + \sum_{i=0}^N \left( \frac{f\left(x+i\Delta x, y+\left(M+\frac{1}{2}\right)\Delta y \right) - f\left( x+i\Delta x, y-\left(\frac{1}{2}\right)\Delta y \right)}{\Delta y} \right).
\end{equation}
Or, rearranging, we find
\begin{equation}
  \sum_{j=0}^M\sum_{i=0}^N \frac{\mathrm{Di}f(x+i\Delta x,y+j\Delta y)}{\mathrm{Di}x} \Delta x \Delta y = \sum_{j=0}^M\left( f\left( x+\left(N+\frac{1}{2}\right)\Delta x,y \right) - f\left( x-\left(\frac{1}{2}\right)\Delta x,y+j\Delta y \right)\right)\Delta y + \sum_{i=0}^N \left( f\left(x+i\Delta x, y+\left(M+\frac{1}{2}\right)\Delta y \right) - f\left( x+i\Delta x, y-\left(\frac{1}{2}\right)\Delta y \right) \right)\Delta x.
\end{equation}
This corresponds to a proof for the divergence theorem for a square area $(x,y)\in[(x_0,x_1)\times(y_0,y_1)]$,
\begin{align}
  \iint_V \nabla \cdot \mathbf{F} \mathrm{d} x \mathrm{d} y & = \int_{y_0}^{y_1} f(x_1,y)\mathrm{d}y - \int_{y_0}^{y_1} f(x_0,y)\mathrm{d}y + \int_{x_0}^{x_1} f(x,y_1) \mathrm{d}x - \int_{x_0}^{x_1} f(x,y_0) \mathrm{d}x,\\
  & = \oint_{\partial V} \mathbf{F}\cdot\mathbf{n}\mathrm{d}(\partial V).
\end{align}

Okay, the math isn't pretty. But once you realize you're dealing with a telescoping sum of finite-difference terms, the rest is conceptually straightforward!
