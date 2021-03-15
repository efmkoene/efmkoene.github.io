---
layout: post
title: "Five introductions to the finite-difference method"
subtitle: ""
tags: [finite-difference,FDM,algorithms]
---

In my PhD thesis, I worked on finite-difference (FD) methods that can be used in physics simulations. I was playing with the idea to include a chapter in my thesis that wraps up different theories concerning FD operators.
In the end, I decided against this idea. I couldn't decide what part was relevant and what part wasn't. In the end, I would just be wasting the reader's time.
However, I still think that some of these observations may be useful. For example, for someone learning about FD methods, and wanting to understand them more clearly. Or, for someone studying FD methods in a research setting.

In the following, I lay out five different ways of approaching, studying, conceptualizing, or understanding finite-difference operators. Typically, when learning about FD methods, you will learn just one of these possible ways. But I think that some understanding of all five ways is useful.

Note that I use "finite-difference operators" in the sense of linear combinations of multiple (sampled) values, in order to approximate the underlying function or its derivative.

# 1. The finite-difference method as approximation of a derivative
The most obvious way to approach the FD method is through an approximation of a derivative. In high-school you probably defined a derivative with some kind of 'limit' operation,

$$ f'(x) \equiv \lim_{h\to 0} \frac{f(x+h)-f(x)}{h}. $$

*(Limits have a [strict formal definition](https://en.wikipedia.org/wiki/(ε,_δ)-definition_of_limit), but the long and short of it is the simple [rise over run](https://www.onlinemath4all.com/rise-over-run-formula.html) formula.)*

Limits are typically computed with *symbolic manipulation*, i.e., with a pen and a piece of paper. We, instead, want to use a computer. They are good only at $+$, $-$, $\times$ and $\div$ operations. The logical thing, for a computer, is thus to compute the limit with direct assignment, skipping the limit operation.

$$ \frac{f(x+h)-f(x)}{h}\bigg|_{h=0} = \frac{0}{0} = \text{NaN}, $$

where $\text{NaN}$ means 'not a number', as we're doing [division by zero](https://en.wikipedia.org/wiki/Division_by_zero). Not very useful.

So, what's the closest thing we can do instead? We can define a parameterized function $f'_h$ as

$$ f_h'(x) \equiv \frac{f(x+h)-f(x)}{h} $$

for which it holds that

$$ \lim_{h\to 0} f_h'(x) = f'(x). $$

Now, if we evaluate $f_h'$ for a *very small* value of $h$, we will get a good approximation of $f'(x)$. This is what you see in the gif, below, for some arbitrary function.

![derivative approximation](../assets/img/FD_as_limit.gif)


We rarely use the FD method in isolation; typically, they're used in the setting of an equation. For example, we have some equation of the form (typically, $x$ may be replace with $t$ to describe problems that are dynamic in time)

$$ \frac{\text{d}f(x)}{\text{d}x} = g(x,f(x)) $$

Now, we write down its approximation $(f(x+h)-f(x))/h$, and manipulate it into the following form

$$ f(x+h) - f(x) = h\cdot g(x,f(x)) $$

This equation can then be solved for future values, $f(x+h) = f(x) + h\cdot g(x,f(x))$. This is a truly magical feature of the FD methods. It turns differential equations into simple $+$, $-$, $\times$ and $\div$ operations.

For large $h$, the derivative approximation is poor. For small $h$, the approximation gets really good! It is intuitively clear that for smaller and smaller (but non-zero!) values of $h$, our approximation to the derivative should get better. 
However, this introduction to FD approximations does not really help us in quantifying the errors we've made thus far, or to help us understand where the FD method can fail.

# 2. The finite-difference method as exact evaluation of integral
With the previous definition still in mind (that the FD method is an approximation), this second definition will blow your mind. Starting again with the evolution equation $f'(x)=g(x,f(x))$, we integrate from $x$ to $x+h$,

$$ \int_{x}^{x+h} f'(x) \text{d}\! x = f(x+h) - f(x) = \int_x^{x+h} g(x,f(x)) \text{d}\! x. $$

The magical thing? With the previous definition, we derived the relation $f(x+h)-f(x)=h\cdot g(x,f(x))$ in a hand-wavy way. Now, we have derived the very same left-hand side in a *precise* form, using the fundamental theorem of calculus. There is no approximation involved!

However, our issues haven't gone away. Unless we can nicely integrate the right-hand side (by knowing $G=\int g\text{d}\! x$), we have to now choose how we want to integrate $g$ over the interval. The simplest thing to do is to make the approximation that $g$ does not change much over the interval $[x,x+h]$, such that we can approximate

$$ f(x+h) - f(x) = \int_x^{x+h} g(x,f(x)) \text{d}\! x \approx g(x,f(x)) \int_x^{x+h} \text{d}\! x = g(x,f(x))\cdot h. $$

and this has brought us back right to where we were, the same FD equation as before! However, we start to get a better grip on the problem. It is clear that we use a poor approximation of $g$, and it is clear that we can use other (better) numerical integration methods. This is where [Runge-Kutta methods](https://en.wikipedia.org/wiki/Runge–Kutta_methods) and other techniques come into form, which are essentially ways to better appproximate the integral of $g$ by using intermediate evaluations (or approximations) of $g$ in the interval $[x,x+h]$. 

Hence, considering the FD method not in isolation, but in the context of a differential equation, we could give an *exact* FD step. However, we're left with another integral which must still be found numerically.

# 3. The finite-difference method as polynomial interpolation
