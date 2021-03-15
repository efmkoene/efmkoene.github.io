---
layout: post
title: "Four introductions to the finite-difference method"
subtitle: ""
tags: [finite-difference,FDM,algorithms]
---

In my PhD thesis, I worked on finite-difference (FD) methods that can be used in physics simulations. I was playing with the idea to include a chapter in my thesis that wraps up different theories concerning FD operators.
In the end, I decided against this idea. I couldn't decide what part was relevant and what part wasn't. In the end, I would just be wasting the reader's time.
However, I still think that some of these observations may be useful. For example, for someone learning about FD methods, and wanting to understand them more clearly. Or, for someone studying FD methods in a research setting.

In the following, I lay out four different ways of approaching, studying, conceptualizing, or understanding finite-difference operators. Typically, when learning about FD methods, you will learn just one of these possible ways. But I think that some understanding of all four ways is useful.

Note that I use "finite-difference operators" in the sense of linear combinations of multiple (sampled) values, in order to approximate the underlying function or its derivative.

# 1. The finite-difference method as approximation of a derivative
The most obvious way to approach the FD method is through an approximation of a derivative. In high-school or university, you probably learned that a derivative is defined through a limiting procedure,

$$ f'(x) = \equiv \lim_{h\to 0} \frac{f(x+h)-f(x)}{h}. $$

*(In mathematics, the limit has a [long history of being sufficiently formalized](https://en.wikipedia.org/wiki/(ε,_δ)-definition_of_limit), but the long and short of it is the simple [rise over run](https://www.onlinemath4all.com/rise-over-run-formula.html).
If you could graph the function, it would be the vertical difference over the horizontal difference $\updownarrow/\leftrightarrow$.)*

To compute these limits, we typically need *symbolic manipulation*, i.e., pen and paper. The computers that we consider are typically not so good at pen and paper operations, and prefer to spend their time doing $+$, $-$, $\times$ and $\div$ operations.
Thus, a computer may try to compute the derivative operation by trying to replace the limit operation with direct assignment,

$$ \frac{f(x+h)-f(x)}{h}\bigg|_{h=0} = \frac{0}{0} = \text{NaN}, $$

where $\text{NaN}$ means 'not a number'. Okay, that's not very useful.

So, what's the closest thing we can do? We can define a parameterized function $f'_h$ as

$$ f_h'(x) = \frac{f(x+h)-f(x)}{h} $$

for which it holds that

$$ \lim_{h\to 0} f_h'(x) = f'(x). $$

Now, if we evaluate $f_h'$ for a *very small* value of $h$, we will get a good approximation of $f'(x)$. 

![derivative approximation](./assets/img/FD_as_limit.gif)


# 2. The finite-difference method as exact evaluation of integral

# 3. The finite-difference method as polynomial interpolation
