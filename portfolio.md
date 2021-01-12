---
layout: page
title: Portfolio
subtitle: Projects I've worked on
---

## Time dispersion corrections
In the modeling of seismic waves, we typically integrate over time using a finite-difference approximation. When we use small steps in time, these approximations are more-or-less correct, but these simulations can be computationally intensive. When we use a large time-step for these finite-difference approximations, the simulation is cheaper but can also contain seeming errors, as shown on the gif below.

![timedispersion](/assets/img/timedispersion.gif)

We have explored a way to use the large time-step for the finite-difference approximations (i.e., compute results in a cheap way) and filter out the errors using a pre- and post-propgation filter of negligible cost. We have shown that this even holds for viscoelastic wave simulations.

![viscousexample](/assets/img/timedispersionexample.png)
