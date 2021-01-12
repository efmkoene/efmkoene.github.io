---
layout: page
title: Portfolio
subtitle: Projects I've worked on
---

## Time dispersion corrections
In the modeling of seismic waves, we typically integrate over time using a finite-difference approximation. When we use small steps in time, these approximations are more-or-less correct, but these simulations can be computationally intensive. When we use a large time-step for these finite-difference approximations, the simulation is cheaper but can also contain seeming errors, as shown on the gif below.

![timedispersion](/assets/img/time_dispersion.gif)

We have explored a way to use the large time-step for the finite-difference approximations (i.e., compute results in a cheap way) and filter out the errors using a pre- and post-propgation filter of negligible cost. We have shown that this even holds for viscoelastic wave simulations. For example, on the figure below you see in red a computationally expensive simulation (2+ hours) with a small step in time. In blue you see a computationally cheap simulation (80 seconds) using large steps in time -- but the blue solution is clearly erroneous when compared to the red solution. The dotted solutions uses the filters for the cheap simulation (coming out to 80.1 seconds, 0.1 seconds extra for the pre- and post-propagation filters), and overlaps with the correct solution excellently. Hence, we can obtain accurate results at a very low cost.

![viscousexample](/assets/img/timedispersionexample.png)
