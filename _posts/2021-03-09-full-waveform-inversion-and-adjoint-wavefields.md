---
layout: post
title: "Full waveform inversion and adjoint wavefields"
subtitle: ""
tags: [FWI,geophysics,adjoint methods]
---

# Reverse-time imaging
Being an exploration geophysicist by training, I am quite familiar with the ideas underlying reverse-time migration. It's a really rather simple method. We can say that, roughly speaking, a wavefield is composed of two parts: a source wavefield (which travels into the earth and is transmitted on each rock interface) and a receiver wavefield (which is generated at places where the source wavefield hits the interface between two rocks and which creates a reflection off those interfaces back up towards the surface). If you can get your hands on those two, separate, wavefields, then it is a simple observation to make that only at those locations where the source & reflected (or 'receiver') wavefield coincide in space and time, a rock interface was present. How do you get those two wavefields? Well, simply by modeling them in a smooth background model!

We can understand this easily when we take it slowly. Consider a pressure wavefield in a simple model, where receivers are placed at (or close to) the surface. This shows you the idea of seismics in a simple way: waves propagate outwards from a source, and they interact with the different rock interfaces (in two different shades of brown).

![easywavefield](/assets/img/easywavefield.gif)



$$ x = y^2 $$

