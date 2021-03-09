---
layout: post
title: "Full waveform inversion and adjoint wavefields"
subtitle: ""
tags: [FWI,geophysics,adjoint methods]
---

# Reverse-time imaging
Being an exploration geophysicist by training, I am quite familiar with the ideas underlying reverse-time migration. It's a really rather simple method. 

## A simple modeling example
Consider a pressure wavefield in a simple model, where receivers are placed at (or close to) the surface. This shows you the idea of seismics in a simple way: waves propagate outwards from a source, and they interact with the different rock interfaces (in two different shades of brown). And in geophysics, we have access to the recordings of the wavefield over time. This is the image that you see at the top. It's as if you're dragging a piece of paper upwards while recording the properties of the wavefield at the surface (positive or negative pressure fluctuations, or particle displacements, or ...).

![easywavefield](/assets/img/easywavefield.gif)

In reality, wavefields can be quite a bit more complicated. A famous model for seismic testing is the Marmousi model, and here we see a wavefield propagating inside the Marmousi model. Every shade of brown corresponds to a rock with slightly different properties in terms of the compressional wave velocity (i.e., how fast a wave propagates through that medium in units of "kilometers per hour")! Again, in seismology you only have access to the wavefield over time, as measured at the recording states, thus: the gray part of the image on the top is all we have.

![easymarmousi](/assets/img/easymarmousi.gif)

## Splitting the wavefield into two parts
Now, we're going to make a very obvious and very fundamental observation. We can say that, roughly speaking, a wavefield is composed of two parts: a source wavefield (which travels into the earth and is transmitted on each rock interface) and a receiver wavefield (which is generated at places where the source wavefield hits the interface between two rocks and which creates a reflection off those interfaces back up towards the surface). If you can get your hands on those two separate portions, then it is a simple observation to make that only at those locations where the source & reflected (or 'receiver') wavefield coincide *in space and time*, a rock interface was present. In the image below, each green line corresponds to a a wavefront at a given time.
If you focus in particular on the (fictitious) event at 1 seconds, you see that the source and receiver wavefield only overlap at two points; the location of the rock interface.

![sourcereceiverwavefields](/assets/img/sourcereceiverwavefields.png)

It's a simple observation to make that, if you would point-wise multiply snapshots of the source and receiver wavefield for each time-step, you would only have non-zero entries at the location of the interfaces. This is what we will do for reverse time migration.

## Creating the separated wavefields
How do you get those two wavefields separately, as we only have recorded data? Well, simply by modeling them in a smooth background model! Such a smooth background model may typically be available from other geophysical methods (such as traveltime tomography, full-waveform inversion, ...). The use of a smooth model is absolutely essential, it means that as we model our source wavefield, it creates no additional reflections, and simply "purely" propagates into our domain.

![forwardmarmousi](/assets/img/forwardmarmousi.gif)

And then what about the reflected/receiver wavefield? Well, this one we model by using the actually recorded data in the following way: we model the wave equation in *reverse time*, and use our recordings as *sources*, in the way that we inject the recorded trace as a source into our model.

![reflectedmarmousi](/assets/img/reflectedmarmousi.gif)

(this image is in gray-scale, such that you can still see the faint events on the recorded data, propagating into the model. I chose to make the gifs I use in these examples relatively small, so the quality isn't great, but I want to keep the page loading time under control).

The creation of a reverse-time wave simulation is entirely straightforward, because the wave equation itself is fully time-reversible! Freely speaking, under zero initial conditions, that means that a solution satisfying the forward time wave equation also satisfies the reverse-time wave equation. All we have to do to create the reverse-time simulation, thus, is to appropriately inject our recorded wavefield as sources into the secondary simulation.

### Time-reversing the source wavefield simulation
We require one additional step. Remember that I said that only where wavefields coincide in space *and* time, we have interfaces present? Well, right now, we have access to one wavefield simulation running in forward time, and another one running in reverse time. It's hard to multiply snapshots at identical times. We do this either by storing all snapshots of the forward wave simulation in forward time (which can get quite expensive in terms of storage and I/O operations)...or by recomputing the wavefield in reverse-time. Similar to how we recompute the reflected wavefield in reverse-time, we can place virtual recording stations all around our model, and reconstruct the source wavefield in reverse time!

![reversemarmousi](/assets/img/reversemarmousi.gif)

## Creating the reverse time migrated image
Now, we have access to the source and receiver wavefield in reverse time. All we have to do, is to multiply each snapshot. In a language like MATLAB, we have a for-loop running from `nt` back to time `1`, and we simply compute `image = image + source .* reflected`, assuming `source` and `reflected` correspond to the current snapshots of the wavefield. If we do that, we end up with a simple reverse time migrated image!

![RTMimage](/assets/img/RTMimage.gif)

If you repeat this for a lot of shots, the noise will stack out, and you'll be left with a good image. Other techniques are possible too, to further clean up the image. If you look at the final RTM image, for a large number of shots, you'll see that it reflects the structural qualities of the original Marmousi model quite well! Just from running two wavefield simulations in parallel (and running three wavefield simulations in total), we can see a lot of details popping up! Quite a magical procedure.

![RTMimage](/assets/img/RTM_final.gif)



# Full waveform inversion
In full waveform inversion (FWI), our goal is more ambitious. The idea is to find the velocity model that can explain our data. It turns out to work in a way that is extremely similar to a reverse time migration. A good way to get your hands on it, is by using the [FWI lab](https://csim.kaust.edu.sa/files/ErSE328.2013/LAB/Chapter.FWIa/index.html) from KAUST, which runs the following FWI example in MATLAB. It's missing one file, which you must create after downloading the MATLAB files, which takes care of writing out the results to a binary short format.
```
function [] = write_bin(destination,seis)
fileID = fopen(destination,'w');
fwrite(fileID,seis,'float');
fclose(fileID);
end
```

## The theory of full waveform inversion
The theory behind FWI is relatively complicated. Typically, the introductions take you pretty deep into the realms of Lagrangians, operator theory, Born approximations, and all that. I think that for a first introduction, this complexity is not actually needed, and we can get away with a LOT less mathematics. This is what I intend to do here in this blog post: de-mystify FWI for a first-time user.

### The adjoint trick
A large part of the efficiency of FWI is based on the following linear algebra 'trick': $\langle\mathbf{y},A\mathbf{x}\rangle \equiv \mathbf{y}^T(A\mathbf{x})=(\mathbf{y}A)\mathbf{x}=(A^T\mathbf{y})^T\mathbf{x}\equiv \langle A^T\mathbf{y},\mathbf{x} \rangle$.








$$ x = y^2 $$

