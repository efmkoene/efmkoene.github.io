---
layout: post
title: "Implementing the Cagniard-de Hoop Method for anisotropic 2-D media"
subtitle: ""
tags: [geophysics, Julia]
---

I need an exercise to learn the programming language [Julia](https://juliacomputing.com/products/juliapro). I have used in in the past as a toy program to implement a finite-difference solver for seismic waves; but despite the reputation of the language (somewhat along the lines of *as easy as MATLAB, as fast as C*), I then found it to be rather hard, and not so straightforward to port my codes from MATLAB to Julia. It requires quite a bit more work to initialize appropriate dataformats, and I find that it's not as easy to slice data in Julia as it is in MATLAB.

Anyhow, I decided to implement the [Cagniard-de Hoop method following the theory of Jos van der Hijden](http://resolver.tudelft.nl/uuid:9923aa31-eb5b-4c35-8acc-078f53fa01fd). I have done this before in MATLAB, but then implemented the method in the general 3-D sense, to then only use it on 2-D wave propagation anyhow. The MATLAB code is alright, but can require up to 7 hours to compute appropriately accurate solutions, which is of course ridiculously long. I want to get this computational time down significantly, but the 3-D nature of the code actually prevents me from doing this (this has to do with the fact that it involves complex analysis along different branches, and the selection procedure to catch the appropriate branch in 3-D is prone to failures). I hope that a 2-D code can actually simplify matters a little bit.

On this blog, I'll try and catalog my code as I write it. Any reader should refer to the dissertation of Van der Hijden for more details.

#### Bond transform matrix
This matrix allows us to rotate a stiffness matrix. Typically, the Bond transform is defined for the (6x6) stiffness *matrix* that applies in 3-D media. However, I simplify it for a (3x3) stiffness matrix, as applicable in 2-D media. It takes an input variable `d`, which is the angle (in degrees) with which the medium is rotated in the counterclockwise direction, assuming x pointing to the right and z pointing up. *(Note that we typically assume z pointing down, in which case this transform thus results in a clockwise direction of the medium.)*

<!-- {% raw %} -->
```Julia
# --- Bond transform matrix
function BondMatrix2D(d)
    b = [cosd(d)     cosd(90 - d);
         cosd(90 + d) cosd(d)     ];
# --- Let's be clever about this. Submatrix (left-top) is just the b^2
    M1 = b.^2;
# --- Each column is the product of the 'other 2' colums of b
    M2 = 2 * b[:, 2] .* b[:, 1];
# --- These are just determinants (but with + rather than - the crossed term)
    M4 = b[1,1]*b[2,2] + b[1,2]*b[2,1]
# --- Output the 2D Bond transform matrix
    M = [M1    M2;
        -M2'/2 M4]
    return M
end
```
<!-- {% endraw %} -->

#### System matrix 'A' and source vector 'F'
Matrix 'A' contains the primary information about a medium; it must be constructed for each medium within our domain. The source vector follows similarly, and is computed here as the building blocks based on the matrix C are already present. Simplified for 2-D media (where C is 3x3, rather than 6x6), it may be computed in the following way, where I have explicitly written out all matrix-matrix operations, in the hope that this is faster and retains the symmetry of the system better... I used *Wolfram Mathematica* to simplify the 2-D expressions,

<!-- {% raw %} -->
```Mathematica
Inverse[{{C55, C53}, {C53, C33}}].Transpose[{{C15, C13}, {C55, C53}}] // MatrixForm // FullSimplify
{{C11, C15}, {C15,C55}} - {{C15, C13}, {C55, C53}}.Inverse[{{C55, C53}, {C53, C33}}].Transpose[{{C15, C13}, {C55, C53}}] // MatrixForm // FullSimplify
```
<!-- {% endraw %} -->


and then the matrices are computed as follows.

<!-- {% raw %} -->
```Julia
function AFmatrix(C,s1,rho,f,h)
    DD = (C[2,3]^2-C[2,2]*C[3,3]);
    CzziCzx = [(C[1,3]*C[2,2]-C[1,2]*C[2,3])/DD 1;
               (C[1,3]*C[2,3]-C[1,2]*C[3,3])/DD 0];
    CxxmCxzCzziCzx = [(C[1,3]^2*C[2,2]-2*C[1,2]*C[1,3]*C[2,3]+C[1,1]*C[2,3]^2+C[1,2]^2*C[3,3]-C[1,1]*C[2,2]*C[3,3])/DD 0;
                       0                    0];
    A1 = -s1 * CzziCzx;
    A2 = inv([C[3,3] C[2,3];
              C[2,3] C[2,2]]);
    A3 = rho*I - s1^2*CxxmCxzCzziCzx;
    A = [ A1 A2           ;
          A3 transpose(A1) ];
    F = [ CzziCzx*h[:,1]+h[:,2] ; f + s1 * CxxmCxzCzziCzx * h[:,1] ];
    return A,F
end
```
<!-- {% endraw %} -->

#### Eigenvalue decomposition
Next, we require the (repeated) eigenvalue decomposition of matrix `A`. Ideally, I can plot the results. This took me quite a while in Julia, as the plotting parameters are not in the slightest as comfortable as those in MATLAB. Eventually, I settled on using the PyPlot back-end. In this way, we can plot the slowness surface of a given stiffness matrix. For example, the following code that puts the ideas from above into practice may be run (ending on `gcf()`, to display the result within `juno`; the `close()` command that follows it is required to restart drawing, as the `plot()` command just writes on top of the previous result.).

```
M = BondMatrix2D(22.5)
C0 = [64 18 -4;
     18 46 -2;
     -4 -2 7 ]*3.1e8;
C = M * C0 * M'
rho = 2000;
s1=0;
f = [1;0]
h = [0 0; 0 0];
# === Get and sort D matrix (0th iteration, with s1=0)
using PyPlot
for s in 0:100
    s1  = s/100000;
    A,F = AFmatrix(C,s1,rho,f,h);
    L,D = eigen(A);
    plot([s1;s1;s1;s1],real(L),marker=".",linestyle="none",color="black")
end
gcf()
savefig("slowness_pyplot.png") # Saves the CURRENT_PLOT as a .png
close()
```

The resulting figure is as shown below, which is the correct solution (the horizontal axis of the plot corresponds to the slowness along the "X" axis; the vertical axis corresponds to the slowness along the "Z" axis.

![slownesses](/assets/img/slowness_pyplot.png)
