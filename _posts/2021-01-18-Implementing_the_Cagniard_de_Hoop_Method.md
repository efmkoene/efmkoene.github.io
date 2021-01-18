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
This matrix allows us to rotate a stiffness matrix. Typically, the Bond transform is defined for the (6x6) stiffness *matrix* that applies in 3-D media. However, I simplify it for a (3x3) stiffness matrix, as applicable in 2-D media. It takes an input variable `d`, which is the angle (in degrees) with which the medium is rotated in the clockwise direction (!).

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

#### System matrix 'A'
This matrix contains the primary information about a medium (per medium). Simplified for 2-D media, it may be computed in the following way. I'm thinking about whether simplifications may be possible here -- I know for certain that some 

```Julia
function Amatrix(C,s1,rho)
    Cxx = [C[1,1] C[1,3];
           C[1,3] C[3,3]];
    Czz = [C[3,3] C[2,3];
           C[2,3] C[2,2]];
    Cxz = [C[1,3] C[1,2];
           C[3,3] C[2,3]];
    Czzi = inv(Czz);
    CxzCzziCzx = Cxz * inv(Czz) * (Cxz');
    CxzCzziCzx = (CxzCzziCzx+CxzCzziCzx')/2; # ---> Ensure symmetry is preserved
    A1 = -s1 * CxzCzziCzx;
    A2 = Czzi;
    A3 = rho*I -s1^2* ([C[1,1] C[1,3];C[1,3] C[3,3]] - CxzCzziCzx)
    A = [ A1 A2           ;
          A3 tranpose(A1) ]; # ---> There may be room for improvement here; e.g., A3 is quite simple!
    return A
end
```
