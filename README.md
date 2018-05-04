# ADUPROP (Algorithmic Differentiation for Uncertainty Propagation)

## Overview

ADUPROP is a scalable implementation for uncertainty propagation using the
method of moments. The Taylor coefficients for the moments first (mean) and second
moment (covariance matrix) require the efficient computation of derivatives.
ADUPROP leverages the technique of Automatic Differentiation (AD) to generate
the implementation of the derivative computations. 

ADUPROP solves a differential equation

x' = f(x),

where f(x) is a nonlinear function and the flow phi(x) of the the differential
equation is obtained using the usual discretization techniques:

x_k+1 = phi(x_k)

The user has only to provide f(x) and its Jacobian. For a given input mean and
covariance matrix ADUPROP propagates this uncertainty information according to
the method of moments.

## Installation

* Install AD tool CoDiPack. Technically, ADUPROP supports any AD tool.
* Copy Makefile.sample to Makefile.inc and edit accordingly.
* Run 'make lib' to build the library
* Run 'make test' to verify that everything correctly and no error is returned.
* Look for examples in the examples folder and get inspired.


