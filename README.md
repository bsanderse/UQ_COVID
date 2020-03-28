# UQ_COVID
Uncertainty quantification (UQ) on COVID ODE model

This Matlab code provides the possibility to study the effect of parametric uncertainty in a simple SEIR-model that are typically used to study the spread of a virus outbreak like COVID-19.
The SEIR model is taken from the website of Peter Forsyth (https://cs.uwaterloo.ca/~paforsyt/SEIR.html).
The UQ framework is based on UQLab (https://www.uqlab.com/).

In order to run the code, you need:
- install UQLab from https://www.uqlab.com/
- execute UQ_corona.m

The type of UQ analysis can be changed in UQ_corona.m. In this file one should also indicate the case-file to be loaded (e.g. Wuhan.m). The uncertain parameters and their distributions, which are case specific, are specified in this case-file. The assumed distributions should follow the UQLab specification.
