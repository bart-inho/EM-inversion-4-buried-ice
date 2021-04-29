# Bachelor Thesis - Anhorn Barthélémy

The goal of the project is to build a program to invert geophysical data collected in the pro-glacial domain with EM methods, in order to identify resistive layers that may correspond to permafrost or buried ice. Using MATLAB, I developed a simple program that need only a few functions, which inverts the data at each measurement point. A large weight matrix connects the different measurement points and their physical characteristics allowing to apply a horizontal and vertical regularization. To test my MATLAB code, I will quickly programme in python with the basic functions of simPEG and PyGimli some one-dimensional inversions. 

## Theory

To build our sub-surface model we need to use some tools. All the equations used in this programm has been defined in a [Geonics publication](http://www.geonics.com/pdfs/technicalnotes/tn6.pdf). In this work we are mainly talking about apparent conductivity. But the measurement of this value is indirect. We use a formula which combines the magnetic field emitted by the measuring device (known) and the magnetic field induced by the subsurface (measured). In the Geonics documentation, the apparent conductivity is defined as :

![sigma_a equation](https://latex.codecogs.com/gif.latex?%5Csigma_a%20%3D%20%5Cfrac%7B4%7D%7B%5Comega%20%5Cmu_0%20s%5E2%7D%20%5Cfrac%7BH_s%7D%7BH_p%7D)

$$E = mc^2$$

Tools : PyGimli, SimPEG, MATLAB
