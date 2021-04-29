# Bachelor Thesis - Anhorn Barthélémy

The goal of the project is to build a program to invert geophysical data collected in the pro-glacial domain with EM methods, in order to identify resistive layers that may correspond to permafrost or buried ice. Using MATLAB, I developed a simple program that need only a few functions, which inverts the data at each measurement point. A large weight matrix connects the different measurement points and their physical characteristics allowing to apply a horizontal and vertical regularization. To test my MATLAB code, I will quickly programme in python with the basic functions of simPEG and PyGimli some one-dimensional inversions. 

## Theory

To build our sub-surface model we need to use some tools. All the equations used in this programm has been defined in a [Geonics publication](http://www.geonics.com/pdfs/technicalnotes/tn6.pdf). In this work we are mainly talking about apparent conductivity. But the measurement of this value is indirect. We use a formula which combines the magnetic field emitted by the measuring device (known) and the magnetic field induced by the subsurface (measured). In the Geonics documentation, the apparent conductivity is defined as :

![sigma_a equation](https://latex.codecogs.com/gif.latex?%5Csigma_a%20%3D%20%5Cfrac%7B4%7D%7B%5Comega%20%5Cmu_0%20s%5E2%7D%20%5Cfrac%7BH_s%7D%7BH_p%7D)

![sigma_ a equation description](https://latex.codecogs.com/gif.latex?%5C%5C%20H_s%20%3D%20%5Ctextrm%7Bsecondary%20magnetic%20fiedl%20at%20the%20reciever%20coil%7D%20%5C%5C%20H_b%20%3D%20%5Ctextrm%7Bprimary%20magnetic%20filed%20at%20the%20reciever%20coil%7D%20%5C%5C%20%5Comega%3D%202%5Cpi%20f%20%5C%5C%20f%20%3D%20%5Ctextrm%7Bfrequency%20%5BHz%5D%7D%5C%5C%20%5Cmu_0%20%3D%20%5Ctextrm%7Bpermeability%20of%20free%20space%7D%5C%5C%20s%20%3D%20%5Ctextrm%7Bintercoil%20spacing%20%5Bm%5D%7D)

Tools : PyGimli, SimPEG, MATLAB
