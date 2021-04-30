# Bachelor Thesis - Anhorn Barthélémy

The goal of the project is to build a program to invert geophysical data collected in the pro-glacial domain with EM methods, in order to identify resistive layers that may correspond to permafrost or buried ice. Using MATLAB, I developed a simple program that need only a few functions, which inverts the data at each measurement point. A large weight matrix connects the different measurement points and their physical characteristics allowing to apply a horizontal and vertical regularization. To test my MATLAB code, I will quickly programme in python with the basic functions of simPEG and PyGimli some one-dimensional inversions. 

## Theory

### Geonics physics

To build our sub-surface model from apparents conductivity we need to use equations. All the equations used in this programm has been defined in a [Geonics publication](http://www.geonics.com/pdfs/technicalnotes/tn6.pdf). In this paper, Geonics presents a simplified physical model that generates correct approximations. However, one must be all the more attentive to the framework in which the data were collected. 

In this work we are mainly talking about apparent conductivity. But the measurement of this value is indirect. We use a formula which combines the magnetic field emitted by the measuring device (known) and the magnetic field induced by the subsurface (measured). In the Geonics documentation, the apparent conductivity is defined as :

![sigma_a equation](https://latex.codecogs.com/gif.latex?%5Csigma_a%20%3D%20%5Cfrac%7B4%7D%7B%5Comega%20%5Cmu_0%20s%5E2%7D%20%5Cfrac%7BH_s%7D%7BH_p%7D)

![sigma_ a equation description](https://latex.codecogs.com/gif.latex?%5C%5C%20H_s%20%3D%20%5Ctextrm%7Bsecondary%20magnetic%20fiedl%20at%20the%20reciever%20coil%7D%20%5C%5C%20H_b%20%3D%20%5Ctextrm%7Bprimary%20magnetic%20filed%20at%20the%20reciever%20coil%7D%20%5C%5C%20%5Comega%3D%202%5Cpi%20f%20%5C%5C%20f%20%3D%20%5Ctextrm%7Bfrequency%20%5BHz%5D%7D%5C%5C%20%5Cmu_0%20%3D%20%5Ctextrm%7Bpermeability%20of%20free%20space%7D%5C%5C%20s%20%3D%20%5Ctextrm%7Bintercoil%20spacing%20%5Bm%5D%7D)

In our case, we synthetically create our data from a *True Model* that we arbitrarily define to fit the context we want to study. To do that we use an another physical relationship that depends of the physical caracteristic of the subsurface to calculate the weight of each layer :

![weight function](https://latex.codecogs.com/gif.latex?%5C%5C%20R_V%28z%29%20%3D%20%5Cfrac%7B1%7D%7B%284z%5E2&plus;1%29%5E%7B%5Cfrac%7B1%7D%7B2%7D%7D%7D%5C%5C%20R_H%28z%29%20%3D%20%5Csqrt%7B4z%5E2&plus;1%7D%20-%202z%5C%5C%20%5C%5C%20z%20%3D%20%5Cfrac%7Bz_%7Btop%7D%7D%7Bs%7D%5C%5C%20z_%7Btop%7D%20%3D%20%5Ctextrm%7Bvertical%20coordinate%20of%20the%20top%20of%20the%20layer%7D)

Once the weights of each layer are defined, we can weight the different conductivity values of the real model such as for a three layer model :

![weight tot](https://latex.codecogs.com/gif.latex?%5Csigma_a%20%3D%20%5Csigma_1%5B1-R%28z_1%29%5D%20&plus;%20%5Csigma_2%5BR%28z_1%29-R%28z_2%29%5D%20&plus;%20%5Csigma_3%20R%28z_2%29)

### Forward model
In a real-world setting, the apparent data is obtained directly from the field. In our case, we need to create synthetic data to be inverted. We do this using the weight system defined above by Geonics. This is called the *"forward model"*. First, we need to define the parameters of our forward model. We will need the conductivity values of the layers, the layer boundary coordinates, the coil spacing, the coil orientations and the horizontal discretization.

Exemple :
```matlab
xlog = 0:0.1:20; %[m] horizontal discretization
nmeasure = length(xlog); % number of horizontal measurments
ztop = repmat([0; 1; 4; 7], 1, nmeasure); % top layer vertical coordinate
sig = repmat([20e-3; 1e-3; 20e-3; 10e-3], 1, nmeasure); % true model map
coilsep = repmat(0.1:0.1:10, nmeasure, 1)'; % coilseparations
ori = repmat([0 1], nmeasure, size(coilsep, 1)/2)'; % orientation of the dipole (0 = vertical, 1 = horizontal)
```

Then we can use the weighting formula defined by Geonics to generate the synthetic apparent conductivities. Then we can use the weighting formula defined by Geonics to generate the synthetic apparent conductivities. The data are stored in a row by row matrix containing the geometry of the model (apparent sigma, coilspacing, orientation, x-coordinate).

Exemple :
```matlab
for i = 1:nmeasure
    % generate datas in a matrix that contains physical properties
    data = [data; forwardEM2D(sig(:, i), ztop(:, i), coilsep(:, i), ori(:, i), xlog(i))];
end
```

### Inversion method



in progress...


### Focus on regularization
in progress...
## Applications

in progress...



Tools : PyGimli, SimPEG, MATLAB
