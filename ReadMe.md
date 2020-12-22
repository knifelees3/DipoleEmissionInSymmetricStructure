# Dipole's emission in symmetric structure

[TOC]

## Introduction

This project is to calculate the emission properties when the dipole is put near symmetric structure.

Dipole is a very useful physical model which can be used in different fields. In my own research, dipole's emission in different media can be used to describe the properties of Atom-photon interactions in a way.  

Symmetry can be translation symmetry or rotation symmetry.

Dipole's emission in arbitrary media can't be described theoretically but when the structures have symmetries we can use some special functions to describe its properties.

## Different Structure

I will calculate the following case:

1. Dipole in free space .

2. Dipole near sphere (Notes can be found here [Dipole’s emission near sphere](https://knifelees3.github.io/2020/06/20/A_En_DipoleEmissionNearSphere/#Introduction)).

   ![](https://raw.githubusercontent.com/knifelees3/my_pictures/master/picgoup/20200107214942441_16601.png)

3. Dipole in multi layered structures.

   ![](https://raw.githubusercontent.com/knifelees3/my_pictures/master/picgoup/20200223215829792_1734.jpg)

1.   Dipole near waveguide.

for dipole near waveguide, only when the cross section has rotation symmetry then the mode can be described theoretically.



## Goal of this project

* Develop a full package that can simulate the dipole's emission in different structures with some symmetries. 
* The package should better based on open source package rather than business software.
* Try to make use of the power of theoretical expressions so that the simulation time would be fast.
* The derivations and use should be complete. 



## Update Log

### 2020 09 20

For the "02-Dipole-In-Multilayer"

Create the structure visualization function named `Show_Structure`. This function was written in the file "Fun_BFP_Image.py" and added in the end of class file "class_BFP_Image_QD.py". The corresponding test note book also updated.

[RunAllFunctions.ipynb](https://nbviewer.jupyter.org/github/knifelees3/DipoleEmissionInSymmetricStructure/blob/master/02-Dipole-In-Multilayer/PythonProgram/RunAllFunctions.ipynb)

You can see the updated note book here

[RunAllFunctions.ipynb](https://nbviewer.jupyter.org/github/knifelees3/DipoleEmissionInSymmetricStructure/blob/master/02-Dipole-In-Multilayer/PythonProgram/RunAllFunctions.ipynb)

Now we can see the structure via a simple plot.