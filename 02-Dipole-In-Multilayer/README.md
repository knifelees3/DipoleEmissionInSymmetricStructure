# Dipole-Emission-In-Multilayer

## Motivation

To calculate the dipole's emission rate and emission pattern in multilayered structure. Since the multilayer structure has translation symmetry in two directions, the dipole's emission properties can be described theoretically. 

The original motivation is to write a program to fit the experimental data and then obtain a 2D-dipole's orientation.

I think this program can also do other things so I will continue develop this program and may use machine learning to fit the angle.

## Jupyter File List

I list the corresponding jupyer file here

[RunAllFunctions.ipynb](https://nbviewer.jupyter.org/github/knifelees3/DipoleEmissionInSymmetricStructure/blob/master/02-Dipole-In-Multilayer/PythonProgram/RunAllFunctions.ipynb)

[Verify_ChenNatPhoton.ipynb](https://nbviewer.jupyter.org/github/knifelees3/DipoleEmissionInSymmetricStructure/blob/master/02-Dipole-In-Multilayer/PythonProgram/Verify_ChenNatPhoton.ipynb)

# Update Log

## 2020 03 23 

Create the repository first and upload the program first.

## 2020 06 22

The theoretical derivations for this code is updated in my personal blog:

[Dipole's Emission In Multi-Layered Structure](https://knifelees3.github.io/2020/06/22/A_En_DipoleInMultiLayerCartesian/#Relations-of-amplitudes-in-different-layers)

## 2020 07 19

* A verification of this program with the reference is added into this folder.

The corresponding file name is 

"Verify_Chen_NatPhoton.ipynb"

You can read the Jupyter files here

[Verify_ChenNatPhoton.ipynb](https://nbviewer.jupyter.org/github/knifelees3/DipoleEmissionInSymmetricStructure/blob/master/02-Dipole-In-Multilayer/PythonProgram/Verify_ChenNatPhoton.ipynb)

* The corresponding MATLAB code is also generated

## 2020 09 20

Create the structure visualization function named `Show_Structure`. This function was written in the file "Fun_BFP_Image.py" and added in the end of class file "class_BFP_Image_QD.py". The corresponding test note book also updated.

[RunAllFunctions.ipynb](https://nbviewer.jupyter.org/github/knifelees3/DipoleEmissionInSymmetricStructure/blob/master/02-Dipole-In-Multilayer/PythonProgram/RunAllFunctions.ipynb)

You can see the updated note book here

[RunAllFunctions.ipynb](https://nbviewer.jupyter.org/github/knifelees3/DipoleEmissionInSymmetricStructure/blob/master/02-Dipole-In-Multilayer/PythonProgram/RunAllFunctions.ipynb)

Now we can see the structure via a simple plot.