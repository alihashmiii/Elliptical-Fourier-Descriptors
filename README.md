# Elliptical-Fourier-Analysis

`@Author: Ali Hashmi`

`Date: 02-03-2017`

This script is a Wolfram Language (WL) implementation (of a Python code by Henrik Blidh) for Elliptical Fourier Descriptors for 2D shape classification. The Lobe Contribution Elliptical Fourier Analysis has been implemented in the package as well. The lobe contribution can be used to classify complex shapes with lobes, for instance cells.

please check the other branch for the python implementation of the code

Note: Lobe-Contribution EFA is based on the following paper: by Yara Sanchez-Corrales et al https://www.biorxiv.org/content/early/2017/06/30/157842


`Package usage:`


![Alt text](https://user-images.githubusercontent.com/10793580/34066620-a7260d84-e211-11e7-9972-ca6f7c0f0272.png)

![Alt text](https://user-images.githubusercontent.com/10793580/34066621-a7534768-e211-11e7-84b2-4cb1cf271915.png)

![Alt text](https://user-images.githubusercontent.com/10793580/34066622-a78f21ac-e211-11e7-8980-5eb7ae64bb88.png)


`shapes of cells and complex objects can be quantified as well`

![Alt text](https://user-images.githubusercontent.com/10793580/34066829-ec903170-e215-11e7-9386-658b86eead64.png)


**`Note: when the lobeContribution option is set to True we obtain a list with three members: the first member being the elliptical fourier coefficients, the second member is the object's outline superimposed with the contour determined from the coefficients and the third member is the information from Lobe-Contribution.
With the option turned off, we obtain a a list without Lobe-Contribution information`**
