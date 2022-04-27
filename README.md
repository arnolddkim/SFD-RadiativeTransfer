# SFD-RadiativeTransfer
Matlab codes to solve the RTE in the spatial frequency domain. This code is used to generate the numerical results contained in the manuscript entitled, "Radiance backscattered by a strongly scattering medium in the high spatial frequency limit" by Boaz Ilan, Arnold D. Kim, and Vasan Venugopalan.

The main driver is BackscatteredRadianceSFRTE.m. This code solves the 3D radiative transfer equation (RTE) in the spatial frequency domain for a fixed spatial frequency.

A spatially modulated plane wave incident obliquely on a half space composed of a uniform absorbing and scattering medium. Scattering is governed by the Henyey-Greenstein scattering model. Therefore, the optical properties of the medium are the albedo, anisotropy factor, and relative refractive index. The user must also specify the angle of incidence (on the xz-plane) and the spatial frequency.

This code implements the discrete ordinate method (DOM) using a product quadrature rule. This product quadrature rule makes use of the product Gauss quadrature rule by Atkinson (Atkinson, ANZIAM J. 23, 332–347 (1982)) and the spherical harmonics.

To interpolate the radiance with respect to angles, this code implements a source integration interpolation method contained in the file SourceIntegrationMethod.m.

NOTE: The source is prescribed just inside the medium, so it does not take into account specular reflection, nor the transmission coefficient from outside of the medium to inside of it.
