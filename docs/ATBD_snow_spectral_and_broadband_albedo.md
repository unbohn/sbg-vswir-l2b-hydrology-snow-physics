**SBG VSWIR Snow Physics Products – Algorithm Theoretical Basis Document (ATBD)**

**N. Bohn**<sup>1</sup>

<sup>1</sup>Jet Propulsion Laboratory, California Institute of Technology

Corresponding author: Niklas Bohn (urs.n.bohn@jpl.nasa.gov)

**Key Points:**

- Collection of snow retrieval algorithms to support the SBG VSWIR snow physics suite
- Includes estimation of snow spectral and broadband albedo
- Presents candidate EMIT scenes as highlights for potential snow applications

**Version:** 1.0

**Release Date:** TBD

**DOI:** TBD

# Abstract

# Plain Language Summary

## Keywords: snow spectral albedo, snow broadband albedo

# 1 Version Description

This is Version 1.0 of the SBG VSWIR snow physics algorithms - Part II: spectral and broadband albedo.

# 2 Introduction

# 3 Context/Background

## 3.1 Historical Perspective

IS-SnARF has been developed by Painter et al. (2013).

## 3.2 Additional Information

# 4 Algorithm Description

![img.png](figs/img.png)

**Figure 1.** _Map of snow physics SRR core product algorithms._

## 4.1 Imaging Spectrometer-Snow Albedo and Radiative Forcing (IS-SnARF)

### 4.1.1 Scientific theory

Airborne and spaceborne imaging spectrometers, such as EMIT and the future SBG VSWIR instruments, do not provide the measurements to allow a direct inversion of snow spectral albedo, as they measure reflected photons on their directional optical path and not from the entire hemisphere.

Retrieved L2A surface reflectance, defined as the hemispherical-directional reflectance factor (HDRF), must therefore be converted to spectral albedo. IS-SnARF utilizes the relationship between HDRF and spectral albedo, referred to as the spectral anisotropy factor $c$, all of which are a function of wavelength $\lambda$:

$$
c_{\theta_r,\phi_0-\phi_v;r;\lambda} = \frac{\alpha(r;\lambda)}{HDRF(\theta_0,\theta_v,\phi_0-\phi_v;r;\lambda)} \tag{1}
$$ 

where $\theta$ and $\phi$ are zenith and azimuth angles, and the subscripts 0, r, and v signify incident, reflected, and view angles. $r$ is the snow grain radius and $\alpha$ is the snow spectral albedo as a function of $r$.

The scalar $c$ is a function of the local geometry between solar illumination and observation angles as well as snow grain size, and rearranging Equation (1) provides an estimate of the snow spectral albedo:

$$
\alpha(r;\lambda) = HDRF(\theta_0,\theta_v,\phi_0-\phi_v;r;\lambda) * c_{\theta_r,\phi_0-\phi_v;r;\lambda} \tag{2}
$$

IS-SnARF uses a comprehensive look-up-table (LUT) of pre-calculated $c$ scalars, modeled by running the multistream DISORT code (Stamnes et al., 1988; Painter et al., 2003), to obtain the respective anisotropy factor on runtime.

### 4.1.1.1 Scientific theory assumptions

### 4.1.2 Mathematical theory

### 4.1.3 Algorithm Input Variables

### 4.1.4 Algorithm Output Variables

### 4.2.1 Scientific theory

###

### 4.2.1.1 Scientific theory assumptions

### 4.2.2 Mathematical theory

### 4.2.2.1 Mathematical theory assumptions

### 4.2.3 Algorithm Input Variables

### 4.2.4 Algorithm Output Variables

# 5 Algorithm Usage Constraints

# 6 Performance Assessment

## 6.1 Validation Methods

## 6.2 Uncertainties

## 6.3 Validation Errors

# 7 Algorithm Implementation

## 7.1 Algorithm Availability

## 7.2 Input Data Access

Candidate EMIT snow scenes:

- EMIT\_L2A\_RFL\_001\_20230216T211823\_2304714\_007
- EMIT\_L2A\_RFL\_001\_20230220T102026\_2305107\_004
- EMIT\_L2A\_RFL\_001\_20230220T194126\_2305113\_011
- EMIT\_L2A\_RFL\_001\_20230323T144706\_2308209\_038

## 7.3 Output Data Access

## 7.4 Important Related URLs

# 8 Significance Discussion

# 9 Open Research

# 10 Acknowledgements

# 11 Contact Details

# References
