**SBG VSWIR Snow Physics Products – Algorithm Theoretical Basis Document (ATBD)**

**N. Bohn**** 1**

1Jet Propulsion Laboratory, California Institute of Technology

Corresponding author: Niklas Bohn (urs.n.bohn@jpl.nasa.gov)

**Key Points:**

- Collection of snow retrieval algorithms to support the SBG VSWIR snow physics suite
- Includes snow fractional cover and grain size estimation plus calculation of snow albedo
- Presents candidate EMIT scenes as highlights for potential snow applications

**Version:** 1.0

**Release Date:** TBD

**DOI:** TBD

[**Comment to Template Users:** All elements above should be the only content on the title page. In the remainder of the document, blue text should be replaced with content. You can structure subheadings as needed to best present your information. Do not change section order or section titles. When your ATBD is completed, remove any remaining blue content (such as this comment). We recommend APA style. If you are unsure of any formatting, it is best to check the [AGU Grammar and Style](https://www.agu.org/Publish-with-AGU/Publish/Author-Resources/Grammar-Style-Guide) guide or for detailed information see the Publication Manual of the American Psychological Association, Sixth Edition. If submitting to a journal use 1.5 spacing. If submitting to APT for distribution, use 1.15 spacing for sentences. This template has been set at 1.5 spacing.]

# Abstract

The abstract is limited to 250 words or less. The abstract should state the general question or problem that is being presented, provide background on the question or problem, briefly describe analyses, and describe the key results.

# Plain Language Summary

The plain language summary should be written for a broad audience that includes journalists and the scientific-interested public. It should state the general problem, describe what research was conducted, the result, and why the findings are important. For tips, learn how to write a [plain language summary](https://www.agu.org/Share-and-Advocate/Share/Community/Plain-language-summary), and check out [this good example](https://agupubs.onlinelibrary.wiley.com/doi/pdfdirect/10.1029/2020MS002301).

## Keywords: snow fractional cover, snow grain size, snow albedo

# 1 Version Description

This is Version 1.0 of the SBG VSWIR snow physics algorithms.

# 2 Introduction

The main text should start with an introduction. Provide a brief description of the ATBD which concisely gives users the information needed to understand the relevance and usefulness of the ATBD.

APT ATBDs have standardized sections and subsections (denoted in this template), but user defined subsections can be added where appropriate. Subsections should be numbered 1.1, 1.2; 1.1.1, 1.2.1, and so on.

# 3 Context/Background

## 3.1 Historical Perspective

MEMSCAG has been developed by Painter et al. (2003) and is an extension of the MESMA algorithm (Roberts et al., 1998).

Provide a brief description of the algorithm development history or lineage. Include a description of how the algorithm came to be in its current form. You may use a table to identify changes between versions, if appropriate.

## 3.2 Additional Information

Provide additional information such as statistical model description or other preparatory work that led to the algorithm described in this document, as needed.

# 4 Algorithm Description

![](RackMultipart20240207-1-mu1w4m_html_d7aa822bba02b1d5.png)

**Figure 1.** _Map of snow physics SRR core product algorithms._

## 4.1 Multiple Endmember Snow Cover and Grain Size (MEMSCAG)

### 4.1.1 Scientific theory

MEMSCAG features a joint estimation of snow grain size and fractional cover by coupling a spectral unmixing approach with a radiative transfer model. Both the number of endmembers and the endmembers themselves are flexible on a per-pixel basis to address spatial heterogeneity. In its default version, MEMSCAG uses a spectral library containing three types of surface endmembers: snow, vegetation, and rock. These can be extended, if desired, by soil and lake ice spectral endmembers.

Snow endmembers are simulated by combining Mie scattering and the discrete-ordinates radiative transfer model (DISORT) (Stamnes et al., 1988) for grain radii of 10 – 1100 µm, with steps of 10 µm. The simulations include variations with respect to differing solar geometry and diffuse and direct components of irradiance, and represent the hemispherical-directional reflectance factor (HDRF):

,(1)

where and are zenith and azimuth angles, and the subscripts 0 and r signify incident and reflected. is reflected radiance, is the direct, and is the diffuse irradiance illuminating the surface. Endmembers for all other surface types are derived from ASD spectral measurements in the field, which are then transformed into HDRF using Equation 1.

MEMSCAG analyzes linear spectral mixtures for all possible combinations of two or more endmembers by fitting a set of linear equations to the HDRF measured by the instrument. The linear spectral mixture model is expressed as:

,(2)

where is the measured HDRF, is the fraction of endmember _i_, is the HDRF of endmember _i_, and is the residual error at wavelength . The system of equations is then solved by modified Gram-Schmidt orthogonalization (see Section 4.1.2). The residual error is expressed accordingly:

.(3)

As goodness-of-fit criterion, MEMSCAG uses the root mean squared error (RMSE) as suggested by Painter et al. (1998) and Roberts et al. (1998):

,(4)

where _M_ is the number of instrument bands. As final step, MEMSCAG normalizes the estimated snow fractional cover by the additive complement of the shade fraction to account for topographic effects on irradiance:

.(5)

For the selection of valid mixture models, MEMSCAG applies specific constraints:

1. No more than one endmember from a surface class is present.
2. The spectral fractions sum to 1.0.
3. Spectral fractions are in the range [-0.01, 1.01].
4. Overall RMSE \< 2.5%.
5. No seven consecutive residuals (i.e., for consecutive wavelengths) exceed 2.5%.

From all _n_-endmember collections that meet the constraints, MEMSCAG selects the models with the lowest RMSE and then assigns snow fraction and grain size from the model with the fewest endmembers.

### 4.1.1.1 Scientific theory assumptions

MEMSCAG applies the following assumptions and limitations:

1. At-sensor radiance is a linear combination of radiances reflected from individual surfaces, i.e., linear spectral unmixing is a valid approach.
2. is the additive complement to the sum of the physical spectral fractions (i.e., )
3. The surface is flat and horizontal, i.e., Lambertian.
4. Variability in HDRF for a given solar geometry is negligible within the range of possible reflectance angles that occur under realistic acquisition conditions.
5. Light-absorbing particles (LAP) and liquid water in snow as well as thin snow do not impact the retrieval of snow fraction and grain size.
6. Effects of sub-pixel surface roughness on the HDRF of snow are negligible.
7. Vegetation canopy is snow-free.

### 4.1.2 Mathematical theory

In the following, all uppercase letters represent matrices, while lowercase letters typify vectors. MEMSCAG solves the set of linear equations in the mixture model by modified Gram-Schmidt orthogonalization. First, it finds the orthogonal and upper-triangular matrices to the transpose of the spectral library dataset matrix by applying QR-factorization, such that

,(6)

where . The orthogonal matrix _Q_ is then used to compute the regression coefficients _b_ based on a given measurement _y_:

,(7)

where _y_ = . Finally, the endmember fractional cover values _p_ are calculated using both _Q_ and the upper-triangular matrix _R_:

.(8)

### 4.1.3 Algorithm Input Variables

**Table 2.** _Input variables to MEMSCAG_

| Name | Description | Unit | Required |
| --- | --- | --- | --- |
| HDRF | hemispherical-directional reflectance factor per wavelength | unitless | true |
| Endmember Spectral Library | collection of spectral surface endmembers (HDRF per wavelength) | unitless | true |
| Observation Geometry | solar zenith angle, view zenith angle, relative azimuth angle | degree | false |
| Solar Irradiance | direct and diffuse downwelling solar irradiance per wavelength |
 | false |

### 4.1.4 Algorithm Output Variables

**Table 2.** _Output variables from MEMSCAG_

| Name | Description | Unit |
| --- | --- | --- |
| Snow Fractional Cover | sub-pixel snow fractional cover | % |
| Snow Grain Size | effective radius of snow grains |
 |
| Spectral Residuals | deviation between measured and modeled HDRF per wavelength | % |
| RMSE | root mean square error for the spectral residuals | % |

**4.2**  **Nolin/Dozier Model**

4.2.1 Scientific theory

The Nolin/Dozier model is an alternative to MEMSCAG in terms of estimating snow grain size. However, it does not include the calculation of snow fractional cover. The method uses the scaled ice absorption band area around 1030 nm and relates it to snow grain radius obtained from a look-up-table (LUT). The scaled absorption band area is a dimensionless quantity and calculated by integrating the scaled absorption band depth over the wavelengths of the ice absorption feature:

,(9)

### where is the continuum reflectance, and is the measured reflectance spectrum. The continuum end points are determined by averaging pairs of reflectance values at both 950 and 960 nm, and 1080 and 1090 nm. The integral of Equation 9 is then calculated using the trapezoidal rule.

The LUT is calculated by using the same models and input variables as for the MEMSCAG algorithm: a combination of Mie scattering and DISORT for varying grain radii (50 – 1000 µm), solar illumination angles, as well as direct and diffuse irradiances.

###

### 4.2.1.1 Scientific theory assumptions

### 4.2.2 Mathematical theory

### 4.2.2.1 Mathematical theory assumptions

### 4.2.3 Algorithm Input Variables

### 4.2.4 Algorithm Output Variables

**4.3** **Imaging Spectrometer-Snow Albedo and Radiative Forcing (IS-SnARF)**

### 4.3.1 Scientific theory

### 4.3.1.1 Scientific theory assumptions

### 4.3.2 Mathematical theory

### 4.3.2.1 Mathematical theory assumptions

### 4.3.3 Algorithm Input Variables

### 4.3.4 Algorithm Output Variables

## 4.1 Scientific Theory

Provide a description of the scientific principles used in product retrieval including the physics and associated geophysical phenomena observed. Use subheadings as needed to organize the information.

Note: Equations must be centered and numbered.

(1)

(2)

### 4.1.2 Scientific theory assumptions

Describe the assumptions or limitations inherent in this algorithm. Any tables must be numbered in chronological order, must include a descriptive caption centered above the table, and should be referenced in the document text.

**Table 1.** _Provide a caption that describes the information displayed in the table_

| Header 1 | Header 2 |
| --- | --- |
| Content 1.1 | Content 2.1 |
| Content 1.2 | Content 2.2 |
| Content 1.3 | Content 2.3 |

## 4.2 Mathematical Theory

Provide a description of the mathematical theory essential to the algorithm development. Figures must be numbered in chronological order, include a descriptive caption centered below the figure, and be referenced in the document text.

![](RackMultipart20240207-1-mu1w4m_html_8d28932808a4c31f.png)

**Figure 1.** _Provide a caption that describes the figure_

### 4.2.1 Mathematical theory assumptions

Provide a description of any mathematical assumptions, simplifications and approximations made when deriving the algorithm. All mathematical theories have built in assumptions which need to be communicated to users.

## 4.3 Algorithm Input Variables

The data that is used as input values must be listed in this section. Data sources are provided in Section 8.

**Table 2.** _Provide an appropriate input variable table title caption_

| Name | Long name | Unit |
| --- | --- | --- |
| In Var 1 | Input Variable 1 | e.g., m2 |
| In Var 2 | Input Variable 2 | e.g., s-2 |

## 4.4 Algorithm Output Variables

The data resulting from running the algorithm. These would be variables found in the data product released to the public.

**Table 3.** _Provide an appropriate output variable table title caption_

| Name | Long name | Unit |
| --- | --- | --- |
| Out Var 1 | Output Variable 1 | e.g., km2 |
| Out Var 2 | Output Variable 2 | e.g., degrees |

# 5 Algorithm Usage Constraints

Describe the constraints or limitations on output data use based on the algorithm and the various assumptions made in producing the data product(s).

# 6 Performance Assessment

## 6.1 Validation Methods

Describe the details of the scientific methods utilized for algorithm performance assessment validation. Details provided should match the current algorithm maturity.

## 6.2 Uncertainties

Describe the uncertainties applicable to the validation methods and data used. This may include uncertainty in scientific or mathematical methods and/or errors associated with observation retrievals used for comparisons.

## 6.3 Validation Errors

Provide estimated errors of the output data product variables. Each variable included in the data product should have an associated error value. Use tables and/or figures as appropriate.

# 7 Algorithm Implementation

## 7.1 Algorithm Availability

Provide the web address for the algorithm open source code and any relevant information needed to execute the algorithm source code. For example, this may include execution instructions, memory requirements, programming languages, dependencies, who to contact, software rights. If this information is already provided at a permanent location, then provide the link to that location with a description.

## 7.2 Input Data Access

Candidate EMIT snow scenes:

- EMIT\_L2A\_RFL\_001\_20230216T211823\_2304714\_007
- EMIT\_L2A\_RFL\_001\_20230220T102026\_2305107\_004
- EMIT\_L2A\_RFL\_001\_20230220T194126\_2305113\_011
- EMIT\_L2A\_RFL\_001\_20230323T144706\_2308209\_038

Provide the URL or location of the input data product(s) used in the algorithm. Also, provide the relevant information needed to execute the algorithm source code using the input data. For example, this may include execution instructions, memory requirements, programming languages and dependencies. If this information is already provided at a permanent location, then provide the link. Use a table if needed to identify the URL and the description for each data product.

## 7.3 Output Data Access

Provide the algorithm output data persistent identifier, such as the DOI provided that the data products are already publicly available. If none is available leave blank, but instead provide any details needed to access the output data, including the name of the point of contact or any special instructions. If the data products are not yet publicly available but will be at a NASA DAAC, then note which DAAC and state the identifier information is TBD.

## 7.4 Important Related URLs

This section is for providing any alternative URLs that may be needed, including links to machine services, ordering services and DAAC websites. Provide details needed to access input or output data through alternative methods, if needed. Ensure readers have the information needed to understand what the related URL is for. If there is no related URLs, then state "No related URLs are needed".

# 8 Significance Discussion

Describe the significance of the work or the implications for or contributions to the relevant science fields.

# 9 Open Research

The data availability statements are required for journal submission and optional for inclusion in APT. Provide an Availability Statement for the underlying data needed to understand, evaluate, and build upon the reported research at the time of peer review and publication. Authors must include an availability statement for the software that has a significant impact on the research.

# 10 Acknowledgements

If you intend to submit this ATBD to an AGU journal, the acknowledgements section is required. This section is optional for inclusion of the ATBD in the APT. List all needed acknowledgements which may include people, organizations, funding programs, etc.

# 11 Contact Details

Provide in this section the contact information and respective role(s) for all contributors to this ATBD. Roles should match the CRediT taxonomy. The following information should be included:

First Name, Last Name (required)

UUID (optional) - universally unique identifier - such as ORCID

URL (optional) - relevant URL related to contacting the person - such as LinkedIn URL

Contact mechanism (required) - phone number, email, etc.

Role(s) related to this ATBD (required) - use CRediT - such as writing - original draft, data curation, methodology, software, etc.

Affiliation - name of associated organization

# References

The list should contain the full citation of all references cited within the text. References for any data or software used should also be included. Information related to reference formatting can be found in the [AGU Grammar and Style](https://www.agu.org/Publish-with-AGU/Publish/Author-Resources/Grammar-Style-Guide#referenceformat) guide.

The following is an example of a reference in the proper style:

Deng, A., & Stauffer, D. R. (2006), On improving 4-km mesoscale model simulations. _Journal of Applied Meteorology and Climatology_, _45_(3), 361–381. doi:10.1175/JAM2341.1

In-text citation style is (last name of author, YYYY) up to two authors. For more than 2 authors, use ,et al. after the first author. Example in-text references:

(Smith, 2020) - one author

(Doe & Smith, 2021) - two authors

(Doe et al., 2021) - three or more authors
