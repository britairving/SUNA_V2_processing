# Introduction
This processing package was developed to process data collected with a SUNA V2 deployed on a mooring. The steps described in this document are similar those described in the SUNA v2 manual and available through Sea-Bird’s UCI software. This toolbox has been developed to allow more control and better understanding throughout the process, and produce consistent data products year to year while maintaining flexibility.

All processing is done using MATLAB routines built in-house and includes data handling, integration, and reprocessing of nitrate data from the SUNA V2 sensor. Goals of the processing package are detailed documentation, flexibility and portability so others may use it, and incorporation of best practices for data reproducibility, quality control, and uncertainty estimation, etc. 

Processing scripts were written in MATLAB R2017b using Microsoft Windows 10 Pro and tested in Matlab R2019a on MacOS Sierra and was written with the hope it will port to other MATLAB versions and operating systems.

Questions, feedback, and comments are welcomed and encouraged!

Thanks to [Seth Danielson](https://www.uaf.edu/cfos/people/faculty/detail/seth-danielson.php) and [Tyler Hennon](https://uaf.edu/cfos/people/research-staff-and-post-docs/detail/tyler-hennon.php) for their feedback and debugging help on this toolbox!

Please review the following wiki pages. 
  
[1. Assumptions, notes, and recommendations](https://github.com/britairving/SUNA_V2_processing/wiki/1.-Assumptions,-notes,-and-recommendations)

[2. Processing background and workflow](https://github.com/britairving/SUNA_V2_processing/wiki/2.-Processing-background-and-workflow)

[3. Automatic QC procedures](https://github.com/britairving/SUNA_V2_processing/wiki/3.-Automatic-QC-procedures)

[4. SUNA calibration corrections](https://github.com/britairving/SUNA_V2_processing/wiki/4.-SUNA-calibration-corrections)

# Citation and acknowledgement
If you use this toolbox, please cite as follows: Brita Irving, (2021), SUNA_V2_processing, GitHub repository, https://github.com/britairving/SUNA_V2_processing/

# References
Daniel, A., Laës-Huon, A., Barus, C., Beaton, A.D., Blandfort, D., Guigues, N., Knockaert, M., Muraron, D., Salter, I., Woodward, E.M.S., Greenwood, N., Achterberg, E.P. (2020). Toward a harmonization for using in situ nutrient sensors in the marine environment. Front. Mar. Sci., 6,  p. 773,doi:10.3389/fmars.2019.00773.

Frank, C., Meier, D., Voß, D., and Zielinski, O. (2014). Computation of nitrate concentrations in coastal waters using an in-situ ultraviolet spectrophotometer: Behavior of different computation methods in a case study a steep salinity gradient in the southern North Sea. Methods Oceanogr. 9, 34–43. doi: 10.1016/ j.mio.2014.09.002.

Intergovernmental Oceanographic Commission IOC (2013). Ocean Data Standards in Recommendation for a Quality Flag Scheme for the Exchange of Oceanographic and Marine Meteorological Data, (Paris: UNESCO-IOC), doi: 10.25607/OBP-6.

Johnson, K. S., and L. J. Coletti (2002). In-situ ultraviolet spectrophotometry for high resolution and long-termmonitoring of nitrate, bromide and bisulfide in the ocean, Deep Sea Res. Part I, 49(7), 1291–1305, doi:10.1016/s0967-0637(02)00020-1.

Mordy et al. 2020. Seasonal and interannual variability of nitrate in the eastern Chukchi Sea: Transport and winter replenishment. https://doi.org/10.1016/j.dsr2.2020.104807

Pellerin, B. A., Bergamaschi, B. A., Downing, B. D., Saraceno, J. F., Garrett, J. D., and Olsen, L. D. (2013). Optical techniques for the determination of nitrate in environmental waters: Guidelines for instrument selection, operation, deployment, maintenance, quality assurance, and data reporting. U.S. Geological Survey Techniques and Methods 1-D5, 37. doi: 10.3133/t m1D5.

Randelhoff, A., Fer, I., Sundfjord, A., Tremblay, J.-É., & Reigstad, M. (2016). Vertical Fluxes of Nitrate in the Seasonal Nitracline of the Atlantic Sector of the Arctic Ocean. Journal of Geophysical Research: Oceans, 121(7), 5282–5295. https://doi.org/10.1002/2016JC011779.

Sakamoto, C. M., Johnson, K. S., Coletti, L. J., and Jannasch, H. W. (2017). Pressure correction for the computation of nitrate concentrations in seawater using an in-situ ultraviolet spectrophotometer. Limnol. Oceanogr. Methods 15, 897–902. doi: 10.1002/lom3.10209.

Sea-Bird Scientific (2018). User manual Submersible Ultraviolet Nitrate Analyzer SUNA V2.


