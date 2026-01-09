##Readme.txt file for compiling the WRF Urban Model  [WRF-UACM Code (v3.8)]
##Prepared by: Utkarsh P. Bhautmage [PhD, ENVR, HKUST(Hong Kong) and UNC(USA)]
##Contact: upbhautmage@connect.ust.hk 
##Last edited on: 29-SEP-2025.
===================================================================================
===================================================================================
Begin
==================

NOTE1:- For running the WRF-UACM Code (v3.8), copy following files into the respective directories and then recompile the WRF code!!!

/WRF/frame/
module_domain_type.f90


/WRF/phys/
get_env.F
Makefile
module_bl_acm.F
module_pbl_driver.F
module_physics_init.F
module_sf_pxlsm.F
module_surface_driver.F
module_sf_pxlsm_data.F


/WRF/dyn_em/
module_first_rk_step_part1.F
start_em.F


/WRF/Registry
Registry.EM_COMMON


/WRF/run
dhlurb_utk.txt
namelist.input


configure the WRF v3.8 code with option(3): PGI (pgf90/gcc), or option (15): intel before recompiling the WRF code.



NOTE2:- Main time step in namelist.input = 8 sec. (with 1:2 timestep ratio for nested domains)

For running the WRF-UACM model with other time steps, make sure the dtpbl/DTPBL/DT conditions are properly set in:
 module_bl_acm.F
 module_sf_pxlsm.F
 module_surface_driver.F

repectively!


dtpbl/DTPBL/DT should be set in the code in a way to target the innermost domain of interest!




NOTE3:- Set the following variables at line no. 1948 in module_sf_pxlsm.F 

  NOD   = 354.0     !29.0      ! Number of a day in a year   ! *IMPORTANT NOTE*:- Set this value equal to the first day no. of a WRF simulation case (Simulation Day1)
  LATT  = 28.6      !22.3      ! Latitude value              ! *IMPORTANT NOTE*:- Set this value equal to the latitude of the region (domain) of interest
  LONGT = 77.219    !114.1     ! Longitude value             ! *IMPORTANT NOTE*:- Set this value equal to the longitude of the region (domain) of interest
  DTUTC = 5.5       !8.0       ! Time difference             ! *IMPORTANT NOTE*:- Set this value (in absolute) as the difference between the local time (LT)
                                                           !                     and Universal Coordinated Time (UTC) in hours



NOTE4:- Set the times according to the region of interest and delta UTC in the IF conditon at line nos. 1981-1985, 2001-2006 in module_sf_pxlsm.F

16x60x60 = 57600
57600 + 86400 = 144000 




NOTE5:- To add the anthropogenic heat fluxes:

Uncomment the ANTHFX1 and ANTHFX_RF at line nos. 1708 and 1751 respectively in the module_sf_pxlsm.F!!!!
Also Note that the flux distribution ratio is set to 1:9 of the total flux estimation at the ground and roof level resp.!! 




NOTE6:- Check all the lines with 'IMPORTANT NOTE' in below modules to be sure to run the model sucessfully. 
 module_bl_acm.F
 module_sf_pxlsm.F
 module_surface_driver.F

repectively!




NOTE7:- This version of the WRF-UACM model can be used only with the USGS Landuse.





NOTE8:- For creating urban mprphological parameter data file like dhlurb_utk.txt

Prpare the data in the .csv column format text file with columns as:

grid_X | grid_Y | Average_building_height | Plan_area_fraction | Frontal_area_fraction | Canyon_street_orientation | Population_density  

please note - Population density is required only for estimating the anthropogenic fluxes (if want to include the fluxes)

Also, default environment values can be set at below lines in the module_surface_driver.F and module_bl_acm.F for those urban grid cells where real urban morphological parameter values culd not be generated:
         CALL GET_ENV (LF2, 'LAMF', 0.45)
         CALL GET_ENV (LP2, 'LAMP', 0.45)
         CALL GET_ENV (H2, 'HGTBDG', 6.0)
         CALL GET_ENV (STANG2, 'STANG', 0.0)
         CALL GET_ENV (PPLN2, 'PPLN', 100.0)




NOTE9:- Cite the below paper for any work related to the WRF-UACM model implementation.

Bhautmage, U. P., Fung, J. C. H., Pleim, J., & Wong, M. M. F. (2022). Development and evaluation of a new urban parameterization in the weather research and forecasting (WRF) model. Journal of Geophysical Research: Atmospheres, 127(16), e2021JD036338. 




For any queries related to the model implementation please contact on:- upbhautmage@connect.ust.hk OR yukimicky7@gmail.com

Thank you!


==============
End






