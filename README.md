# WRF-UACM Urban Model (v3.8)
### Weather Research and Forecasting – Urban Asymmetric Convetive Model

---

## Overview

The **WRF–Urban Atmospheric Canopy Model (WRF-UACM v3.8)** is an advanced urban parameterization implemented within the **WRF v3.8** modeling framework.  
It enhances the representation of:

- Urban morphology
- Urban surface energy balance
- Boundary-layer processes
- Anthropogenic heat fluxes (optional)

**NOTE**: Use `README.txt` file for the WRF-UACM model setup and configuration setting details.

---

## Author & Contact

**Prepared by:**  
**Utkarsh P. Bhautmage**  
_PhD, ENVR, HKUST (Hong Kong) and UNC (USA)_

**Email:**  
- upbhautmage@connect.ust.hk  
- yukimicky7@gmail.com  

**Last edited:** 29-SEP-2025

---

## Citation

**For any work related to the WRF-UACM model implementation, cite the following paper:**

Bhautmage, U. P., Fung, J. C. H., Pleim, J., & Wong, M. M. F. (2022). Development and evaluation of a new urban parameterization in the weather research and forecasting (WRF) model. Journal of Geophysical Research: Atmospheres, 127(16), e2021JD036338.

---

## Repository Structure

```text
WRF/
├── frame/
│   └── module_domain_type.f90
├── phys/
│   ├── get_env.F
│   ├── Makefile
│   ├── module_bl_acm.F
│   ├── module_pbl_driver.F
│   ├── module_physics_init.F
│   ├── module_sf_pxlsm.F
│   ├── module_sf_pxlsm_data.F
│   └── module_surface_driver.F
├── dyn_em/
│   ├── module_first_rk_step_part1.F
│   └── start_em.F
├── Registry/
│   └── Registry.EM_COMMON
└── run/
    ├── dhlurb_utk.txt
    └── namelist.input



