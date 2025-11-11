# My soft X-ray background tool and database <br> - based on ROSAT all-sky maps -

## Database file
- sxrbg_hpmap_ns64r50_ray.fits <br>
Soft X-ray background spectra of ROSAT PSPC energy bin (7 bands 0.1-2 keV), modelled with Raymond (local hot bubble) + Powerlaw (CXB) function for all healpix of NSIDE = 64 (Number of pixels = 49152).

- xrbkg_rayspec.fits <br>
X-ray background spectra of nominal energy bin (0.1-10 keV / 5 eV step = 1980 bins), with the same model function for the same healpix of NSIDE = 64 (Number of pixels = 49152).

## Examples
- Example to use sxrbg_hpmap_ns64r50_ray.fits <br> 
  - python script [plt_sxbkgspec.py](plt_sxbkgspec.py)
  - juypyter-notebook [plt_sxbkgspec.ipynb](plt_sxbkgspec.ipynb)
- Example to use xrbkg_rayspec.fits <br> 
  - python script [plt_xrbkgspec.py](plt_xrbkgspec.py)
  - juypyter-notebook [plt_xrbkgspec.ipynb](plt_xrbkgspec.ipynb)
