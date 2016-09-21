import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import termcolor
import numpy as np
from astropy.coordinates import SkyCoord as sky
import astropy.units as unit
from astropy.cosmology import LambdaCDM
import astropy.cosmology as cosmo


print(cosmo.z_at_value(LambdaCDM(H0=70,Om0=0.3,Ode0=0.7).distmod,(17+20.89)*unit.mag))