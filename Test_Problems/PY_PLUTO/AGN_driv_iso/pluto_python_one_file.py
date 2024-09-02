#!/usr/bin/env python

import subprocess,sys
import glob
from astropy import constants as c
from astropy import units as u
from scipy.integrate import quad
import pyPLUTO as pp
import numpy as np
import pluto_python_sub as pps
import re

pps.pluto2py_rtheta(4600)
#pps.python_input_file(0,data,cycles=2)
