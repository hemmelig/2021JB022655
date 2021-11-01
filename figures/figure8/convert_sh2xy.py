# -*- coding: utf-8 -*-
"""
This script will convert files containing the vectors representing the maximum
horizontal stress (generated from the stress field modelled in Coulomb) to .xy
files that can be plotted in GMT

Author: Conor Bacon
Date: 26 February 2020
"""

import pandas as pd
import pyproj

# Define a projection
proj = pyproj.Proj(proj="utm", zone="28W", north=True, ellps="WGS84", units="m",
                   no_defs=True)

# Set a reference point
refx, refy = proj(-17.6, 64.5)

vectors = pd.read_csv("SH_files/SH_2km.txt", header=None, delim_whitespace=True)
xy_file = "SH_files/SHmax_2km.xy"
with open(xy_file, "w") as f:
    for i in range(len(vectors)):
        if i % 3 == 1 or i % 3 == 2:
            pass
        else:
            # Convert vector to lat/lon
            p11, p12 = (float(vectors.iloc[i][0])*1000 + refx,
                        float(vectors.iloc[i][1])*1000 + refy)
            p21, p22 = (float(vectors.iloc[i+1][0])*1000 + refx,
                        float(vectors.iloc[i+1][1])*1000 + refy)
            lon1, lat1 = proj(p11, p12, inverse=True)
            lon2, lat2 = proj(p21, p22, inverse=True)
            f.write(f">\n{lon1} {lat1}\n{lon2} {lat2}\n")
