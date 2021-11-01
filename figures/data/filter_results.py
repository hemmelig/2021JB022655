# coding: utf-8
"""
Filters MFAST summary files for inputs used for analysis and figures in:

    On the origin of seismic anisotropy in the shallow crust of the Northern
    Volcanic Zone, Iceland
    Bacon, C.A., Johnson, J., White, R.S., and Rawlinson, N.

which has been submitted to the Journal of Geophysical Research - Solid Earth.

We apply the following filters:

    1. Cluster grade - only use measurements that have A-graded clusters
    2. SNR - only use measurements with an SNR above 4
    3. Error in phi - only use measurements with an error in phi < 10 degrees,
       as determined from the confidence plots
    4. Value of dt - only use measurements with dt < 0.4 * the full grid-search
       range (1.2 s -> 0.48 s)
    5. Error in dt - only use measurements with an error in dt < 0.05 s, as
       determined from the confidence plots
    6. Depth - only use measurements from events shallower than 10 km
    7. Shear-wave window - only use measurements within a SWW of 50 degrees
    8. Remove events during the Holuhraun eruption to minimise impact of stress
       transients
    9. Finally, filter by manual labels (applied to the events that passed up to
       step 8)

"""

# --- Import libraries ---
import pathlib
import pyproj

import pandas as pd
from scipy.stats import circmean, circstd

# --- Read in summary files ---
summary_dir = pathlib.Path.cwd() / "mfast_station_results"
summary_files = summary_dir.glob("*.summ")
summaries = pd.read_csv(next(summary_files))
for summary_file in summary_files:
    summaries = pd.concat([summaries, pd.read_csv(summary_file)])
summaries.reset_index(drop=True, inplace=True)

# --- Filter summary files ---
grade = "ACl"
snr = 4.0
fast_error = 10.0
dt_error = 0.05
dt_max = 0.4 * max(summaries["tlag"].values)
sww = 50.
above = 10

# Grade
summary = summaries[summaries["gradeABCNR"] == grade]

# SNR
summary = summary[summary["20SNR"] >= snr]

# Fast error
summary = summary[summary["Dfast"] <= fast_error]

# Delay time error
summary = summary[summary["Dtlag"] <= dt_error]

# Delay time
summary = summary[summary["tlag"] <= dt_max]

# Filter by min/max depth bounds
summary = summary[(summary["depthkm"] < 30) & (summary["depthkm"] > -0.5)]
summary = summary[summary["depthkm"] <= above]

# Remove data during dike intrusion and eruption in 2014
summary = summary[(summary["year"] != 2014) | (summary["doy_det"] <= 227)]

# Drop data for stations outside of study region
summary = summary[(summary["slat"] <= 65.2520) & (summary["slat"] >= 64.9470)
                  & (summary["slon"] <= -16.1148) & (summary["slon"] >= -16.8808)]
        
# Order the DataFrame by depth
summary.sort_values(by=["depthkm"], inplace=True)

summary.reset_index(drop=True, inplace=True)

# Calculate the midpoint for the event-station pairs
cproj = pyproj.Proj("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
gproj = pyproj.Proj("+proj=lcc +lon_0=-16.5 +lat_0=65.1 +lat_1=64.7 +lat_2=65.5 +datum=WGS84 +units=m +no_defs") 

midpoints = pyproj.transform(cproj, gproj,
                             ((summary["slon"] + summary["evlo"]) / 2).values,
                             ((summary["slat"] + summary["evla"]) / 2).values)
summary["midx"], summary["midy"] = midpoints

# Add manual labels (where available)
inspected_splits = pd.read_csv("inspected_splits.csv")
labels = []
for i, event in summary.iterrows():
    evid = ".".join(event["1event"].split(".")[:6])
    label = inspected_splits[inspected_splits["EventID"] == evid]
    if len(label) == 0:
        labels.append("u")
    else:
        labels.append(label["Label"].values[0])
summary["label"] = labels

# Output
summary.to_csv("filtered_shallow_results_unlabelled.summ", index=False)
labelled_results = summary[summary["label"] != 0]
labelled_results.to_csv("filtered_shallow_results.summ", index=False)

print("Number of events before filtering for shear-wave window - "
      f"{len(labelled_results)}")

# Filter by incoming angle
summary = summary[summary["anginc"] <= sww]
summary.to_csv("filtered_shallow_results_sww_unlabelled.summ", index=False)
labelled_results = summary[summary["label"] == 1]
labelled_results.to_csv("filtered_shallow_results_sww.summ", index=False)

fasts = labelled_results["fast"].values
dts = labelled_results["tlag"]
print("Number of events after filtering for shear-wave window - "
      f"{len(labelled_results)}")
print(f"Average delay time       = {circmean(fasts, high=180):5.1f} "
      f"+/- {circstd(fasts, high=180):5.1f}")
print(f"Average fast orientation = {dts.mean():5.2f} +/- {dts.std():5.2f}")
