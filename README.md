# Monitoring Mammals with Remote Camera Deployments

The ABMI monitors mid- to large- sized mammals using remote cameras deployed at terrestrial and wetland sites across Alberta. These data are used to measure habitat use, distribution, and current status of mammals, as well as monitor changes in mammal populations regionally and provincially over time. 

## Background

The ABMI implemented a new field protocol in 2015 to monitor the presence of vertebrates (mainly mammals) using Reconyx PC900 Hyperfire remote camera traps. The aim of this protocol is to compile a comprehensive dataset on species occurrence at ABMI terrestrial and wetland sites while minimizing the time and effort associated with collecting these data. We place Reconyx PC900 Hyperfire remote cameras in a square approximately 600 m apart and centred on an ABMI site to passively monitor and record mid- to large-sized mammals.

## Data Download 

ABMI remote camera mammal data is a raw data file comprised of species occurrence information at surveyed ABMI sites. For each site, occurrence is determined by identifying the species present within each image taken by the remote cameras deployed at a site; all records of detections of each species are provided. To download the raw data, visit the ABMI website [here](https://abmi.ca/home/data-analytics/da-top/da-product-overview/remote-camera-mammal-data/remote-camera-mammal-data-download.html).

---

## Current Code Base

This repository currently hosts the code associated with estimating mammal density using images captured from remote cameras. The ABMI uses a method based on random encounter and staying time (REST), utilizing data collected on the number of individuals, the total time in the camera field of view, the area of the camera field of view, and the total camera operating time in the field. Further details can be found in this [report](https://www.abmi.ca/home/publications/501-550/516).

Available documentation:

1. Code-annotated version of the **Animal Density from Camera Data** report can be viewed [here](https://abbiodiversity.github.io/mammals-camera/01_Process-Cam-Data-Calc-Animal-Density_10-15-2019.html).

Under development:

1. ABMI habitat modeling for mammals using density data.

1. Empirical comparisons to other density methods.

1. Estimating density for user-defined areas of interest. 

## Using the Raw Data

If you would like to use the raw data to estimate density, or be able to run the code available in this repository, please contact Marcus Becker (mabecker@ualberta.ca).
