# Quality_Assurance

This package perfroms the Quality Assurance of CT volumes. It evaluates three different capacities from the CT: CT number calibration, edge detection and the calculation of radiomic features.

## Version:
_Currently under development_

This package uses reference segmentations and images for the CIRS electron density phantom, therefore new images added to the DDBB are registered and resampled to the reference, so that is not necessary dedicated segmentations for every image. The package supplies the reference images and segmentations, although the user can freely change it. 

The required structure for the DDBB is detailed below.

![DDBB structure](images/ddbb_diagram.png)
