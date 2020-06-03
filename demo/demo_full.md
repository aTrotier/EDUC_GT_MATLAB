# Lecture 9 : Protyping at the scanner with MATLAB part 2
Title : Protyping at the scanner with MATLAB part 2

Schedule : June 25, 2020 | 16:00-17:00 

Speakers : Stanislas Rapacchi & Aurélien Trotier

# Full demo

## Sequence and Data

We will use a 3D MP2RAGE sequence with a variable density poisson undersampling mask acquired on a 3T Prisma from Siemens.

Data is available at this [link](https://www.rmsb.u-bordeaux.fr/nextcloud/index.php/s/YyTDmJrw6225D2p)

One dataset (without noise calibration) is available: 

- brain, 0.8mm isotropic, acceleration factor = 20. 

The data has been converted with **siemens_to_ismrmrd**, we will not discuss data conversion here. This will be the object of the following readings.
Put the dataset in the folder **data/**

## Procedure for demo

In order to launch this demo :

* Launch matlab
* Add matlab folder to the pass
* change the path to your bart installation in file **matlab/GT_reco_MP2RAGE_CS** line 64 :
`run('/home/CODE/bart/startup.m');`
* type this command in the matlab command window: `gadgetron.external.listen(18000,@GT_reco_MP2RAGE_CS)`
* open 2 terminals. In the first one, launch gadgetronn server : `gadgetron`
* In the second on launch, launch the command line (change the path to data,output data, xml config) : 
`gadgetron_ismrmrd_client -f /home/atrotier/GITHUB/EDUC_GT_MATLAB/data/MP2RAGE_CS20.h5 -o /home/atrotier/GITHUB/EDUC_GT_MATLAB/data/out_MP2RAGE_CS20.h5 -C /home/atrotier/GITHUB/EDUC_GT_MATLAB/config/MP2RAGE_CS.xml`
* Results will be saved in the **data/** folder and can be visualize using the command : `img = read_image_h5();`


