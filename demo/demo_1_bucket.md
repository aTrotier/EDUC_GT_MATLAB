# Lecture 9 : Protyping at the scanner with MATLAB part 2
Title : Protyping at the scanner with MATLAB part 2

Schedule : June 25, 2020 | 16:00-17:00 

Speakers : Stanislas Rapacchi & Aurélien Trotier

# Demo 1

Objectives understand the input in matlab.



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


