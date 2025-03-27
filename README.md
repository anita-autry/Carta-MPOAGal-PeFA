The three scripts provided, ScriptA, ScriptB and ScriptC are meant to be used in sequence. 
The three scripts were written in Matlab R2018b, changes might be needed if using a different version.
*For the scripts to work properly, the following packages and functions are required:
Signal processing toolbox (Add-on for Matlab)
stdshade.m (available through Matworks)
sem.m (available through Matworks)


ScriptA aligns photometry traces and behavioral data and was specifically designed to analyze data acquired through the Neurophotometrics System hardware and Bonsai software.
We provided 3 files as sample data for script A and B: 
1.	Photometry data with signal evoked by 470nm (calcium-dependent)
2.	Photometry data with signal evoked by 410nm (calcium independent)
3.	Datasheet of behavior scored, in binary format (1=behavior, 0=no behavior) 
The data was collected during an experiment in which a female mouse was exposed to a newborn pup. The duration of the recording was around 10 minutes.
Script A will load the three files, normalize the calcium data and overlay with behavior data.

Immediately after running Script A, do not clear the variables and run Script B, which will find perievent windows around behaviors of interest. 
Scripts A and B are meant to be used on single experiments and will save the results according to experiment name. Once all experiments are analyzed, create a merged .mat file that collects all the experiments, for example all the trials of all females exposed to pups, or other behaviors of interest. 
We provided also a sample data file for Script C:
1.	.mat dataset 
This dataset is from a group of mice that were subjected to tail suspension, and contains deltaF/F traces segmented around behavioral events.
Script C splits the data into female and male data, computes average zscore with standard error of the mean, max zscore and area under the curve. 
