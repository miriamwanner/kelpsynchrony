## Dispersal synchronizes giant kelp forests

This code accompanies the study *Dispersal synchronizes giant kelp forests*. It was created with R Version 4.0.5 (2021-03-31), and the up to date versions of the packages listed in requirements.txt.

Before running any code, the data needs to be downloaded [from this website](https://portal.edirepository.org/nis/mapbrowse?packageid=knb-lter-sbc.162.1) and saved in the ```Data``` folder with the filename ```CAkelpCanopyEnv_2021.nc```. All the code is run in the ```Code``` folder. Running the ```main.R``` file runs all of the experiments, and saves the figures in the ```results``` folder. To recreate specific experiments, see below.

### Recreating Figures and Results
#### Initial Data Processing
The initial data processing is done in ```init_data_processing.R```. The spore survival rate (1 - spore loss rate) is set to 0.1 (90% spore loss rate), but changing line 33 in ```init_data_processing.R``` can be used to recreate experiments with a different survival rate.

#### MRM Results
Once the initial data processing has run, MRM results can be recreated by running ```mrm_results.R```. The results are stored as ```Rds``` files in the ```results``` folder. Note that there is some randomization within the MRM method, and therefore *p*-values will not be exactly the same as in our study, but the same patterns will exist.

#### Spline Correlogram
To recreate the spline correlogram, run ```spline_figure.R```, and the figure will be saved as ```spline_correlogram.pdf``` in the ```results``` folder.

#### Clustering Maps and their Average Time Series
To recreate the clustering figures, run ```clustering_figures.R```, and the figures will be saved as ```all_clustering_map.pdf```, ```mainland_clustering_map.pdf```, and ```island_clustering_map.pdf``` in the ```results``` folder.

To then create the average time series figure for each cluster, run ```avg_ts_figures.R```, and the figures will be saved as ```all_avg_ts.pdf```, ```main_avg_ts.pdf```, and ```island_avg_ts.pdf```. The ```avg_ts_figures.R``` file can only be run after running ```clustering_figures.R```.

#### Models Plotted Against Synchrony
To recreate the models plotted against synchrony, run ```transparent_plots_figure.R```, and the figure will be saved as ```transparent_plots_figure.pdf``` in the ```results``` folder.

#### Synchrony Matrices
To recreate the plotted synchrony matrices, run ```matrices_figure.R```, and the figure will be saved as ```matrices_figure.pdf``` in the ```results``` folder.
