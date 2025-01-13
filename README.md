# GTPase-Patterning

All data analysis was conducted using Python 3.10.13. Numpy 1.24.4 and Pandas 2.1.2 were used for data manipulation and analysis. R version 4.2.2 was utilized for the execution of R scripts. Shell scripts were executed in Zsh version 5.9 for batch processing and file handling.

Image analysis and segmentation were performed using Sci-kit Image 0.22.0, OpenCV 4.7.0, and Scipy 1.11.3. Parallelization was implemented where applicable with the Multiprocessing package to enhance processing efficiency.

Single molecule spots were tracked with TrackMate. Colocalized GTPase and GDI spots were identified using the SpotColocalization.py script. Colocalization statistics were calculated using getStat_TracksColocalized.py and colocalization probabilities at specific positions throughout the lifetime of tracks were determined with calc_ColocalizationProbability_positionSpecific.py. Lipid patterns were segmented using autoSegmentation_dice-N-splice_LipidPatch.py, and subsequent analysis was carried out with GTPase_patterning_analysis.py.


Scripts for analysing colocalizations in dual channel single molecule TIRF microscopy data of GTPase and GDI recruitment/extraction events on lipid bilayers. 
- **SpotColocalization_SingleFrame.py**: Used for identifying spot colocalizations from dual channel single molecule TIRF microscopy data. Requires *All_Spots_statistics.csv* files generated with *TrackMate* in *Fiji* for spot and track detection in both channels.
