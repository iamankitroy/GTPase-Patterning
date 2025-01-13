# GTPase-Patterning

All data analysis was conducted using _Python 3.10.13_. _Numpy 1.24.4_ and _Pandas 2.1.2_ were used for data manipulation and analysis. _R version 4.2.2_ was utilized for the execution of R scripts. Shell scripts were executed in _Zsh version 5.9_ for batch processing and file handling.

Image analysis and segmentation were performed using _Sci-kit Image 0.22.0_, _OpenCV 4.7.0_, and _Scipy 1.11.3_. Parallelization was implemented where applicable with the _Multiprocessing_ package to enhance processing efficiency.

Single molecule spots were tracked with _TrackMate_. Colocalized GTPase and GDI spots were identified using the _SpotColocalization.py_ script. Colocalization statistics were calculated using _getStat_TracksColocalized.py_ and colocalization probabilities at specific positions throughout the lifetime of tracks were determined with _calc_ColocalizationProbability_positionSpecific.py_. Lipid patterns were segmented using _autoSegmentation_dice-N-splice_LipidPatch.py_, and subsequent analysis was carried out with _GTPase_patterning_analysis.py_.


Scripts for analysing colocalizations in dual channel single molecule TIRF microscopy data of GTPase and GDI recruitment/extraction events on lipid bilayers. 
- **SpotColocalization_SingleFrame.py**: Used for identifying spot colocalizations from dual channel single molecule TIRF microscopy data. Requires *All_Spots_statistics.csv* files generated with *TrackMate* in *Fiji* for spot and track detection in both channels.
