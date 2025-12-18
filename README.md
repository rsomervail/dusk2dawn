# dusk2dawn
Dusk2Dawn plugin for the EEGLAB toolbox in MATLAB.

- Dusk2Dawn allows users to easily clean whole-night sleep EEG data using sleep-specific implementations of Artifact Subspace Reconstruction (ASR).   
- The accessible GUI interface also allows users to easily test a range of ASR parameters and visualise the effects on their data (e.g. the effects on Slow-Wave amplitude), before deciding which set of parameters to use.   
- The plugin uses functions from the official 'clean_rawdata' ASR plugin included in recent versions of EEGLAB.  

	Author: Richard Somervail, Istituto Italiano di Tecnologia, 2022  
		r.somervail@gmail.com    
		www.iannettilab.net

  ** NOTE THAT THE GRAPHICAL USER INTERFACE (GUI) IN EEGLAB AND DUSK2DAWN ONLY SEEMS TO WORK ON MATLAB VERSION UP TO 2024b  (2025a+ changed how the GUI works) **
		-> PLEASE USE WITH MATLAB 2022-2024

If you use this plugin, please cite as:
Somervail R, Cataldi J, Stephan AM, Siclari F, Iannetti GD. 2023. 
Dusk2Dawn: an EEGLAB plugin for automatic cleaning of whole-night sleep electroencephalogram using Artifact Subspace Reconstruction. 
Sleep. 1–14.
		  
	History:  
    15/01/2024 ver 3.4.1 Added published paper to cite when using D2D  
    25/10/2023 ver 3.4.0 Added button to GUI to automatically generate a script from current settings. Also misc bugfixes.  
	27/09/2023 ver 3.3.1 GUI bugfixes & misc improvements.  
    07/09/2023 ver 3.3.0 Crucial bugfix to allow to computing of frequency spectrum for datasets containing no events  
	18/08/2023 ver 3.2.0 Crucial bugfix to allow datasets containing no events   
	09/05/2023 ver 3.1.0 Crucial bugfix to sliding window ASR  
	22/03/2023 ver 3.0.0 Crucial bugfix to sliding window ASR, major efficiency upgrades, all ASR parameters now variable from GUI + additional D2D-specific features.  
	09/03/2023 ver 2.1.0 Crucial bugfixes. Also added most ASR parameters + various quality of life improvements.  
	06/03/2023 ver 2.0.1 Various quality of life upgrades (e.g. printing time elapsed for each dataset).  
	24/02/2023 ver 2.0.0 Added group-level/batch processing of multiple datasets.  
	30/01/2023 ver 1.0.1 Patch to fix initial bugs and solidify basic functionality.  
	19/01/2023 ver 1.0.0 Created  
	
