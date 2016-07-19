Work-Flow: by Zhihong Ye, written in 07/19/2016

1, Use "FASPD.C" to read in each background root-file, 
    evaluate the energy deposit of each backgorund particle in FASPD.
    Save to a new root-file for the energy-deposition together with 
    other important quantities related to FASPD.


2, Use "GetHisto_EM.C" to create background distribution into histograms.
   For example, EDep vs. radius, so we can randomly read the histrogran
   in the next step to roughly simualte the random energy-depositions with
   fluaction. 

3, Use "FASPD_Segment.C" to read in the root-file for high-energy photons from
   pi0 decay; study the energy-deposition when photons passing through FASPD,
   like via electron-pair productions. Then combine other background source, 
   generated in step 1. and 2., so we get a total energy deposition for a high
   energy photon which is acampanied with many other background particles.

4, Use "GetRject.C" to define a time-window and a threshold cut, so we can study
   the photon-rejection vs. how many segmentations. 

The complete files are located at: /w/work6501/yez/solid/spd_study

Note that: 
1, I also have some other scripts used to evaluate rates, e.g. GetRate.C
2, You need to run SPD simualtion with GEMC and different generated background
   (or ask Zhiwen to run for you, or use my old files)
3, You also need the trigger files which can be found in my work-directory.


