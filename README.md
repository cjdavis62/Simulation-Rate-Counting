# Simulation-Rate-Counting
Rate counting script for calibration simulations

To run Th232 peaks, just run plot_AllString_Th232Peaks.cc


To run Co56 peaks, cd into Co56_PeakFits/ and then run
>python Fitting.py -i inputfile.root

This will generate plots in Co56_PeakFits/output and a file with the efficiency (signal / background) of each peak

Then run cd up a directory back to Simulation-Rate-Counting and run
> root plot_AllString_Co56Peaks.cc

This will generate the root file with the rates/times for each peak
