# Ganymede TakeHome Project: Liquid Chromatography

This repo contains a take home project for my interview process with Ganymede for the Scientific Software Engineer role.  The code takes data from a liquid chromatography machine, parses it, and performs calculations to find the peaks and the peak volumes of the data.  All of this data is stored in the __ChromatogramRun__ class, which includes 3 data classes to store metadata, arrays to store raw data and calculations for peaks and integrations, plus methods to calculate these.

# Parsing
The code parses metadata and raw data from regular txt files, storing the metadata in dataclasses and the raw data in an array.  
* __Injection Info:__ Info about the date, injection parameters, and instruments used
* __Chromatogram Data Info:__ Info about the data collection and some base statistics about the run
* __Signal Parameter Info:__ Type of signal parameter
* __Raw Data:__ raw data from the machine including time (min) - step (s) - value (EU)

The code assumes well formatted data, error handling can be added for future iterations

# Peak Finding
Peak finding might seem trivial to the naked eye but is anything but.  The method chosen in this program is a form of z-score smoothing, based off of this StackOverflow post (https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data).  The idea is to keep a rolling average of your data based on a lag parameter that you can vary.  Then if you encounter data that is a certain number of standard deviations away from the rolling average (the threshold), then you have detected a peak.  You then record the start of the peak, the max of the peak, and the end of the peak, based on how long you stay above the threshold.  You also have an influence parameter that controls how much new data influences the rolling average during a peak.

The peak data is stored in a list of lists where each element includes the peak max location, the start location, end location, and peak height.

For a given set of data you normally have to adjust the parameters a little bit to get good results.  In this data, I've found good results with a threshold of 3.5, 0 influence, and a lag of 100.  

# Peak Integration
There are many different methods for numerical integration.  For this task I chose Simpson's rule, which uses quadratic interpolation to calculate areas.  The error for this method is proportional to -(h^5)/90 which should be sufficient to keep the error low.  The calculations are made between the thresholds for each peak using the raw data to calculate the heights to provide and accurate area.  The area for each peak is stored in the peakVolumes array and then totaled to give a final result.

# Visualization
The visualize method plots the raw data and peak data if it is present.  This helps to adjust parameters to ensure accurate results

# How to Run
There is one file to run, GanymedeLCTakeHome.py  

* __ChromatogramRun(file_name)__: initializes chromatogram class; requires filename of the chromatogram data
* __.find_peaks(self, threshold, influence, lag)__: method to find peaks in the data.  Threshold: number of standard deviations away from mean to detect peak. Influence - influence of new data on rolling average (0-1) (0 for stationary data).  Lag: How far back you look for the rolling average
* __.elutionVolume(self)__: calculates the total elution volume of the peaks
* __.visualize(self)__: creates a visual of the raw data and any found peaks

By default, it will cycle through the example data, creating a new ChromatogramRun instance for each loop.  It will then call the methods to find the peaks and peak volumes and then visualize the data.  There is also a loop that is commented out where you can input the files you want to read on the command line and they will be inserted into the loop.  External libraries used include, sys, dataclasses, numpy, and matplotlib

# Limitations
There are some limitations with the current state of the code to address.
* __Peak Finding:__ The method used does a good job of smoothing out the noise in the data, but it is very sensitive to the parameters and struggles to resolve double peaks, instead lumping them into one finding.  There are methods that you can use which looks at the slope of the data and starts finding a peak once that reaches a certain baseline.  The program would then find when the slope went +0- (top of the peak) and then cut off once the slope levels off again.
* __Integration:__ The same issue that affects the peak finding method also affects the integration scheme.  Where one might want two different areas calculated, this method will only be able to resolve one peak, and combine both of the areas together.  Other methods address this issue by setting a baseline for the double beak, then dropping a line down from the trough to the baseline and separating the areas at that point.
* __Continuous Data:__ This current method uses a rolling average to find peaks, but only for a given dataset.  You can adjust the method to use continuous new data by ensuring the new data is properly integrated into the rolling averages and stds.
* __Testing:__ I wasn't able to find sources of truth for the given data, so my testing has been primarily heuristic.  It would be good to find some test data to confirm the accuracy of the methods being used.
