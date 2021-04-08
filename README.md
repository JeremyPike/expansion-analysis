# Analysis scripts for analysis of nuclear nano structure from expansion microscopy datatsets. 

1. expansion_nuclear_analysis.groovy is a [Fiji](https://imagej.net/Fiji/Downloads) script . It makes use of the [TrackMate](https://imagej.net/TrackMate) plugin to detect spots in both the site and satellite channels. Site spots are then clustered to determine the location of individual sites. Cropped images of sites are presented to the user who is asked to classify them into predefined categories for downstream processing. Site crops along with various spot statistics are saved to disk.
2. particle_averaging.m is a Matlab script which takes the cropped sites outputted by the Fiji script and averages them. Binned radial plots are produced and the average site for each class is saved to disk.

## Citation

If these scipts are useful please consider citing:

Imaging Nanoscale Nuclear Structures with Expansion Microscopy\  
Emma L. Faulkner, Jeremy A. Pike, Ruth M. Densham, Evelyn Garlick, Steven G.Thomas, Robert K. Neely1 and Joanna R. Morris. \
In preperation.