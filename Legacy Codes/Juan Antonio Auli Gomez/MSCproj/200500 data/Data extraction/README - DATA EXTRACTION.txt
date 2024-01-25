Codes for data extraction from STAR-CCM+



Data-extraction_DNS.m is the main code which extracts the data from internal XYZ tables. 
There are two sets of CSV files. The first contains the velocity (or other values of interest) along the presentation grid.
the second set contains the values along the fish wall (arbitrary grid). For the specific case of this work, the tables generated a double layer data 
along the fish, where the upper layer was not of interest. This is why in the code the indecies for removing the data are present. 
Beware that the order for the XYZ tables from star-ccm+ is not consistent. Use this code only for a guide on how to write
your own code to generate the .mat file required for further analysis.



fishdata.m extracts the values from the fish wall CSV files. (as describes above)

get_data.m extracts the data from the presentation grid. (Grid values)


