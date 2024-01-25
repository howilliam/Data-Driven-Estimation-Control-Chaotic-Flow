README FOR SYSTEM ID AND THE KALMAN FILTER


ONLY ONE OF THESE CODES IS NEEDED TO START. 

all of these codes can be used to perform system ID and kalman filter. 

START WITH:

SID_V5_wQR.m is the code to 


SID_V4_iter.m is an iterative sensor search which replaces the sensor index (sensor_point) to evaluate the FIT values of a

SID_V5_fullcode.m also performs balanced trucnation, which is not discussed in the report. The method was not carried out. Discard unless requested. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE ON THE Q MATRIX


Beware while doing the system ID and extimating few modes. Sometimes the noise variance is not calculated by
n4sid, and the result is an IDENTITY matrix. This may be because not enough training dataset is used. 

The original definition by flavio Q=sysID.K*sysID.NoiseVariance*sysID.K'; works for systems IDs wich do not encounter
this error. This definition is encouraged. All system IDs in the report worked with this definition. #


Alternatiely the codes present a further definiton of the Q matrix bades on the prediction errors. It should 
be an equivalent deficiniton, however in practice there is a small difference due to numerical calculations. 
However, reconstruction results are exactly the same when using both definitions for Q. 


Again, Q=sysID.K*sysID.NoiseVariance*sysID.K'; is encouraged. 

%%%%%%%%%%%%%%%


