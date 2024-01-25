Codes for data extraction from STAR-CCM+

This folder contains two important aspects of the controller:
	- the MATLAB code which deisigns the controller from the excited system
	- the JAVA code that is used to connect with STAR-CCM+ to do the closed loop control

All files included should be sufficient

	- MATLAB code saves matrices into CSV file, Java code reads the CSV files to do calculation within the STAR-CCM+ macro
	- linear algebra was done useng Apache Commons version 3, all appropriate JAR files are included in the 'lib' folder 
	- You will need to connect the directory of the lib folder in STAR-CCM+, and then simply run the macro


	ADVICE:
	- all macro syntax can be found in the STAR-CCM+ user guide, however I found all relevant syntax by playing with the macro itself
	- CFD Online forum and Siemens forums were especially helpful for anything STAR-CCM+ related 
	
	



Check this link for: shows the solution for some sort of closed loop control using a PID controller:

https://community.sw.siemens.com/s/question/0D54O00007mbeDKSAY/i-would-like-to-have-a-pressure-outlet-condition-change-based-on-a-target-inlet-pressure-the-inlet-boundary-condition-will-be-mass-flow-and-the-outlet-pressure-out-is-this-possible-with-code-or-design-manager-or-some-other-method


