# AcousticRoutinesPSU

Step 1: Write a script pXXXX_sa similar to p4640_sa to synchronize the acoustic data with the biax data. "sa" stands for "synchronize acoustics". Display the Sync channel to pick the triggers. The output of the sa file is a mat called "p4640_sync_run1.mat" (one mat file per acoustic run).

Step 2: Write a script pXXXX_a similar to p4640_a to display waveforms (ShowMeWFs function) and retrieves TimeShift, Amplitudes, etc (ProcessAc_Tomo). 


