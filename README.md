# AcousticRoutinesPSU

Step 1: Write a script pXXXX_sa similar to p4640_sa to synchronize the acoustic data with the biax data. "sa" stands for "synchronize acoustics". Display the Sync channel to pick the triggers. The output of the sa file is a mat called "p4640_sync_run1.mat" (one mat file per acoustic run).

Step 2: Write a script pXXXX_a similar to p4640_a to display waveforms (ShowMeWFs function) and retrieves TimeShift, Amplitudes, etc (ProcessAc_Tomo). 

Step 3: Depending on your experiment, you might want to use different options to compute the cross-correlation and retrieve TimeShift. Look at p4640_comparison_abs_rel_mix.m as an example.
