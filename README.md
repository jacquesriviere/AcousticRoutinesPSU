# AcousticRoutinesPSU

The main program is p4473_a. 
It uses SyncAcData to synchronize acoustic data from Verasonics with mechanical data from the biax recorder.
Then, it runs ProcessAcData which computes the change in time of flight (timeshift = TS), rms amplitude (RMS) and maximum of intercorrelation (MI) using the first 50 waveforms as a reference (first 50 waveforms are stacked and used as a pattern).

After running the main program, you'll get a time vector for acoustic data with corresponding TS, RMS and MI.
The next step will be to "merge" the two time vectors (mechanical and acoustical) to be able to plot, for instance TS versus slip (I started to write a routine for that...coming soon).
