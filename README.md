# AcousticRoutinesPSU

# The main program is p4473_a. 
# It uses SyncAcData to synchronize acoustic data from Verasonics with mechanical data from the biax recorder.
# Then, it runs ProcessAcData which compute the change in time of flight, rms amplitude and maximum of intercorrelation using the first 50 waveforms as a reference (first 50 waveforms are stacked and used as a pattern).
