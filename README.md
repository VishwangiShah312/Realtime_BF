The RTB.cpp computes the hour angle and phases for the realtime beamformer. The unix time and RA and DEC of the source wrt the mean equinox is manually
input to the code. The RA and DEC of the source wrt the mean equinox are computed by the j2000tonow code. The unix time is taken directly from the 
'date_unix' returned by the Off_BF_phase_gen code. 

The j2000tonow code has functions taken directly from the Tracking Beam Scheduler which give the RA and DEC of the source wrt the equinox of the time of 
observation. I compute the RA and DEC from here and input it manually into the RTB.cpp.

The Off_BF_phase_gen code contains parts of the code from baseband_analysis which compute the hour angle and phases. The 'date_unix' returned in this code 
is maually input into the RTB.cpp to find the hour angle for the same time. This code needs to run in a baseband-analysis docker container on frb-analysis.
