# COVID Project

Adaptation of clustering COVID time-series trajectories using Dynamic Time Warping (https://www.r-bloggers.com/coronadash-app-use-case-clustering-countries-covid-19-active-cases-trajectories/). This method splits UK UA trajectories into 14 day subtrajectories and then clusters the differential of estimated current cases in that period.

The second project uses the same 14 day trajectories and applies a LSTM CNN to predict the 15th day of activity. 
