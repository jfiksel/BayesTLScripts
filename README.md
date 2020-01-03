### Recreating results

* Figure 1 can be recreated by running  `phmrc_eda/error-matrix-all-observations.R`

* All simulation results can be obtained by first running `fakeCOD.R` in `SimulationStudy`. Note that this file was run as an SGE array job using the `fakeCOD.sh` script, and thus would have to be modified accordingly (i.e. using a for loop) if you do not have access to a HPC environment. Simulation results are collected using `collectResults.R` and then plotted with `plotResults.R`

* All files for recreating figures and results for the PHMRC data analysis is in the `CSMFPredictionsWeightedResampling` directory. First, run `trainModelsChild.R`. Then run `generateCalibrationIndicesChild.R` via the shell script `generateCalibrationIndicesChild.sh`. The main simulations are done in `runSimulationChild.R` (again run using the associated shell script). Simulations for Section S7 are in `SimulationChildExpandedCauses.R`. All simulation results are then collected using `collectResults.R` and plotted in `plotResults.R`. The file `plotSympDifferences.R` is used for creating Figure S8. 