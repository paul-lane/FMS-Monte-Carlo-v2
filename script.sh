cp Batch001/Inputs/binParameters.inp Inputs/
cp Batch001/Inputs/chamberDimensions.inp Inputs/
cp Batch001/Inputs/directionInputs.inp Inputs/
cp Batch001/Inputs/fmInputs.inp Inputs/
cp Batch001/Inputs/herriottCellGeom.inp Inputs/
cp Batch001/Inputs/mirrorSpots.inp Inputs/
cp Batch001/Inputs/overwrite.inp Inputs/
cp Batch001/Inputs/timingInputs.inp Inputs/
mkdir Batch001/Outputs
mkdir Batch001/Outputs/Processing1
mkdir Batch001/Outputs/Speed
./MCScattering_FM.out
./TimeSliceBinning.out
./StepSizeBinning.out
./Absorption.out
mv Outputs/AppProfile.txt Batch001/Outputs
mv Outputs/IngoingCellProfile.txt Batch001/Outputs
mv Outputs/IngoingCollimatorProfile.txt Batch001/Outputs
mv Outputs/IngoingProbed.txt Batch001/Outputs
mv Outputs/MissedWheel.txt Batch001/Outputs
mv Outputs/ScatteredCellProfile.txt Batch001/Outputs
mv Outputs/ScatteredProbed.txt Batch001/Outputs
mv Outputs/Statistics.txt Batch001/Outputs
mv Outputs/WheelProfile.txt Batch001/Outputs
mv Outputs/Abs_001.txt Batch001/Outputs/
mv Outputs/Abs_002.txt Batch001/Outputs/
mv Outputs/Abs_003.txt Batch001/Outputs/
mv Outputs/Abs_004.txt Batch001/Outputs/
mv Outputs/FMabs_001.txt Batch001/Outputs/
mv Outputs/FMabs_002.txt Batch001/Outputs/
mv Outputs/FMabs_003.txt Batch001/Outputs/
mv Outputs/FMabs_004.txt Batch001/Outputs/
mv Outputs/SpeedDistr_001.txt Batch001/Outputs/
mv Outputs/SpeedDistr_002.txt Batch001/Outputs/
mv Outputs/SpeedDistr_003.txt Batch001/Outputs/
mv Outputs/SpeedDistr_004.txt Batch001/Outputs/
mv Outputs/Processing1/Output_001.txt Batch001/Outputs/Processing1/
mv Outputs/Processing1/Output_002.txt Batch001/Outputs/Processing1/
mv Outputs/Processing1/Output_003.txt Batch001/Outputs/Processing1/
mv Outputs/Processing1/Output_004.txt Batch001/Outputs/Processing1/
mv Outputs/Speed/*.* Batch001/Outputs/Speed/


