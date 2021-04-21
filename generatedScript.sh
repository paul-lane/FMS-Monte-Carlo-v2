cp Batch016/Inputs/*.* Inputs/
mkdir Batch016/Outputs
mkdir Batch016/Outputs/Processing1
mkdir Batch016/Outputs/Speed
 ./MCScattering_FM.out
 ./TimeSliceBinning.out
 ./StepSizeBinning.out
 ./Absorption.out
mv Outputs/AppProfile.txt Batch016/Outputs
mv Outputs/Ingoing*.txt Batch016/Outputs
mv Outputs/Scattered*.txt Batch016/Outputs
mv Outputs/Statistics.txt Batch016/Outputs
mv Outputs/WheelProfile.txt Batch016/Outputs
mv Outputs/Abs*.txt Batch016/Outputs/
mv Outputs/FMabs*.txt Batch016/Outputs/
mv Outputs/Processing1/Output*.txt Batch016/Outputs/Processing1/
mv Outputs/Speed/*.* Batch016/Outputs/Speed
cp Batch017/Inputs/*.* Inputs/
mkdir Batch017/Outputs
mkdir Batch017/Outputs/Processing1
mkdir Batch017/Outputs/Speed
 ./MCScattering_FM.out
 ./TimeSliceBinning.out
 ./StepSizeBinning.out
 ./Absorption.out
mv Outputs/AppProfile.txt Batch017/Outputs
mv Outputs/Ingoing*.txt Batch017/Outputs
mv Outputs/Scattered*.txt Batch017/Outputs
mv Outputs/Statistics.txt Batch017/Outputs
mv Outputs/WheelProfile.txt Batch017/Outputs
mv Outputs/Abs*.txt Batch017/Outputs/
mv Outputs/FMabs*.txt Batch017/Outputs/
mv Outputs/Processing1/Output*.txt Batch017/Outputs/Processing1/
mv Outputs/Speed/*.* Batch017/Outputs/Speed
