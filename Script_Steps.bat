Xcopy /y Batch001\Inputs\binParameters.inp Inputs\
Xcopy /y Batch001\Inputs\chamberDimensions.inp Inputs\
Xcopy /y Batch001\Inputs\directionInputs_v2.inp Inputs\
Xcopy /y Batch001\Inputs\fmInputs.inp Inputs\
Xcopy /y Batch001\Inputs\herriottCellGeom.inp Inputs\
Xcopy /y Batch001\Inputs\mirrorSpots.inp Inputs\
Xcopy /y Batch001\Inputs\overwrite.inp Inputs\
Xcopy /y Batch001\Inputs\timingInputs.inp Inputs\
mkdir Batch001\Outputs
mkdir Batch001\Outputs\Processing1
mkdir Batch001\Outputs\Speed
.\MCScattering_FM_v2.exe
timeout /t 10
.\TimeSliceBinning.exe
timeout /t 10
.\PickingSteps.exe
timeout /t 10
.\Absorption_v2.exe
timeout /t 10
move Outputs\AppProfile.txt Batch001\Outputs
move Outputs\IngoingCellProfile.txt Batch001\Outputs
move Outputs\IngoingCollimatorProfile.txt Batch001\Outputs
move Outputs\IngoingProbed.txt Batch001\Outputs
move Outputs\MissedWheel.txt Batch001\Outputs
move Outputs\ScatteredCellProfile.txt Batch001\Outputs
move Outputs\ScatteredProbed.txt Batch001\Outputs
move Outputs\Statistics.txt Batch001\Outputs
move Outputs\WheelProfile.txt Batch001\Outputs
move Outputs\Abs_001.txt Batch001\Outputs\
move Outputs\Abs_002.txt Batch001\Outputs\
move Outputs\Abs_003.txt Batch001\Outputs\
move Outputs\Abs_004.txt Batch001\Outputs\
move Outputs\FMabs_001.txt Batch001\Outputs\
move Outputs\FMabs_002.txt Batch001\Outputs\
move Outputs\FMabs_003.txt Batch001\Outputs\
move Outputs\FMabs_004.txt Batch001\Outputs\
move Outputs\SpeedDistr_001.txt Batch001\Outputs\
move Outputs\SpeedDistr_002.txt Batch001\Outputs\
move Outputs\SpeedDistr_003.txt Batch001\Outputs\
move Outputs\SpeedDistr_004.txt Batch001\Outputs\
move Outputs\Processing1\Output_001.txt Batch001\Outputs\Processing1\
move Outputs\Processing1\Output_002.txt Batch001\Outputs\Processing1\
move Outputs\Processing1\Output_003.txt Batch001\Outputs\Processing1\
move Outputs\Processing1\Output_004.txt Batch001\Outputs\Processing1\
move Outputs\Speed\*.* Batch001\Outputs\Speed\
