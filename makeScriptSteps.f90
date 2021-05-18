Program makeScript

Implicit none

Integer :: i, start, end

Write(*,*) "Enter first batch number (1 - 999)"
Read(*,*) start
Write(*,*) "Enter last batch number (1 - 999)"
Read(*,*) end

Open(unit=10,file='ScriptSteps.bat')
Do i = start, end
write(10,'(a, I3.3, a)') "Xcopy /y Batch",i,"\Inputs\binParameters.inp Inputs\"
write(10,'(a, I3.3, a)') "Xcopy /y Batch",i,"\Inputs\chamberDimensions.inp Inputs\"
write(10,'(a, I3.3, a)') "Xcopy /y Batch",i,"\Inputs\directionInputs_v2.inp Inputs\"
write(10,'(a, I3.3, a)') "Xcopy /y Batch",i,"\Inputs\fmInputs.inp Inputs\"
write(10,'(a, I3.3, a)') "Xcopy /y Batch",i,"\Inputs\herriottCellGeom.inp Inputs\"
write(10,'(a, I3.3, a)') "Xcopy /y Batch",i,"\Inputs\mirrorSpots.inp Inputs\"
write(10,'(a, I3.3, a)') "Xcopy /y Batch",i,"\Inputs\overwrite.inp Inputs\"
write(10,'(a, I3.3, a)') "Xcopy /y Batch",i,"\Inputs\timingInputs.inp Inputs\"
write(10,'(a, I3.3, a)') "mkdir Batch",i, "\Outputs"
write(10,'(a, I3.3, a)') "mkdir Batch",i, "\Outputs\Processing1"
write(10,'(a, I3.3, a)') "mkdir Batch",i, "\Outputs\Speed"
write(10,*) ".\MCScattering_FM_v2.exe"
write(10,*) "timeout /t 10"
write(10,*) ".\TimeSliceBinning.exe"
write(10,*) "timeout /t 10"
write(10,*) ".\PickingSteps.exe"
write(10,*) "timeout /t 10"
write(10,*) ".\Absorption_v2.exe"
write(10,*) "timeout /t 10"
write(10,'(a, I3.3, a)') "move Outputs\AppProfile.txt Batch",i,"\Outputs"
write(10,'(a, I3.3, a)') "move Outputs\IngoingCellProfile.txt Batch",i,"\Outputs"
write(10,'(a, I3.3, a)') "move Outputs\IngoingCollimatorProfile.txt Batch",i,"\Outputs"
write(10,'(a, I3.3, a)') "move Outputs\IngoingProbed.txt Batch",i,"\Outputs"
write(10,'(a, I3.3, a)') "move Outputs\MissedWheel.txt Batch",i,"\Outputs"
write(10,'(a, I3.3, a)') "move Outputs\ScatteredCellProfile.txt Batch",i,"\Outputs"
write(10,'(a, I3.3, a)') "move Outputs\ScatteredProbed.txt Batch",i,"\Outputs"
write(10,'(a, I3.3, a)') "move Outputs\Statistics.txt Batch",i,"\Outputs"
write(10,'(a, I3.3, a)') "move Outputs\WheelProfile.txt Batch",i,"\Outputs"
write(10,'(a, I3.3, a)') "move Outputs\Abs_001.txt Batch",i,"\Outputs\"
write(10,'(a, I3.3, a)') "move Outputs\Abs_002.txt Batch",i,"\Outputs\"
write(10,'(a, I3.3, a)') "move Outputs\Abs_003.txt Batch",i,"\Outputs\"
write(10,'(a, I3.3, a)') "move Outputs\Abs_004.txt Batch",i,"\Outputs\"
write(10,'(a, I3.3, a)') "move Outputs\FMabs_001.txt Batch",i,"\Outputs\"
write(10,'(a, I3.3, a)') "move Outputs\FMabs_002.txt Batch",i,"\Outputs\"
write(10,'(a, I3.3, a)') "move Outputs\FMabs_003.txt Batch",i,"\Outputs\"
write(10,'(a, I3.3, a)') "move Outputs\FMabs_004.txt Batch",i,"\Outputs\"
write(10,'(a, I3.3, a)') "move Outputs\SpeedDistr_001.txt Batch",i,"\Outputs\"
write(10,'(a, I3.3, a)') "move Outputs\SpeedDistr_002.txt Batch",i,"\Outputs\"
write(10,'(a, I3.3, a)') "move Outputs\SpeedDistr_003.txt Batch",i,"\Outputs\"
write(10,'(a, I3.3, a)') "move Outputs\SpeedDistr_004.txt Batch",i,"\Outputs\"
write(10,'(a, I3.3, a)') "move Outputs\Processing1\Output_001.txt Batch",i,"\Outputs\Processing1\"
write(10,'(a, I3.3, a)') "move Outputs\Processing1\Output_002.txt Batch",i,"\Outputs\Processing1\"
write(10,'(a, I3.3, a)') "move Outputs\Processing1\Output_003.txt Batch",i,"\Outputs\Processing1\"
write(10,'(a, I3.3, a)') "move Outputs\Processing1\Output_004.txt Batch",i,"\Outputs\Processing1\"
write(10,'(a, I3.3, a)') "move Outputs\Speed\*.* Batch",i,"\Outputs\Speed"
End do

Write(*,*) "Script generated"

End Program makeScript
