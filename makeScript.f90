Program makeScript

Implicit none

Integer :: i, start, end

Write(*,*) "Enter first batch number (1 - 999)"
Read(*,*) start
Write(*,*) "Enter last batch number (1 - 999)"
Read(*,*) end

Open(unit=10,file='generatedScript.sh')
Do i = start, end
write(10,'(a, I3.3, a)') "cp Batch",i,"/Inputs/*.* Inputs/"
write(10,'(a, I3.3, a)') "mkdir Batch",i, "/Outputs"
write(10,'(a, I3.3, a)') "mkdir Batch",i, "/Outputs/Processing1"
write(10,'(a, I3.3, a)') "mkdir Batch",i, "/Outputs/Speed"
write(10,*) "./MCScattering_FM.out"
write(10,*) "./TimeSliceBinning.out"
write(10,*) "./StepSizeBinning.out"
write(10,*) "./Absorption.out"
write(10,'(a, I3.3, a)') "mv Outputs/AppProfile.txt Batch",i,"/Outputs"
write(10,'(a, I3.3, a)') "mv Outputs/Ingoing*.txt Batch",i,"/Outputs"
write(10,'(a, I3.3, a)') "mv Outputs/Scattered*.txt Batch",i,"/Outputs"
write(10,'(a, I3.3, a)') "mv Outputs/Statistics.txt Batch",i,"/Outputs"
write(10,'(a, I3.3, a)') "mv Outputs/WheelProfile.txt Batch",i,"/Outputs"
write(10,'(a, I3.3, a)') "mv Outputs/MissedWheel.txt Batch",i,"/Outputs"
write(10,'(a, I3.3, a)') "mv Outputs/Abs*.txt Batch",i,"/Outputs/"
write(10,'(a, I3.3, a)') "mv Outputs/FMabs*.txt Batch",i,"/Outputs/"
write(10,'(a, I3.3, a)') "mv Outputs/SpeedDistr*.txt Batch",i,"/Outputs/"
write(10,'(a, I3.3, a)') "mv Outputs/Processing1/Output*.txt Batch",i,"/Outputs/Processing1/"
write(10,'(a, I3.3, a)') "mv Outputs/Speed/*.* Batch",i,"/Outputs/Speed"
End do

Write(*,*) "Script generated"

End Program makeScript
