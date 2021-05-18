gfortran -ffree-line-length-none -o TimeSliceBinning.out TimeSliceBinning.f90
gfortran -ffree-line-length-none -o PickingSteps.out PickingSteps.f90
gfortran -ffree-line-length-none -o Absorption_v2.out Absorption_v2.f90
./TimeSliceBinning.out
./PickingSteps.out
./Absorption_v2.out
