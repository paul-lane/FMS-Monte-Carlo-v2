gfortran -ffree-line-length-none -o TimeSliceBinning.out TimeSliceBinning.f90
gfortran -ffree-line-length-none -o LorentzianLaser.out LorentzianLaser.f90
gfortran -ffree-line-length-none -o Absorption_v2.out Absorption_v2.f90
./TimeSliceBinning.out
./LorentzianLaser.out
./Absorption_v2.out
