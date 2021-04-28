Program getHerriottCell
! Program to calculate beam positions at centre of H-Cell and Euler angles from mirrorSpots.inp

Implicit none

Double Precision, Allocatable :: x(:), y(:), z(:) 
Double Precision, Allocatable :: my(:), mz(:), th1(:), th2(:), ypos(:), zpos(:)
Double Precision :: pi
Integer :: i, nSpot, nPass
Logical :: flip

pi = 3.14159265359


Open(unit=10,file='mirrorSpots.inp')
Read(10,*) nSpot
Read(10,*) flip  
Read(10,*) 

nPass = nSpot-1

Open(unit=11,file='herriottCellGeom.inp')

Allocate(x(nSpot), y(nSpot), z(nSpot))
Allocate(my(nPass), mz(nPass), th1(nPass), th2(nPass), ypos(nPass), zpos(nPass))

Do i = 1, nSpot					! for number of spots
	Read(10,*) x(i), y(i), z(i)		! read in spot position in x, y & z
	If(flip) then				! if flip = true then 
		if(Mod(i,2) .eq. 0) then	! for even values of i flip the sign of z component
		z(i) = -1.0*z(i)
		end if
	End if
	
	x(i) = 1e-3*x(i)
	y(i) = 1e-3*y(i)
	z(i) = 1e-3*z(i)

End do
Close(10)

Write(11,*) nspot-1
Write(11,*) 0.001				! hardwired laser beam radius
Write(11,*) " yposition, zposition (+ve is closer to valve), theta1 deg (side to side), theta2 deg (up/down)" 
Do i = 1, nspot-1				! number of passes is 1 less than number of spots
	mz(i) = (z(i+1)-z(i))/(x(i+1)-x(i))	! calculate the gradient in y and z directions
	my(i) = (y(i+1)-y(i))/(x(i+1)-x(i))	! to enable calculation of Euler angles
	
	th1(i) = ATAN(mz(i))
	th2(i) = ATAN(my(i))

	th1(i) = (th1(i)/pi)*180.0
	th2(i) = (th2(i)/pi)*180.0

!	ypos(i) = y(i)-my(i)*x(i)
!	zpos(i) = z(i)-mz(i)*x(i)

	ypos(i) = (y(i+1)+y(i))/2.0		! The positions at the centre of the cell are the average of the positions at either end of the cell
	zpos(i) = (z(i+1)+z(i))/2.0

	Write(11,*) ypos(i), zpos(i), th1(i), th2(i)
End do
Close(11)
	
End Program
