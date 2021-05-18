module mathConstants

    ! General use mathematical constants for use in all subroutines

    double precision, parameter :: pi = 3.141592653589793D0			! Pi
    double precision, parameter :: boltzmannConstant = 1.38064852D-23		! Boltzmann Constant
    double precision, parameter :: avagadro = 6.0221409D23			! Avagadro's number
    double precision, parameter :: speedLight =  299792458.0 	                ! Speed of light in ms-1
 

    public :: pi, boltzmannConstant, avagadro, speedLight

end module mathConstants
