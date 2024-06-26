module constantsModule
implicit none
! Programmed by AHMAD EDUARDO GUENNAM
! Departamento de Ciencias de la Computación
! Facultad de Ciencias Exactas y Tecnología (FACET)
! Universidad Nacional de Tucumán - Argentina
! eguennam@herrera.unt.edu.ar
character(32), parameter :: logfile   = 'pusher.log'

real(8), parameter :: dp5=0.5D0, d4p=4.0D0, d2p=2.0D0, d1p=1.0D0, d0p=0.0D0, d5p=5.0D0, d180p=180.0D0

real(8), parameter :: pi = 3.141592653589793D0

real(8), parameter :: AU2m = 150.0D9, m2AU = d1p/150.0D9

real(8), parameter :: c_AUs = 300000.0D0 / 150.0D6 ! UA/s 

real(8), parameter :: c_ms = 300.0D6 ! m/s 

real(8), parameter :: eV2J = 1.602177D-19

real(8), parameter :: deg2rad = 3.141592653589793D0 / 180.0D0


type pusherParams
    character(128) :: title
    character(256) :: outputFileName
    character(32)  :: fieldShape
    real(8)        :: frLength, frRadius
    real(8)        :: dt,tfinal,mass, charge,Energy_eV,alfa,H,dir0(3),xyz0(3)
    real(8)        :: direcc, azimut
    real(8)        :: alpha0, R0,z0
    integer        :: every
    logical        :: addBFldNoise
    real(8)        :: dirBsw(3) !dirBswX, dirBswY, dirBswZ
    real(8)        :: BswMean, BswNoiseStd, BfrMean, BfrNoiseStd, BshMean, BshNoiseStd
    real(8)        :: focalDist
    !Calculated based on input data
    integer        :: ninc,nout
    
end type pusherParams

contains
end module constantsModule