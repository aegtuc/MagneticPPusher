program ppusher
use constantsModule
use subsfuncsModule
implicit none
!----------------------------------------------------
! Programmed by AHMAD EDUARDO GUENNAM
! Departamento de Ciencias de la Computación
! Facultad de Ciencias Exactas y Tecnología (FACET)
! Universidad Nacional de Tucumán - Argentina
! eguennam@herrera.unt.edu.ar
!----------------------------------------------------
type(pusherParams)   :: pp
character(len=12)    :: REAL_CLOCK(6)
INTEGER              :: D_T(8)
character(128)       :: inputfile
integer              :: err
character(len=64)    :: date_time_str
!....................................................
inputfile = 'model.pshr'
CALL DATE_AND_TIME (REAL_CLOCK(1), REAL_CLOCK(2), REAL_CLOCK(3), D_T)
write(date_time_str,'(i2,"/",i2,"/",i4," AT ",i2,":",i2,":",i2)') D_T(3),D_T(2),D_T(1),D_T(5),D_T(6),D_T(7)
call writelog(text="====================================================================", ovrwrt=1)
call writelog(text="**  Relativistic Pusher of Charged Particles in Magnetic Fields   **", ovrwrt=0)
call writelog(text="** Ahmad Eduardo Guennam   -   email: eguennam@herrera.unt.edu.ar **", ovrwrt=0)
call writelog(text="====================================================================", ovrwrt=0)
call WriteLog('',0)
call readInp(inputfile=trim(inputfile),pp=pp)
call calculateOtherPusherParams(pp=pp)
call listParams(pp=pp)
call WriteLog('ANALYSIS STARTS: '//trim(date_time_str),  ovrwrt=0)
call relativisticPusher(pp=pp)       
CALL DATE_AND_TIME (REAL_CLOCK(1), REAL_CLOCK(2), REAL_CLOCK(3), D_T)
write(date_time_str,'(i2,"/",i2,"/",i4," AT ",i2,":",i2,":",i2)') D_T(3),D_T(2),D_T(1),D_T(5),D_T(6),D_T(7)
call WriteLog('ANALYSIS SUCCESSFULLY COMPLETED: '//trim(date_time_str),  ovrwrt=0)
end program ppusher