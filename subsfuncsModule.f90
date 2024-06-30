module subsfuncsModule
use constantsModule
implicit none
! Programmed by AHMAD EDUARDO GUENNAM
! Departamento de Ciencias de la Computación
! Facultad de Ciencias Exactas y Tecnología (FACET)
! Universidad Nacional de Tucumán - Argentina
! eguennam@herrera.unt.edu.ar
contains
    
   


!**************************************************************************************************************
!
!                                          F U N C T I O N:  velocityFromEnergy
!                                          
!**************************************************************************************************************
function velocityFromEnergy(mass, energy)
implicit none
real(8), intent(in) :: mass, energy
real(8) :: velocityFromEnergy
!energy is the total energy:
velocityFromEnergy = c_ms * dsqrt(d1p - ( (mass* c_ms**2) / energy)**2  )
return
end function velocityFromEnergy
!**************************************************************************************************************
!
!                                          F U N C T I O N:  lorentzFactor
!                                           gamma = 1/sqrt( 1 - (v/c)^2 )
!**************************************************************************************************************
function lorentzFactor(v)
implicit none
real(8), intent(in) :: v !Three velovity, u=gamma*v; 
real(8) :: lorentzFactor
lorentzFactor =  d1p / dsqrt(d1p - (v/c_ms)**2)
return
end function lorentzFactor

!**************************************************************************************************************
!
!                                          F U N C T I O N:  cross
!                                          
!**************************************************************************************************************
function cross(a,b)
implicit none
real(8), dimension(3), intent(in) :: a, b
real(8), dimension(3) :: cross
!e1 e2 e3
!a1 a2 a3
!b1 b2 b3
cross(1) =  a(2)*b(3) - a(3)*b(2)
cross(2) = -a(1)*b(3) + a(3)*b(1)
cross(3) =  a(1)*b(2) - a(2)*b(1) 
return
end function cross
!**************************************************************************************************************
!
!                                          S U B R U T I N E: informAdvance
!
!**************************************************************************************************************
logical function informAdvance(ntotal, ncurrent, ns) 
integer, intent(in) :: ntotal, ncurrent, ns
integer :: myntotal
informAdvance = .false.
myntotal = ntotal
if ( mod(ntotal,ns) .ne. 0) myntotal = ntotal - mod(ntotal,ns)
if (myntotal > ns)then
    if (mod(ncurrent,int(myntotal/ns)) .eq. 0) then
        informAdvance = .true. 
    endif    
endif
return
end function informAdvance
!**************************************************************************************************************
!
!                                  S U B R U T I N E: initializePusherParams
!
!**************************************************************************************************************
subroutine initializePusherParams(pp)
implicit none
type(pusherParams),intent(out)    :: pp
pp%title    = ''
pp%frLength = d0p
pp%ninc = 0
pp%dt = d0p
pp%tfinal = d0p
pp%xyz0 = d0p
pp%dir0 = d0p
pp%fieldShape = ''
pp%outputFileName = 'output.csv'
pp%direcc  = d0p
pp%azimut   = d0p
pp%addBFldNoise = .false.
pp%dirBsw = d0p  !dirBswX, dirBswY, dirBswZ
pp%BswMean=d0p
pp%BswNoiseStd=d0p
pp%BfrMean=d0p
pp%BfrNoiseStd=d0p
pp%BshMean=d0p
pp%BshNoiseStd=d0p
return
end subroutine


!**************************************************************************************************************
!
!                                  S U B R U T I N E: calculateOtherPusherParams
!
!**************************************************************************************************************
subroutine calculateOtherPusherParams(pp)
implicit none
type(pusherParams),intent(inout)    :: pp
!------------------------------------------------------------
real(8) :: phi, psi, Rphi(3,3), Rpsi(3,3), v1(3), aux3(3)
!------------------------------------------------------------
pp%ninc = int(pp%tfinal/pp%dt)

pp%nout = int(dfloat(pp%ninc)/dfloat(pp%every))

!Initial Velocity

phi = (pp%alpha0 + d180p + pp%direcc) * deg2rad
psi = -pp%azimut * deg2rad
v1 = (/d1p,d0p,d0p/)
Rphi(1,1) = dcos(phi)
Rphi(1,2) = -dsin(phi)
Rphi(1,3) = d0p
    
Rphi(2,1) = dsin(phi)
Rphi(2,2) = dcos(phi)
Rphi(2,3) = d0p
    
Rphi(3,1) = d0p
Rphi(3,2) = d0p
Rphi(3,3) = d1p
    
Rpsi(1,1) = dcos(psi)
Rpsi(1,2) = d0p
Rpsi(1,3) = dsin(psi)
    
Rpsi(2,1) = d0p
Rpsi(2,2) = d1p
Rpsi(2,3) = d0p
    
Rpsi(3,1) = -dsin(psi)
Rpsi(3,2) = d0p
Rpsi(3,3) = dcos(psi)
   
aux3 = matmul(Rpsi,v1)
    
pp%dir0 = matmul(Rphi,aux3)
  
!Initial Position:
pp%xyz0(1) =  pp%R0*dcos(pp%alpha0 * deg2rad)
pp%xyz0(2) =  pp%R0*dsin(pp%alpha0 * deg2rad)
pp%xyz0(3) =  pp%z0

return
end subroutine calculateOtherPusherParams
!**************************************************************************************************************
!
!                                  S U B R U T I N E: listParams
!
!**************************************************************************************************************
subroutine listParams(pp)
implicit none
type(pusherParams),intent(inout)    :: pp

call writelog(text="--------------------------------------------------", ovrwrt=0)
call writelog(text="           LISTING INPUT PARAMETERS", ovrwrt=0)
call writelog(text="--------------------------------------------------", ovrwrt=0)
call writelog(text="title                : "//trim(pp%title), ovrwrt=0)
call writelog(text="-----------------------------------------", ovrwrt=0)
call writelog(text="fluxRopeLength       : "//d2str(pp%frLength), ovrwrt=0)
call writelog(text="fluxRopeRadius       : "//d2str(pp%frRadius), ovrwrt=0)
call writelog(text="-----------------------------------------", ovrwrt=0)
call writelog(text="Energy               : "//d2str(pp%Energy_eV), ovrwrt=0)
call writelog(text="Mass                 : "//d2str(pp%mass), ovrwrt=0)
call writelog(text="Charge               : "//d2str(pp%charge), ovrwrt=0)
call writelog(text="-----------------------------------------", ovrwrt=0)
call writelog(text="",ovrwrt=0)
call writelog(text="Initial Position: ", ovrwrt=0)
call writelog(text="R0                   : "//d2str(pp%R0), ovrwrt=0)
call writelog(text="alpha0               : "//d2str(pp%alpha0), ovrwrt=0)
call writelog(text="",ovrwrt=0)
call writelog(text="xyz0x                : "//d2str(pp%xyz0(1)), ovrwrt=0)
call writelog(text="xyz0y                : "//d2str(pp%xyz0(2)), ovrwrt=0)
call writelog(text="xyz0z                : "//d2str(pp%xyz0(3)), ovrwrt=0)
call writelog(text="",ovrwrt=0)
call writelog(text="-----------------------------------------", ovrwrt=0)
call writelog(text="",ovrwrt=0)
call writelog(text="Initial Velocity: ", ovrwrt=0)
call writelog(text="",ovrwrt=0)
call writelog(text="direcc               : "//d2str(pp%direcc), ovrwrt=0)
call writelog(text="azimut               : "//d2str(pp%azimut), ovrwrt=0)
call writelog(text="",ovrwrt=0)
call writelog(text="dir0x                : "//d2str(pp%dir0(1)), ovrwrt=0)
call writelog(text="dir0y                : "//d2str(pp%dir0(2)), ovrwrt=0)
call writelog(text="dir0z                : "//d2str(pp%dir0(3)), ovrwrt=0)
call writelog(text="-----------------------------------------", ovrwrt=0)
call writelog(text="BfrMean              : "//d2str(pp%BfrMean), ovrwrt=0)
call writelog(text="BswMean              : "//d2str(pp%BswMean), ovrwrt=0)
call writelog(text="alfa                 : "//d2str(pp%alfa), ovrwrt=0)
call writelog(text="H                    : "//d2str(pp%H), ovrwrt=0)
call writelog(text="Field Shape          : "//trim(pp%fieldShape), ovrwrt=0)
call writelog(text="-----------------------------------------", ovrwrt=0)
call writelog(text="tfinal               : "//d2str(pp%tfinal), ovrwrt=0)
call writelog(text="dt                   : "//d2str(pp%dt), ovrwrt=0)
call writelog(text="nicn                 : "//i2str(pp%ninc), ovrwrt=0)
call writelog(text="every                : "//i2str(pp%every), ovrwrt=0)
call writelog(text="nout                 : "//i2str(pp%nout), ovrwrt=0)
call writelog(text="", ovrwrt=0)
call writelog(text="outputFileName       : "//trim(pp%outputFileName), ovrwrt=0)
call writelog(text="",ovrwrt=0)
call writelog(text="addBFldNoise         : "//b2str(pp%addBFldNoise), ovrwrt=0)


call writelog(text="",ovrwrt=0)
call writelog(text="BswNoiseStd          : "//d2str(pp%BswNoiseStd), ovrwrt=0)
call writelog(text="BshNoiseStd          : "//d2str(pp%BshNoiseStd), ovrwrt=0)
call writelog(text="BfrNoiseStd          : "//d2str(pp%BfrNoiseStd), ovrwrt=0)
call writelog(text="",ovrwrt=0)
call writelog(text="BswMean              : "//d2str(pp%BswMean), ovrwrt=0)
call writelog(text="BshMean              : "//d2str(pp%BshMean), ovrwrt=0)
call writelog(text="BfrMean              : "//d2str(pp%BfrMean), ovrwrt=0)

call writelog(text="",ovrwrt=0)
call writelog(text="dirBswX              : "//d2str(pp%dirBsw(1)), ovrwrt=0)
call writelog(text="dirBswY              : "//d2str(pp%dirBsw(2)), ovrwrt=0)
call writelog(text="dirBswZ              : "//d2str(pp%dirBsw(3)), ovrwrt=0)


call writelog(text="",ovrwrt=0)
call writelog(text="focalDist            : "//d2str(pp%focalDist), ovrwrt=0)
call writelog(text="--------------------------------------------------------", ovrwrt=0)
return
end subroutine listParams


!**************************************************************************************************************
!
!                                          S U B R U T I N E: readInp
!
!**************************************************************************************************************
subroutine readInp(inputfile,pp)
!Read geomparam parameters from input file
character(*),intent(in)       :: inputfile
type(pusherParams),intent(out)    :: pp
!
character(256) :: varname,varval, line
integer :: i,j,ios,varind, i1
!
open(2,file=inputfile,action='read',err=1)
ios=0
read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
pp%title = trim(adjustl(line))

read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(line,*,iostat=ios)  pp%frLength

read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(line,*,iostat=ios)  pp%frRadius

read(2,'(a)',iostat=ios) line

read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
pp%outputFileName = trim(adjustl(line))

read(2,'(a)',iostat=ios) line

read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(line,*,iostat=ios)  pp%Energy_eV

read(2,'(a)',iostat=ios) line

read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(line,*,iostat=ios) pp%dt

read(2,'(a)',iostat=ios) line

read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(line,*,iostat=ios)  pp%tfinal

read(2,'(a)',iostat=ios) line

read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(line,*,iostat=ios)  pp%mass

read(2,'(a)',iostat=ios) line

read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(line,*,iostat=ios)  pp%charge

read(2,'(a)',iostat=ios) line

read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(line,*,iostat=ios)  pp%alfa

read(2,'(a)',iostat=ios) line

read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(line,*,iostat=ios)  pp%H

read(2,'(a)',iostat=ios) line

read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
pp%fieldShape=trim(adjustl(line))

read(2,'(a)',iostat=ios) line

read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(line,*,iostat=ios)  pp%BswMean, pp%BswNoiseStd

read(2,'(a)',iostat=ios) line

read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(line,*,iostat=ios)  pp%BfrMean, pp%BfrNoiseStd

read(2,'(a)',iostat=ios) line

read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(line,*,iostat=ios)  pp%BshMean , pp%BshNoiseStd

read(2,'(a)',iostat=ios) line

read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(line,*,iostat=ios)  pp%dirBsw(1), pp%dirBsw(2), pp%dirBsw(3)

read(2,'(a)',iostat=ios) line

read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(line,*,iostat=ios)  pp%focalDist

read(2,'(a)',iostat=ios) line

read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
pp%addBFldNoise=.false.
if (trim(adjustl(line)) == 'Y')then
  pp%addBFldNoise=.true.
endif

read(2,'(a)',iostat=ios) line

read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(line,*,iostat=ios)  pp%every

read(2,'(a)',iostat=ios) line

read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(line,*,iostat=ios)  pp%R0, pp%z0, pp%alpha0

read(2,'(a)',iostat=ios) line

read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(2,'(a)',iostat=ios) line
read(line,*,iostat=ios)  pp%direcc, pp%azimut

close(2)
return
1 write(*,*) "***ERROR: Unable to read Input file: "//trim(inputfile)//'!'
close(2)
return
end subroutine readInp
!**************************************************************************************************************
!
!                                       S U B R O U T I N E: WriteLog
!Credits: Juan Pablo Marquez
!**************************************************************************************************************
subroutine WriteLog(text,ovrwrt)
!Writes log file
character*(*),intent(in) :: text
integer,      intent(in) :: ovrwrt
!.............................................................
integer :: unit1
unit1=1
select case (ovrwrt)
	case (0)
		open(unit1,access='APPEND',file=logfile,action='write',err=1)
		write(unit1,'(a)') trim(text)
		close(unit1)
	case (1)
		open(unit1,file=logfile,action='write',err=1)
		write(unit1,'(a)') trim(text)
		close(unit1)
end select
return
1 close(unit1)
return
end subroutine WriteLog

!**************************************************************************************************************
!
!                                          F U N C T I O N: d2str
!Credits: Juan Pablo Marquez 
!**************************************************************************************************************
character(24) function d2str(num,fmt)
real(8)     , intent(in)           :: num
character(*), intent(in), optional :: fmt
integer                            :: len1
character(24) :: str1
!.............................................................
len1 = len(fmt)
if (present(fmt)) then      
    write(d2str,"("//fmt(1:len1)//")") num
else
    write(d2str,"(1PG24.15E3)") num 
endif
d2str = trim(adjustl(d2str))
return
end function d2str
!**************************************************************************************************************
!
!                                          F U N C T I O N: i2str
!Credits: Juan Pablo Marquez
!**************************************************************************************************************
character(12) function i2str(int)
integer,intent(in) :: int
!.............................................................
write(i2str,*) int
i2str = trim(adjustl(i2str))
return
end function i2str

!**************************************************************************************************************
!
!                                          F U N C T I O N: b2str
!
!**************************************************************************************************************
character(3) function b2str(boo)
logical,intent(in) :: boo
!.............................................................
b2str = 'NO'
if (boo .eqv. .true.) then
b2str = 'YES'
endif
return
end function b2str
!**************************************************************************************************************
!
!                                S U B R O U T I N E: relativisticPusher
!Reference:
!Ripperda, B., Bacchini, F., Teunissen, J., Xia, C., Porth, O., Sironi, L., Lapenta, G., Keppens, R. (2018, mar). 
!A comprehensive comparison of relativistic particle integrators. The Astrophysical Journal Supplement Series, 
!235 (1), 21.
!**************************************************************************************************************
subroutine relativisticPusher(pp)
implicit none
type(pusherParams),intent(inout)    :: pp
!........................................................................................
real(8)              :: Vel0module, LorFactor, qdto2m, tt, tmp1
integer(8)           :: i1, icont, isave, idim
real(8), allocatable :: r(:,:), v(:,:), time(:), BFldonPath(:,:)
real(8)              :: BFld(3) 
real(8)              :: xyz1(3), uvw1(3)   !time n
real(8)              :: xyz2(3), uvw2(3)   !time n+1
real(8)              :: xnhalf(3)          !time n + 1/2
real(8)              :: t(3)
real(8)              :: s(3) 
real(8)              :: ETIME(2) 
real(8), allocatable :: tempPath(:,:)
integer              :: err
!........................................................................................
Vel0module = velocityFromEnergy(mass=pp%mass, energy=pp%Energy_eV*eV2J)  * m2AU
!Memory Allocation:
allocate(time(pp%nout+1),r(pp%nout+1,3),v(pp%nout+1,3),BFldOnPath(pp%nout+1,3))
r = d0p
v = d0p
BFldOnPath = d0p 

r(1,1) = pp%xyz0(1) !x0
r(1,2) = pp%xyz0(2) !y0
r(1,3) = pp%xyz0(3) !z0

v(1,1) = Vel0module * pp%dir0(1)
v(1,2) = Vel0module * pp%dir0(2)
v(1,3) = Vel0module * pp%dir0(3)

BFld = (/d0p,d0p, d0p/)

call cpu_time(ETIME(1))

time(1) = d0p
xyz1 = r(1,:)   !xn
xyz2 = d0p      !xn+1

uvw1 = v(1,:)   !vn
uvw2 = d0p      !vn+1

icont = 1
isave = 1

write(*,'(a)') '0%--+---+---+---+---+---+---+---+----100%'
write(*,'(a)',advance="no") '|'
LorFactor = lorentzFactor(v=dsqrt(uvw1(1)**2 + uvw1(2)**2 + uvw1(3)**2)*AU2m)
qdto2m = pp%charge / (pp%mass * LorFactor)  * dp5 * pp%dt
do i1=2,pp%ninc

    xnhalf = xyz1 + uvw1 * pp%dt * dp5 
    
    call magneticField(x=xnhalf(1),y=xnhalf(2),z=xnhalf(3),pp=pp,BFld=BFld,err=err)
    
    !Riperda et al, 2018.
    t = qdto2m * BFld
    tt = t(1)**2 + t(2)**2 + t(3)**2 
    tmp1 = d2p  / (d1p + tt) 
    
    s = tmp1 * t
    uvw2 = uvw1 + cross(a=(uvw1+cross(a=uvw1,b=t)),b=s)
    
    xyz2 = xnhalf + uvw2 * pp%dt * dp5   
    
    
    icont = icont + 1
    
    if (icont == pp%every)then
        isave = isave + 1 
        r(isave,:) = xyz2
        v(isave,:) = uvw2
        time(isave) = pp%dt*dfloat(i1)
        
        BFldonPath(isave,1) = BFld(1)
        BFldonPath(isave,2) = BFld(2)
        BFldonPath(isave,3) = BFld(3) 
    
        icont = 0
        
        if ( informAdvance(ntotal=pp%nout, ncurrent=int(isave), ns=10) .eqv. .true.)then
            write(*,'(a)',advance="no") '###|'
        endif    
        
    endif
    xyz1 = xyz2
    uvw1 = uvw2
    
enddo

call cpu_time(ETIME(2))



write(*,'(/,a)') '-----------------------------------------'
write(*,'(a,E15.6,a)') 'Time step                  : ', pp%dt, ' sec.'
write(*,'(a,i16)')     'Number of time steps       : ', pp%ninc
write(*,'(a,F15.6,a)') 'Time to finish the run     : ', (ETIME(2)-ETIME(1)) , ' sec.'
call WriteLog('', ovrwrt=0)
call WriteLog('--------------------------------------------------------', ovrwrt=0)
call WriteLog('Time step                  : '//trim(d2str(pp%dt))//' sec.', ovrwrt=0)
call WriteLog('Number of time steps       : '//trim(i2str(pp%ninc)) ,ovrwrt=0)
call WriteLog('Time to finish the run     : '//trim(d2str(ETIME(2)-ETIME(1)))//' sec.', ovrwrt=0)
call WriteLog('--------------------------------------------------------', ovrwrt=0)
call WriteLog('', ovrwrt=0)
call writeOutputFile(fileName=trim(pp%outputFileName),time=time,r=r,BFld=BFldOnPath)
return
end subroutine relativisticPusher

!**************************************************************************************************************
!
!                                S U B R O U T I N E: writeOutputFile
!
!**************************************************************************************************************
subroutine writeOutputFile(fileName,time,r,BFld)
implicit none
character(64),intent(in)    :: fileName
!........................................................................................
integer(8)           :: i1
real(8), dimension(:,:), intent(in) :: r, BFld
real(8), dimension(:)  , intent(in) :: time

open(1,file=trim(fileName))

do i1=1,size(r,1)
  write(1,'(6(ES15.7E3,","),ES15.7E3)') time(i1),r(i1,1),r(i1,2),r(i1,3),BFld(i1,1),BFld(i1,2),BFld(i1,3) 
enddo
  
close(1)

end subroutine
!**************************************************************************************************************
!
!                                         S U B R O U T I N E:  random_stduniform
!                                           
!**************************************************************************************************************
!Random numbers
!Based on:
!https://masuday.github.io/fortran_tutorial/random.html
subroutine random_stduniform(u)
   implicit none
   real(8),intent(out) :: u
   real(8) :: r
   call random_number(r)
   u = d1p - r
end subroutine random_stduniform    
!**************************************************************************************************************
!
!                                         S U B R O U T I N E:  random_stdnormal
!                                           
!**************************************************************************************************************
subroutine random_stdnormal(x)
    real(8),intent(out) :: x
    real(8) :: u1,u2
    call random_stduniform(u1)
    call random_stduniform(u2)
    x = dsqrt(-d2p*dlog(u1))*dcos(d2p*pi*u2)
end subroutine random_stdnormal
!**************************************************************************************************************
!
!                                         F U N C T I O N:  random_normal_mean_std
!                                           
!**************************************************************************************************************
real(8) function random_normal_mean_std(mean,std)   
real(8),intent(in) :: mean,std
real(8) :: x
call random_stdnormal(x=x)
random_normal_mean_std = std * x + mean
end function random_normal_mean_std

!**************************************************************************************************************
!
!                                         S U B R O U T I N E:  magneticField
!                                           
!**************************************************************************************************************
subroutine magneticField(x,y,z,pp,BFld,err)
implicit none
real(8), intent(in)                       :: x, y, z
type(pusherParams),intent(in)             :: pp
real(8), dimension(3), intent(out)        :: BFld  
integer, intent(out)                      :: err
!....................................................
real(8)                                   :: Bo,Bphi,rho,phi,H,alfa,cosTita,Bsw(3),tempB,xparable,noise(3)
cosTita = x / dsqrt(x**2+y**2)
err = 0
if(pp%addBFldNoise .eqv. .true.)then
  !with noise  
  if (trim(pp%fieldShape) .eq. 'UNIFORM') then
    tempB = random_normal_mean_std(mean=pp%BswMean,std=pp%BswNoiseStd)
    BFld(1) = pp%dirBsw(1) * tempB   
    BFld(2) = pp%dirBsw(2) * tempB
    BFld(3) = pp%dirBsw(3) * tempB
  elseif (trim(pp%fieldShape) .eq. 'BESSEL') then
    if (x**2+y**2 .le. pp%frRadius**2) then
      !Zone 1 (flux rope)
      alfa = pp%alfa
      H=pp%H
      rho = dsqrt(x**2 + y**2)
      Bo = pp%BfrMean 
      noise(1) = random_normal_mean_std(mean=d0p,std=pp%BfrNoiseStd)
      noise(2) = random_normal_mean_std(mean=d0p,std=pp%BfrNoiseStd)
      noise(3) = random_normal_mean_std(mean=d0p,std=pp%BfrNoiseStd)
      
      Bphi = Bo*H*bessel_jn(1,alfa*rho) 
      phi = datan2(y=y, x=x)
      BFld(1) = -Bphi*dsin(phi) + noise(1)
      BFld(2) = Bphi*dcos(phi) + noise(2)
      BFld(3) = Bo*bessel_jn(0,alfa*rho) + noise(3)
    else
      noise(1) = random_normal_mean_std(mean=d0p,std=pp%BswNoiseStd)
      noise(2) = random_normal_mean_std(mean=d0p,std=pp%BswNoiseStd)
      noise(3) = random_normal_mean_std(mean=d0p,std=pp%BswNoiseStd)  
      BFld(1) = pp%dirBsw(1) * pp%BswMean + noise(1)   
      BFld(2) = pp%dirBsw(2) * pp%BswMean + noise(2)
      BFld(3) = pp%dirBsw(3) * pp%BswMean + noise(2)  
    endif
  elseif (trim(pp%fieldShape) .eq. 'BESSEL_SHEATH') then
    xparable = parableSimX(x0=pp%focalDist, y0=d0p, a=-pp%focalDist, y=y)
    if (x**2+y**2 .le. pp%frRadius**2) then
      !Zone 1 (flux rope)
      alfa = pp%alfa
      H=pp%H
      rho = dsqrt(x**2 + y**2)
      Bo = pp%BfrMean 
      noise(1) = random_normal_mean_std(mean=d0p,std=pp%BfrNoiseStd)
      noise(2) = random_normal_mean_std(mean=d0p,std=pp%BfrNoiseStd)
      noise(3) = random_normal_mean_std(mean=d0p,std=pp%BfrNoiseStd)
      
      Bphi = Bo*H*bessel_jn(1,alfa*rho)
      phi = datan2(y=y, x=x)
      BFld(1) = -Bphi*dsin(phi) + noise(1)
      BFld(2) = Bphi*dcos(phi) + noise(2)
      BFld(3) = Bo*bessel_jn(0,alfa*rho) + noise(3)
    elseif(x .gt. xparable) then
      !Zone 2 (before the parable)
      noise(1) = random_normal_mean_std(mean=d0p,std=pp%BswNoiseStd)
      noise(2) = random_normal_mean_std(mean=d0p,std=pp%BswNoiseStd)
      noise(3) = random_normal_mean_std(mean=d0p,std=pp%BswNoiseStd)  
      BFld(1) = pp%dirBsw(1) * pp%BswMean + noise(1)   
      BFld(2) = pp%dirBsw(2) * pp%BswMean + noise(2)
      BFld(3) = pp%dirBsw(3) * pp%BswMean + noise(2)
    elseif(x .le. xparable .and. x .ge. d0p)then
      !Zone 3 (Sheath: between parable and flux rope)
      noise(1) = random_normal_mean_std(mean=d0p,std=pp%BswNoiseStd)
      noise(2) = random_normal_mean_std(mean=d0p,std=pp%BswNoiseStd)
      noise(3) = random_normal_mean_std(mean=d0p,std=pp%BswNoiseStd)  
      Bsw(1) = pp%dirBsw(1) * pp%BswMean + noise(1)   
      Bsw(2) = pp%dirBsw(2) * pp%BswMean + noise(2)
      Bsw(3) = pp%dirBsw(3) * pp%BswMean + noise(2)
      
      noise(1) = random_normal_mean_std(mean=d0p,std=pp%BshNoiseStd)
      noise(2) = random_normal_mean_std(mean=d0p,std=pp%BshNoiseStd)
      noise(3) = random_normal_mean_std(mean=d0p,std=pp%BshNoiseStd)
      
      BFld(1) = (pp%dirBsw(1) * (pp%BshMean-pp%BswMean) + noise(1))*cosTita +  Bsw(1)     
      BFld(2) = (pp%dirBsw(2) * (pp%BshMean-pp%BswMean) + noise(2))*cosTita +  Bsw(2)
      BFld(3) = (pp%dirBsw(3) * (pp%BshMean-pp%BswMean) + noise(3))*cosTita +  Bsw(3)
      
    else
      !Zone 4 (from flux rope center backwars) 
      noise(1) = random_normal_mean_std(mean=d0p,std=pp%BswNoiseStd)
      noise(2) = random_normal_mean_std(mean=d0p,std=pp%BswNoiseStd)
      noise(3) = random_normal_mean_std(mean=d0p,std=pp%BswNoiseStd)  
      BFld(1) = pp%dirBsw(1) * pp%BswMean + noise(1)   
      BFld(2) = pp%dirBsw(2) * pp%BswMean + noise(2)
      BFld(3) = pp%dirBsw(3) * pp%BswMean + noise(2)
      
    endif
  else
    !Error message
    
  endif
else
  !Without Noise
  if (trim(pp%fieldShape) .eq. 'UNIFORM') then
    tempB = pp%BswMean
    BFld(1) = pp%dirBsw(1) * tempB   
    BFld(2) = pp%dirBsw(2) * tempB
    BFld(3) = pp%dirBsw(3) * tempB
  elseif (trim(pp%fieldShape) .eq. 'BESSEL')then
    if (x**2+y**2 .le. pp%frRadius**2) then
      !Zone 2 (flux rope)
      alfa = pp%alfa
      H=pp%H
      rho = dsqrt(x**2 + y**2)
      Bo = pp%BfrMean
      Bphi = Bo*H*bessel_jn(1,alfa*rho)
      phi = datan2(y=y, x=x)
      BFld(1) = -Bphi*dsin(phi)
      BFld(2) = Bphi*dcos(phi)
      BFld(3) = Bo*bessel_jn(0,alfa*rho)  
    else
      !Zone 4 (from flux rope center backwards)
      tempB = pp%BswMean
      BFld(1) = pp%dirBsw(1) * tempB
      BFld(2) = pp%dirBsw(2) * tempB
      BFld(3) = pp%dirBsw(3) * tempB
    endif  
  elseif (trim(pp%fieldShape) .eq. 'BESSEL_SHEATH')then
    xparable = parableSimX(x0=pp%focalDist, y0=d0p, a=-pp%focalDist, y=y)
    if (x**2+y**2 .le. pp%frRadius**2) then
      !Zone 1 (flux rope)
      alfa = pp%alfa
      H=pp%H
      rho = dsqrt(x**2 + y**2)
      Bo = pp%BfrMean
      Bphi = Bo*H*bessel_jn(1,alfa*rho)
      phi = datan2(y=y, x=x)
      BFld(1) = -Bphi*dsin(phi)
      BFld(2) = Bphi*dcos(phi)
      BFld(3) = Bo*bessel_jn(0,alfa*rho)
    elseif(x .gt. xparable) then
      !Zone 2 (before the parable)
      tempB = pp%BswMean
      BFld(1) = pp%dirBsw(1) * tempB   
      BFld(2) = pp%dirBsw(2) * tempB
      BFld(3) = pp%dirBsw(3) * tempB  
    elseif(x .le. xparable .and. x .ge. d0p)then
      !Zone 3 (between parable and flux rope) (sheath)
      tempB = pp%BswMean
      Bsw(1) = pp%dirBsw(1) * tempB   
      Bsw(2) = pp%dirBsw(2) * tempB
      Bsw(3) = pp%dirBsw(3) * tempB
      tempB = pp%BshMean-pp%BswMean
      BFld(1) = pp%dirBsw(1) * tempB *cosTita +  Bsw(1)     
      BFld(2) = pp%dirBsw(2) * tempB *cosTita +  Bsw(2)
      BFld(3) = pp%dirBsw(3) * tempB *cosTita +  Bsw(3)  
    else
      !Zone 4 (desde el centro del flux rope hacia atras)
      tempB = pp%BswMean
      BFld(1) = pp%dirBsw(1) * tempB
      BFld(2) = pp%dirBsw(2) * tempB
      BFld(3) = pp%dirBsw(3) * tempB
    endif
  else
    !Mensaje de error  
  endif    
endif

return
1 call WriteLog(text="ERROR. Magnetic field shape not specified.",ovrwrt=0)
err = err + 1
return
end subroutine magneticField
!**************************************************************************************************************
!
!                                          F U N C T I O N:  parableSimX
!                                          
!**************************************************************************************************************
real(8) function parableSimX(x0,y0,a,y)
implicit none
real(8), intent(in) :: x0,y0,a,y
parableSimX = (y-y0)**2 / (d4p*a) + x0
end function parableSimX

end module subsfuncsModule
