!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module cooling
!
! Gas cooling
!  Current options:
!     0 = none
!     1 = 'explicit' cooling [implicitly calculated, ironically] [default]
!         (not sure what this actually means, but requires a lot of non-default inputs)
!     2 = Townsend (2009) cooling tables [implicitly calculated]
!     3 = Gammie cooling [explicitly calculated]
!     5 = Koyama & Inutuska (2002) [explicitly calculated]
!     6 = Koyama & Inutuska (2002) [implicitly calculated]
!     7 = Joung & Mac Low (2006) in the framework of Koyama & Inutuska (2002) [implicitly calculated]
!         (suitable for high-temp simulations up to 10^8K)
!
! :References:
!   Koyama & Inutsuka (2002), ApJL 564, 97-100
!   Vazquez-Semadeni, et.al (2007), ApJ 657, 870-883
!   Townsend (2009), ApJS 181, 391-397
!   Gail & Sedlmayr textbook Physics and chemistry of Circumstellar dust shells
!   Joung & Mac Low (2006), ApJ 653, 1266-1279
!   De Rijcke,et.al (2013), MNRAS 433, 3005-3016
!   SPHNG coolcurve.F & thermeq.f modified by Ian Bonnell
!
! :Owner: Lionel Siess
!
! :Runtime parameters:
!   - C_cool               : *factor controlling cooling timestep*
!   - Tfloor               : *temperature floor (K); on if > 0*
!   - beta_cool            : *beta factor in Gammie (2001) cooling*
!   - bowen_Cprime         : *radiative cooling rate (g.s/cm³)*
!   - cooltable            : *data file containing cooling function*
!   - habund               : *Hydrogen abundance assumed in cooling function*
!   - icool_dust_collision : *dust collision on/off*
!   - icool_radiation_H0   : *H0 cooling on/off*
!   - icool_relax_bowen    : *Bowen (diffusive) relaxation on/off*
!   - icool_relax_stefan   : *radiative relaxation on/off*
!   - icooling             : *cooling function (0=off, 1=explicit, 2=Townsend table, 3=Gammie, 5=KI02)*
!   - temp_floor           : *Minimum allowed temperature in K for Townsend cooling table*
!
! :Dependencies: chem, datafiles, dim, eos, h2cooling, infile_utils, io,
!   options, part, physcon, timestep, units
!

 use options,  only:icooling
 use timestep, only:C_cool

 implicit none
 character(len=*), parameter :: label = 'cooling'

 public :: init_cooling,calc_cooling_rate,energ_cooling
 public :: write_options_cooling, read_options_cooling
 public :: find_in_table, implicit_cooling, exact_cooling
 public :: lambdacoolJML,lambdacoolDR
 logical, public :: calc_Teq

 logical, public :: cooling_implicit
 logical, public :: cooling_explicit
 real,    public :: bowen_Cprime = 3.000d-5
 real,    public :: GammaKI_cgs = 2.d-26 ! [erg/s] heating rate for Koyama & Inutuska cooling

 integer, public, parameter :: maxcoolingoff = 20
 logical, public :: snecoolingoff = .false.      ! Temporarily switch off cooling around sne
 integer, public :: ncoolingoff = 0              ! Number of coords to switch off cooling at current time
 real,    public :: xyzh_coolingoff(4,maxcoolingoff)
 real,    public :: dt_cooloff = 2E-4            ! Time period to disable cooling
 real,    public :: range_cooloff = 5E-2         ! Radius around sn location to disable cooling (code units)
 logical, public :: range_cooloff_useh = .false. ! or use resolution length of progenitor

 integer, public, parameter :: maxpart_nocooling = 1E8
 logical, public :: part_nocooling = .true.            ! Disable cooling for some particles
 integer, public :: ipart_nocooling(maxpart_nocooling) ! Indicies of particles to stop cooling
 integer, public :: npartnocool                        ! Current total number of particles to stop cooling

 private
 integer, parameter :: nTg = 64
 integer, parameter :: maxt = 5000
 real,    parameter :: Tref = 1.d5, T_floor = 10.   ! required for exact_cooling
 integer :: nt
 real    :: temper(maxt),lambda(maxt),slope(maxt),yfunc(maxt)
 real    :: beta_cool  = 3.
 real    :: habund     = 0.7
 real    :: temp_floor = 1.e4                       ! required for exact_cooling_table
 real    :: Tgrid(nTg)

 real    :: LambdaKI_coef,GammaKI
 real    :: KI02_rho_min_cgs = 1.0d-30  ! minimum density of the KI02 cooling curve
 real    :: KI02_rho_max_cgs = 1.0d-14  ! maximum density of the KI02 cooling curve
 real    :: KI02_rho_min,KI02_rho_max
 real    :: rhov4_KI02(2,maxt)

 real    :: rhominJML_cgs = 6.35d-29   ! density range of which JML06 cooling curve has equilibium solution(s)
 real    :: rhomaxJML_cgs = 1.d-14    !  -Note: rhominJML_cgs is hard limit; rhomaxJML_cgs can be increased by lowering TminJML
 real    :: TminJML = 1E0            ! temperature range of JML06 cooling curve
 real    :: TmaxJML = 1E9
 real    :: rhoueqJML_table(5,maxt)

 integer :: icool_radiation_H0 = 0, icool_relax_Bowen = 0, icool_dust_collision = 0, icool_relax_Stefan = 0
 character(len=120) :: cooltable = 'cooltable.dat'
 !--Minimum temperature (failsafe to prevent u < 0); optional for ALL cooling options
 real,    public :: Tfloor = 0. ! [K]; set in .in file.  On if Tfloor > 0.
 real,    public :: ufloor = 0. ! [code units]; set in init_cooling

contains

!-----------------------------------------------------------------------
!+
!  Initialise cooling
!+
!-----------------------------------------------------------------------
subroutine init_cooling(id,master,iprint,ierr)
 use dim,       only:maxvxyzu
 use units,     only:utime,umass,udist,unit_ergg
 use physcon,   only:mass_proton_cgs,kboltz
 use io,        only:fatal
 use eos,       only:gamma,gmw
 use part,      only:h2chemistry
 use h2cooling, only:init_h2cooling
 use chem,      only:init_chem
 integer, intent(in)  :: id,master,iprint
 integer, intent(out) :: ierr

 if (h2chemistry) then
    if (id==master) write(iprint,*) 'initialising cooling function...'
    call init_chem()  ! nothing in init_chem
    call init_h2cooling()
 else
    !you can't have cool_relaxation_Stefan and cool_relaxation_Bowen at the same time
    if (icool_relax_bowen == 1 .and. icool_relax_stefan == 1) then
       call fatal(label,'you can"t have bowen and stefan cooling at the same time')
    endif

#ifdef KROME
    !krome calculates its own cooling rate
    icool_radiation_H0 = 0
    icool_dust_collision = 0
#else
    !if no cooling flag activated, disable cooling
    if (icooling == 1 .and. (icool_radiation_H0+icool_relax_Bowen+icool_dust_collision+&
          icool_relax_Stefan == 0)) then
       icooling = 0
       calc_Teq = .false.
       return
    endif
#endif
    calc_Teq = (icool_relax_Bowen == 1) .or. (icool_relax_Stefan == 1) .or. (icool_dust_collision == 1)

    !--initialise remaining variables
    if (icooling == 2) then
       call init_cooltable(ierr)
    elseif (icooling == 5 .or. icooling == 6 .or. icooling == 7) then
       LambdaKI_coef = GammaKI_cgs*umass*utime**3/(mass_proton_cgs**2 * udist**5)
       GammaKI       = GammaKI_cgs*utime**3/(mass_proton_cgs*udist**2)
       if (icooling == 7) then
          call init_rhoutableJML(ierr)
          if (ierr > 0) call fatal('init_cooling','Failed to create JML06 cooling table')
       else
          call init_hv4table(ierr)
          if (ierr > 0) call fatal('init_cooling','Failed to create KI02 cooling table')
       endif
    elseif (icooling > 0) then
       call set_Tgrid
    endif
 endif

 !--Determine if this is implicit or explicit cooling
 cooling_implicit = .false.
 cooling_explicit = .false.
 if (h2chemistry) then
    if (icooling > 0) cooling_implicit = .true.    ! cooling is calculated implicitly in step
 elseif (icooling > 0) then
    if (icooling == 3 .or. icooling == 5) then
       cooling_explicit = .true.                   ! cooling is calculated explicitly in force
    else
       cooling_implicit = .true.                   ! cooling is calculated implicitly in step
    endif
 endif

 !--calculate the energy floor in code units
 if (Tfloor > 0.) then
    if (gamma > 1.) then
       ufloor = kboltz*Tfloor/((gamma-1.)*gmw*mass_proton_cgs)/unit_ergg
    else
       ufloor = 3.0*kboltz*Tfloor/(2.0*gmw*mass_proton_cgs)/unit_ergg
    endif
    if (maxvxyzu < 4) ierr = 1
 else
    ufloor = 0.
 endif

end subroutine init_cooling

!-----------------------------------------------------------------------
!+
!  read cooling table from file and initialise arrays
!+
!-----------------------------------------------------------------------
subroutine init_cooltable(ierr)
 use io,        only:fatal
 use datafiles, only:find_phantom_datafile
 integer, intent(out) :: ierr
 integer, parameter :: iu = 127
 integer :: i
 character(len=120) :: filepath

 !
 ! read the cooling table from file
 !
 filepath=find_phantom_datafile(cooltable,'cooling')
 open(unit=iu,file=filepath,status='old',iostat=ierr)
 if (ierr /= 0) call fatal('cooling','error opening cooling table')
 i = 0
 do while(ierr==0 .and. i < maxt)
    i = i + 1
    read(iu,*,iostat=ierr) temper(i),lambda(i)
 enddo
 nt = i-1
 if (nt==maxt) call fatal('cooling','size of cooling table exceeds array size')
 if (nt < 2) call fatal('cooling','size of cooling table is too small',ival=nt,var='nt')
 !
 ! calculate the slope of the cooling function
 !
 do i=1,nt-1
    slope(i) = log(lambda(i+1)/lambda(i))/log(temper(i+1)/temper(i))
 enddo
 slope(nt) = slope(nt-1)

 !
 ! initialise the functions required for Townsend method
 !
 yfunc(nt) = 0.
 do i=nt-1,1,-1
    !Lionel Siess : I think there is an error. in yfunc slope(nt) should be replaced by lambda(nt)
    ! Eq A6
    if (abs(slope(i)-1.) < tiny(0.)) then
       !ori yfunc(i) = yfunc(i+1) - slope(nt)*temper(i)/(lambda(i)*temper(nt))*log(temper(i)/temper(i+1))
       yfunc(i) = yfunc(i+1) - lambda(nt)*temper(i)/(lambda(i)*temper(nt))*log(temper(i)/temper(i+1))
    else
       !ori yfunc(i) = yfunc(i+1) - slope(nt)*temper(i)/((1. - slope(i))*lambda(i)*temper(nt))&
       yfunc(i) = yfunc(i+1) - lambda(nt)*temper(i)/((1. - slope(i))*lambda(i)*temper(nt))&
                 *(1.- (temper(i)/temper(i+1))**(slope(i) - 1.))
    endif
 enddo
end subroutine init_cooltable

!-----------------------------------------------------------------------
!+
!  create a h-v4 table based upon the cooling curve of KI02
!+
!-----------------------------------------------------------------------
subroutine init_hv4table(ierr)
 use part,    only:hrho,igas
 use physcon, only:mass_proton_cgs,kboltz
 use units,   only:unit_density,unit_velocity
 use eos,     only:gmw,gamma
 integer, intent(out) :: ierr
 integer              :: i,ctr
 real                 :: nrho0_min,nrho0_max,nrho,dnrho,dGammaKI,Lambda,dLambda
 real                 :: T,Tnew,Trat,fatT,faTdT
 logical              :: iterate
 logical              :: print_cc = .true. ! Print the cooling curve (for testing)

 !--Initialise densities
 KI02_rho_min = KI02_rho_min_cgs/unit_density
 KI02_rho_max = KI02_rho_max_cgs/unit_density
 nrho0_min    = KI02_rho_min_cgs/mass_proton_cgs
 nrho0_max    = KI02_rho_max_cgs/mass_proton_cgs
 dnrho        = (log10(nrho0_max) - log10(nrho0_min))/maxt
 !--Initialise additional variables
 dGammaKI     = 0.0
 ierr         = 0

 if (print_cc) open(unit=1031,file='coolingcurve.dat')
 if (print_cc) write(1031,'(3A15)') 'n','T','Lambda(T)'

 !--Iterate (in cgs units)!
 T = 20000.
 do i = 1,maxt
    ctr     = 0
    iterate = .true.
    nrho    = 10**(log10(nrho0_min) + (i-1)*dnrho)
    do while ( iterate )
       Lambda  = 1.d7*exp(-1.184d5/(T+1.d3)) + 0.014*sqrt(T)*exp(-92./T) ! This is actually Lamda / Gamma
       dLambda = 0.007*exp(-92./T)*(T+184.)*T**(-1.5) + 1.184d12*exp(-1.184d5/(T+1.d3))*(T+1.d3)**(-2)
       fatT    = Lambda*GammaKI_cgs*nrho - GammaKI_cgs
       faTdT   = dLambda*GammaKI_cgs*nrho - dGammaKI
       Tnew    = abs(T - fatT/faTdT)
       Trat    = abs( 1.0 - T/Tnew )
       T       = Tnew
       ctr     = ctr + 1
       !--converged
       if (Trat < 1.0d-6) iterate = .false.
       !--failed to converge
       if (T < 0. .or. ctr > 2000) then
          iterate = .false.
          ierr    = 1
       endif
    enddo
    if (print_cc) then
       if (nrho > 0) write(1031,'(E15.5E2,F15.5,E15.5E2)') nrho,T,Lambda*GammaKI_cgs
    endif
    rhov4_KI02(1,i) = nrho
    rhov4_KI02(2,i) = T
 enddo
 if (print_cc) close(1031)

 !--Convert to useful values
 do i = 1,maxt
    rhov4_KI02(1,i) = rhov4_KI02(1,i)*mass_proton_cgs/unit_density                               ! number density (cm^-3) -> mass density (code units)
    rhov4_KI02(2,i) = kboltz*rhov4_KI02(2,i)/(gmw*mass_proton_cgs*(gamma-1.0))/unit_velocity**2  ! T -> internal energy (code units)
 enddo

end subroutine init_hv4table

!-----------------------------------------------------------------------
!
!  Procedures for initiating rho-u table using the cooling curve of JML06;
!  For each rho, there can be up to 3 thermal equilibium solutions
!
!-----------------------------------------------------------------------
subroutine init_rhoutableJML(ierr)
 use io,      only:fatal
 use physcon, only:mass_proton_cgs,kboltz
 use units,   only:unit_density,unit_ergg
 use eos,     only:gmw,gamma
 integer, intent(out) :: ierr
 integer :: i,irterr1,irterr2,irterr3,numroots
 real    :: nrhomin_cgs,nrhomax_cgs,dnrho_cgs,nrho_cgs,rho
 real    :: T01,T02,Teq1,Teq2,Teq3,ueq1,ueq2,ueq3
 logical :: write_table = .true. 

 nrhomin_cgs = rhominJML_cgs/mass_proton_cgs
 nrhomax_cgs = rhomaxJML_cgs/mass_proton_cgs
 dnrho_cgs = (log10(nrhomax_cgs)-log10(nrhomin_cgs))/maxt
 ierr = 0

 if (write_table) open(unit=3020,file='rho_Teq_table.dat',status='replace')

 ! Temp of local max and min of JML06 cooling curve (for root-search brackets)
 T01 = 198609.4917357372
 T02 = 31622776.6016837933

 ! Solve for Teq1, Teq2 and Teq3 for each value of n
 do i = 1,maxt
    nrho_cgs = nrhomin_cgs * 10**((i-1)*dnrho_cgs)
    call root_bisection(nrho_cgs,TminJML,T01,Teq1,irterr1)
    call root_bisection(nrho_cgs,T01,T02,Teq2,irterr2)
    call root_bisection(nrho_cgs,T02,TmaxJML,Teq3,irterr3)

    ! Get total number of roots
    numroots = 0
    if (irterr3 == 1) then
       if (irterr2 == 1) then
          if (irterr1 == 1) then
             ierr = 1
             call fatal('cooling','No roots found for one of the densities in table')
          else
             numroots = 1
          endif
       else
          numroots = 2
       endif
    else
       numroots = 3
    endif

    ! Convert nrho_cgs and T to rho and u in code units
    rho = nrho_cgs * mass_proton_cgs / unit_density
    ueq1 = kboltz * Teq1 / (gmw*mass_proton_cgs*(gamma-1.)) /unit_ergg
    ueq2 = kboltz * Teq2 / (gmw*mass_proton_cgs*(gamma-1.)) /unit_ergg
    ueq3 = kboltz * Teq3 / (gmw*mass_proton_cgs*(gamma-1.)) /unit_ergg

    rhoueqJML_table(1,i) = rho
    rhoueqJML_table(2,i) = numroots
    rhoueqJML_table(3,i) = ueq1
    rhoueqJML_table(4,i) = ueq2
    rhoueqJML_table(5,i) = ueq3

    if (write_table) write(3020,*) rho*unit_density,numroots,Teq1,Teq2,Teq3 
 enddo

 if (write_table) close(3020)

end subroutine init_rhoutableJML

!
! Cooling function mimicking the cooling curve of Joung & Mac Low (2006) Fig. 1;
! created by modifying the cooling function of Koyama & Inutsuka (2002)
!
real function lambdacoolJML(temp)
 real, intent(in) :: temp
 real  :: lambdagamma1,lambdagamma2,lambdagamma

 ! First term of KI02
 if (temp > 10**4.15) then
    lambdagamma1 = 4.69414E-4 * 1E7 * exp(-1.184E5 * 1.15983E6 / ((temp*10**(-0.08))**2.68935 + 1000))
 else
    lambdagamma1 = 1E7 * exp(-1.184E5 / ((temp*10**(-0.16)) + 1000))* temp**0.18
 endif
 ! Second term of KI02
 lambdagamma2 = 0.215 * 0.014 * (temp*10**(-0.2))**0.66 * exp(-92/(temp*10**(-0.12))) *10**0.2

 ! Shift up
 lambdagamma = (lambdagamma1 + lambdagamma2) * 10**0.75

 ! Modify low temp part
 if (temp < 10**3.7) then
    lambdagamma = lambdagamma * 10**1.665 * temp**(-0.45)
 endif
 ! Modify high temp part
 if (temp > 10**5.3) then
    if (temp < 10**6.5) then
       lambdagamma = lambdagamma * 10**5.3 * temp**(-1.)
    else
       lambdagamma = lambdagamma * 10**0.425 * temp**(-0.25)
    endif
 endif
 ! Modify highest temp part
 if (temp > 10**7.5) then
    lambdagamma = lambdagamma * 10**(-5.25) * temp**0.7
 endif

 lambdacoolJML = lambdagamma * GammaKI_cgs

end function lambdacoolJML



real function lambdacoolDR(temp)
 real, intent(in) :: temp
 real  :: lambdagamma1,lambdagamma2,lambdagamma

 if (temp > 10**4.15) then
    lambdagamma1 = 4.69414E-4 * 1E7 * exp(-1.184E5 * 1.15983E6 / ((temp*10**(-0.08))**2.68935 + 1000))
 else
    lambdagamma1 = 1E7 * exp(-1.184E5 / ((temp*10**(-0.16)) + 1000))* temp**0.18
 endif

 lambdagamma2 = 0.115 * 0.014 * (temp*10**(-0.75))**0.60 * exp(-72/(temp*10**(-0.12)))
 lambdagamma2 = lambdagamma2 * temp**(0.15)

 lambdagamma = (lambdagamma1 + lambdagamma2) * 10**(0.95)

 if (temp < 10**3.7) then
    lambdagamma = lambdagamma * 10**1.665 * temp**(-0.45)
 endif

 if (temp > 10**5.3) then
    if (temp < 10**6.5) then
        lambdagamma = lambdagamma * 10**5.3 * temp**(-1.)
    else
        lambdagamma = lambdagamma * 10**0.425 * temp**(-0.25)
    endif
 endif

 if (temp > 10**7.5) then
    lambdagamma = lambdagamma * 10**(-5.25) * temp**(0.7)
 endif

 lambdagamma = lambdagamma / 10**(0.35)

 lambdacoolDR = lambdagamma * 2E-26

end function lambdacoolDR


!
! Thermal equilibium: f_T = n*lambda(T) - gamma
!
real function equifunc(nrho_cgs,temp)
 real, intent(in) :: nrho_cgs,temp

 equifunc = nrho_cgs*GammaKI_cgs - nrho_cgs**2*lambdacoolJML(temp)

end function equifunc

!
! Brackets the equilibium solution Teq with given Tmin and Tmax;
! returns Teq = 0. and irterr = 1 if no roots found
!
subroutine root_bisection(nrho,Tmin0,Tmax0,Teq,irterr)
 use io,  only:fatal
 real,    intent(in)  :: nrho,Tmin0,Tmax0
 integer, intent(out) :: irterr
 real,    intent(out) :: Teq
 integer :: niter
 real    :: Tmin,Tmax,Tmid,Tminsign,Tmaxsign,func_min,func_max,func_mid
 real    :: tol = 0.05   ! tolerance in root-search
 logical :: converged

 converged = .false.  ! root flag
 irterr = 0           ! no-root indicator
 niter = 0

 Tminsign = sign(1.,equifunc(nrho,Tmin0))
 Tmaxsign = sign(1.,equifunc(nrho,Tmax0))
 if (Tminsign == Tmaxsign) then
    Teq = 0.
    irterr = 1
 endif

 Tmin = Tmin0
 Tmax = Tmax0
 do while (.not.converged .and. irterr == 0)

    func_max = equifunc(nrho,Tmax)
    func_min = equifunc(nrho,Tmin)
    if (func_max*func_min > 0) call fatal('heating_cooling_cmi','not bracketing Teq root')

    Tmid = (Tmax+Tmin)/2.
    func_mid = equifunc(nrho,Tmid)

    if (func_mid*func_min < 0) then
       Tmax = Tmid
    elseif (func_mid*func_max < 0) then
       Tmin = Tmid
    endif

    if ((Tmax-Tmin) < tol) then
       converged = .true.
       Teq = (Tmax+Tmin)/2.
    endif

    niter = niter + 1
    if (niter > 1000) then
       Teq = 0.
       irterr = 1
    endif
 enddo

end subroutine root_bisection

!-----------------------------------------------------------------------
!
!  calculate cooling rates
!
!-----------------------------------------------------------------------
subroutine calc_cooling_rate(Q, dlnQ_dlnT, rho, T, Teq, mu, K2, kappa)
 use units,   only:unit_ergg,unit_density
 real, intent(in) :: rho, T, Teq !rho in code units
 real, intent(in), optional :: mu, K2, kappa !cgs
 real, intent(out) :: Q, dlnQ_dlnT !code units
 real :: Q_cgs,Q_H0, Q_relax_Bowen, Q_col_dust, Q_relax_Stefan, rho_cgs
 real :: dlnQ_H0, dlnQ_relax_Bowen, dlnQ_col_dust, dlnQ_relax_Stefan

 rho_cgs = rho*unit_density
 Q_H0 = 0.
 Q_relax_Bowen = 0.
 Q_col_dust = 0.
 Q_relax_Stefan = 0.
 dlnQ_H0 = 0.
 dlnQ_relax_Bowen = 0.
 dlnQ_col_dust = 0.
 dlnQ_relax_Stefan = 0.
 if (icool_radiation_H0 == 1)   call cooling_neutral_hydrogen(T, rho_cgs, Q_H0, dlnQ_H0)
 if (icool_relax_Bowen == 1)    call cooling_Bowen_relaxation(T, Teq, rho_cgs, mu, Q_relax_Bowen, dlnQ_relax_Bowen)
 if (icool_dust_collision == 1) call cooling_dust_collision(T, Teq, rho_cgs, K2, mu, Q_col_dust, dlnQ_col_dust)
 if (icool_relax_Stefan == 1)   call cooling_radiative_relaxation(T, Teq, kappa, Q_relax_Stefan, dlnQ_relax_Stefan)
 Q_cgs = Q_H0 + Q_relax_Bowen+ Q_col_dust+ Q_relax_Stefan
 dlnQ_dlnT = (Q_H0*dlnQ_H0 + Q_relax_Bowen*dlnQ_relax_Bowen+ Q_col_dust*dlnQ_col_dust+ Q_relax_Stefan*dlnQ_relax_Stefan)/Q_cgs
 !limit exponent to prevent overflow
 dlnQ_dlnT = sign(min(50.,abs(dlnQ_dlnT)),dlnQ_dlnT)
 Q = Q_cgs/unit_ergg
end subroutine calc_cooling_rate


!-----------------------------------------------------------------------
!+
!  Bowen 1988 cooling term
!+
!-----------------------------------------------------------------------
subroutine cooling_Bowen_relaxation(T, Teq, rho, mu, Q, dlnQ_dlnT)
! all quantities in cgs
 use eos,     only:gamma
 use physcon, only:Rg
 real, intent(in) :: T, Teq, rho, mu
 real, intent(out) :: Q,dlnQ_dlnT

 Q = Rg/((gamma-1.)*mu)*rho*(Teq-T)/bowen_Cprime
 dlnQ_dlnT = -T/(Teq-T+1.d-10)

end subroutine cooling_Bowen_relaxation

!-----------------------------------------------------------------------
!+
!  collisionnal cooling
!+
!-----------------------------------------------------------------------
subroutine cooling_dust_collision(T, Teq, rho, K2, mu, Q, dlnQ_dlnT)
! all quantities in cgs
 use physcon, only: kboltz, mass_proton_cgs, pi
 real, intent(in) :: T, Teq, rho, K2, mu
 real, intent(out) :: Q,dlnQ_dlnT

 real, parameter :: f = 0.15, a0 = 1.28e-8
 real :: A

 A = 2. * f * kboltz * a0**2/(mass_proton_cgs**2*mu) &
         * (1.05/1.54) * sqrt(2.*pi*kboltz/mass_proton_cgs) * 2.*K2 * rho
 Q = A * sqrt(T) * (Teq-T)
 if (Q  >  1.d6) then
    print *, f, kboltz, a0, mass_proton_cgs, mu
    print *, mu, K2, rho, T, Teq, A, Q
    stop 'cooling'
 else
    dlnQ_dlnT = 0.5+T/(Teq-T+1.d-10)
 endif
end subroutine cooling_dust_collision

!-----------------------------------------------------------------------
!+
!  Woitke (2006 A&A) cooling term
!+
!-----------------------------------------------------------------------
subroutine cooling_radiative_relaxation(T, Teq, kappa, Q, dlnQ_dlnT)
 use physcon, only: steboltz
 real, intent(in) :: T, Teq, kappa
 real, intent(out) :: Q,dlnQ_dlnT

 Q = 4.*steboltz*(Teq**4-T**4)*kappa
 dlnQ_dlnT = -4.*T**4/(Teq**4-T**4+1.d-10)

end subroutine cooling_radiative_relaxation

!-----------------------------------------------------------------------
!+
!  Cooling due to neutral H (Spitzer)
!+
!-----------------------------------------------------------------------
subroutine cooling_neutral_hydrogen(T, rho, Q, dlnQ_dlnT)
 use physcon, only: mass_proton_cgs, pi
 real, intent(in) :: T, rho
 real, intent(out) :: Q,dlnQ_dlnT

 real, parameter :: f = 0.2! 1.d0 !0.2
 real :: eps_e

 if (T > 3000.) then
    eps_e = calc_eps_e(T)
    !Q = -f*7.3d-19*eps_e*exp(-118400./T)*rho/(mass_per_H)**2
    Q = -f*7.3d-19*eps_e*exp(-118400./T)*rho/(1.4*mass_proton_cgs)**2
    dlnQ_dlnT = 118400.d0/T+log(calc_eps_e(1.001*T)/eps_e)/log(1.001)
 else
    Q = 0.
    dlnQ_dlnT = 0.
 endif
end subroutine cooling_neutral_hydrogen

!-----------------------------------------------------------------------
!+
!  compute electron equilibrium abundance (Palla et al 1983)
!+
!-----------------------------------------------------------------------
real function calc_eps_e(T)
 real, intent(in) :: T
 real :: k1, k2, k3, k8, k9, p, q

 k1 = 1.88d-10 / T**6.44e-1
 k2 = 1.83d-18 * T
 k3 = 1.35d-9
 k8 = 5.80d-11 * sqrt(T) * exp(-1.58d5/T)
 k9 = 1.7d-4 * k8
 p = .5*k8/k9
 q = k1*(k2+k3)/(k3*k9)
 calc_eps_e = (p + sqrt(q+p**2))/q
end function calc_eps_e

subroutine set_Tgrid
 integer :: i
 real :: dlnT
 dlnT = log(Tref)/(nTg-1)

 do i = 1,nTg
    Tgrid(i) = exp((i-1)*dlnT)
 enddo
end subroutine set_Tgrid

!-----------------------------------------------------------------------
!+
!   Gammie (2001) cooling
!+
!-----------------------------------------------------------------------
subroutine cooling_Gammie(xi,yi,zi,ui,dudti)
 real, intent(in)    :: ui,xi,yi,zi
 real, intent(inout) :: dudti

 real :: omegai,r2,tcool1

 r2     = xi*xi + yi*yi + zi*zi
 Omegai = r2**(-0.75)
 tcool1 = Omegai/beta_cool
 dudti  = dudti - ui*tcool1

end subroutine cooling_Gammie

!-----------------------------------------------------------------------
!+
!   Cooling rate as per Koyama & Inutuska (2002; eqns 4 & 5);
!   typos corrected as per Vazquez-Semadeni+ (2007)
!   This is for the explicit calculation
!   In equilibrium, n*LambdaKI = (rho/mp)*LambdaKI = GammaKI
!+
!-----------------------------------------------------------------------
subroutine cooling_KoyamaInutuska_explicit(rhoi,Tgas,dudti)
 real, intent(in)    :: rhoi,Tgas
 real, intent(inout) :: dudti
 real                :: LambdaKI

 ! Derivation to obtain correct units; used Koyama & Inutuska (2002) as the reference
 !LambdaKI = GammaKI_cgs * (1.d7*exp(-118400./(Tgas+1000))+0.014*sqrt(Tgas)*exp(-92./Tgas)) ! The cooling rate in erg cm^3/s = g cm^5/s^3
 !LambdaKI = LambdaKI/mass_proton_cgs**2                                                    ! units are now cm^5/(g s^3) ! since [u] = erg/g = cm^2/s^2
 !LambdaKI = LambdaKI*umass*utime**3/udist**5                                               ! convert to from cm^5/(g s^3) to code units
 !dudti    = dudti - LambdaKI*rhoi*fac                                                      ! multiply by rho (code) to get l^5/(m t^3) * m/l^3 = l^2/s^3 = [u]
 !
 !GammaKI = GammaKI_cgs                                                                     ! The heating rate in erg /s = g cm^2/s^3
 !GammaKI = GammaKI/mass_proton_cgs                                                         ! divide by proton mass.  Units are now g cm^2 / s^3 / g = cm^2/s^3
 !GammaKI = GammaKI*utime**3/udist**2                                                       ! convert from cm^2/s^3 to code units
 !dudti   = dudti + GammaKI                                                                 ! units and dependencies are correct

 LambdaKI = LambdaKI_coef*(1.d7*exp(-118400./(Tgas+1000.))+0.014*sqrt(Tgas)*exp(-92./Tgas))
 dudti    = dudti - LambdaKI*rhoi + GammaKI

end subroutine cooling_KoyamaInutuska_explicit

!-----------------------------------------------------------------------
!+
!   Cooling rate as per Koyama & Inutuska (2002; eqns 4 & 5);
!   typos corrected as per Vazquez-Semadeni+ (2007)
!   This is the implicit method given by (5)-(6) in Vazquez-Semadeni+ (2007)
!+
!-----------------------------------------------------------------------
subroutine cooling_KoyamaInutuska_implicit(eni,rhoi,dt,dudti)
 use eos, only:gamma,temperature_coef,gmw
 real, intent(in)    :: rhoi,eni,dt
 real, intent(out)   :: dudti
 integer             :: i,j,jm1
 real                :: ponrhoi,tempi,eni_equil,eni_final,deni,tau1,LambdaKI

 !--Determine the indicies surrounding the input h
 i = minloc(abs(rhov4_KI02(1,1:maxt)-rhoi), 1)
 if (i==1) then
    !print*, 'min density too large! extrapolating using two smallest densities'
    j = 2
 elseif (i==maxt) then
    !print*, 'max density too small! extrapolating using two largest densities'
    j = maxt
 elseif (rhov4_KI02(1,i-1) <= rhoi .and. rhoi <= rhov4_KI02(1,i  )) then
    j = i
 elseif (rhov4_KI02(1,i  ) <= rhoi .and. rhoi <= rhov4_KI02(1,i+1)) then
    j = i+1
 else
    print*, rhoi,rhov4_KI02(1,i-1:i+1)
    print*, 'this should not happen'
    stop
 endif

 !--Calculate the equilibrium energy by linear interpolation
 jm1       = j - 1
 eni_equil = rhov4_KI02(2,j) + (rhov4_KI02(2,jm1)-rhov4_KI02(2,j))/(rhov4_KI02(1,jm1)-rhov4_KI02(1,j))*(rhoi-rhov4_KI02(1,j))

 !--Determine the inverse time require to radiate/acquire excess/deficit energy & Update energy
 ponrhoi  = (gamma-1.)*eni
 tempi    = temperature_coef*gmw*ponrhoi
 LambdaKI = LambdaKI_coef*(1.d7*exp(-118400./(tempi+1000.))+0.014*sqrt(tempi)*exp(-92./tempi))
 dudti    = LambdaKI*rhoi - GammaKI
 deni     = eni - eni_equil

 if (abs(deni) > 0.) then
    ! in both limits, this will approach the correct value
    tau1      = abs(dudti/deni)
    eni_final = eni_equil + deni*exp(-dt*tau1)
    dudti     = -(eni - eni_final)/dt
 else
    ! in the unlikly chance deni = 0
    dudti = -dudti
 endif

end subroutine cooling_KoyamaInutuska_implicit

!-----------------------------------------------------------------------
!
! Implicit cooling based on Koyama & Inutuska (2002), adapted to deal with cooling
! curve of Joung & Mac Low (2006) which has more than one equilibium solution
!
!-----------------------------------------------------------------------
subroutine cooling_JoungMacLow_implicit(eni,rhoi,dt,dudti)
 use io,      only:fatal
 use physcon, only:kboltz,mass_proton_cgs
 use eos,     only:gamma,gmw
 use units,   only:unit_ergg,unit_density 
 real, intent(in)  :: rhoi,eni,dt
 real, intent(out) :: dudti
 integer :: i,j,r,numroots
 real    :: ueqs(3),ueq_final,tempi,LambdaKI,deni,tau,eni_new
 real    :: rhotable(maxt),numroottable(maxt),utable(maxt)

 rhotable = rhoueqJML_table(1,:)
 numroottable = rhoueqJML_table(2,:)

 !- First entry in rhoueq_table is the minimum density of which there exists at least one root;
 !- Directly use first entry if rhoi is smaller than this limit
 if (rhoi < rhotable(1)) then
    numroots = int(numroottable(1))
    ueqs = (/ rhoueqJML_table(1,3), rhoueqJML_table(1,4), rhoueqJML_table(1,5) /)
 else

    !- Interpolate between table entries to give the final ueqs
    ! ([i] closest index; [j] lower bound; [j+1] upper bound around rhoi)
    i = minloc(abs(rhotable(:)-rhoi),1)
    j = 0
    if (i == 1) then
       j = 1
    elseif (i == maxt) then
       j = i-1
    elseif (rhotable(i) >= rhoi .and. rhotable(i-1) < rhoi) then
       j = i-1
    elseif (rhotable(i) < rhoi .and. rhotable(i+1) >= rhoi) then
       j = i
    endif

    !- Interpolate only if both ueq[j] and ueq[j+1] exist,
    !- otherwise, use the nearest root
    numroots = int(numroottable(i))
    if (numroots == 0) call fatal('cooling_JoungMacLow_implicit','no equilibrium solution found')
    ueqs = (/ 0., 0., 0. /)
    each_root: do r = 1,numroots
       utable = rhoueqJML_table(2+r,:)
       if (utable(j) > 0. .and. utable(j+1) > 0.) then
          ueqs(r) = utable(j) + (rhoi-rhotable(j))*(utable(j+1) - utable(j))/(rhotable(j+1)-rhotable(j))
       else
          if (utable(j) == 0. .and. utable(j+1) > 0.) then
             ueqs(r) = utable(j+1)
          elseif (utable(j) > 0. .and. utable(j+1) == 0.) then
             ueqs(r) = utable(j)
          else
             call fatal('cooling_JoungMacLow_implicit','erroneous equilibrium solutions')
          endif
       endif
    enddo each_root
 endif

 ! Determine the right ueq with the particle's eni
 ! Note: ueqs(1) and ueqs(3) are stable; ueqs(2) is unstable
 if (numroots == 1) then
    ueq_final = ueqs(1)
 else
    if (eni <= ueqs(2)) then
       ueq_final = ueqs(1)
    else
       if (numroots == 2) then
          ueq_final = kboltz*TmaxJML / (gmw*mass_proton_cgs*(gamma-1)) / unit_ergg
!          print*,'going to Tmax',rhoi*unit_density 
       elseif (numroots == 3) then
          ueq_final = ueqs(3)
!          print*,'going to Teq3',rhoi*unit_density 
       endif
    endif
 endif

 ! Isothermal temp of particle
 tempi = gmw*mass_proton_cgs/kboltz*(gamma-1.)*eni * unit_ergg

 ! Compute the new internal energy
 LambdaKI = LambdaKI_coef * lambdacoolJML(tempi)/GammaKI_cgs
 dudti = rhoi*LambdaKI - GammaKI
 deni = eni - ueq_final

 if (abs(deni) > 0.) then
     tau = abs(deni/dudti)
     eni_new = ueq_final + deni*exp(-dt/tau)
     dudti = -(eni-eni_new)/dt
 else
     dudti = -dudti
 endif

end subroutine cooling_JoungMacLow_implicit

!-----------------------------------------------------------------------
!
!   explicit cooling
!
!-----------------------------------------------------------------------
subroutine explicit_cooling (ui, dudt, rho, dt, Trad, mu_in, K2, kappa)
 use eos,     only:gamma,gmw
 use physcon, only:Rg
 use units,   only:unit_ergg
 real, intent(in) :: ui, rho, dt, Trad !code units
 real, intent(in), optional :: mu_in, K2, kappa
 real, intent(out) :: dudt !code units

 real :: u,Q,dlnQ_dlnT,T,mu,T_on_u

 if (.not.present(mu_in)) then
    mu = gmw
 else
    mu = mu_in
 endif
 T_on_u = (gamma-1.)*mu*unit_ergg/Rg
 T = T_on_u*ui
 call calc_cooling_rate(Q, dlnQ_dlnT, rho, T, Trad, mu, K2, kappa)
 if (-Q*dt  > ui) then   ! assume thermal equilibrium
    u = Trad/T_on_u
    dudt = (u-ui)/dt
 else
    dudt = Q
 endif
 !print *,T,Teq,T_on_u*u,'dT=',T_on_u*Q*dt,u,Q*dt

end subroutine explicit_cooling

!-----------------------------------------------------------------------
!
!   implicit cooling
!
!-----------------------------------------------------------------------
subroutine implicit_cooling (ui, dudt, rho, dt, Trad, mu_in, K2, kappa)
 use eos,     only:gamma,gmw
 use physcon, only:Rg
 use units,   only:unit_ergg
 real, intent(in) :: ui, rho, dt
 real, intent(in), optional :: Trad, mu_in, K2, kappa
 real, intent(out) :: dudt

 real, parameter :: tol = 1.d-4 ! to be adjusted
 integer, parameter :: iter_max = 200
 real :: u,Q,dlnQ_dlnT,T,mu,T_on_u,delta_u,term1,term2,term3
 integer :: iter

 if (.not.present(mu_in)) then
    mu = gmw
 else
    mu = mu_in
 endif
 u = ui
 T_on_u = (gamma-1.)*mu*unit_ergg/Rg
 delta_u = 1.d-3
 iter = 0
 !The pdv_work also depends on the internal energy and could also be included
 !in this loop provided this contribution was not accounted for in Force.F90
 ! see PP flag : IMPLICIT COOLING - pb: we need div(v) and it is only real*4
 !term2 = 1.-(gamma-1.)*dt*divcurlv !pdv=(gamma-1.)*vxyzu(4,i)*divcurlv(1,i)*dt
 term2 = 1.
 term1 = u !initial internal energy without cooling contributions
 do while (abs(delta_u) > tol .and. iter < iter_max)
    T = u*T_on_u
    call calc_cooling_rate(Q,dlnQ_dlnT, rho, T, Trad, mu, K2, kappa)
    term3 = u*term2-Q*dt
    delta_u = (term1-term3)/(term2-Q*dlnQ_dlnT*dt/u)
    u = u+delta_u
    iter = iter + 1
 enddo
 dudt =(u-term1)/dt
 if (u < 0. .or. isnan(u)) then
    print *,u
    stop ' u<0'
 endif

end subroutine implicit_cooling

!-----------------------------------------------------------------------
!
!   this routine returns the effective cooling rate du/dt
!
!-----------------------------------------------------------------------
subroutine energ_cooling(xi,yi,zi,ui,dudt,rho,dt,Trad,mu_in,K2,kappa,Tgas)
 use io,   only: fatal
 real, intent(in)           :: xi,yi,zi,ui,rho,dt         ! in code units
 real, intent(in), optional :: Tgas,Trad,mu_in,K2,kappa   ! in cgs units
 real, intent(inout)        :: dudt                       ! in code units
 integer :: icf
 real    :: dist(3),magdist,old_dudt

 select case (icooling)
 case (3)
    call cooling_Gammie(xi,yi,zi,ui,dudt)
 case (2)
    call exact_cooling_table(ui,rho,dt,dudt)
 case (5)
    if (present(Tgas)) then
       call cooling_KoyamaInutuska_explicit(rho,Tgas,dudt)
    else
       call fatal('energ_cooling','Koyama & Inutuska cooling requires gas temperature')
    endif
 case (6)
    call cooling_KoyamaInutuska_implicit(ui,rho,dt,dudt)
 case (7)
    call cooling_JoungMacLow_implicit(ui,rho,dt,dudt)
 case default
    !call exact_cooling(u, dudt, rho, dt, Trad, mu_in, K2, kappa)
    !call implicit_cooling(u, dudt, rho, dt, Trad, mu_in, K2, kappa)
    if (present(Trad) .and. present(mu_in) .and. present(K2) .and. present(kappa)) then
       call explicit_cooling(ui, dudt, rho, dt, Trad, mu_in, K2, kappa)
    else
       call fatal('energ_cooling','default requires optional arguments; change icooling or ask D Price or L Siess to patch')
    endif
 end select
 !
 ! Temporarily reduce cooling for particles near SNe
 !
 if (snecoolingoff .and. ncoolingoff >= 1) then
!    if (xyzh_coolingoff(1,1)>0.) print*, 'xyzh_coolingoff in cooling' , xyzh_coolingoff
    each_coord: do icf = 1,ncoolingoff
       dist(1:3) = xyzh_coolingoff(1:3,icf) - (/ xi,yi,zi /)
       magdist = sqrt(dist(1)**2 + dist(2)**2 + dist(3)**2)
       if (range_cooloff_useh) range_cooloff = 2*xyzh_coolingoff(4,icf)  ! full radius = 2h
       if (magdist <= range_cooloff) then
          old_dudt = dudt
          ! Scale down dudt by a fraction depending on distance from sn
          dudt = dudt - dudt*cooloff_frac(magdist,range_cooloff)
          print*,'reducing cooling: old dudt, new dudt', old_dudt,dudt
       endif
    enddo each_coord
 endif

end subroutine energ_cooling


!
! Compute the fraction of dudtcool to reduce around supernova
! Smoothed using cubic spline to avoid discontinuities caused by switching off cooling
!
real function cooloff_frac(dist,totrange)
 real, intent(in) :: dist,totrange
 real :: r

 r = dist/totrange*2.  ! Scale to range of cubic kernel
 if (r >= 0. .and. r < 1.) then
    cooloff_frac = 1 - 3./2.*r**2 + 3./4.*r**3
 elseif (r >= 1. .and. r < 2.) then
    cooloff_frac = 1./4.*(2-r)**3
 elseif (r >= 2) then
    cooloff_frac = 0.
 endif

end function cooloff_frac

!-----------------------------------------------------------------------
!
!   cooling using Townsend (2009), ApJS 181, 391-397 method with
!   analytical cooling rate prescriptions
!
!-----------------------------------------------------------------------
subroutine exact_cooling (u, dudt, rho, dt, Trad, mu_in, K2, kappa)
 use eos,     only:gamma,gmw
 use physcon, only:Rg
 use units,   only:unit_ergg
 real, intent(in) :: u, rho, dt, Trad
 real, intent(in), optional :: mu_in, K2, kappa
 real, intent(out) :: dudt

 real, parameter :: tol = 1.d-12
 real :: Qref,dlnQref_dlnT,Q,dlnQ_dlnT,Y,Yk,Yinv,Temp,dy,T,mu,T_on_u
 integer :: k

 if (.not.present(mu_in)) then
    mu = gmw
 else
    mu = mu_in
 endif
 T_on_u = (gamma-1.)*mu*unit_ergg/Rg
 T = T_on_u*u

 if (T < T_floor) then
    Temp = T_floor
 elseif (T > Tref) then
    call calc_cooling_rate(Q, dlnQ_dlnT, rho, T, Trad, mu, K2, kappa)
    Temp = T+T_on_u*Q*dt
 else
    call calc_cooling_rate(Qref,dlnQref_dlnT, rho, Tref, Trad, mu, K2, kappa)
    Y = 0.
    k = nTg
    Q = Qref                  ! default value if Tgrid < T for all k
    dlnQ_dlnT = dlnQref_dlnT  ! default value if Tgrid < T for all k
    do while (Tgrid(k) > T)
       k = k-1
       call calc_cooling_rate(Q, dlnQ_dlnT, rho, Tgrid(k), Trad, mu, K2, kappa)
       ! eqs A6
       if (abs(dlnQ_dlnT-1.) < tol) then
          y = y - Qref*Tgrid(k)/(Q*Tref)*log(Tgrid(k)/Tgrid(k+1))
       else
          y = y - Qref*Tgrid(k)/(Q*Tref*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/Tgrid(k+1))**(dlnQ_dlnT-1.))
       endif
    enddo
    !eqs A5
    yk = y
    if (abs(dlnQ_dlnT-1.) < tol) then
       y = yk + Qref*Tgrid(k)/(Q*Tref)*log(Tgrid(k)/T)
    else
       y = yk + Qref*Tgrid(k)/((Q*Tref)*(1.-dlnQ_dlnT))*(1.-(Tgrid(k)/T)**(dlnQ_dlnT-1))
    endif
    !eq 26
    dy = Qref*dt*T_on_u/Tref
    y = y + dy
    !compute Yinv (eqs A7)
    if (abs(dlnQ_dlnT-1.) < tol) then
       Temp = max(Tgrid(k)*exp(-Q*Tref*(y-yk)/(Qref*Tgrid(k))),T_floor)
    else
       Yinv = 1.-(1.-dlnQ_dlnT)*Q*Tref/(Qref*Tgrid(k))*(y-yk)
       if (Yinv > 0.) then
          Temp = Tgrid(k)*(Yinv**(1./(1.-dlnQ_dlnT)))
       else
          Temp = T_floor
       endif
    endif
 endif

 dudt = (Temp-T)/T_on_u/dt
 !note that u = Temp/T_on_u

end subroutine exact_cooling

!-----------------------------------------------------------------------
!+
!  cooling using Townsend (2009) method with tabulated rate
!
!   Implements cooling defined using a tabulated cooling table
!   produced e.g. by CLOUDY.
!+
!-----------------------------------------------------------------------
subroutine exact_cooling_table(uu,rho,dt,dudt)
 use eos,     only:gamma,gmw
 use physcon, only:atomic_mass_unit,kboltz,Rg
 use units,   only:unit_density,unit_ergg,utime
 real, intent(in)  :: uu, rho,dt
 real, intent(out) :: dudt
 real    :: gam1,density_cgs,dt_cgs,amue,amuh,dtemp
 real    :: sloperef,slopek,temp,temp1,tref,yfunx,yinv0
 integer :: k

 gam1 = gamma - 1.
 temp = gam1*uu/Rg*gmw*unit_ergg

 tref     = temper(nt)
 sloperef = slope(nt)

 if (temp < temp_floor) then
    temp1 = temp_floor
 else
    amue = 2.*atomic_mass_unit/(1. + habund)
    amuh = atomic_mass_unit/habund
    density_cgs = rho*unit_density
    dt_cgs      = dt*utime

    !Lionel Siess : I think there is an error. in dtemp sloperef should be replaced by lambda(nt)
    !original dtemp = gam1*density_cgs*(atomic_mass_unit*gmw/(amue*amuh*kboltz))* &
    !     sloperef/tref*dt_cgs
    ! Eq 26
    dtemp = gam1*density_cgs*(atomic_mass_unit*gmw/(amue*amuh*kboltz))*lambda(nt)/tref*dt_cgs

    k = find_in_table(nt,temper,temp)
    slopek = slope(k)
    ! Eq A5
    if (abs(slopek - 1.) < tiny(0.)) then
       yfunx = yfunc(k) + lambda(nt)*temper(k)/(lambda(k)*temper(nt))*log(temper(k)/temp)
    else
       yfunx = yfunc(k) + lambda(nt)*temper(k)/(lambda(k)*temper(nt)*(1. - slopek)) &
                          *(1. - (temper(k)/temp)**(slopek-1.))
    endif
    yfunx = yfunx + dtemp
    ! Eq A7
    if (abs(slopek - 1.) < tiny(0.)) then
       temp1 = max(temper(k)*exp(-lambda(k)*temper(nt)/(lambda(nt)*temper(k))*(yfunx-yfunc(k))),temp_floor)
    else
       yinv0 = 1. - (1. - slopek)*lambda(k)*temper(nt)/(lambda(nt)*temper(k))*(yfunx-yfunc(k))
       if (yinv0 > 0.) then
          temp1 = max(temper(k)*yinv0**(1./(1. - slopek)),temp_floor)
       else
          temp1 = temp_floor
       endif
    endif
 endif

 dudt = (temp1 - temp)*Rg/(gam1*gmw*unit_ergg)/dt

end subroutine exact_cooling_table

!-----------------------------------------------------------------------
!+
!  utility to find the index of closest value in a table
!+
!-----------------------------------------------------------------------
pure integer function find_in_table(n,table,val) result(i)
 integer, intent(in) :: n
 real,    intent(in) :: table(n), val
 integer :: i0,i1

 i0 = 0
 i1 = n + 1
 do while (i1 - i0 > 1)
    i = (i0 + i1)/2
    if ((table(n) >= table(1)).eqv.(val >= table(i))) then
       i0 = i
    else
       i1 = i
    endif
 enddo
 if (abs(val-table(1)) < tiny(0.)) then
    i = 1
 elseif (abs(val-table(n)) < tiny(0.)) then
    i = n-1
 else
    i = i0
 endif

end function find_in_table

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_cooling(iunit)
 use infile_utils, only:write_inopt
 use h2cooling,    only:write_options_h2cooling
 use part,         only:h2chemistry
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling cooling'
 call write_inopt(C_cool,'C_cool','factor controlling cooling timestep',iunit)
 if (h2chemistry) then
    call write_inopt(icooling,'icooling','cooling function (0=off, 1=on)',iunit)
    if (icooling > 0) then
       call write_options_h2cooling(iunit)
    endif
 else
    call write_inopt(icooling,'icooling','cooling function (0=off, 1=explicit, 2=Townsend table, 3=Gammie,&
                   & 5=KI02 explicit, 6=KI02 implicit, 7=JML06 implicit)',iunit)
    select case(icooling)
    case(1)
       call write_inopt(icool_radiation_H0,'icool_radiation_H0','H0 cooling on/off',iunit)
       call write_inopt(icool_relax_bowen,'icool_relax_bowen','Bowen (diffusive) relaxation on/off',iunit)
       call write_inopt(icool_relax_stefan,'icool_relax_stefan','radiative relaxation on/off',iunit)
       call write_inopt(icool_dust_collision,'icool_dust_collision','dust collision on/off',iunit)
       call write_inopt(bowen_Cprime,'bowen_Cprime','radiative cooling rate (g.s/cm³)',iunit)
    case(2)
       call write_inopt(cooltable,'cooltable','data file containing cooling function',iunit)
       call write_inopt(habund,'habund','Hydrogen abundance assumed in cooling function',iunit)
       call write_inopt(temp_floor,'temp_floor','Minimum allowed temperature in K for Townsend cooling table',iunit)
    case(3)
       call write_inopt(beta_cool,'beta_cool','beta factor in Gammie (2001) cooling',iunit)
    end select
 endif
 if (icooling > 0) call write_inopt(Tfloor,'Tfloor','temperature floor (K); on if > 0',iunit)
 if (ufloor > 0.) call write_inopt(ufloor,'ufloor','internal energy floor; on if > 0',iunit)

end subroutine write_options_cooling

!-----------------------------------------------------------------------
!+
!  reads sink particle options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_cooling(name,valstring,imatch,igotall,ierr)
 use part,         only:h2chemistry
 use h2cooling,    only:read_options_h2cooling
 use io,           only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 logical :: igotallh2,igotallcf

 imatch  = .true.
 igotall = .false.  ! cooling options are compulsory
 igotallh2 = .true.
 igotallcf = .true.
 select case(trim(name))
 case('icooling')
    read(valstring,*,iostat=ierr) icooling
    ngot = ngot + 1
 case('icool_radiation_H0')
    read(valstring,*,iostat=ierr) icool_radiation_H0
    ngot = ngot + 1
 case('icool_relax_bowen')
    read(valstring,*,iostat=ierr) icool_relax_bowen
    ngot = ngot + 1
 case('icool_relax_stefan')
    read(valstring,*,iostat=ierr) icool_relax_stefan
    ngot = ngot + 1
 case('icool_dust_collision')
    read(valstring,*,iostat=ierr) icool_dust_collision
    ngot = ngot + 1
 case('C_cool')
    read(valstring,*,iostat=ierr) C_cool
    ngot = ngot + 1
 case('cooltable')
    read(valstring,*,iostat=ierr) cooltable
    ngot = ngot + 1
 case('habund')
    read(valstring,*,iostat=ierr) habund
    ngot = ngot + 1
 case('temp_floor')
    read(valstring,*,iostat=ierr) temp_floor
    ngot = ngot + 1
 case('bowen_Cprime')
    read(valstring,*,iostat=ierr) bowen_Cprime
    ngot = ngot + 1
 case('beta_cool')
    read(valstring,*,iostat=ierr) beta_cool
    ngot = ngot + 1
    if (beta_cool < 1.) call fatal('read_options','beta_cool must be >= 1')
 case('Tfloor')
    ! not compulsory to read in
    read(valstring,*,iostat=ierr) Tfloor
case('ufloor')
    read(valstring,*,iostat=ierr) ufloor
 case default
    imatch = .false.
    if (h2chemistry) then
       call read_options_h2cooling(name,valstring,imatch,igotallh2,ierr)
    endif
 end select
 if (icooling == 3 .and. ngot >= 1) igotall = .true.
 if (icooling == 2 .and. ngot >= 3) igotall = .true.
 if (icooling == 1 .and. ngot >= 5) igotall = .true.
 if (igotallh2 .and. ngot >= 1) igotall = .true.

end subroutine read_options_cooling

end module cooling
