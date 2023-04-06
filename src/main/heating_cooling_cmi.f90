!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module heatcool_cmi
!
! The CMI suite: photoionize_cmi.f90 kdtree_cmi.f90 hnode_cmi.f90 *heating_cooling_cmi.f90*
! This module contains all the subroutines necessary for doing the heating and cooling
! processes involved in photoionization.
!
! :References: Koyama & Inutsuka (2002), ApJL 564, 97-100
!              Vazquez-Semadeni,et.al (2007), ApJ 657, 870-883
!              Joung & Mac Low (2006), ApJ 653, 1266-1279
!              Osterbrock (1974), Astrophysics of Gaseous Nebulae, W.H. Freeman and Company
!
! :Owner: Cheryl Lau
!
! :Runtime parameters:
!   - rhomin/max_cgs       : *range of densities of which equifunc has at least one solution*
!   - gammamin/max_cgs     : *range of photoionization heating rate*
!   - Tmin/max             : *temperature range of cooling curve*
!   - use_const_recombcoef : *adopt a constant recombination coef, otherwise it varies with temp*
!   - get_ueq_by_interp    : *perform 2D interpolation, otherwise simply grabs closest ueq*
!   - incl_recomb_cooling  : *include recombination cooling in heating term gamma*
!
! :Dependencies: io,eos,physcon,units,parts,cooling
!
! :Note: The variable gamma in this mod refers to heating term, not to be confused with gamma
!        (polytropic constant) from eos. In subroutines where both are used, gamma (heating term)
!        will be referred to as gammaheat.
!
!
 use cooling, only:lambdacoolJML,lambdacoolDR  !- cooling function lambda(T) [erg cm^3 s^-1]

 implicit none

 public :: check_to_stop_cooling,init_ueq_table,precompute_uterms
 public :: compute_heating_term,compute_cooling_term,get_ueq,compute_du,compute_dudt

 logical, public :: stop_cooling  !- flag to stop cooling at current step

 !- Options for controlling heating/cooling physics
 logical, public :: use_const_alphaA  = .false.
 logical, public :: get_ueq_by_interp = .true.

 !- Heating from cosmic rays, X-rays, H2 formation and destruction etc. (as of KI02)
 real,    public :: gamma_background_cgs = 2E-26

 private

 !- Pre-computed table of equilibiurms ueq(rho,gamma)
 integer, parameter :: maxrho   = 1000
 integer, parameter :: maxgamma = 1000
 real   :: rho_gamma_ueq_table(maxrho,maxgamma,6)  ! stores: rho,gamma,numroots,ueq1,ueq2,ueq3

 real   :: rhomin_cgs   = 7E-29
 real   :: rhomax_cgs   = 1E-19
 real   :: gammamin_cgs = 2E-26   !- min background heating
 real   :: gammamax_cgs = 1E-18
 real   :: Tmin  = 1E1
 real   :: Tmax  = 1E8

 !- Pre-computed vars
 real   :: one_over_mH,one_over_mH2,one_over_mH_cgs,one_over_mH2_cgs
 real   :: one_over_unit_gamma,one_over_unit_lambda
 real   :: mass_proton,pmass,pmass_cgs

contains

!
! This module already handles all heating and cooling processes when photoionization
! is being switched on; flag to temporarily disable orginal cooling routines in code
!
subroutine check_to_stop_cooling(nphotosrc)
 integer, intent(in) :: nphotosrc

 stop_cooling = .false.
 if (nphotosrc >= 1) stop_cooling = .true.

end subroutine check_to_stop_cooling

!
! Pre-compute some terms involved in du/dudt-update for speeding up
!
subroutine precompute_uterms
 use physcon, only:mass_proton_cgs
 use part,    only:massoftype,igas
 use units,   only:umass,utime,udist,unit_energ

 one_over_mH_cgs  = 1./mass_proton_cgs
 one_over_mH2_cgs = 1./mass_proton_cgs**2

 mass_proton  = mass_proton_cgs/umass
 one_over_mH  = 1./mass_proton
 one_over_mH2 = 1./mass_proton**2

 pmass     = massoftype(igas)
 pmass_cgs = pmass*umass

 one_over_unit_gamma  = 1./(unit_energ/utime)
 one_over_unit_lambda = 1./(unit_energ*udist**3/utime)

end subroutine precompute_uterms

!
! Pre-solves the equilibriums ueq as function of both rho and gamma
! and organize into a 3D table: rho_gamma_ueq_table(irho,igamma,(rho,gamma,numroot,ueq123))
!
subroutine init_ueq_table
 use physcon, only:mass_proton_cgs,kboltz
 use units,   only:utime,unit_density,unit_ergg
 use eos,     only:gmw,gamma
 use io,      only:fatal
 integer :: igamma,irho,numroots,irterr1,irterr2,irterr3
 real    :: Tlocalmax,Tlocalmin,gamma_cgs,dgamma_cgs,drho_cgs,rho_cgs,nrho_cgs
 real    :: rho,gammaheat
 real    :: Teq1,Teq2,Teq3,ueq1,ueq2,ueq3,gmw0
 logical :: write_ueq = .true.
 logical :: write_Teq = .true.

 print*,'Pre-computing thermal equilibrium solutions'

 gmw0 = gmw
 gmw = 0.5   !- temporarily changing mean molecular weight

 ! Temp of local max and min of JML06 cooling curve (for root-search brackets)
 Tlocalmax = 198609.4917357372
 Tlocalmin = 31622776.6016837933

 dgamma_cgs = (log10(gammamax_cgs)-log10(gammamin_cgs))/maxgamma
 drho_cgs   = (log10(rhomax_cgs)-log10(rhomin_cgs))/maxrho

 if (write_ueq) then
    open(2035,file='rho_gamma_ueqs.txt')
    write(2035,'(6a20)') 'rho','gamma','numroot','ueq1','ueq2','ueq3'
 endif
 if (write_Teq) then
    open(2036,file='rho_gamma_Teqs.txt')
    write(2036,'(6a20)') 'rho','gamma','numroot','Teq1','Teq2','Teq3'
 endif

 over_gamma: do igamma = 1,maxgamma
    gamma_cgs = gammamin_cgs * 10**((igamma-1)*dgamma_cgs)

    over_rho: do irho = 1,maxrho
       rho_cgs  = rhomin_cgs * 10**((irho-1)*drho_cgs)
       nrho_cgs = rho_cgs/mass_proton_cgs

       call root_bisection(nrho_cgs,Tmin,Tlocalmax,gamma_cgs,Teq1,irterr1)
       call root_bisection(nrho_cgs,Tlocalmax,Tlocalmin,gamma_cgs,Teq2,irterr2)
       call root_bisection(nrho_cgs,Tlocalmin,Tmax,gamma_cgs,Teq3,irterr3)

       ! Get total number of roots (irterr=1 means no roots found)
       if (irterr3 == 1) then
          if (irterr2 == 1) then
             if (irterr1 /= 1) then
                numroots = 1
             endif
          else
             numroots = 2
          endif
       else
          numroots = 3
       endif

       !- convert to code units and store
       rho   = rho_cgs/unit_density
       gammaheat = gamma_cgs*one_over_unit_gamma
       ueq1  = kboltz * Teq1 / (gmw*mass_proton_cgs*(gamma-1.)) /unit_ergg  ! gamma here is polytrope const
       ueq2  = kboltz * Teq2 / (gmw*mass_proton_cgs*(gamma-1.)) /unit_ergg
       ueq3  = kboltz * Teq3 / (gmw*mass_proton_cgs*(gamma-1.)) /unit_ergg

       rho_gamma_ueq_table(irho,igamma,1) = rho
       rho_gamma_ueq_table(irho,igamma,2) = gammaheat
       rho_gamma_ueq_table(irho,igamma,3) = numroots
       rho_gamma_ueq_table(irho,igamma,4) = ueq1
       rho_gamma_ueq_table(irho,igamma,5) = ueq2
       rho_gamma_ueq_table(irho,igamma,6) = ueq3

       if (write_ueq) write(2035,*) rho_cgs,gamma_cgs,numroots,ueq1*unit_ergg,ueq2*unit_ergg,ueq3*unit_ergg
       if (write_Teq) write(2036,*) rho_cgs,gamma_cgs,numroots,Teq1,Teq2,Teq3

    enddo over_rho
 enddo over_gamma

 gmw = gmw0

 if (write_ueq) close(2035)
 if (write_Teq) close(2036)

end subroutine init_ueq_table

!
! Thermal equilibium: f_T = n*lambda(T) - gamma
!
real function equifunc(nrho_cgs,temp,gamma_cgs)
 real, intent(in) :: nrho_cgs,temp,gamma_cgs

 equifunc = nrho_cgs*lambdacoolDR(temp) - gamma_cgs

end function equifunc

!
! Brackets the equilibium solution Teq with given Tmin and Tmax;
! returns Teq = 0. and irterr = 1 if no roots found
!
subroutine root_bisection(nrho,Tmin0,Tmax0,gamma,Teq,irterr)
 use io,  only:fatal
 real,    intent(in)  :: nrho,Tmin0,Tmax0,gamma
 integer, intent(out) :: irterr
 real,    intent(out) :: Teq
 integer :: niter
 real    :: Tmin,Tmax,Tmid,Tminsign,Tmaxsign,func_min,func_max,func_mid
 real    :: tol = 0.05   ! tolerance in root-search
 logical :: converged

 converged = .false.  ! root flag
 irterr = 0           ! no-root indicator
 niter = 0

 Tminsign = sign(1.,equifunc(nrho,Tmin0,gamma))
 Tmaxsign = sign(1.,equifunc(nrho,Tmax0,gamma))
 if (Tminsign == Tmaxsign) then
    Teq = 0.
    irterr = 1
 endif

 Tmin = Tmin0
 Tmax = Tmax0
 do while (.not.converged .and. irterr == 0)

    func_max = equifunc(nrho,Tmax,gamma)
    func_min = equifunc(nrho,Tmin,gamma)
    if (func_max*func_min > 0) call fatal('heating_cooling_cmi','not bracketing Teq root')

    Tmid = (Tmax+Tmin)/2.
    func_mid = equifunc(nrho,Tmid,gamma)

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

!
! Compute photoionization heating rate gamma in code units using nH from CMI
! Note that Osterbrock74 G(H) = n*gamma
!
subroutine compute_heating_term(nH,rho,u,totlumin,nphotosrc,freq_photon,gammaheat)
 use physcon, only:mass_proton_cgs,mass_electron_cgs,kboltz,steboltz,pi,planckh
 use units,   only:umass,utime,unit_density,unit_ergg,unit_energ
 use eos,     only:gamma,gmw
 use io,      only:fatal
 integer, intent(in)  :: nphotosrc
 real,    intent(in)  :: nH,rho,u
 real,    intent(in)  :: totlumin,freq_photon
 real,    intent(out) :: gammaheat
 real    :: num_neu,num_ion,num_e,num_p
 real    :: rho_cgs,nrho_cgs,vpart_cgs,temp,alphaA,betaA,Ne,Np
 real    :: energ_photon_cgs,lumin_cgs,temp_star
 real    :: heating_rate_cgs,recomb_cooling_rate_cgs,gamma_cgs
 real    :: rstar_cgs = 6.957E11 !- Radius of a typical O-star

 rho_cgs  = rho*unit_density
 nrho_cgs = rho_cgs/mass_proton_cgs

 !- number density of electrons and protons
 Ne = (1.-nH)*nrho_cgs
 Np = (1.-nH)*nrho_cgs

 !- current temp from u
 temp = u/kboltz*(gmw*mass_proton_cgs*(gamma-1.))*unit_ergg

 !- compute approx alphaA
 alphaA = get_alphaA(temp)

 !- temp of newly created photoelectrons (~ temp of star)
 energ_photon_cgs = planckh*freq_photon
 lumin_cgs = totlumin/nphotosrc*energ_photon_cgs   ! [erg s^-1]
 temp_star = (lumin_cgs/(4.*pi*rstar_cgs**2*steboltz))**(1./4.)

 !- Photoionization heating G(H) as of Osterbrock74 eqn 3.2
 heating_rate_cgs = Ne*Np*alphaA*(3./2.)*kboltz*temp_star   ! [erg cm^-3 s-1]

 !- gamma
 nrho_cgs  = rho_cgs*one_over_mH_cgs
 gamma_cgs = heating_rate_cgs/nrho_cgs  ! [erg s-1]

 !- Include background heating
 gamma_cgs = gamma_cgs + gamma_background_cgs

 if (gamma_cgs < gammamin_cgs .or. gamma_cgs > gammamax_cgs) then
    print*,'gamma_cgs',gamma_cgs
    call fatal('heating_cooling_cmi','gamma_cgs exceeded range')
 endif

 !- convert to code units
 gammaheat = gamma_cgs*one_over_unit_gamma

end subroutine compute_heating_term

!
! Recombination coef alphaA as func of temp obtained from fitting table 2.1
! of Osterbrock74 [cm^3 s^-1]
!
real function get_alphaA(temp)
 real, intent(in) :: temp
 real :: a = 6.113723867E-11
 real :: b = -1.85763022E-13

 if (use_const_alphaA) then
    get_alphaA = 2.7e-13
 else
    get_alphaA = a* temp**(-1./2.) + b  !- P.16: recomb coef ~ T^(-1/2)
 endif

end function get_alphaA


!
! Routine to read cooling rate off temp directly from cooling curve
!
subroutine compute_cooling_term(u,lambda)
 use physcon, only:mass_proton_cgs,kboltz
 use eos,     only:gmw,gamma
 use units,   only:unit_ergg
 real, intent(in)  :: u
 real, intent(out) :: lambda
 real :: temp,lambda_cgs

 temp = gmw*mass_proton_cgs/kboltz*(gamma-1.)*u *unit_ergg
 lambda_cgs = lambdacoolDR(temp)  ! [erg cm^3 s^-1]

 !- convert to code units
 lambda = lambda_cgs*one_over_unit_lambda

end subroutine compute_cooling_term

!
! Extract ueq(s) from pre-computed table with the given rho and gamma
! Selects the right ueq using the input vxyzu
!
subroutine get_ueq(rho,gammaheat,u,ueq_final)
 use physcon, only:kboltz,mass_proton_cgs
 use units,   only:unit_ergg,unit_energ,utime
 use eos,     only:gmw,gamma
 use io,      only:fatal
 real, intent(in)  :: rho,gammaheat,u
 real, intent(out) :: ueq_final
 integer :: irho,igamma,jrho,jgamma,numroots,r
 real    :: rhotable(maxrho),gammatable(maxgamma)  !- unpacked from rho_gamma_ueq_table
 real    :: ueqs(3)                                !- holds the 3 possible ueq solutions
 real    :: x1,x2,y1,y2,q11,q12,q21,q22            !- for 2D interpolation

 rhotable   = rho_gamma_ueq_table(:,1,1)
 gammatable = rho_gamma_ueq_table(1,:,2)

 !print*,'gamma (code units)',gammaheat

 !- Closest index
 irho   = minloc(abs(rhotable(:)-rho),1)
 igamma = minloc(abs(gammatable(:)-gammaheat),1)

 numroots = int(rho_gamma_ueq_table(irho,igamma,3))
 if (numroots == 0) call fatal('cooling_heating_cmi','no equilibrium solution found')
 ueqs = (/ 0.,0.,0. /)

 if (get_ueq_by_interp) then
    ! Note: [i] closest index; [j] lower bound; [j+1] upper bound bracketing input var

    ! get bracketing indices [j] and [j+1] for rho
    if (irho == 1) then
       jrho = 1
    elseif (irho == maxrho) then
       jrho = irho-1
    elseif (rhotable(irho) >= rho .and. rhotable(irho-1) < rho) then
       jrho = irho-1
    elseif (rhotable(irho) < rho .and. rhotable(irho+1) >= rho) then
       jrho = irho
    endif
    ! get bracketing indices [j] and [j+1] for gamma
    if (igamma == 1) then
       jgamma = 1
    elseif (igamma == maxgamma) then
       jgamma = igamma-1
    elseif (gammatable(igamma) >= gammaheat .and. gammatable(igamma-1) < gammaheat) then
       jgamma = igamma-1
    elseif (gammatable(igamma) < gammaheat .and. gammatable(igamma+1) >= gammaheat) then
       jgamma = igamma
    endif

    ! Bilinear interpolation
    !- Interpolate only if all four bracketing points has a solution; else use nearest root
    x1 = rhotable(jrho)
    x2 = rhotable(jrho+1)
    y1 = gammatable(jgamma)
    y2 = gammatable(jgamma+1)
    each_root: do r = 1,numroots
       q11 = rho_gamma_ueq_table(jrho,jgamma,3+r)
       q12 = rho_gamma_ueq_table(jrho,jgamma+1,3+r)
       q21 = rho_gamma_ueq_table(jrho+1,jgamma,3+r)
       q22 = rho_gamma_ueq_table(jrho+1,jgamma+1,3+r)
       if (q11 > 0. .and. q12 > 0. .and. q21 > 0. .and. q22 > 0.) then
          ueqs(r) = 1./((x2-x1)*(y2-y1)) * ( q11*(x2-rho)*(y2-gammaheat) + q21*(rho-x1)*(y2-gammaheat) + &
                    q12*(x2-rho)*(gammaheat-y1) + q22*(rho-x1)*(gammaheat-y1) )
       else
          ueqs(r) = rho_gamma_ueq_table(irho,igamma,3+r)
       endif
    enddo each_root

 else !- simply use the closest solution to save comp time
    do r = 1,numroots
       ueqs(r) = rho_gamma_ueq_table(irho,igamma,3+r)
    enddo
 endif

 ! Select the right ueq using the input u
 ! Note: ueqs(1) and ueqs(3) are stable; ueqs(2) is unstable
 if (numroots == 1) then
    ueq_final = ueqs(1)
 else
    if (u <= ueqs(2)) then
       ueq_final = ueqs(1)
    else
       if (numroots == 2) then
          ueq_final = kboltz*Tmax / (gmw*mass_proton_cgs*(gamma-1)) / unit_ergg
       elseif (numroots == 3) then
          ueq_final = ueqs(3)
       endif
    endif
 endif

end subroutine get_ueq

!
! With the ueq, compute new u after dt to give du in code units
! method based on Vazquez-Semadeni07
!
subroutine compute_du(dt,rho,u,ueq,gamma,lambda,du)
 real, intent(in)  :: dt,rho,u,ueq,gamma,lambda
 real, intent(out) :: du
 real :: unew,tau

 !- timescale to radiate thermal energy excess
 tau  = abs((u-ueq)/(rho*one_over_mH2*lambda - one_over_mH*gamma))

 unew = ueq + (u-ueq)*exp(-dt/tau)
 du   = unew - u

end subroutine compute_du

!
! Balance heating and cooling term to give dudt in code units
!
subroutine compute_dudt(rho,gamma,lambda,dudt)
 real, intent(in)  :: rho,gamma,lambda
 real, intent(out) :: dudt

 dudt = one_over_mH*gamma - rho*one_over_mH2*lambda

end subroutine compute_dudt

! Note:
! According to KI02 and VS07, it should be nrho*gamma - nrho^2*lambda, where nrho = rho/mass_proton.
! But since u is in erg/g, we further divide nrho by rho to make the equations dimensionally correct,
! giving dudt = (1/mass_proton)*gamma - (rho/mass_proton**2)*lambda; [erg/g/s].

end module heatcool_cmi
