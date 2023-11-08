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
! :References: Osterbrock (1974), Astrophysics of Gaseous Nebulae, W.H. Freeman and Company
!              Vazquez-Semadeni,et.al (2007), ApJ 657, 870-883
!              Koyama & Inutsuka (2002), ApJL 564, 97-100
!              De Rijcke,et.al (2013), MNRAS 433, 3005-3016
!
! :Owner: Cheryl Lau
!
! :Runtime parameters:
!   - rhomin/max_cgs       : *range of densities of which equifunc has at least one solution*
!   - gammamin/max_cgs     : *range of photoionization heating rate*
!   - Tmin/max             : *temperature range of cooling curve*
!   - use_const_alpha      : *adopt a constant recombination coef, otherwise it varies with temp*
!   - get_ueq_by_interp    : *perform 2D interpolation, otherwise simply grabs closest ueq*
!   - gamma_background_cgs : *background heating rate*
!
! :Dependencies: io, eos, physcon, units, parts, kdtree, linklist
!
! :Note: The variable gamma in this mod refers to heating term, not to be confused with gamma
!        (polytropic constant) from eos. In subroutines where both are used, gamma (heating term)
!        will be referred to as gammaheat.
!
!
 implicit none

 public :: check_to_stop_cooling,init_ueq_table,precompute_uterms
 public :: heating_term,cooling_term,get_ueq,compute_du,compute_dudt
 public :: compute_Rtype_time

 logical, public :: stop_cooling  !- flag to stop cooling at current step

 ! Options for controlling heating/cooling physics
 logical, public :: use_const_alpha   = .true.
 logical, public :: get_ueq_by_interp = .true.

 ! Heating from cosmic rays, X-rays, H2 formation and destruction etc. (as of KI02)
 real,    public :: gamma_background_cgs = 2E-26

 private

 ! Pre-computed table of equilibrium ueq(rho,gamma)
 integer, parameter :: maxrho   = 10000
 integer, parameter :: maxgamma = 1000
 real   :: rho_gamma_ueq_table(maxrho,maxgamma,6)  ! stores: rho,gamma,numroots,ueq1,ueq2,ueq3

 real   :: rhomin_cgs   = 1E-26
 real   :: rhomax_cgs   = 1E-13
 real   :: gammamin_cgs = 1E-27
 real   :: gammamax_cgs = 1E-11
 real   :: Tmin  = 1E0
 real   :: Tmax  = 1E9

 ! Pre-computed vars
 real   :: one_over_mH,one_over_mH2,one_over_mH_cgs,one_over_mH2_cgs
 real   :: one_over_unit_gamma,one_over_unit_lambda
 real   :: mass_proton,pmass,pmass_cgs

 ! Switches for plotting/debugging
 logical :: write_ueq = .false.
 logical :: write_Teq = .true.

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
! Creates the 3D rho_gamma_ueq_table(irho,igamma,(rho,gamma,numroot,ueq1,ueq2,ueq3))
!
subroutine init_ueq_table
 use physcon, only:mass_proton_cgs,kboltz
 use units,   only:unit_density,unit_ergg
 use eos,     only:gmw,gamma
 use io,      only:fatal
 integer :: igamma,irho,numroots,irterr1,irterr2,irterr3
 real    :: Tlocalmax,Tlocalmin,gamma_cgs,dgamma_cgs,drho_cgs,rho_cgs,nrho_cgs
 real    :: rho,gammaheat
 real    :: Teq1,Teq2,Teq3,ueq1,ueq2,ueq3,gmw0

 print*,'Pre-computing thermal equilibrium solutions'

 gmw0 = gmw
 gmw = 0.5   !- temporarily changing mean molecular weight

 ! Temp of local max and min of DeRijcke13 cooling curve (for root-search brackets)
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

       call solve_temp_equil(nrho_cgs,Tmin,Tlocalmax,gamma_cgs,Teq1,irterr1)
       call solve_temp_equil(nrho_cgs,Tlocalmax,Tlocalmin,gamma_cgs,Teq2,irterr2)
       call solve_temp_equil(nrho_cgs,Tlocalmin,Tmax,gamma_cgs,Teq3,irterr3)

       ! Get total number of roots (irterr=1 means no roots found)
       if (irterr3 == 1) then
          if (irterr2 == 1) then
             if (irterr1 == 1) then
                numroots = 0
             else
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
! Thermal equilibium: f_T = n*gamma - n^2*lambda(T)
!
real function equifunc(nrho_cgs,temp,gamma_cgs)
 real, intent(in) :: nrho_cgs,temp,gamma_cgs

 equifunc = nrho_cgs*gamma_cgs - nrho_cgs**2*lambdacoolDR(temp)

end function equifunc

!
! Brackets the equilibium solution Teq with given Tmin and Tmax;
! returns Teq = 0. and irterr = 1 if no roots found
!
subroutine solve_temp_equil(nrho,Tmin0,Tmax0,gamma,Teq,irterr)
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
 bisection: do while (.not.converged .and. irterr == 0)

    func_max = equifunc(nrho,Tmax,gamma)
    func_min = equifunc(nrho,Tmin,gamma)

    Tmid = (Tmax+Tmin)/2.
    func_mid = equifunc(nrho,Tmid,gamma)

    if (func_mid*func_min < 0) then
       Tmax = Tmid
    elseif (func_mid*func_max < 0) then
       Tmin = Tmid
    endif

    if ((Tmax-Tmin) < tol) then
       converged = .true.
       Teq = (Tmax+Tmin)/2.  !- get midpoint
    endif

    niter = niter + 1
    if (niter > 1000) then   !- assume no roots found
       Teq = 0.
       irterr = 1
    endif
 enddo bisection

end subroutine solve_temp_equil

!
! Compute photoionization heating rate gamma in code units using nH from CMI
! Note that Osterbrock74 G(H) = nrho*gamma
!
subroutine heating_term(nH,rho,u,temp_star,gammaheat)
 use physcon, only:mass_proton_cgs,kboltz
 use units,   only:unit_density,unit_ergg
 use eos,     only:gamma,gmw
 use io,      only:fatal
 real, intent(in)  :: nH,rho,u,temp_star
 real, intent(out) :: gammaheat
 real :: rho_cgs,nrho_cgs,temp,alphaA,Ne,Np
 real :: heating_rate_cgs,gamma_cgs

 if (nH < 1.) then

    !- number density of H atoms
    rho_cgs = rho*unit_density
    if (rho_cgs < rhomin_cgs .or. rho_cgs > rhomax_cgs) then
       print*,'rho_cgs',rho_cgs
       call fatal('heating_cooling_cmi','rho_cgs exceeded range')
    endif
    nrho_cgs = rho_cgs*one_over_mH_cgs

    !- number density of free electrons and protons
    Ne = (1.-nH)*nrho_cgs
    Np = (1.-nH)*nrho_cgs

    !- compute approx alphaA
    temp = u/kboltz*(gmw*mass_proton_cgs*(gamma-1.))*unit_ergg
    alphaA = get_alphaA(temp)

    !- Photoionization heating G(H) as of Osterbrock74 eqn 3.2
    heating_rate_cgs = Ne*Np*alphaA*3./2.*kboltz*temp_star   ! [erg cm^-3 s-1]

    !- gamma
    gamma_cgs = heating_rate_cgs/nrho_cgs  ! [erg s-1]

 else
    gamma_cgs = 0.
 endif

 !- Include background heating
 gamma_cgs = gamma_cgs + gamma_background_cgs

 if (gamma_cgs < gammamin_cgs .or. gamma_cgs > gammamax_cgs) then
    print*,'nH; gamma_cgs',nH,gamma_cgs
    call fatal('heating_cooling_cmi','gamma_cgs exceeded range')
 endif

 !- convert to code units
 gammaheat = gamma_cgs*one_over_unit_gamma

end subroutine heating_term

!
! Routine to read cooling rate off temp directly from cooling curve
!
subroutine cooling_term(u,lambda)
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

end subroutine cooling_term

!
! Extract ueq(s) from pre-computed table with the given rho and gamma
! Selects the right ueq using the current u
!
subroutine get_ueq(rho,gammaheat,u,numroots,ueq_final)
 use physcon, only:kboltz,mass_proton_cgs
 use units,   only:unit_ergg
 use eos,     only:gmw,gamma
 use io,      only:fatal,warning
 real,    intent(in)  :: rho,gammaheat,u
 integer, intent(out) :: numroots
 real,    intent(out) :: ueq_final
 integer :: irho,igamma,jrho,jgamma,r
 real    :: rhotable(maxrho),gammatable(maxgamma)  !- unpacked from rho_gamma_ueq_table
 real    :: ueqs(3)                                !- holds the 3 possible ueq solutions
 real    :: x1,x2,y1,y2,q11,q12,q21,q22            !- for 2D interpolation

 rhotable   = rho_gamma_ueq_table(:,1,1)
 gammatable = rho_gamma_ueq_table(1,:,2)

 !- Closest index
 irho   = minloc(abs(rhotable(:)-rho),1)
 igamma = minloc(abs(gammatable(:)-gammaheat),1)
 if (irho == 1 .or. irho == maxrho) call warning('heating_cooling_cmi','rho range too small')
 if (igamma == 1 .or. igamma == maxgamma) call warning('heating_cooling_cmi','gamma range too small')

 numroots = int(rho_gamma_ueq_table(irho,igamma,3))
 if (numroots == 0) then
    call warning('heating_cooling_cmi','no equilibrium solution found - drift to Tmax')
    ueq_final = kboltz*Tmax / (gmw*mass_proton_cgs*(gamma-1)) / unit_ergg
    return
 endif

 ueqs = 0.
 if (get_ueq_by_interp) then
    ! Note: [i] closest index; [j] lower bound; [j+1] upper bound bracketing input var

    ! get bracketing indices [j] and [j+1] for rho
    jrho = 0
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
    jgamma = 0
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
!       print*,'u is greater than ueq2! Teq2 is ',ueqs(2)/kboltz*(gmw*mass_proton_cgs*(gamma-1))*unit_ergg
       if (numroots == 2) then
          ueq_final = kboltz*Tmax / (gmw*mass_proton_cgs*(gamma-1)) / unit_ergg
       elseif (numroots == 3) then
          ueq_final = ueqs(3)
       endif
    endif
 endif

end subroutine get_ueq

!
! With the ueq, compute new u after dt to give du in code units as of VS07
!
subroutine compute_du(dt,rho,u,ueq,gamma,lambda,tau,du)
 real, intent(in)  :: dt,rho,u,ueq,gamma,lambda
 real, intent(out) :: tau,du
 real  :: unew

 !- time required to radiate energy excess / gain energy
 tau  = abs((u-ueq)/(rho*one_over_mH2*lambda - one_over_mH*gamma))

 unew = ueq + (u-ueq)*exp(-dt/tau)
 du   = unew - u

end subroutine compute_du

!
! Balance heating and cooling term to give dudt in code units
!
subroutine compute_dudt(dt,rho,u,gamma,lambda,dudt)
 real, intent(in)    :: dt,rho,u
 real, intent(inout) :: gamma,lambda
 real, intent(out)   :: dudt

 dudt = one_over_mH*gamma - rho*one_over_mH2*lambda

end subroutine compute_dudt

! Note:
! According to KI02 and VS07, it should be nrho*gamma - nrho^2*lambda, where nrho = rho/mass_proton.
! But since u is in erg/g, we further divide nrho by rho to make the equations dimensionally correct,
! giving dudt = (1/mass_proton)*gamma - (rho/mass_proton**2)*lambda; [erg/g/s].


!
! Calculates timescale of R-type phase i.e. time taken to reach Stromgren radius
!
subroutine compute_Rtype_time(nphotosrc,xyz_photosrc,xyzh,t_recomb)
 use units,    only:udist,utime,unit_density
 use physcon,  only:mass_proton_cgs
 use part,     only:rhoh,massoftype,igas
 use kdtree,   only:getneigh
 use linklist, only:node,ifirstincell,listneigh
 use io,       only:fatal
 integer, intent(in)  :: nphotosrc
 real,    intent(in)  :: xyz_photosrc(3,nphotosrc)
 real,    intent(in)  :: xyzh(:,:)
 real,    intent(out) :: t_recomb
 integer, parameter   :: neighcachesize = 1E5
 integer :: isrc,nneigh,n,ip
 real    :: xyzcache(3,neighcachesize)
 real    :: pos_src(3),rcut,hmean,pmass,rhomean,rhomean_cgs,alphaB
 real    :: t_recomb_src
 real    :: rcut_cgs = 3.1E18    !- sample within 1pc around source

 rcut = rcut_cgs/udist
 alphaB = get_alphaB(1E4)  !- assumes constant temp of 1E4 K in HII region

 t_recomb = epsilon(t_recomb)

 each_source: do isrc = 1,nphotosrc
    pos_src = xyz_photosrc(3,isrc)

    !- Get mean density within a radius of rcut
    call getneigh(node,pos_src,0.,rcut,3,listneigh,nneigh,xyzh,xyzcache,neighcachesize,ifirstincell,.false.)
    if (nneigh < 100) call fatal('heating_cooling_cmi','not sampling enough particles for estimating t_recomb')

    hmean = 0.
    over_neigh: do n = 1,nneigh
       ip = listneigh(n)
       hmean = hmean + xyzh(4,ip)
    enddo over_neigh
    hmean = hmean/nneigh

    pmass = massoftype(igas)
    rhomean = rhoh(hmean,pmass)
    rhomean_cgs = rhomean *unit_density

    t_recomb_src = mass_proton_cgs/(alphaB*rhomean_cgs) /utime
    t_recomb = max(t_recomb,t_recomb_src)
 enddo each_source

end subroutine compute_Rtype_time

!
! Case-A recombination coeff alphaA as func of temp obtained from fitting table 2.1
! of Osterbrock74 [cm^3 s^-1]
!
real function get_alphaA(temp)
 real, intent(in) :: temp
 real :: a = 6.113723867E-11
 real :: b = -1.85763022E-13

 if (use_const_alpha) then
    get_alphaA = 2.7E-13
 else
    get_alphaA = a* temp**(-1./2.) + b  !- P.16: recomb coef ~ T^(-1/2)
 endif
 get_alphaA = max(get_alphaA,1.E-13)

end function get_alphaA

!
! Case-B recombination coeff alphaB as func of temp obtained from fitting table 2.1
! of Osterbrock74 [cm^3 s^-1]
!
real function get_alphaB(temp)
 real, intent(in) :: temp
 real :: a = 4.417139142E-11
 real :: b = -1.73910208E-13

 if (use_const_alpha) then
    get_alphaB = 2.7E-13
 else
    get_alphaB = a* temp**(-1./2.) + b
 endif
 get_alphaB = max(get_alphaB,1.E-13)

end function get_alphaB

!
! Analytic function which mimics the cooling curve of DeRijcke et al. (2013)
! Note: This function is a smoothed-out best fit of the tabulated curve which removes all
!       the 'wriggles' to limit the number of equil solutions to a max of 3.
!
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

end module heatcool_cmi
