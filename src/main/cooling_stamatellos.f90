!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module cooling_stamatellos
!
! Cooling method of Stamatellos et al. 2007
!
! :References: Stamatellos et al. 2007
!
! :Owner: Alison Young
!
! :Runtime parameters: None
!
! :Dependencies: eos_stamatellos, physcon,io,units,part
!
 
 implicit none
 real, public :: Lstar ! in units of L_sun
 integer :: isink_star ! index of sink to use as illuminating star
 integer :: od_method
 public :: cooling_S07,write_options_cooling_stamatellos,read_options_cooling_stamatellos
 public :: init_star
 
 contains

subroutine init_star()
  use part,    only:nptmass,xyzmh_ptmass
  integer :: i,imin
  real :: rsink2,rsink2min
  
  rsink2min = 0d0
  if (nptmass == 0 .or. Lstar == 0.0) then
     isink_star = 0 ! no stellar heating
     print *, "No stellar heating."
  elseif (nptmass == 1) then
     isink_star = 1
  else
     do i=1,nptmass
        rsink2 = xyzmh_ptmass(1,i)**2 + xyzmh_ptmass(2,i)**2 + xyzmh_ptmass(3,i)**2
        if (i==1 .or. (rsink2 < rsink2min) ) then
           rsink2min = rsink2
           imin = i
        endif
     enddo
     isink_star = imin
  endif
  if (isink_star > 0)  print *, "Using sink no. ", isink_star, "as illuminating star."
end subroutine init_star
  
!
! Do cooling calculation
!
   subroutine cooling_S07(rhoi,ui,dudti_cool,xi,yi,zi,Tfloor,dudti_sph,dt,i)
     use io,       only:warning
     use physcon,  only:steboltz,pi,solarl
     use units,    only:umass,udist,unit_density,unit_ergg,utime
     use eos_stamatellos, only:getopac_opdep,getintenerg_opdep,gradP_cool,Gpot_cool, &
          iunitst
     use part,       only:eos_vars,igasP,xyzmh_ptmass
     real,intent(in) :: rhoi,ui,dudti_sph,xi,yi,zi,Tfloor,dt
     integer,intent(in) :: i
     real,intent(out) :: dudti_cool
     real            :: coldensi,kappaBari,kappaParti,ri2
     real            :: gammai,gmwi,Tmini4,Ti,dudt_rad,Teqi
     real            :: tcool,ueqi,umini,tthermi,poti,presi
     
     poti = Gpot_cool(i)
     presi = eos_vars(igasP,i)
!    Tfloor is from input parameters and is background heating
!    Stellar heating
     if (isink_star > 0 .and. Lstar > 0) then
        ri2 = (xi-xyzmh_ptmass(1,isink_star))**2d0 &
             + (yi-xyzmh_ptmass(2,isink_star))**2d0 &
             + (zi-xyzmh_ptmass(3,isink_star))**2d0
        ri2 = ri2 *udist*udist
! Tfloor + stellar heating
        Tmini4 = Tfloor**4d0 + (Lstar*solarl/(16d0*pi*steboltz*ri2))
     else
        Tmini4 = Tfloor**4d0
     endif

! get opacities & Ti for ui
     call getopac_opdep(ui*unit_ergg,rhoi*unit_density,kappaBari,kappaParti,&
           Ti,gmwi,gammai)
     select case (od_method)
     case (1)
        coldensi = sqrt(abs(poti*rhoi)/4.d0/pi)
        coldensi = 0.368d0*coldensi ! n=2 in polytrope formalism Forgan+ 2009
        coldensi = coldensi*umass/udist/udist ! physical units
     case(2)
! testing Lombardi+ method of estimating the mean column density
        coldensi = 1.014d0 * presi / abs(-gradP_cool(i))! 1.014d0 * P/(-gradP/rho) Lombardi+ 2015
        coldensi = coldensi *umass/udist/udist ! physical units
     end select
!     write(iunitst,'(5E12.5)') coldensi,presi,gradP_cool(i)
     
     tcool = (coldensi**2d0)*kappaBari +(1.d0/kappaParti) ! physical units
     dudt_rad = 4.d0*steboltz*(Tmini4 - Ti**4.d0)/tcool/unit_ergg*utime! code units
! calculate Teqi
     Teqi = dudti_sph*(coldensi**2.d0*kappaBari + (1.d0/kappaParti))*unit_ergg/utime
     Teqi = Teqi/4.d0/steboltz
     Teqi = Teqi + Tmini4
     if (Teqi < Tmini4) then
        Teqi = Tmini4**(1.0/4.0)
     else
        Teqi = Teqi**(1.0/4.0)
     endif
     call getintenerg_opdep(Teqi,rhoi*unit_density,ueqi)
     ueqi = ueqi/unit_ergg
     call getintenerg_opdep(Tmini4**(1.0/4.0),rhoi*unit_density,umini)
     umini = umini/unit_ergg
! calculate thermalization timescale
     if ((dudti_sph + dudt_rad) == 0.d0) then
        tthermi = 0d0
        write(iunitst,'(A)') "ttherm=0"
     else
        tthermi = abs((ueqi - ui)/(dudti_sph + dudt_rad))
     endif
     
! internal energy update -> put in form where it'll work as dudtcool
     if (tthermi == 0d0) then
        dudti_cool = 0.d0 ! condition if denominator above is zero
     else
        dudti_cool = (ui*exp(-dt/tthermi) + ueqi*(1.d0-exp(-dt/tthermi)) -ui)/dt !code units
     endif
     
     if (isnan(dudti_cool)) then
        print *, "kappaBari=",kappaBari, "kappaParti=",kappaParti
        print *, "poti=",poti, "rhoi=",rhoi, "Ti=", Ti
        print *, "tcool=",tcool,"coldensi=",coldensi,"dudti_sph",dudti_sph
        print *, "Teqi=",Teqi, "dt=",dt,"tthermi=", tthermi,"ueqi=", ueqi
        call warning("In Stamatellos cooling","dudticool=NaN. ui",val=ui)
        stop
     else if (dudti_cool < 0.d0 .and. abs(dudti_cool) > ui/dt) then
        dudti_cool = (umini - ui)/dt
     endif
     
   end subroutine cooling_S07


 subroutine write_options_cooling_stamatellos(iunit)
 use infile_utils, only:write_inopt
 use eos_stamatellos, only: eos_file
 integer, intent(in) :: iunit

 !N.B. Tfloor handled in cooling.F90
 call write_inopt(eos_file,'EOS_file','File containing tabulated EOS values',iunit)
 call write_inopt(od_method,'OD method','Method for estimating optical depth: (1) potential (2) pressure',iunit)
 call write_inopt(Lstar,'Lstar','Luminosity of host star for calculating Tmin (Lsun)',iunit)

end subroutine write_options_cooling_stamatellos

 subroutine read_options_cooling_stamatellos(name,valstring,imatch,igotallstam,ierr)
 use io, only:warning,fatal
use eos_stamatellos, only: eos_file
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotallstam
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0


 imatch  = .true.
 igotallstam = .false. ! cooling options are compulsory
 select case(trim(name))
 case('Lstar')
    read(valstring,*,iostat=ierr) Lstar
    ngot = ngot + 1
 case('OD method')
    read(valstring,*,iostat=ierr) od_method
    if (od_method < 1 .or. od_method > 2) then
       call fatal('cooling options','od_method must be 1 or 2',var='od_method',ival=od_method)
    endif
    ngot = ngot + 1
 case('EOS_file')
	read(valstring,*,iostat=ierr) eos_file
	ngot = ngot + 1
 case default
    imatch = .false.
 end select
 if (od_method .ne. 1 .and. od_method .ne. 2) then
    call warning('cooling_stamatellos','optical depth method unknown')
 endif

 if (ngot >= 3) igotallstam = .true.

end subroutine read_options_cooling_stamatellos

 end module cooling_stamatellos
