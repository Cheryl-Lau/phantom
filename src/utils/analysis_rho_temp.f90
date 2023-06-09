!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine for plotting T [K] vs log rho [g cm^-3]
!
! :References: None
!
! :Owner: Cheryl Lau
!
! :Runtime parameters: None
!
! :Dependencies:
!

 implicit none
 character(len=20), parameter, public :: analysistype = 'rho_temp'

 public :: do_analysis

 private
 integer, parameter :: num_rhobin = 100
 integer :: binned_temp(num_rhobin)
 real    :: rhobin_min_cgs = 1E-25
 real    :: rhobin_max_cgs = 1E-18
 real    :: logrhobin_min,logrhobin_max

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part,    only:hfact
 use units,   only:unit_density,unit_ergg
 use physcon, only:kboltz,mass_proton_cgs
 use eos,     only:gamma,gmw
 use io,      only:fatal
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: i,irhobin,nbinnedpart(num_rhobin)
 real    :: rho,drho,rhobin(num_rhobin),temp
 character(len=70) :: filename

 filename = 'rho_temp_parts_'//TRIM(dumpfile)//'.dat'
 open(2032,file=filename)

 !
 ! rho in log-scale
 !
 logrhobin_min = log10(rhobin_min_cgs/unit_density)
 logrhobin_max = log10(rhobin_max_cgs/unit_density)
 drho = (logrhobin_max-logrhobin_min)/num_rhobin
 !
 ! Init binned array
 !
 do irhobin = 1,num_rhobin
    binned_temp(irhobin) = 0.
    nbinnedpart(irhobin) = 0
 enddo
 !
 ! Bin all particles by density
 !
 do i = 1,npart
    !
    ! Convert smoothing length h to density rho
    !
    rho = particlemass*(hfact/abs(xyzh(4,i)))**3
    !
    ! Calculate particle temperature
    !
    temp = vxyzu(4,i)/kboltz*(gmw*mass_proton_cgs*(gamma-1.))*unit_ergg

    write(2032,*) rho*unit_density, temp
    !
    ! Bin temp by rho
    !
    irhobin = (log10(rho)-logrhobin_min+drho)/drho
    if (irhobin < 0) call fatal('analysis_bindensity','require smaller logrhobin_min')
    if (irhobin > num_rhobin) call fatal('analysis_bindensity','require larger logrhobin_max')
    binned_temp(irhobin) = binned_temp(irhobin) + temp
    nbinnedpart(irhobin) = nbinnedpart(irhobin) + 1
 enddo
 !
 ! Convert to physical units / Take mean and store results
 !
 do irhobin = 1,num_rhobin
    rhobin(irhobin) = (10**(logrhobin_min + (irhobin-1)*drho))*unit_density
    if (nbinnedpart(irhobin) > 0) then
       binned_temp(irhobin) = binned_temp(irhobin)/nbinnedpart(irhobin)
    else
       binned_temp(irhobin) = 0.
    endif
 enddo

 filename = 'rho_temp_binned_'//TRIM(dumpfile)//'.dat'
 open(unit=2021,file=filename)
 do irhobin = 1,num_rhobin
    write(2021,*) rhobin(irhobin), binned_temp(irhobin)
 enddo

 close(2021)
 close(2032)

end subroutine do_analysis
!--------------------------------------------------------------------------
end module analysis
