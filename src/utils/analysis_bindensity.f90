!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine which bins the particles by density
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies:
!

 implicit none
 character(len=20), parameter, public :: analysistype = 'bin_density'

 public :: do_analysis

 private
 integer, parameter :: num_rhobin = 100
 integer :: binned_rho(num_rhobin)
 real    :: rhobin_min_cgs = 5E-23
 real    :: rhobin_max_cgs = 4E-18
 real    :: logrhobin_min,logrhobin_max

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part,  only:hfact
 use units, only:unit_density
 use io,    only:fatal
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: i,irhobin
 real    :: rho,drho,rhobin(num_rhobin)
 character(len=70) :: filename

 logrhobin_min = log10(rhobin_min_cgs/unit_density)
 logrhobin_max = log10(rhobin_max_cgs/unit_density)

 !
 ! Separation of rho in log-scale
 !
 drho = (logrhobin_max-logrhobin_min)/num_rhobin
 !
 ! Init binned array
 !
 do irhobin = 1,num_rhobin
    binned_rho(irhobin) = 0.
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
    ! Bin rho
    !
    irhobin = (log10(rho)-logrhobin_min+drho)/drho
    if (irhobin < 0) call fatal('analysis_bindensity','require smaller logrhobin_min')
    if (irhobin > num_rhobin) call fatal('analysis_bindensity','require larger logrhobin_max')
    binned_rho(irhobin) = binned_rho(irhobin) + 1
 enddo
 !
 ! Convert to physical units and store results
 !
 do irhobin = 1,num_rhobin
    rhobin(irhobin) = (10**(logrhobin_min + (irhobin-1)*drho))*unit_density
 enddo

 filename = 'density_binned_'//TRIM(dumpfile)//'.dat'
 open(unit=2021,file=filename)
 do irhobin = 1,num_rhobin
    write(2021,*) rhobin(irhobin), binned_rho(irhobin)
 enddo
 close(2021)

end subroutine do_analysis
!--------------------------------------------------------------------------
end module analysis
