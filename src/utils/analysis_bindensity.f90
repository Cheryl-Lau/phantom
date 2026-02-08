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
! :Owner: Cheryl Lau
!
! :Runtime parameters: None
!
! :Dependencies:
!

 implicit none
 character(len=20), parameter, public :: analysistype = 'bin_density'

 public :: do_analysis

 private
 integer, parameter :: num_rhobin = 500
 integer :: binned_rho(num_rhobin)
 real    :: rhobin_min_cgs = 1E-28
 real    :: rhobin_max_cgs = 1E-8
 real    :: logrhobin_min,logrhobin_max
 real    :: centre(3) = (/ 0.,0.,0. /)
 real    :: radius = 70. 
 logical :: box_only = .true. 

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part,  only:hfact
 use units, only:unit_density
 use io,    only:fatal
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: i,irhobin,npart_box
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
 npart_box = 0
 do i = 1,npart
    !
    ! Filter 
    !
    if (box_only) then 
       if (xyzh(1,i) < centre(1)-radius .or. xyzh(1,i) > centre(1)+radius .or. &
           xyzh(2,i) < centre(2)-radius .or. xyzh(2,i) > centre(2)+radius .or. &
           xyzh(3,i) < centre(3)-radius .or. xyzh(3,i) > centre(3)+radius) cycle 
    endif 
    !
    ! Convert smoothing length h to density rho
    !
    rho = particlemass*(hfact/abs(xyzh(4,i)))**3
    !
    ! Bin rho
    !
    irhobin = (log10(rho)-logrhobin_min+drho)/drho
    if (irhobin < 0) then 
       print*,'rho_cgs',rho*unit_density
       call fatal('analysis_bindensity','require smaller logrhobin_min')
    elseif (irhobin > num_rhobin) then 
       print*,'rho_cgs',rho*unit_density
       call fatal('analysis_bindensity','require larger logrhobin_max')
    endif 
    binned_rho(irhobin) = binned_rho(irhobin) + 1

    npart_box = npart_box + 1 
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

 print*,'Number of particles in box',npart_box

end subroutine do_analysis
!--------------------------------------------------------------------------
end module analysis
