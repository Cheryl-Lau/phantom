!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine which measures the cloud turbulence timescale 
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
 character(len=20), parameter, public :: analysistype = 't_turb'

 public :: do_analysis

 private

 real    :: rholimit_cgs     = 1e-20
 logical :: only_highdenpart = .true. 

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use units,    only:udist,utime,umass,unit_velocity,unit_density,unit_pressure,unit_ergg
 use io,       only:fatal,warning 
 use part,     only:hfact,rhoh,massoftype,igas
 use physcon,  only:pc 
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer :: ip,n
 real    :: rho,pmass,mean_v,sigma_v,t_turb_cgs,t_turb_Myr,rholimit
 real    :: xmin,xmax,ymin,ymax,zmin,zmax,dx,dy,dz,rmax
 character(len=70) :: filename

 !- Particle mass
 pmass = massoftype(igas)

 rholimit = rholimit_cgs/unit_density 

 !
 ! Find mean velocity
 !
 n = 0
 mean_v = 0. 
 !$omp parallel do default(none) shared(npart,xyzh,vxyzu,pmass,hfact,only_highdenpart,rholimit) &
 !$omp private(ip,rho) &
 !$omp reduction(+:mean_v,n) &
 !$omp schedule(runtime)
 get_mean: do ip = 1,npart 
    rho = pmass*(hfact/abs(xyzh(4,ip)))**3 
    if (only_highdenpart .and. rho < rholimit) cycle get_mean 
    n = n + 1 
    mean_v = mean_v + sqrt(mag2(vxyzu(1:3,ip)))
 enddo get_mean
 !$omp end parallel do
 mean_v = mean_v/real(n) 

 !
 ! Find velocity dispersion
 ! 
 sigma_v = 0. 
 !$omp parallel do default(none) shared(npart,xyzh,vxyzu,pmass,hfact,mean_v,only_highdenpart,rholimit) &
 !$omp private(ip,rho) &
 !$omp reduction(+:sigma_v) &
 !$omp schedule(runtime)
 get_sigma: do ip = 1,npart 
    rho = pmass*(hfact/abs(xyzh(4,ip)))**3 
    if (only_highdenpart .and. rho < rholimit) cycle get_sigma 
    sigma_v = sigma_v + (sqrt(mag2(vxyzu(1:3,ip))) - mean_v)**2
 enddo get_sigma 
 !$omp end parallel do
 sigma_v = sqrt(sigma_v/real(n)) 

 !
 ! Find extent of cloud 
 !
 xmax = -huge(xmax)
 ymax = -huge(ymax)
 zmax = -huge(zmax)
 xmin =  huge(xmin)
 ymin =  huge(ymin)
 zmin =  huge(zmin)
 !$omp parallel do default(none) shared(npart,xyzh,pmass,hfact,only_highdenpart,rholimit) &
 !$omp private(ip,rho) &
 !$omp reduction(min:xmin,ymin,zmin) &
 !$omp reduction(max:xmax,ymax,zmax) &
 !$omp schedule(runtime)
 get_sizescale: do ip = 1,npart
    rho = pmass*(hfact/abs(xyzh(4,ip)))**3 
    if (only_highdenpart .and. rho < rholimit) cycle get_sizescale
    xmin = min(xmin,xyzh(1,ip))
    ymin = min(ymin,xyzh(2,ip))
    zmin = min(zmin,xyzh(3,ip))
    xmax = max(xmax,xyzh(1,ip))
    ymax = max(ymax,xyzh(2,ip))
    zmax = max(zmax,xyzh(3,ip))
 enddo get_sizescale
 !$omp end parallel do

 dx = abs(xmax - xmin)
 dy = abs(ymax - ymin)
 dz = abs(zmax - zmin)
 rmax = (dx+dy+dz)/3.d0

 !
 ! Turbulent crossing timescale 
 !
 t_turb_cgs = (rmax*udist)/(sigma_v*unit_velocity) 
 t_turb_Myr = t_turb_cgs/(1d6*365*24*60*60)

 print*,'sizescale [pc] = ',rmax*udist/pc 
 print*,'sigma_v [cm/s] = ',sigma_v*unit_velocity 
 print*,'t_turb [Myr] = ',t_turb_Myr


end subroutine do_analysis


real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2

!--------------------------------------------------------------------------
end module analysis
