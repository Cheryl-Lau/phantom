!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine which measures the velocity dispersion of all/high-density particles
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
 character(len=20), parameter, public :: analysistype = 'ionized_vel_dispersion'

 public :: do_analysis

 private

 integer, parameter :: nrad  = 10
 real    :: rad_thresh(nrad) = (/ 6e-2,8e-2,1e-1,2e-1,4e-1,6e-1,8e-1,1e0,2e0,4e0 /)

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use linklist, only:node,ifirstincell,listneigh,set_linklist
 use kdtree,   only:getneigh
 use dim,      only:maxneigh
 use units,    only:udist,utime,umass,unit_velocity,unit_density,unit_pressure,unit_ergg
 use io,       only:fatal,warning 
 use part,     only:hfact,rhoh,massoftype,igas
 use eos,      only:gamma,gmw
 use physcon,  only:kboltz,mass_proton_cgs,gg,pi
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer, parameter :: neighcachesize = 1E5
 integer :: ip,ineigh,nneigh,ixyzcachesize,n,irad,ip_neigh
 real    :: xyzcache(3,neighcachesize)
 real    :: vel,temp 
 real    :: pmass,maxrad,dist2,rho,rad2_limit,mean_v,sigma_v,mean_rho
 real    :: sigma_v_allrad(nrad),virial_term_allrad(nrad),rho_avg_allrad(nrad)
 real,   allocatable :: dumxyzh(:,:)
 character(len=70) :: filename

 !- Particle mass
 pmass = massoftype(igas)


 !- Build tree 
 allocate(dumxyzh(4,npart))
 dumxyzh = xyzh
 call set_linklist(npart,npart,dumxyzh,vxyzu)


 open(unit=2028,file='vel_temp_'//TRIM(dumpfile)//'.dat',status='replace')
 open(unit=2029,file='velocity_dispersion_ionized_'//TRIM(dumpfile)//'.dat',status='replace')
 open(unit=2030,file='virial_term_ionized_'//TRIM(dumpfile)//'.dat',status='replace')
 open(unit=2031,file='density_sizescale_ionized_'//TRIM(dumpfile)//'.dat',status='replace')


 do ip = 1,npart
    vel = sqrt(mag2(vxyzu(1:3,ip)))
    temp = vxyzu(4,ip)/kboltz*(gmw*mass_proton_cgs*(gamma-1.))*unit_ergg
    write(2028,'(2f30.10)') vel*unit_velocity, temp

    ionized: if (temp > 9e3) then 
       ! Get all neigh 
       call getneigh(node,xyzh(1:3,ip),0.,maxrad,3,listneigh,nneigh,xyzh,xyzcache,neighcachesize,ifirstincell,.false.)

       over_rad: do irad = 1,nrad
         rad2_limit = (rad_thresh(irad))**2 

         !- Find mean velocity 
         mean_v = 0. 
         mean_rho = 0.
         n = 0
         over_neigh_mean: do ineigh = 1,nneigh
            ip_neigh = listneigh(ineigh)
            dist2 = mag2(xyzh(1:3,ip)-xyzh(1:3,ip_neigh))
            if (dist2 < rad2_limit) then 
               n = n + 1 
               mean_v = mean_v + sqrt(mag2(vxyzu(1:3,ip_neigh)))
               mean_rho = mean_rho + pmass*(hfact/abs(xyzh(4,ip_neigh)))**3
            endif 
         enddo over_neigh_mean
         if (n == 0) cycle over_rad
         mean_v = mean_v/real(n) 
         mean_rho = mean_rho/real(n) 

         !- Find dispersion (std)
         sigma_v = 0.
         over_neigh_sigma: do ineigh = 1,nneigh
            ip_neigh = listneigh(ineigh)
            dist2 = mag2(xyzh(1:3,ip)-xyzh(1:3,ip_neigh)) 
            if (dist2 < rad2_limit) then 
               sigma_v = sigma_v + (sqrt(mag2(vxyzu(1:3,ip_neigh))) - mean_v)**2 
            endif 
         enddo over_neigh_sigma
         sigma_v = sqrt(sigma_v/real(n)) 

         sigma_v_allrad(irad) = sigma_v*unit_velocity
         if (sigma_v > 0) then   
            virial_term_allrad(irad) = 2.d0*gg*(n*pmass*umass)/((sigma_v*unit_velocity)**2*rad_thresh(irad)*udist)
         else 
            virial_term_allrad(irad) = -1 ! to be removed during analysis 
         endif 

         !- Mean density 
         rho_avg_allrad(irad) = (n*pmass*umass)/(4./3.*pi*(rad_thresh(irad)*udist)**3)

       enddo over_rad

       write(2029,'(1i15,10f30.10)') ip, sigma_v_allrad(1:nrad)
       write(2030,'(1i15,10f30.10)') ip, virial_term_allrad(1:nrad)
       write(2031,'(1i15,10es30.10)') ip, rho_avg_allrad(1:nrad)

    endif ionized
 enddo

 close(2028)
 close(2029)
 close(2030) 
 close(2031) 
 deallocate(dumxyzh)


end subroutine do_analysis


real function mag2(vec)
 real,   intent(in) :: vec(3)

 mag2 = dot_product(vec,vec)

end function mag2

!--------------------------------------------------------------------------
end module analysis
