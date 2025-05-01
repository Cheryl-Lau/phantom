!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module moddump
!
! To impose a turbulence velocity field onto the particles 
!
! :References: 
!
! :Owner: Cheryl Lau
!
! :Runtime parameters: None
!
! :Dependencies: part
!

 implicit none

 real    :: rms_mach = 1.d0 

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 use part,         only:igas,xyzmh_ptmass,vxyz_ptmass,nptmass,ihsoft,ihacc 
 use units,        only:unit_velocity,unit_ergg
 use velfield,     only:set_velfield_from_cubes
 use datafiles,    only:find_phantom_datafile
 use eos,          only:gmw,gamma 
 use physcon,      only:mass_proton_cgs,kboltz
 use io,           only:fatal
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i,ierr
 real    :: xmini,xmaxi,ymini,ymaxi,zmini,zmaxi,boxlength,turbboxsize
 real    :: rmsmach,v2i,turbfac,cs0,cs0_cgs,mean_v,sigma_v,u_mean,temp 
 character(len=20), parameter :: filevx = 'cube_v1.dat'
 character(len=20), parameter :: filevy = 'cube_v2.dat'
 character(len=20), parameter :: filevz = 'cube_v3.dat'
 character(len=120) :: filex,filey,filez

 !
 ! Find box size 
 ! 
 xmini = huge(xmini); xmaxi = tiny(xmaxi)
 ymini = huge(ymini); ymaxi = tiny(ymaxi)
 zmini = huge(zmini); zmaxi = tiny(zmaxi)
 do i = 1,npart
    xmini = min(xmini,xyzh(1,i)); xmaxi = max(xmaxi,xyzh(1,i))
    ymini = min(ymini,xyzh(2,i)); ymaxi = max(ymaxi,xyzh(2,i))
    zmini = min(zmini,xyzh(3,i)); zmaxi = max(zmaxi,xyzh(3,i))
 enddo 
 boxlength = max(xmaxi-xmini,ymaxi-ymini,zmaxi-zmini)

 !
 ! Import cubes 
 !
 filex  = find_phantom_datafile(filevx,'velfield_sphng_small')
 filey  = find_phantom_datafile(filevy,'velfield_sphng_small')
 filez  = find_phantom_datafile(filevz,'velfield_sphng_small')
 !- [To convert endian for different vfield files: setenv GFORTRAN_CONVERT_UNIT big/small_endian]

 !
 ! Calculate sound speed 
 !
 u_mean = 0.d0
 do i = 1,npart
    u_mean = u_mean + vxyzu(4,i)
 enddo 
 u_mean = u_mean/npart
 temp = u_mean/kboltz*(gmw*mass_proton_cgs*(gamma-1.))*unit_ergg
 cs0_cgs = sqrt(temp * (gamma*kboltz) / (gmw*mass_proton_cgs))
 cs0 = cs0_cgs/unit_velocity 

 !
 ! Set turbulence 
 !
 turbboxsize = boxlength
 call set_velfield_from_cubes(xyzh(:,1:npart),vxyzu(:,1:npart),npart, &
                           filex,filey,filez,1.,turbboxsize,.false.,ierr)
 if (ierr /= 0) call fatal('setup','error setting up velocity field on clouds')

 rmsmach = 0.0
 do i = 1,npart
    v2i     = mag(vxyzu(1:3,i))**2
    rmsmach = rmsmach + v2i/cs0**2
 enddo
 rmsmach = sqrt(rmsmach/npart)
 if (rmsmach > 0.) then
    turbfac = rms_mach/rmsmach ! normalise the energy to the desired mach number
 else
    turbfac = 0.
 endif

 !- Set velocities and get mean 
 mean_v = 0.
 do i = 1,npart
    vxyzu(1:3,i) = turbfac*vxyzu(1:3,i)
    mean_v = mean_v + mag(vxyzu(1:3,i))
 enddo
 mean_v = mean_v/real(npart) 

 !- Get initial dispersion 
 sigma_v = 0. 
 do i = 1,npart
    sigma_v = sigma_v + (mag(vxyzu(1:3,i)) - mean_v)**2 
 enddo 
 sigma_v = sqrt(sigma_v/real(npart))


 open(2040,file='turbbox_velocities.dat',status='replace')
 write(2040,*) 'Mean velocity [cm s^-1]: ',mean_v*unit_velocity 
 write(2040,*) 'Dispersion [cm s^-1]: ',sigma_v*unit_velocity 
 close(2040)

 return
end subroutine modify_dump


real function mag(vec)
 real,   intent(in) :: vec(3)

 mag = sqrt(dot_product(vec,vec))

end function mag



end module moddump
