!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module ghostpart
!
! This module contains all the subroutines necessary for adding ghost particles -
! reference points following a grid-like structure that tracks flows in SPH particles;
! used for computing densities and other SPH fluid properties at regions which are not
! sufficently resolved by SPH particles. Suitable for MCRT grid construction.
!
!- Ghost particles will gather at regions where an astrophysical flow is detected.
!- Motion of ghost particles can be similar to that of longitudinal waves, hence can be
!- controlled using curves that define their separations. We superpose Gaussian functions
!- onto the curve at the positions of the detected flows to achieve this.
!
! :References: none
!
! :Owner: Cheryl Lau
!
! :Runtime parameters:
!   - sphpart_to_ghpart_ratio  : *Approx fraction of ghost particles to SPH particles*
!   - velmag_diff_tol          : *Tolerance in mag(velocity) difference between ghost particles
!                                 and sph particles beyond which we will re-setup the grid*
!   - velvec_diff_tol          : *Tolerence in angle differnce of velocities between ghost particles
!                                 and sph particles beyond which we will re-setup the grid*
!   - sphpart_vel_ghpart_sep   : *Free param for controlling the proportionality between flow
!                                 velocity and closeness of ghost particles*
!
! :Dependencies: infile_utils, allocutils, io,
!
!
 implicit none

 public :: init_ghostpart,set_ghostpart,allocate_ghostpart,deallocate_ghostpart
 public :: read_options_ghostpart,write_options_ghostpart

 integer, public :: sphpart_to_ghpart_ratio = 50
 real   , public :: boxsize_ghpart = (/-10.,-10.,-10.,10.,10.,10./) ! xyzmin xyzmax
 logical, public :: auto_boxsize   = .true.

 real,    public :: velmag_diff_tol = 1E-1      ! code units
 real,    public :: velvec_diff_tol = 0.78      ! radians
 real,    public :: sphpart_vel_ghpart_sep = 1E-4

 integer, public :: nghpart
 real   , public :: xyz_ghpart(:,:),vxyz_ghpart(:,:),rcut_ghpart(:),rho_ghpart(:)

 private

 type gausswave
  real :: centre
  real :: amplitude
  real :: sigma
  real :: velocity
 end type gausswave

 real,    allocatable :: edgepoint(:,:),veli_edgepoint(:,:)      !- pos and vel of edge lines
 real,    allocatable :: deltax_rel_curve(:,:,:)                 !- curve defining separation of points
 real,    allocatable :: width_stripe(:,:),veli_stripe(:,:)      !- wdith and max sphpart vel in between points
 logical, allocatable :: flowflag_stripe(:,:),isinglestripe(:,:) !- detecting flows in between points
 type(gausswave), allocatable :: curvewave(:,:)                  !- stores waves added onto deltax-curve

 real,    allocatable :: velmag_parts_stripe(:,:,:)     !- mag(vel) of individual sphpart in each stipe
 real,    allocatable :: veli_parts_stripe(:,:,:)       !- vel_idim of individual sphpart in each stipe

 integer :: npoint(3),nwave(3)
 real    :: edgemin(3),edgemax(3),dx_curve(3),hmean
 logical :: first_call


contains

!----------------------------------------------------------------
!+
! Initialize constant params and storage arrays for setting up ghost particles
!+
!----------------------------------------------------------------
subroutine init_ghostpart(npart,xyzh)
 use io,  only:warning
 integer, intent(in) :: npart
 real,    intent(in) :: xyzh(:,:)

 !
 ! Set box size containing all ghost particles (fixed throughout sim)
 !
 if (auto_boxsize) then
    xmax = -huge(xmax)
    ymax = -huge(ymax)
    zmax = -huge(zmax)
    xmin =  huge(xmin)
    ymin =  huge(ymin)
    zmin =  huge(zmin)
    !$omp parallel do default(none) shared(npart,xyzh) private(ip) &
    !$omp reduction(min:xmin,ymin,zmin) &
    !$omp reduction(max:xmax,ymax,zmax)
    do ip = 1,npart
      xmin = min(xmin,xyzh(1,ip))
      ymin = min(ymin,xyzh(2,ip))
      zmin = min(zmin,xyzh(3,ip))
      xmax = max(xmax,xyzh(1,ip))
      ymax = max(ymax,xyzh(2,ip))
      zmax = max(zmax,xyzh(3,ip))
    enddo
    !$omp end parallel do
    margin = 1E-5*(xmax-xmin)
    edgemin(:) = (/xmin-margin,ymin-margin,zmin-margin/)
    edgemax(:) = (/xmax+margin,ymax+margin,zmax+margin/)
 else
    edgemin(:) = boxsize_ghpart(1:3)
    edgemax(:) = boxsize_ghpart(4:6)
 endif
 !
 ! Set number of ghost particles on each edge
 !
 volbox    = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)
 lsphpart  = (npart/volbox)**(1./3.)   !- avg number of sph particles per unit length
 lghpart   = lsphpart/sphpart_to_ghpart_ratio
 do idim = 1,3
    npoint(idim) = nint((edgemax(idim)-edgemin(idim))*lghpart)
 enddo
 nghpart   = npoint(1)*npoint(2)*npoint(3)

 npartdiff = abs(nghpart - npart/sphpart_to_ghpart_ratio)
 if (npartdiff > epsilon(npartdiff)) print*,'Warning - number of ghost particles is off &
    the required sphpart:ghpart ratio by ',nint(npartdiff/(npart/sphpart_to_ghpart_ratio)*100.),'%'
 print*,'Placing ',nghpart,' ghost particles'
 !
 ! Allocate cubical storages for all dimensions
 !
 maxpoint = max(npoint(1),npoint(2),npoint(3))
 maxghpart = maxpoint**3
 call allocate_ghostpart(maxghpart,maxpoint,npart)
 !
 ! Init deltax-rel curve
 !
 do idim = 1,3
    dx_curve(idim) = (edgemax(idim)-edgemin(idim))/npoint(idim)  !- curve sampling interval
    call reset_deltaxcurve(idim,edgemin(idim),edgemax(idim),npoint(idim),dx_curve(idim))
 enddo
 first_call = .true.
 !
 ! Mean resolution length of sph particles - for setting lower limit in deltax
 !
 hmean = 0.
 !$omp parallel do default(none) shared(npart,xyzh) private(ip) &
 !$omp reduction(+:hmean)
 do ip = 1,npart
    hmean = hmean + xyzh(4,ip)
 enddo
 !$omp end parallel do
 hmean = hmean/npart

end subroutine init_ghostpart


subroutine allocate_ghostpart(maxghpart,maxpoint,npart)
 use allocutils, only:allocate_array
 integer, intent(in) :: maxghpart,maxpoint,npart

 call allocate_array('xyz_ghpart',xyz_ghpart,3,maxghpart)
 call allocate_array('vxyz_ghpart',vxyz_ghpart,3,maxghpart)
 call allocate_array('rho_ghpart',rho_ghpart,maxghpart)
 call allocate_array('rcut_ghpart',rcut_ghpart,maxghpart)
 call allocate_array('edgepoint',edgepoint,maxpoint,3)
 call allocate_array('deltax_rel_curve',deltax_rel_curve,2,maxpoint-1,3)
 call allocate_array('curvewave',curvewave,maxpoint-1,3)
 call allocate_array('flowflag_stripe',flowflag_stripe,maxpoint-1,3)
 call allocate_array('isinglestripe',isinglestripe,maxpoint-1,3)
 call allocate_array('width_stripe',width_stripe,maxpoint-1,3)
 call allocate_array('veli_stripe',veli_stripe,maxpoint-1,3)
 call allocate_array('veli_edgepoint',veli_edgepoint,maxpoint,3)

 call allocate_array('velmag_parts_stripe',velmag_parts_stripe,npart,maxpoint-1,3)
 call allocate_array('veli_parts_stripe',veli_parts_stripe,npart,maxpoint-1,3)

end subroutine allocate_ghostpart


subroutine deallocate_ghostpart()

 if (allocated(xyz_ghpart))       deallocate(xyz_ghpart)
 if (allocated(vxyz_ghpart))      deallocate(vxyz_ghpart)
 if (allocated(rho_ghpart))       deallocate(rho_ghpart)
 if (allocated(rcut_ghpart))      deallocate(rcut_ghpart)
 if (allocated(edgepoint))        deallocate(edgepoint)
 if (allocated(deltax_rel_curve)) deallocate(deltax_rel_curve)
 if (allocated(curvewave))        deallocate(curvewave)
 if (allocated(flowflag_stripe))  deallocate(flowflag_stripe)
 if (allocated(isinglestripe))    deallocate(isinglestripe)
 if (allocated(width_stripe))     deallocate(width_stripe)
 if (allocated(veli_stripe))      deallocate(veli_stripe)
 if (allocated(veli_edgepoint))   deallocate(veli_edgepoint)

 if (allocated(velmag_parts_stripe))  deallocate(velmag_parts_stripe)
 if (allocated(veli_parts_stripe))    deallocate(veli_parts_stripe)

end subroutine deallocate_ghostpart

!----------------------------------------------------------------
!+
! Sets up ghost particles based on motion of SPH particles at current timestep
!+
!----------------------------------------------------------------
subroutine set_ghostpart(npart,xyzh,vxyzu,dt)
 use io,   only:fatal
 integer, intent(in) :: npart
 real,    intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,    intent(in) :: dt
 real   :: sphp_velvec(3)

 !
 ! Check if ghost particles are still coupled to SPH particles
 !- If yes, keep wave velocities from previous iteration and continue to evolve,
 !- otherwise, redo the ghost particle setup process
 !
 resetup_ghostpart = .false.
 if (first_call) then
    resetup_ghostpart = .true.
    first_call = .false.
 else
    check_coupling: do igp = 1,nghpart
       ghp_pos  = xyz_ghpart(1:3,igp)
       ghp_velmag = mag(vxyz_ghpart(1:3,igp))
       ghp_velvec = vxyz_ghpart(1:3,igp)/ghp_velmag
       !- do brute-force for now; will use tree to get neighbouring sph particles,
       !- also need a better method to do the vel-checking
       nsphp = 0
       sphp_velmag    = 0.
       sphp_velvec(:) = (/ 0.,0.,0. /)
       neighbours: do ip = 1,npart
          sphp_pos = xyzh(1:3,ip)
          if (mag(sphp_pos-ghp_pos) < 1E-1) then
             nsphp = nsphp + 1
             sphp_velmag = sphp_velmag + mag(vxyzu(1:3,ip))
             sphp_velvec = sphp_velvec + vxyzu(1:3,ip)/mag(vxyzu(1:3,ip))
          endif
       enddo neighbours
       sphp_velmag = sphp_velmag/nsphp
       angle_diff = anglediff(ghp_velvec,sphp_velvec)
       ! if differed by too much, flag resetup_ghostpart
       if (abs(sphp_velmag-ghp_velmag) > velmag_diff_tol .or. &
           angle_diff > velvec_diff_tol) then
          resetup_ghostpart = .true.
          exit check_coupling
       endif
    enddo check_coupling
 endif
 !
 ! Detect SPH particle velocities to set up ghost particles
 !
 each_dimension: do idim = 1,3

    start    = edgemin(idim)
    end      = edgemax(idim)
    numpoint = npoint(idim)
    dx       = dx_curve(idim)

    redo_grid: if (resetup_ghostpart) then

       ! construct edgepoint array from deltax_rel_curve
       call build_edgepoint(idim,start,end,numpoint,dx)

       ! Search for flows in between each edgepoint and store velocities
       call check_sphpart_in_stripes(idim,start,end,numpoint,xyzh,vxyzu,npart)

       ! Check whether the flagged stripes are consecutive (i.e. resolved)
       call check_flowflag_consecutiveness(idim,numpoint,nsingstripe)

       ! Drag down curve around stripes where flow is detected
       iwave = 0
       do istripe = 1,numpoint-1
          if (flowflag_stripe(istripe,idim) == .true.) then
             cen_addwave = edgepoint(istripe,idim)
             ! drag down by approx 1/4 of dx plus a small amount proportional to velocity detected
             amp_addwave  = 1./4.*dx - sphpart_vel_ghpart_sep*abs(veli_stripe(istripe,idim))
             if (sphpart_vel_ghpart_sep*abs(veli_stripe(istripe,idim)) > 1E-1*dx) call fatal('ghostpart',&
                'insensible sphpart_vel_ghpart_sep; require a lower value ')
             sig_addwave  = width_stripe(istripe,idim)
             veli_addwave = veli_stripe(istripe,idim)
             call addwave_deltaxcurve(idim,numpoint,dx,iwave,cen_addwave,amp_addwave,sig_addwave,veli_addwave)
          endif
       enddo

       ! Redo grid until all flows are resolved
       if (nsingstripe > 0) then
          niter = 0
          ghpart_setup_done = .false.

          resolve_iter: do while (.not.ghpart_setup_done)
             niter = niter + 1
             if (niter > 10) call fatal('ghostpart','something wrong with resolve_iter')

             ! For the flows which are not yet resolved, drag down the curve even further
             do isingstripe = 1,nsingstripe
                istripe = isinglestripe(isingstripe,idim)
                cen_addwave  = edgepoint(istripe,idim)
                ! refine onto the small-scale flows by slowly bringing up wave amplitude
                amp_addwave  = (1.-1./(2.**niter))*dx - 1./niter*sphpart_vel_ghpart_sep*abs(veli_stripe(istripe,idim))
                sig_addwave  = 2.*width_stripe(istripe,idim) ! avoid sudden change in deltax
                veli_addwave = veli_stripe(istripe,idim)
                call addwave_deltaxcurve(idim,numpoint,dx,iwave,cen_addwave,amp_addwave,sig_addwave,veli_addwave)
             enddo

             ! do not resolve further than sphp mean resolution length
             if (dx - amp_addwave < hmean) ghpart_setup_done = .true.

             ! reconstruct edgepoint array from the updated deltax_rel_curve
             call build_edgepoint(idim,start,end,numpoint,dx)

            ! detect velocities with the new points and check again
             call check_sphpart_in_stripes(idim,start,end,numpoint,xyzh,vxyzu,npart)
             call check_flowflag_consecutiveness(idim,numpoint,nsingstripe)

             ! flag to stop if all flows are resolved
             if (nsingstripe == 0) ghpart_setup_done = .true.

          enddo resolve_iter
       endif

       ! build the final edgepoint and veli_edgepoint arrays
       call build_edgepoint(idim,start,end,numpoint,dx)
       call set_veli_edgepoint_from_stripe(idim,numpoint)

    else
       ! Using the new deltax_rel_curve updated from previous iteration
       ! reconstruct edgepoint array
       call build_edgepoint(idim,start,end,numpoint,dx)
       ! reconstruct veli_edgepoint array
       call set_veli_edgepoint_from_wave(idim,numpoint)

    endif redo_grid

    !
    ! Evolve the curve to prepare for the next timestep
    !
    call evolve_curvewave(idim,nwave,dt,start,end,numpoint,dx)

 enddo each_dimension

 !
 ! Output ghost particles positions and velocities
 !
 igp = 1
 do ix = 1,npoint(1)
    do iy = 1,npoint(2)
       do iz = 1,npoint(3)
          xyz_ghpart(:,igp)  = (/ edgepoint(ix,1),edgepoint(iy,2),edgepoint(iz,3) /)
          vxyz_ghpart(:,igp) = (/ veli_edgepoint(ix,1),veli_edgepoint(iy,2),veli_edgepoint(iz,3) /)
          igp = igp + 1
       enddo
    enddo
 enddo

 !
 ! Compute densities of ghost particles
 ! by summing masses from neighbouring sph particles (use tree?)
 !

 !- need to make sure that ghost particles gets mass from ALL sph particles
 !- to ensure mass conservation (?)


end subroutine set_ghostpart

!----------------------------------------------------------------
!+
! Procedures to operate on the ghost particle arrays for a given dimension
!+
!----------------------------------------------------------------
!
! Routine to reconstruct the edge array from the updated deltax-rel curve
!
subroutine build_edgepoint(idim,start,end,numpoint,dx)
 use io,    only:fatal
 integer, intent(in)  :: idim,numpoint
 real,    intent(in)  :: start,end,dx

 ! re-build edpoint array
 deltax_rel_sum = 0.
 edgepoint(1,idim) = start
 do ipoint = 1,numpoint-1
    edgepoint(ipoint+1,idim) = edgepoint(ipoint,idim) + deltax_rel_curve(2,ipoint,idim)
    deltax_rel_sum = deltax_rel_sum + deltax_rel_curve(2,ipoint,idim)
 enddo
 ! scale everything back up to match total length
 do ipoint = 1,numpoint
    edgepoint(ipoint) = edgepoint(ipoint) * (end-start)/deltax_rel_sum
 enddo
 if (abs(edgepoint(numpoint)-end) > epsilon(end)) call fatal('ghostpart','Error in setting up edgepoints')

 ! update width_stripe array
 do istripe = 1,numpoint-1
    width_stripe(istripe,idim) = edgepoint(ipoint+1,idim) - edgepoint(ipoint,idim)
 enddo

end subroutine build_edgepoint

!
! Routine to superpose a given gaussian wave onto the current deltax-rel curve
! i.e. dragging it down to reduce the separations
!
subroutine addwave_deltaxcurve(idim,numpoint,dx,iwave,centre,amplitude,sigma,veli_stripe)
 use io,    only:warning
 integer, intent(in)    :: idim,numpoint
 real,    intent(in)    :: dx,centre,amplitude,sigma,veli_stripe
 integer, intent(inout) :: iwave

 wave_has_effect = .false.

 do istripe = 1,numpoint-1
    x_val = deltax_rel_curve(1,istripe,idim)
    deltax_rel_waveval = dx - amplitude*exp(-(x_val-centre)**2/(2*sigma**2))
    !- Add to curve only if it is being dragged further down i.e. smaller than current value
    if (deltax_rel_waveval < deltax_rel_curve(2,istripe,idim)) then
       deltax_rel_curve(2,istripe,idim) = deltax_rel_waveval
       wave_has_effect = .true.
    endif
 enddo

 if (wave_has_effect) then
    ! store wave
    iwave = iwave + 1
    curvewave(iwave,idim)%amplitude = amplitude
    curvewave(iwave,idim)%centre    = centre
    curvewave(iwave,idim)%sigma     = sigma
    curvewave(iwave,idim)%velocity  = veli_stripe
 else
    call warning('ghostpart','flow detected but no changes are applied to deltax-curve')
 endif

end subroutine addwave_deltaxcurve

!
! Routine to reset deltax-rel curve to initial uniform separation
!
subroutine reset_deltaxcurve(idim,start,end,numpoint,dx)
 use io,    only:fatal
 integer, intent(in)    :: idim,numpoint
 real,    intent(in)    :: start,end,dx

 x_val = start
 do istripe = 1,numpoint-1
    deltax_rel_curve(1,istripe,idim) = x_val
    deltax_rel_curve(2,istripe,idim) = dx
    x_val = x_val + dx
 enddo
 if (abs(x_val-end) > epsilon(x_val)) call fatal('ghostpart','Error in resetting deltax-curve')

end subroutine reset_deltaxcurve

!
! Routine to bin sph particles according to the current edgepoints array
! and check for flows in each stripe.
! Sets the velocity veli for each stripe for evolving ghost particles
!
subroutine check_sphpart_in_stripes(idim,start,end,numpoint,xyzh,vxyzu,npart)
 integer, intent(in) :: idim,numpoint,npart
 real,    intent(in) :: start,end
 real,    intent(in) :: xyzh(:,:),vxyzu(:,:)
 integer :: pcount(numpoint-1)

 ! init particle counters for each stripe
 do istripe = 1,numpoint-1
    pcount(istripe) = 0
 enddo

 ! collect velocities of all particles falling within each stripe
 do ip = 1,npart
    xi_part = xyzh(idim,ip)
    iclosest = minloc(abs(edgepoint(1:numpoint,idim)-xi_part),1)
    if (xi_part < edgepoint(iclosest)) then
      istripe_part = iclosest-1
    else
      istripe_part = iclosest
    endif
    !- store velocity of this particle
    pcount(istripe_part) = pcount(istripe_part) + 1
    velmag_parts_stripe(pcount(istripe_part),istripe_part,idim) = mag(vxyzu(1:3,ip))
    veli_parts_stripe(pcount(istripe_part),istripe_part,idim)   = vxyzu(idim,ip)
 enddo

 ! search for astrophysical flows [require a better method]
 over_stripes: do istripe = 1,numpoint-1
    numpart_here = pcount(istripe)

    mean_velmag = 0.
    max_velmag  = epsilon(max_velmag)
    max_veli = epsilon(max_veli)
    over_parts_in_stripe: do ip = 1,numpart_here
       mean_velmag = mean_velmag + velmag_parts_stripe(ip,istripe,idim)
       max_velmag  = max(max_velmag,velmag_parts_stripe(ip,istripe,idim))
       !- estimate velocity of the flow by taking the largest particle vi
       veli_ip = veli_parts_stripe(ip,istripe,idim)
       if (abs(veli_ip) > abs(max_veli)) then
          max_veli = veli_ip
       endif
    enddo over_parts_in_stripe
    mean_velmag = mean_velmag/numpart_here

    !- we will roughly say there is a flow if max(vel) >> mean(vel)
    if (abs(max_velmag - mean_velmag) > 100.*mean_velmag) then
       flowflag_stripe(istripe,idim) = .true.
    else
       flowflag_stripe(istripe,idim) = .false.
    endif
    ! set veli_stripe to be the estimated velocity of the detected flow
    veli_stripe(istripe,idim) = max_veli
 enddo over_stripes

end subroutine check_sphpart_in_stripes

!
! Check whether or not all stripes that detected a flow are consecutive (i.e. resolved)
! returns icount; if icount=0 means all resolved
!
subroutine check_flowflag_consecutiveness(idim,numpoint,icount)
 integer, intent(in)  :: idim,numpoint
 integer, intent(out) :: icount

 ! Clear memory from prev runs
 do i = 1,numpoint-1
    isinglestripe(i,idim) = 0.
 enddo

 ! Find stripes with flows detected which are single
 icount = 0
 do istripe = 1,numpoint-1
    if (flowflag_stripe(istripe,idim) == .true.) then
       if (istripe == 1) then
          !- only check right-hand side
          if (flowflag_stripe(2,idim) == .false.) then
             icount = icount + 1
             isinglestripe(icount,idim) = istripe
          endif
       elseif (istripe == numpoint-1) then
          !- only check left-hand side
          if (flowflag_stripe(numpoint-2,idim) == .false.) then
             icount = icount + 1
             isinglestripe(icount,idim) = istripe
          endif
       else
          !- check both sides
          if (flowflag_stripe(istripe-1,idim) == .false. .and. &
              flowflag_stripe(istripe+1,idim) == .false.) then
             icount = icount + 1
             isinglestripe(icount,idim) = istripe
          endif
       endif
    endif
 enddo

end subroutine check_flowflag_consecutiveness

!
! Set edge point velocities with veli_stripe
!- to be then used for setting up ghost particle velocities
!  when re-doing the grid
!
subroutine set_veli_edgepoint_from_stripe(idim,numpoint)
 integer, intent(in) :: idim,numpoint

 do ipoint = 1,numpoint
    if (ipoint == 1 .or. ipoint == numpoint) then
       !- ghost particles at corners are stationary
       veli_edgepoint(ipoint,idim) = 0.
    else
       ! set velocity of point to be the mean between the two stripes around it
       vel_low  = veli_stripe(ipoint-1,idim)
       vel_high = veli_stripe(ipoint,idim)
       veli_edgepoint(ipoint,idim) = (vel_low+vel_high)/2.
    endif
 enddo

end subroutine set_veli_edgepoint_from_stripe

!
! Set edge point velocities with waves on the updated deltax-curve
!- to be then used for setting up ghost particle velocities
!  when NOT re-doing the grid
!
!- Note: waves on the evolved curve are the only thing that carries information
!  regarding the velocity of the previously detected flow, NOT the velocity of ghost
!  particles nor the stripes (since they do not actually move with the flow).
!  Hence veli_edgepoint will be re-built at every timestep regardless of whether or
!  not the grid is being re-constructed.
!
subroutine set_veli_edgepoint_from_wave(idim,numpoint)
 integer, intent(in) :: idim,numpoint

 do ipoint = 1,numpoint
    if (ipoint == 1 .or. ipoint == numpoint) then
      veli_edgepoint(ipoint,idim) = 0.
    else
      xi_point = edgepoint(ipoint,idim)
      do iwave = 1,numwave
          wavepos = curvewave(iwave,idim)%centre
          wavesig = curvewave(iwave,idim)%sigma
          if (xi_point >= wavepos-2.*wavesig .and. xi_point <= wavepos+2.*wavesig) then
             veli_edgepoint(ipoint,idim) = curvewave(iwave,idim)%velocity
          else
             veli_edgepoint(ipoint,idim) = 0. ! assume regions outside of waves are stationary
          endif
      enddo
    endif
 enddo

end subroutine set_veli_edgepoint_from_wave

!
! Evolves the waves on the curve for the next timestep
!
subroutine evolve_curvewave(idim,nwave,dt,start,end,numpoint,dx)
 integer, intent(in)    :: idim,numpoint
 integer, intent(inout) :: nwave
 real,    intent(in)    :: start,end,dx
 real,    intent(in)    :: dt

 ! reset curve back to uniform
 call reset_deltaxcurve(idim,start,end,numpoint,dx)

 ! For each wave stored, use its v, compute new pos on curve
 ! and add it back to the curve
 nwave_new = 0
 each_wave: do iwave = 1,nwave
    waveobj   = curvewave(iwave,idim)
    curr_pos  = waveobj%centre
    wave_veli = waveobj%velocity
    wave_amp  = waveobj%amplitude
    wave_sig  = waveobj%sigma
    !- clear previously stored wave
    call clear_wave(iwave,idim)
    ! update position of wave centre, keep the rest of the properties the same
    new_pos = curr_pos + wave_veli*dt
    if (new_pos >= start .and. new_pos <= end) then
       call addwave_deltaxcurve(idim,numpoint,dx,iwave,new_pos,wave_amp,wave_sig,wave_veli)
       nwave_new = nwave_new + 1
    endif
 enddo each_wave
 nwave = nwave_new

end subroutine evolve_curvewave


subroutine clear_wave(iwave,idim)
 integer, intent(in) :: idim,iwave

 curvewave(iwave,idim)%centre    = 0.
 curvewave(iwave,idim)%amplitude = 0.
 curvewave(iwave,idim)%sigma     = 0.
 curvewave(iwave,idim)%velocity  = 0.

end subroutine clear_wave


real function mag(vec)
 real, intent(in) :: vec(3)

 mag = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)

end function mag

real function anglediff(vec1,vec2)
 real, intent(in) :: vec1(3),vec2(3)

 anglediff = acos(dot_product(vec1,vec2)/mag(vec1)*mag(vec2))

end function anglediff

!-----------------------------------------------------------------------
!+
!  writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_ghostpart(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling ghost particles'
 call write_inopt(,'','',iunit)

end subroutine write_options_ghostpart

!-----------------------------------------------------------------------
!+
!  reads input options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_ghostpart(name,valstring,imatch,igotall,ierr)
 use io,    only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer,          save        :: ngot = 0
 character(len=30), parameter  :: label = 'read_options_ghostpart'

 imatch = .false.
 igotall = .true.
 ierr = 0

 select case(trim(name))
 case('')
    read(valstring,*,iostat=ierr)
    ngot = ngot + 1
    if ( <= 0.) call fatal(label,'invalid setting for  (<=0)')

 case default
    imatch = .false.
 end select
 igotall = ( ngot >=  )

end subroutine read_options_ghostpart

end module ghostpart
