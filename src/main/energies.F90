!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2021 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.bitbucket.io/                                          !
!--------------------------------------------------------------------------!
module energies
!
! Calculates global quantities on the particles
!   To Developer: See instructions in evwrite.F90 about adding values
!                 to the .ev file
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: dim, dust, eos, externalforces, fastmath, gravwaveutils,
!   io, metric_tools, mpiutils, nicil, options, part, ptmass, units,
!   utils_gr, vectorutils, viscosity
!
 use dim, only: maxdusttypes,maxdustsmall
 use units, only:utime
 implicit none

 logical,         public    :: gas_only,track_mass,track_lum
 real,            public    :: ekin,etherm,emag,epot,etot,totmom,angtot,mtot,xyzcom(3)
 real,            public    :: vrms,rmsmach,accretedmass,mdust(maxdusttypes),mgas
 real,            public    :: xmom,ymom,zmom
 real,            public    :: totlum
 integer,         public    :: iquantities
 integer(kind=8), public    :: ndead,npartall,np_cs_eq_0,np_e_eq_0
 integer,         public    :: iev_time,iev_ekin,iev_etherm,iev_emag,iev_epot,iev_etot,iev_totmom,iev_com(3),&
                               iev_angmom,iev_rho,iev_dt,iev_dtx,iev_entrop,iev_rmsmach,iev_vrms,iev_rhop(6),&
                               iev_alpha,iev_B,iev_divB,iev_hdivB,iev_beta,iev_temp,iev_etao,iev_etah(2),&
                               iev_etaa,iev_vel,iev_vhall,iev_vion,iev_n(7),&
                               iev_dtg,iev_ts,iev_dm(maxdusttypes),iev_momall,iev_angall,iev_maccsink(2),&
                               iev_macc,iev_eacc,iev_totlum,iev_erot(4),iev_viscrat,iev_gws(4)
 integer,         public    :: iev_erad
 real,            public    :: erad
 integer,         parameter :: inumev  = 150  ! maximum number of quantities to be printed in .ev
 integer,         parameter :: iev_sum = 1    ! array index of the sum of the quantity
 integer,         parameter :: iev_max = 2    ! array index of the maximum of the quantity
 integer,         parameter :: iev_min = 3    ! array index of the minimum of the quantity
 integer,         parameter :: iev_ave = 4    ! array index of the average of the quantity
 ! Subroutines
 public  :: compute_energies,ev_data_update
 private :: get_erot,initialise_ev_data,collate_ev_data,finalise_ev_data
 ! Arrays
 real,             public :: ev_data(4,0:inumev),erot_com(6)

contains

!----------------------------------------------------------------
!+
!  Subroutine to compute global quantities on the particles
!+
!----------------------------------------------------------------
subroutine compute_energies(t)
 use dim,            only:maxp,maxvxyzu,maxalpha,maxtypes,mhd_nonideal,&
                          lightcurve,use_dust,store_temperature,&
                          maxdusttypes,gws,do_radiation
 use part,           only:rhoh,xyzh,vxyzu,massoftype,npart,maxphase,iphase,&
                          npartoftype,alphaind,Bxyz,Bevol,divcurlB,iamtype,&
                          igas,idust,iboundary,istar,idarkmatter,ibulge,&
                          nptmass,xyzmh_ptmass,vxyz_ptmass,isdeadh,&
                          isdead_or_accreted,epot_sinksink,imacc,ispinx,ispiny,&
                          ispinz,mhd,gravity,poten,dustfrac,eos_vars,itemp,&
                          nden_nimhd,eta_nimhd,iion,ndustsmall,graindens,grainsize,&
                          iamdust,ndusttypes,rad,iradxi
 use part,           only:pxyzu,fxyzu,fext
 use gravwaveutils,  only:calculate_strain
 use eos,            only:polyk,utherm,gamma,equationofstate,gamma_pwp
 use io,             only:id,fatal,master
 use externalforces, only:externalforce,externalforce_vdependent,was_accreted,accradius1
 use options,        only:iexternalforce,calc_erot,alpha,alphaB,ieos,use_dustfrac
 use mpiutils,       only:reduceall_mpi
 use ptmass,         only:get_accel_sink_gas
 use viscosity,      only:irealvisc,shearfunc
 use nicil,          only:nicil_update_nimhd,nicil_get_halldrift,nicil_get_ambidrift, &
                     use_ohm,use_hall,use_ambi,n_data_out,n_warn
#ifdef GR
 use part,           only:metrics,metricderivs
 use metric_tools,   only:unpack_metric
 use utils_gr,       only:dot_product_gr,get_geodesic_accel
 use vectorutils,    only:cross_product3D
#ifdef FINVSQRT
 use fastmath,       only:finvsqrt
#endif
#endif
#ifdef LIGHTCURVE
 use part,           only:luminosity
#endif
#ifdef KROME
 use part, only: gamma_chem
#endif
#ifdef DUST
 use dust,           only:get_ts,idrag
 integer :: iregime,idusttype
 real    :: tsi(maxdustsmall)
#endif
 real, intent(in) :: t
 real    :: ev_data_thread(4,0:inumev)
 real    :: xi,yi,zi,hi,vxi,vyi,vzi,v2i,Bxi,Byi,Bzi,Bi,B2i,rhoi,angx,angy,angz
 real    :: xmomacc,ymomacc,zmomacc,angaccx,angaccy,angaccz,xcom,ycom,zcom,dm
 real    :: epoti,pmassi,dnptot,dnpgas
 real    :: xmomall,ymomall,zmomall,angxall,angyall,angzall,rho1i,vsigi
 real    :: ponrhoi,spsoundi,dumx,dumy,dumz,divBi,hdivBonBi,alphai,valfven2i,betai
 real    :: n_total,n_total1,n_ion,shearparam_art,shearparam_phys,ratio_phys_to_av
 real    :: gasfrac,rhogasi,dustfracisum,dustfraci(maxdusttypes),dust_to_gas(maxdusttypes)
 real    :: etaohm,etahall,etaambi,vhall,vion
 real    :: curlBi(3),vhalli(3),vioni(3),data_out(n_data_out)
 real    :: erotxi,erotyi,erotzi,fdum(3)
 real    :: ethermi
#ifdef GR
 real    :: pdotv,bigvi(1:3),alpha_gr,beta_gr_UP(1:3),lorentzi,pxi,pyi,pzi
 real    :: gammaijdown(1:3,1:3),angi(1:3),fourvel_space(3)
#endif
 integer :: i,j,itype,iu
 integer :: ierrlist(n_warn)
 integer(kind=8) :: np,npgas,nptot,np_rho(maxtypes),np_rho_thread(maxtypes)

 real    :: hx,hp,hxx,hpp
 real, allocatable :: axyz(:,:)

 ! initialise values
 itype  = igas
 pmassi = massoftype(igas)
 ekin   = 0.
 etherm = 0.
 if (maxvxyzu < 4 .and. gamma < 1.0001 .and. ieos/=9) etherm = 1.5*polyk
 epot = 0.
 emag = 0.
 etot = 0.
 erad = 0.
 xcom = 0.
 ycom = 0.
 zcom = 0.
 mtot = 0.
 dm   = 0.
 xmom = 0.
 ymom = 0.
 zmom = 0.
 angx = 0.
 angy = 0.
 angz = 0.
 iu   = 4
 np   = 0
 vrms = 0.
 rmsmach = 0.
 npgas   = 0
 xmomacc = 0.
 ymomacc = 0.
 zmomacc = 0.
 angaccx = 0.
 angaccy = 0.
 angaccz = 0.
 mgas    = 0.
 mdust   = 0.
 mgas    = 0.
 np_cs_eq_0 = 0
 np_e_eq_0  = 0
 ierrlist = 0
 if (maxalpha==maxp) then
    alphai = 0.
 else
    alphai = alpha
 endif
 np_rho      = 0
 call initialise_ev_data(ev_data)

!$omp parallel default(none) &
!$omp shared(maxp,maxphase,maxalpha) &
!$omp shared(xyzh,vxyzu,pxyzu,rad,iexternalforce,npart,t,id,npartoftype) &
!$omp shared(alphaind,massoftype,irealvisc,iu) &
!$omp shared(ieos,gamma,nptmass,xyzmh_ptmass,vxyz_ptmass,xyzcom) &
!$omp shared(Bxyz,Bevol,divcurlB,alphaB,iphase,poten,dustfrac,use_dustfrac) &
!$omp shared(use_ohm,use_hall,use_ambi,nden_nimhd,eta_nimhd) &
!$omp shared(ev_data,np_rho,erot_com,calc_erot,gas_only,track_mass) &
!$omp shared(iev_erad,iev_rho,iev_dt,iev_entrop,iev_rhop,iev_alpha) &
!$omp shared(iev_B,iev_divB,iev_hdivB,iev_beta,iev_temp,iev_etao,iev_etah) &
!$omp shared(iev_etaa,iev_vel,iev_vhall,iev_vion,iev_n) &
!$omp shared(iev_dtg,iev_ts,iev_macc,iev_totlum,iev_erot,iev_viscrat) &
!$omp shared(eos_vars,grainsize,graindens,ndustsmall) &
#ifdef KROME
!$omp shared(gamma_chem) &
#endif
!$omp private(i,j,xi,yi,zi,hi,rhoi,vxi,vyi,vzi,Bxi,Byi,Bzi,Bi,B2i,epoti,vsigi,v2i) &
!$omp private(ponrhoi,spsoundi,ethermi,dumx,dumy,dumz,valfven2i,divBi,hdivBonBi,curlBi) &
!$omp private(rho1i,shearparam_art,shearparam_phys,ratio_phys_to_av,betai) &
!$omp private(gasfrac,rhogasi,dustfracisum,dustfraci,dust_to_gas,n_total,n_total1,n_ion) &
!$omp private(etaohm,etahall,etaambi,vhalli,vhall,vioni,vion,data_out) &
!$omp private(erotxi,erotyi,erotzi,fdum) &
!$omp private(ev_data_thread,np_rho_thread) &
!$omp firstprivate(alphai,itype,pmassi) &
#ifdef GR
!$omp shared(metrics) &
!$omp private(pxi,pyi,pzi,gammaijdown,alpha_gr,beta_gr_UP,bigvi,lorentzi,pdotv,angi,fourvel_space) &
#endif
#ifdef DUST
!$omp shared(idrag) &
!$omp private(tsi,iregime,idusttype) &
#endif
#ifdef LIGHTCURVE
!$omp shared(luminosity,track_lum) &
#endif
!$omp reduction(+:np,npgas,np_cs_eq_0,np_e_eq_0) &
!$omp reduction(+:xcom,ycom,zcom,mtot,xmom,ymom,zmom,angx,angy,angz,mdust,mgas) &
!$omp reduction(+:xmomacc,ymomacc,zmomacc,angaccx,angaccy,angaccz) &
!$omp reduction(+:ekin,etherm,emag,epot,erad,vrms,rmsmach,ierrlist)
 call initialise_ev_data(ev_data_thread)
 np_rho_thread = 0
!$omp do
 do i=1,npart
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    hi = xyzh(4,i)
    if (.not.isdead_or_accreted(hi)) then
       if (maxphase==maxp) then
          itype = iamtype(iphase(i))
          if (itype <= 0) call fatal('energies','particle type <= 0')
          pmassi = massoftype(itype)
       endif

       rhoi = rhoh(hi,pmassi)
       call ev_data_update(ev_data_thread,iev_rho,rhoi)
       if (.not.gas_only) then
          select case(itype)
          case(igas)
             call ev_data_update(ev_data_thread,iev_rhop(1), rhoi)
             np_rho_thread(igas) =  np_rho_thread(igas) + 1
          case(idust)
             call ev_data_update(ev_data_thread,iev_rhop(2),rhoi)
             np_rho_thread(idust) =  np_rho_thread(idust) + 1
          case(iboundary)
             call ev_data_update(ev_data_thread,iev_rhop(3), rhoi)
             np_rho_thread(iboundary) =  np_rho_thread(iboundary) + 1
          case(istar)
             call ev_data_update(ev_data_thread,iev_rhop(4),rhoi)
             np_rho_thread(istar) =  np_rho_thread(istar) + 1
          case(idarkmatter)
             call ev_data_update(ev_data_thread,iev_rhop(5),  rhoi)
             np_rho_thread(idarkmatter) =  np_rho_thread(idarkmatter) + 1
          case(ibulge)
             call ev_data_update(ev_data_thread,iev_rhop(6), rhoi)
             np_rho_thread(ibulge) =  np_rho_thread(ibulge) + 1
          end select
       endif

       np   = np + 1

       vxi  = vxyzu(1,i)
       vyi  = vxyzu(2,i)
       vzi  = vxyzu(3,i)

#ifdef GR
       pxi  = pxyzu(1,i)
       pyi  = pxyzu(2,i)
       pzi  = pxyzu(3,i)

       !  linear momentum
       xmom = xmom + pmassi*pxi
       ymom = ymom + pmassi*pyi
       zmom = zmom + pmassi*pzi

       call unpack_metric(metrics(:,:,:,i),betaUP=beta_gr_UP,alpha=alpha_gr,gammaijdown=gammaijdown)
       bigvi    = (vxyzu(1:3,i)+beta_gr_UP)/alpha_gr
       v2i      = dot_product_gr(bigvi,bigvi,gammaijdown)
#ifdef FINVSQRT
       lorentzi = finvsqrt(1.-v2i)
#else
       lorentzi = 1./sqrt(1.-v2i)
#endif
       pdotv    = pxi*vxi + pyi*vyi + pzi*vzi

       ! angular momentum
       fourvel_space = (lorentzi/alpha_gr)*vxyzu(1:3,i)
       call cross_product3D(xyzh(1:3,i),fourvel_space,angi) ! position cross with four-velocity
       angx = angx + pmassi*angi(1)
       angy = angy + pmassi*angi(2)
       angz = angz + pmassi*angi(3)

       ! kinetic energy
       ekin     = ekin + pmassi*(pdotv + alpha_gr/lorentzi - 1.) ! The 'kinetic term' in total specific energy, minus rest mass
#else
       ! centre of mass
       xcom = xcom + pmassi*xi
       ycom = ycom + pmassi*yi
       zcom = zcom + pmassi*zi
       mtot = mtot + pmassi

       ! linear momentum
       xmom = xmom + pmassi*vxi
       ymom = ymom + pmassi*vyi
       zmom = zmom + pmassi*vzi

       ! angular momentum
       angx = angx + pmassi*(yi*vzi - zi*vyi)
       angy = angy + pmassi*(zi*vxi - xi*vzi)
       angz = angz + pmassi*(xi*vyi - yi*vxi)

       ! kinetic energy & rms velocity
       v2i  = vxi*vxi + vyi*vyi + vzi*vzi
       ekin = ekin + pmassi*v2i
#endif

       vrms = vrms + v2i

       ! rotational energy around each axis through the Centre of mass
       ! note: for efficiency, centre of mass is from the previous time energies was called
       if (calc_erot) then
          call get_erot(xi,yi,zi,vxi,vyi,vzi,xyzcom,pmassi,erotxi,erotyi,erotzi)
          call ev_data_update(ev_data_thread,iev_erot(1),erotxi)
          call ev_data_update(ev_data_thread,iev_erot(2),erotyi)
          call ev_data_update(ev_data_thread,iev_erot(3),erotzi)
       endif

       if (iexternalforce > 0) then
          dumx = 0.
          dumy = 0.
          dumz = 0.
#ifdef GR
          epoti = 0.
#else
          call externalforce(iexternalforce,xi,yi,zi,hi,t,dumx,dumy,dumz,epoti,ii=i)
          call externalforce_vdependent(iexternalforce,xyzh(1:3,i),vxyzu(1:3,i),fdum,epoti)
#endif
          epot = epot + pmassi*epoti
       endif
       if (nptmass > 0) then
          dumx = 0.
          dumy = 0.
          dumz = 0.
          call get_accel_sink_gas(nptmass,xi,yi,zi,hi,xyzmh_ptmass,dumx,dumy,dumz,epoti)
          epot = epot + pmassi*epoti
       endif
       if (gravity) epot = epot + poten(i)
#ifdef DUST
       if (iamdust(iphase(i))) then
          idusttype = ndustsmall + itype - idust + 1
          mdust(idusttype) = mdust(idusttype) + pmassi
       endif
#endif
       if (do_radiation) erad = erad + pmassi*rad(iradxi,i)
       !
       ! the following apply ONLY to gas particles
       !
       isgas: if (itype==igas) then

          npgas = npgas + 1
          if (use_dustfrac) then
             dustfraci    = dustfrac(:,i)
             dustfracisum = sum(dustfraci)
             gasfrac      = 1. - dustfracisum
             dust_to_gas  = dustfraci(:)/gasfrac
             do j=1,ndustsmall
                call ev_data_update(ev_data_thread,iev_dtg,dust_to_gas(j))
             enddo
             mdust(1:ndustsmall) = mdust(1:ndustsmall) + pmassi*dustfraci(1:ndustsmall)
          else
             dustfraci    = 0.
             dustfracisum = 0.
             gasfrac      = 1.
          endif
          mgas = mgas + pmassi*gasfrac

          ! thermal energy
          if (maxvxyzu >= 4) then
             ethermi = pmassi*utherm(vxyzu(iu,i),rhoi)*gasfrac
#ifdef GR
             ethermi = (alpha_gr/lorentzi)*ethermi
#endif
             etherm = etherm + ethermi
#ifdef KROME
             call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xi,yi,zi,eni=vxyzu(iu,i),&
                                   gamma_local=gamma_chem(i))
#else
             if (store_temperature) then
                call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xi,yi,zi,vxyzu(iu,i),tempi=eos_vars(itemp,i))
             else
                call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xi,yi,zi,vxyzu(iu,i))
             endif
#endif
             if (vxyzu(iu,i) < tiny(vxyzu(iu,i))) np_e_eq_0 = np_e_eq_0 + 1
             if (spsoundi < tiny(spsoundi) .and. vxyzu(iu,i) > 0. ) np_cs_eq_0 = np_cs_eq_0 + 1
          else
             call equationofstate(ieos,ponrhoi,spsoundi,rhoi,xi,yi,zi)
             if (ieos==2 .and. gamma > 1.001) then
                !--thermal energy using polytropic equation of state
                etherm = etherm + pmassi*ponrhoi/(gamma-1.)*gasfrac
             elseif (ieos==9) then
                !--thermal energy using piecewise polytropic equation of state
                etherm = etherm + pmassi*ponrhoi/(gamma_pwp(rhoi)-1.)*gasfrac
             endif
             if (spsoundi < tiny(spsoundi)) np_cs_eq_0 = np_cs_eq_0 + 1
          endif
          vsigi = spsoundi
          ! entropy

#ifdef KROME
          call ev_data_update(ev_data_thread,iev_entrop,pmassi*ponrhoi*rhoi**(1.-gamma_chem(i)))
#else
          call ev_data_update(ev_data_thread,iev_entrop,pmassi*ponrhoi*rhoi**(1.-gamma))
#endif

#ifdef DUST
          ! min and mean stopping time
          if (use_dustfrac) then
             rhogasi = rhoi*gasfrac
             do j=1,ndustsmall
                call get_ts(idrag,j,grainsize(j),graindens(j),rhogasi,rhoi*dustfracisum,spsoundi,0.,tsi(j),iregime)
                call ev_data_update(ev_data_thread,iev_ts,tsi(j))
             enddo
          endif
#endif

#ifdef LIGHTCURVE
          if (track_lum) call ev_data_update(ev_data_thread,iev_totlum,real(luminosity(i)))
#endif

          ! rms mach number
          if (spsoundi > 0.) rmsmach = rmsmach + v2i/spsoundi**2

          ! max of dissipation parameters
          if (maxalpha==maxp) then
             alphai = alphaind(1,i)
             call ev_data_update(ev_data_thread,iev_alpha,alphai)
          endif

          ! physical viscosity
          if (irealvisc /= 0) then
             shearparam_art  = 0.1*alphai*hi*vsigi
             shearparam_phys = shearfunc(xi,yi,zi,spsoundi)
             if (shearparam_art > 0.) then
                ratio_phys_to_av = shearparam_phys/shearparam_art
             else
                ratio_phys_to_av = 0.
             endif
             call ev_data_update(ev_data_thread,iev_viscrat,ratio_phys_to_av)
          endif

          ! mhd parameters
          if (mhd) then
             Bxi       = Bevol(1,i)*rhoi
             Byi       = Bevol(2,i)*rhoi
             Bzi       = Bevol(3,i)*rhoi
             B2i       = Bxi*Bxi + Byi*Byi + Bzi*Bzi
             Bi        = sqrt(B2i)
             rho1i     = 1./rhoi
             valfven2i = B2i*rho1i
             vsigi     = sqrt(valfven2i + spsoundi*spsoundi)
             emag      = emag + pmassi*B2i*rho1i

             divBi     = abs(divcurlB(1,i))
             if (B2i > 0.) then
                hdivBonBi = hi*divBi/Bi
                betai     = 2.0*ponrhoi*rhoi/B2i ! plasma beta
             else
                hdivBonBi = 0.
                betai     = 0.
             endif
             call ev_data_update(ev_data_thread,iev_B,    Bi       )
             call ev_data_update(ev_data_thread,iev_divB, divBi    )
             call ev_data_update(ev_data_thread,iev_hdivB,hdivBonBi)
             call ev_data_update(ev_data_thread,iev_beta, betai    )

             if ( mhd_nonideal ) then
                call nicil_update_nimhd(0,etaohm,etahall,etaambi,Bi,rhoi, &
                                        eos_vars(itemp,i),nden_nimhd(:,i),ierrlist,data_out)
                curlBi = divcurlB(2:4,i)
                call ev_data_update(ev_data_thread,iev_temp,eos_vars(itemp,i))
                if (use_ohm) then
                   call ev_data_update(ev_data_thread,iev_etao,   etaohm      )
                endif
                if (use_hall) then
                   call nicil_get_halldrift(etahall,Bxi,Byi,Bzi,curlBi,vhalli)
                   vhall = sqrt( dot_product(vhalli,vhalli) )
                   call ev_data_update(ev_data_thread,iev_etah(1),etahall     )
                   call ev_data_update(ev_data_thread,iev_etah(2),abs(etahall))
                   call ev_data_update(ev_data_thread,iev_vhall  ,vhall       )
                endif
                if (use_ambi) then
                   call nicil_get_ambidrift(etaambi,Bxi,Byi,Bzi,curlBi,vioni)
                   vion   = sqrt( dot_product(vioni,  vioni  ) )
                   call ev_data_update(ev_data_thread,iev_etaa,   etaambi     )
                   call ev_data_update(ev_data_thread,iev_vel,    sqrt(v2i)   )
                   call ev_data_update(ev_data_thread,iev_vion,   vion        )
                endif
                n_ion = 0
                do j = 9,21
                   n_ion = n_ion + data_out(j)
                enddo
                n_total = n_ion + data_out(4)
                if (n_total > 0.) then
                   n_total1 = 1.0/n_total
                else
                   n_total1 = 0.0         ! only possible if eta_constant = .true.
                endif
                eta_nimhd(iion,i) = n_ion*n_total1    ! Save ionisation fraction for the dump file
                call ev_data_update(ev_data_thread,iev_n(1),n_ion*n_total1)
                call ev_data_update(ev_data_thread,iev_n(2),data_out( 8)*n_total1)
                call ev_data_update(ev_data_thread,iev_n(3),data_out( 8))
                call ev_data_update(ev_data_thread,iev_n(4),data_out( 4))
                call ev_data_update(ev_data_thread,iev_n(5),data_out(24))
                call ev_data_update(ev_data_thread,iev_n(6),data_out(23))
                call ev_data_update(ev_data_thread,iev_n(7),data_out(22))
             endif
          endif
       endif isgas

    elseif (was_accreted(iexternalforce,hi)) then
!
!--count accretion onto fixed potentials (external forces) separately
!
       vxi = vxyzu(1,i)
       vyi = vxyzu(2,i)
       vzi = vxyzu(3,i)
       if (maxphase==maxp) then
          pmassi = massoftype(iamtype(iphase(i)))
       else
          pmassi = massoftype(igas)
       endif
       xmomacc = xmomacc + pmassi*vxi
       ymomacc = ymomacc + pmassi*vyi
       zmomacc = zmomacc + pmassi*vzi

       angaccx = angaccx + pmassi*(yi*vzi - zi*vyi)
       angaccy = angaccy + pmassi*(zi*vxi - xi*vzi)
       angaccz = angaccz + pmassi*(xi*vyi - yi*vxi)

       call ev_data_update(ev_data_thread,iev_macc,pmassi)

    endif
 enddo
!$omp enddo
!
!--add contribution from sink particles
!
 if (id==master) then
    !$omp do
    do i=1,nptmass
       xi     = xyzmh_ptmass(1,i)
       yi     = xyzmh_ptmass(2,i)
       zi     = xyzmh_ptmass(3,i)
       pmassi = xyzmh_ptmass(4,i)
       if (pmassi < 0.) cycle

       vxi    = vxyz_ptmass(1,i)
       vyi    = vxyz_ptmass(2,i)
       vzi    = vxyz_ptmass(3,i)

       !phii   = fxyz_ptmass(4,i)

       xcom = xcom + pmassi*xi
       ycom = ycom + pmassi*yi
       zcom = zcom + pmassi*zi
       mtot = mtot + pmassi

       xmom   = xmom + pmassi*vxi
       ymom   = ymom + pmassi*vyi
       zmom   = zmom + pmassi*vzi

       angx   = angx + pmassi*(yi*vzi - zi*vyi)
       angy   = angy + pmassi*(zi*vxi - xi*vzi)
       angz   = angz + pmassi*(xi*vyi - yi*vxi)

       angx   = angx + xyzmh_ptmass(ispinx,i)
       angy   = angy + xyzmh_ptmass(ispiny,i)
       angz   = angz + xyzmh_ptmass(ispinz,i)

       v2i    = vxi*vxi + vyi*vyi + vzi*vzi
       ekin   = ekin + pmassi*v2i

       ! rotational energy around each axis through the origin
       if (calc_erot) then
          call get_erot(xi,yi,zi,vxi,vyi,vzi,xyzcom,pmassi,erotxi,erotyi,erotzi)
          call ev_data_update(ev_data_thread,iev_erot(1),erotxi)
          call ev_data_update(ev_data_thread,iev_erot(2),erotyi)
          call ev_data_update(ev_data_thread,iev_erot(3),erotzi)
       endif
    enddo
    !$omp enddo
 endif

!$omp critical(collatedata)
 call collate_ev_data(ev_data_thread,ev_data)
 if (.not.gas_only) then
    do i = 1,maxtypes
       np_rho(i) = np_rho(i) + np_rho_thread(i)
    enddo
 endif
!$omp end critical(collatedata)
!$omp end parallel

 !--Determing the number of active gas particles
 nptot    = reduceall_mpi('+',np)
 npgas    = reduceall_mpi('+',npgas)
 npartall = reduceall_mpi('+',npart)
 ndead    = npartall - nptot
 if (nptot > 0) then
    dnptot = 1./real(nptot)
 else
    dnptot = 0.
 endif
 if (npgas > 0) then
    dnpgas = 1./real(npgas)
 else
    dnpgas = 0.
 endif
 !--Number of gas particles without a sound speed or energy
 np_cs_eq_0 = reduceall_mpi('+',np_cs_eq_0)
 np_e_eq_0  = reduceall_mpi('+',np_e_eq_0)
 !--Finalise the arrays & correct as necessary;
 !  Almost all of the average quantities are over gas particles only
 call finalise_ev_data(ev_data,dnpgas)
#ifndef GR
 ekin = 0.5*ekin
#endif
 emag = 0.5*emag
 ekin = reduceall_mpi('+',ekin)
 if (maxvxyzu >= 4 .or. gamma >= 1.0001) etherm = reduceall_mpi('+',etherm)
 emag = reduceall_mpi('+',emag)
 epot = reduceall_mpi('+',epot)
 erad = reduceall_mpi('+',erad)
 if (nptmass > 1) epot = epot + epot_sinksink

 etot = ekin + etherm + emag + epot + erad

 xcom = reduceall_mpi('+',xcom)
 ycom = reduceall_mpi('+',ycom)
 zcom = reduceall_mpi('+',zcom)
 mtot = reduceall_mpi('+',mtot)
 if (mtot > 0.0) dm = 1.0 / mtot
 xcom = xcom * dm
 ycom = ycom * dm
 zcom = zcom * dm

 xmom = reduceall_mpi('+',xmom)
 ymom = reduceall_mpi('+',ymom)
 zmom = reduceall_mpi('+',zmom)
 totmom = sqrt(xmom*xmom + ymom*ymom + zmom*zmom)

 angx = reduceall_mpi('+',angx)
 angy = reduceall_mpi('+',angy)
 angz = reduceall_mpi('+',angz)
 angtot = sqrt(angx*angx + angy*angy + angz*angz)

 vrms    = reduceall_mpi('+',vrms)
 rmsmach = reduceall_mpi('+',rmsmach)
 vrms    = sqrt(vrms*dnptot)
 rmsmach = sqrt(rmsmach*dnpgas)

 !--fill in the relevant array elements for energy & momentum
 ev_data(iev_sum,iev_time  ) = t
 ev_data(iev_sum,iev_ekin  ) = ekin
 ev_data(iev_sum,iev_etherm) = etherm
 ev_data(iev_sum,iev_emag  ) = emag
 ev_data(iev_sum,iev_epot  ) = epot
 ev_data(iev_sum,iev_etot  ) = etot
 ev_data(iev_sum,iev_erad  ) = erad
 ev_data(iev_sum,iev_totmom) = totmom
 ev_data(iev_sum,iev_angmom) = angtot
 ev_data(iev_sum,iev_com(1)) = xcom
 ev_data(iev_sum,iev_com(2)) = ycom
 ev_data(iev_sum,iev_com(3)) = zcom
 do i=1,ndusttypes
    ev_data(iev_sum,iev_dm(i)) = mdust(i)
 enddo
 ev_data(iev_sum,iev_vrms   ) = vrms
 ev_data(iev_sum,iev_rmsmach) = rmsmach
 xyzcom(1) = xcom
 xyzcom(2) = ycom
 xyzcom(3) = zcom

 if (calc_erot) then
    ev_data(iev_sum,iev_erot(1)) = 0.5*ev_data(iev_sum,iev_erot(1))
    ev_data(iev_sum,iev_erot(2)) = 0.5*ev_data(iev_sum,iev_erot(2))
    ev_data(iev_sum,iev_erot(3)) = 0.5*ev_data(iev_sum,iev_erot(3))
    ev_data(iev_sum,iev_erot(4)) = sqrt(ev_data(iev_sum,iev_erot(1))**2 &
                                 +      ev_data(iev_sum,iev_erot(2))**2 &
                                 +      ev_data(iev_sum,iev_erot(3))**2)
 endif

 if (use_dust) then
    mgas  = reduceall_mpi('+',mgas)
    mdust = reduceall_mpi('+',mdust)
 endif

 if (.not. gas_only) then
    do i = 1,maxtypes
       np_rho(i) = reduceall_mpi('+',np_rho(i))
    enddo
    ! correct the average densities so that division is by n_p and not n_gas
    ev_data(iev_ave,iev_rho) = ev_data(iev_ave,iev_rho)*real(npgas)*dnptot
    if (np_rho(idust)       > 0) ev_data(iev_ave,iev_rhop(2)) = ev_data(iev_ave,iev_rhop(2))*real(npgas)/real(np_rho(idust))
    if (np_rho(iboundary)   > 0) ev_data(iev_ave,iev_rhop(3)) = ev_data(iev_ave,iev_rhop(3))*real(npgas)/real(np_rho(iboundary))
    if (np_rho(istar)       > 0) ev_data(iev_ave,iev_rhop(4)) = ev_data(iev_ave,iev_rhop(4))*real(npgas)/real(np_rho(istar))
    if (np_rho(idarkmatter) > 0) ev_data(iev_ave,iev_rhop(5)) = ev_data(iev_ave,iev_rhop(5))*real(npgas)/real(np_rho(idarkmatter))
    if (np_rho(ibulge)      > 0) ev_data(iev_ave,iev_rhop(6)) = ev_data(iev_ave,iev_rhop(6))*real(npgas)/real(np_rho(ibulge))
 endif

 if (iexternalforce > 0) then
    xmomacc   = reduceall_mpi('+',xmomacc)
    ymomacc   = reduceall_mpi('+',ymomacc)
    zmomacc   = reduceall_mpi('+',zmomacc)

    xmomall   = xmom + xmomacc
    ymomall   = ymom + ymomacc
    zmomall   = zmom + zmomacc
    ev_data(iev_sum,iev_momall) = sqrt(xmomall*xmomall + ymomall*ymomall + zmomall*zmomall)

    angaccx = reduceall_mpi('+',angaccx)
    angaccy = reduceall_mpi('+',angaccy)
    angaccz = reduceall_mpi('+',angaccz)
    angxall = angx + angaccx
    angyall = angy + angaccy
    angzall = angz + angaccz
    ev_data(iev_sum,iev_angall) = sqrt(angxall*angxall + angyall*angyall + angzall*angzall)
 endif

 if (track_mass) then
    accretedmass = ev_data(iev_sum,iev_macc)
    ev_data(iev_sum,iev_eacc) = accretedmass/accradius1 ! total accretion energy
 endif
 if (track_lum) totlum = ev_data(iev_sum,iev_totlum)

 if (gws) then
    allocate(axyz(3,npart))
#ifdef GR
    call get_geodesic_accel(axyz,npart,vxyzu(1:3,:),metrics,metricderivs)
#else
    axyz   = fxyzu(1:3,1:npart) + fext(1:3,1:npart)
#endif
    pmassi = massoftype(igas)
    call calculate_strain(hx,hp,hxx,hpp,xyzh,vxyzu(1:3,:),axyz,pmassi,npart)
    deallocate(axyz)
    hx  = reduceall_mpi('+',hx)
    hp  = reduceall_mpi('+',hp)
    hxx = reduceall_mpi('+',hxx)
    hpp = reduceall_mpi('+',hpp)
    ev_data(iev_sum,iev_gws(1)) = hx
    ev_data(iev_sum,iev_gws(2)) = hp
    ev_data(iev_sum,iev_gws(3)) = hxx
    ev_data(iev_sum,iev_gws(4)) = hpp
 endif

 return
end subroutine compute_energies
!----------------------------------------------------------------
!+
!  calculates rotational energy
!+
!----------------------------------------------------------------
subroutine get_erot(xi,yi,zi,vxi,vyi,vzi,xyzcom,pmassi,erotxi,erotyi,erotzi)
 real, intent(in)  :: xi,yi,zi,vxi,vyi,vzi,pmassi,xyzcom(3)
 real, intent(out) :: erotxi,erotyi,erotzi
 real              :: dx,dy,dz,dvx,dvy,dvz
 real              :: rcrossvx,rcrossvy,rcrossvz,radxy2,radyz2,radxz2
 !
 erotxi = 0.0
 erotyi = 0.0
 erotzi = 0.0

 dx  = xi  - xyzcom(1)
 dy  = yi  - xyzcom(2)
 dz  = zi  - xyzcom(3)
 dvx = vxi              ! results are less reliable if subtracting vcom
 dvy = vyi
 dvz = vzi

 rcrossvx = (dy*dvz - dz*dvy)
 rcrossvy = (dz*dvx - dx*dvz)
 rcrossvz = (dx*dvy - dy*dvx)

 radxy2 = dx*dx + dy*dy
 radyz2 = dy*dy + dz*dz
 radxz2 = dx*dx + dz*dz

 if (radyz2 > 0.) erotxi = pmassi*rcrossvx*rcrossvx/radyz2
 if (radxz2 > 0.) erotyi = pmassi*rcrossvy*rcrossvy/radxz2
 if (radxy2 > 0.) erotzi = pmassi*rcrossvz*rcrossvz/radxy2

end subroutine get_erot
!----------------------------------------------------------------
!+
!  initiallised the ev_data array
!+
!----------------------------------------------------------------
subroutine initialise_ev_data(evdata)
 real,    intent(inout) :: evdata(4,0:inumev)

 evdata            = 0.0
 evdata(iev_max,:) = -huge(evdata(iev_max,:))
 evdata(iev_min,:) =  huge(evdata(iev_min,:))

end subroutine initialise_ev_data
!----------------------------------------------------------------
!+
!  update the ev_data_array
!+
!----------------------------------------------------------------
subroutine ev_data_update(evdata,itag,val)
 integer, intent(in)    :: itag
 real,    intent(in)    :: val
 real,    intent(inout) :: evdata(4,0:inumev)

 evdata(iev_sum,itag) =     evdata(iev_sum,itag)+val
 evdata(iev_max,itag) = max(evdata(iev_max,itag),val)
 evdata(iev_min,itag) = min(evdata(iev_min,itag),val)

end subroutine ev_data_update
!----------------------------------------------------------------
!+
!  combines the ev_data from the various threads
!+
!----------------------------------------------------------------
subroutine collate_ev_data(evdata_thread,evdata)
 real,            intent(in)    :: evdata_thread(4,0:inumev)
 real,            intent(inout) :: evdata(4,0:inumev)
 integer                        :: i

 do i = 1,iquantities
    evdata(iev_sum,i) =     evdata(iev_sum,i)+evdata_thread(iev_sum,i)
    evdata(iev_max,i) = max(evdata(iev_max,i),evdata_thread(iev_max,i))
    evdata(iev_min,i) = min(evdata(iev_min,i),evdata_thread(iev_min,i))
 enddo

end subroutine collate_ev_data
!----------------------------------------------------------------
!+
!  Performs final generic housekeeping on the ev_data array
!+
!----------------------------------------------------------------
subroutine finalise_ev_data(evdata,dnptot)
 use mpiutils, only:reduceall_mpi
 real,            intent(inout) :: evdata(4,0:inumev)
 real,            intent(in)    :: dnptot
 integer                        :: i

 do i = 1,iquantities
    evdata(iev_sum,i) = reduceall_mpi('+',  evdata(iev_sum,i))
    evdata(iev_max,i) = reduceall_mpi('max',evdata(iev_max,i))
    evdata(iev_min,i) = reduceall_mpi('min',evdata(iev_min,i))
    evdata(iev_ave,i) = evdata(iev_sum,i)*dnptot
 enddo

end subroutine finalise_ev_data
!----------------------------------------------------------------
end module energies
