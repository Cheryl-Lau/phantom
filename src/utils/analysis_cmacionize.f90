!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: analysis
!
!  DESCRIPTION:
!  Analysis routine for discs
!
!  REFERENCES: None
!
!  OWNER: Chris Nixon
!
!  $Id$
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: infile_utils, io, physcon
!+
!--------------------------------------------------------------------------
module analysis
 implicit none
 character(len=20), parameter, public :: analysistype = 'CMacIonize'
 public :: do_analysis, cmacionize_write_param

 private

contains

subroutine do_analysis(dumpfile,numfile,xyzh,vxyzu,pmass,npart,time,iunit)
! use io,      only:fatal
 use physcon, only:pi
 use units,          only:umass,utime,udist
 use eos, only : temperature_coef,gmw,gamma,polyk,equationofstate,ieos,get_temperature_from_ponrho

use cmi_fortran_library

 character(len=*), intent(in) :: dumpfile
 real,             intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(inout) :: pmass,time
 integer,          intent(in) :: npart,iunit,numfile
 character(len=19) :: output
 character(len=19) :: outstats
 character(len=30) :: snapshot
 integer           :: i, j, idist, idistm, idistp, ienergy_input, npart_zoom
 real              :: num_dens, volume, neutral_frac, factor, energy
 real ::  box_anchor(3), box_sides(3)
 real ::  pc, Myr, Rcmi, Rcmim, Rcmip, Rstr, Rsp, Rhi, Rhi1, dist, dist2, distm, distp, alphaH, csound
 real :: x(npart), y(npart), z(npart), h(npart), m(npart), nH(npart)
 real :: ponrhoi,spsoundi,rhoi, temperature, rad, con
 real :: xmin, xmax, ymin, ymax, zmin, zmax
 logical :: lloyd = .true.

! Print the analysis being done
 write(*,'("Performing analysis type ",A)') analysistype
 write(*,'("Input file name is ",A)') dumpfile

 !write(output,"(a10,i5.5,a4)") 'voro_sites',numfile,'.txt'
 write(output,"(a19)") 'voro_sites00000.txt'
 write(*,'("Output file name is ",A)') output
 !write(snapshot,"(a10,i5.5,a15)") 'voro_sites',numfile,'snapshot020.txt'
 write(snapshot,"(a30)") 'voro_sites00000snapshot020.txt'
 write(outstats,"(a10,i5.5,a4)") 'voro_stats',numfile,'.txt'

 factor = 1.0/(temperature_coef*gmw*(gamma-1))
 ienergy_input = 0

  ! pc in m
  pc = 3.086d16
  ! Myr in s
  Myr = 60.d0*60.d0*24.d0*365.25d0*1.d6

  idist = 0
  idistm = 0
  idistp = 0
  dist = 0.d0
  dist2 = 0.d0
  distm = 0.d0
  distp = 0.d0

  ! Alpha parameter from the paper
  alphaH = 2.7d-19
  ! Alpha parameter that Bert uses
  !alphaH = 4.91452e-19 ! m^3 s^-1

  ! Sound speed in ionized gas in m/s
  csound = 12850.d0

  Rstr = (0.75*1.d49/pi*1.67d-27/5.19d-18*1.67d-27/5.19d-18/alphaH)**real(1./3.)
  Rsp = Rstr * (1.d0 + 7.d0/4.d0*csound*(time-0)*utime/Rstr) ** real(4./7.)
  Rhi = Rstr * (1.d0 + sqrt(4.d0/3.d0)*7.d0/4.d0*csound*(time-0)*utime/Rstr) ** real(4./7.)

  if(time < 0.001) then
     Rhi1 = 1.5 * Rhi
  else
     Rhi1 = 1.2 * Rhi
  end if

  xmin = -Rhi1
  xmax = Rhi1
  ymin = -Rhi1
  ymax = Rhi1
  zmin = -Rhi1
  zmax = Rhi1

  npart_zoom = 0

  do i=1, npart
    rad = sqrt(xyzh(1,i)*xyzh(1,i) + xyzh(2,i)*xyzh(2,i) + xyzh(3,i)*xyzh(3,i))*udist/100.
    if(rad < Rhi1) then
!     if(abs(xyzh(1,i)*udist/100.)<Rhi1 .and. abs(xyzh(2,i)*udist/100.)<Rhi1 .and. &
!          & abs(xyzh(3,i)*udist/100.)<Rhi1) then
       npart_zoom = npart_zoom + 1
       x(npart_zoom) = xyzh(1,i)*udist/100.
       y(npart_zoom) = xyzh(2,i)*udist/100.
       z(npart_zoom) = xyzh(3,i)*udist/100.
       h(npart_zoom) = xyzh(4,i)*udist/100. *2.
       m(npart_zoom) = pmass*umass/1000.
       nH(npart_zoom) = 0.
    end if
  end do

  ! simulation box dimensions (these are the same as in test_CMI_library.param)
  box_anchor(1) = 1.1 * xmin
  box_anchor(2) = 1.1 * ymin
  box_anchor(3) = 1.1 * zmin
  box_sides(1) = 1.1 * (xmax - xmin)
  box_sides(2) = 1.1 * (ymax - ymin)
  box_sides(3) = 1.1 * (zmax - zmin)


 open(iunit,file=output)
 do i=1, npart_zoom
   write(iunit,'(3(es18.10,1X))') x(i), y(i), z(i)
 end do

 close(iunit)

 call cmacionize_write_param(output, npart_zoom, box_sides, box_anchor, numfile, lloyd)

! call cmi_init_periodic_dp("phantomtest.param", 1, 1.d0, 1.d0, box_anchor, &
!                            box_sides)
 call cmi_init("phantomtest.param", 16, 1.d0, 1.d0)


  ! run the simulation
  call cmi_compute_neutral_fraction_dp(x(:npart_zoom), y(:npart_zoom), z(:npart_zoom), h(:npart_zoom), &
       & m(:npart_zoom), nH(:npart_zoom), int8(npart_zoom))


 open(iunit,file=outstats)
 open(666,file='ionization-radius.txt',status='old',action='write',form='formatted',position="append")

 j = 1
 do i=1, npart
   rad = sqrt(xyzh(1,i)*xyzh(1,i) + xyzh(2,i)*xyzh(2,i) + xyzh(3,i)*xyzh(3,i))*udist/100.
   if(rad < Rhi1) then
      neutral_frac = nH(j)
      if (neutral_frac < 0.5) then
         energy = 10000. * 2. * factor
         if(energy > vxyzu(4,i)) then
            vxyzu(4,i) = energy
            ienergy_input = ienergy_input + 1
         end if
      endif

      if ((neutral_frac > 0.2).and.(neutral_frac < 0.8)) then
         rad =  sqrt(xyzh(1,i)*xyzh(1,i) + xyzh(2,i)*xyzh(2,i) + xyzh(3,i)*xyzh(3,i)) !+ xyzh(4,i)
         dist = dist + rad
         dist2 = dist2 + rad*rad
         idist = idist + 1
      end if

      write(iunit,'(3(es18.10,1X))') sqrt(xyzh(1,i)**2+xyzh(2,i)**2+xyzh(3,i)**2)*udist/100./pc, neutral_frac, vxyzu(4,i)
      j = j + 1
   end if
 end do

 print *, "Number of particles with modified energies: ", ienergy_input

 if(idist > 0) then
   Rcmi = dist/real(idist)
   con = udist/100./pc

   write(666,'(4(es18.10,1X))') time*utime/Myr, Rcmi*con, Rsp /pc, Rhi /pc
 else
   print *,"Poor resolution of the ionization front!"
 end if

 close(iunit)
 close(666)

  ! clean up the library
  call cmi_destroy()

end subroutine do_analysis

subroutine cmacionize_write_param(output, npart, box_sides, box_anchor, numfile, lloyd)
 character(len=*), intent(in) :: output
 integer,          intent(in) :: npart, numfile
 real,             intent(in) :: box_sides(3), box_anchor(3)
 logical,          intent(in) :: lloyd

 open(666,file='phantomtest.param')

 write(666,'(a11)') "Abundances:"
 write(666,'(2x,a11)') "carbon: 0.0"
 write(666,'(2x,a11)') "helium: 0.0"
 write(666,'(2x,a9)') "neon: 0.0"
 write(666,'(2x,a13)') "nitrogen: 0.0"
 write(666,'(2x,a11)') "oxygen: 0.0"
 write(666,'(2x,a12)') "sulphur: 0.0"

 write(666,'(a22)') "TemperatureCalculator:"
 write(666,'(2x,a33)') "do temperature calculation: false"

 write(666,'(a13)') "PhotonSource:"
 write(666,'(2x,a20)') "diffuse field: false"

 write(666,'(a14)') "SimulationBox:"
 write(666,'(2x,a9,e12.3,a4,e12.3,a4,e12.3,a3)') "anchor: [",box_anchor(1)," m, ",box_anchor(2)," m, ",box_anchor(3)," m]"
 write(666,'(2x,a8,e12.3,a4,e12.3,a4,e12.3,a3)') "sides: [",box_sides(1)," m, ",box_sides(2)," m, ",box_sides(3)," m]"
 write(666,'(2x,a34)') "periodicity: [false, false, false]"

 write(666,'(a12)') "DensityGrid:"
 write(666,'(2x,a13)') "type: Voronoi"
 write(666,'(2x,a14)') "grid type: Old"

 if(lloyd .eqv. .true.) then
 write(666,'(2x,a29)') "number of Lloyd iterations: 5"
 end if

 write(666,'(2x,a29)') "VoronoiGeneratorDistribution:"
 write(666,'(4x,a21,i7)') "number of positions: ",npart
 write(666,'(4x,a9)') "type: SPH"
 write(666,'(4x,a10,a19)') "filename: ", output

 write(666,'(a16)') "DensityFunction:"
 write(666,'(2x,a17)') "type: Homogeneous"

 write(666,'(a25)') "PhotonSourceDistribution:"
 write(666,'(2x,a16)') "type: SingleStar"
 write(666,'(2x,a28)') "position: [0. m, 0. m, 0. m]"
 write(666,'(2x,a23)') "luminosity: 1.0e49 s^-1"
 write(666,'(a21)') "PhotonSourceSpectrum:"
 write(666,'(2x,a19)') "type: Monochromatic"
 write(666,'(2x,a18)') "frequency: 13.6 eV"

 write(666,'(a21)') "IonizationSimulation:"
 write(666,'(2x,a24)') "number of iterations: 10"
 write(666,'(2x,a26)') "number of photons: 1000000"

 write(666,'(a18)') "DensityGridWriter:"
 write(666,'(2x,a15)') "type: AsciiFile"
 write(666,'(2x,a8,a15,a8,i2)') "prefix: ",output,"snapshot",numfile
 write(666,'(2x,a10)') "padding: 3"

 write(666,'(a14)') "CrossSections:"
 write(666,'(2x,a16)') "type: FixedValue"
 write(666,'(2x,a24)') "hydrogen_0: 6.3e-18 cm^2"
 write(666,'(2x,a16)') "helium_0: 0. m^2"
 write(666,'(2x,a16)') "carbon_1: 0. m^2"
 write(666,'(2x,a16)') "carbon_2: 0. m^2"
 write(666,'(2x,a18)') "nitrogen_0: 0. m^2"
 write(666,'(2x,a18)') "nitrogen_1: 0. m^2"
 write(666,'(2x,a18)') "nitrogen_2: 0. m^2"
 write(666,'(2x,a16)') "oxygen_0: 0. m^2"
 write(666,'(2x,a16)') "oxygen_1: 0. m^2"
 write(666,'(2x,a14)') "neon_0: 0. m^2"
 write(666,'(2x,a14)') "neon_1: 0. m^2"
 write(666,'(2x,a17)') "sulphur_1: 0. m^2"
 write(666,'(2x,a17)') "sulphur_2: 0. m^2"
 write(666,'(2x,a17)') "sulphur_3: 0. m^2"

 write(666,'(a19)') "RecombinationRates:"
 write(666,'(2x,a16)') "type: FixedValue"
 write(666,'(2x,a29)') "hydrogen_1: 2.7e-13 cm^3 s^-1"
 write(666,'(2x,a21)') "helium_1: 0. m^3 s^-1"
 write(666,'(2x,a21)') "carbon_2: 0. m^3 s^-1"
 write(666,'(2x,a21)') "carbon_3: 0. m^3 s^-1"
 write(666,'(2x,a23)') "nitrogen_1: 0. m^3 s^-1"
 write(666,'(2x,a23)') "nitrogen_2: 0. m^3 s^-1"
 write(666,'(2x,a23)') "nitrogen_3: 0. m^3 s^-1"
 write(666,'(2x,a21)') "oxygen_1: 0. m^3 s^-1"
 write(666,'(2x,a21)') "oxygen_2: 0. m^3 s^-1"
 write(666,'(2x,a19)') "neon_1: 0. m^3 s^-1"
 write(666,'(2x,a19)') "neon_2: 0. m^3 s^-1"
 write(666,'(2x,a22)') "sulphur_2: 0. m^3 s^-1"
 write(666,'(2x,a22)') "sulphur_3: 0. m^3 s^-1"
 write(666,'(2x,a22)') "sulphur_4: 0. m^3 s^-1"

 close(666)

end subroutine cmacionize_write_param

end module analysis
