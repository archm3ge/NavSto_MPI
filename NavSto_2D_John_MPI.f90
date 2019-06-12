!  Solves 2D periodic Navier Stokes equation using spectral methods
!  FFT using fftw library from REAL physical space to complex spectral space
!  Written by Hossein A Kafiabad fall 2017
!  time-stepping using forward Euler, 3-step Adams-Bashforth, and Crank-Nicolson (for the moment)

!       NOTES:
!       the equation that is solved:
!       zz : vorticity
!       Dzz/Dt = \nu.\nabla^n(zz) (n = 2 for the moment, newtonian dissipation)
!       psi stream function ---> \nabla^2(psi) = -zz
!       u =   psi_y zonal velocity (along x-axis)
!       v = - psi_x meridional velocity (along y-axis)

!       ADDITIONAL NOTES:
!       tint = time interval over which numerical method runs
!       tstop = final time used
!       tstart = initial time

!       ICREAD = UNUSED
!       ICRK = initial condition in REAL SPACE

!       VISC = molecular viscosity in the 2D Navier-Stokes equation

!       ikx = index from 1 to (size of mesh in x-direction)/2 + 1 (see param_John for size of mesh)
!       iky = index from 1 to size of mesh in y-direction (see param_John for size of mesh)
!       jj = variable used in calculating y-component of wave number vector
!       i = index variable only used in subroutines
!       j = index variable only used in subroutines
!       k = index variable for output file

!       kx = variable storing wave numbers created in kxa
!       ky = variable storing wave numbers created in kya
!       kh = square root of dissipation term k^2
!       k2 = dissipation term k^2
!       Lino = stand in variable for the term in Crank-Nicolson nu*k2/2
!       xx = index variable only used in subroutines
!       yy = index variable only used in subroutines

!       ur = real array of "u" component of velocity going into fft
!       uk = complex array of "u" component of velocity coming out of fft
!       vr = real array of "v" component of velocity going into fft
!       vk = complex array of "v" component of velocity coming out of fft
!       zzr = real array of vorticities going into fft
!       zzk = complex array of vorticities coming out of fft
!       nur = real array of updated velocity "u" component
!       nuk = complex array of updated velocity "u" component
!       nvr = real array of updated velocity "v" component
!       nvk = complex array of updated velocity "v" component

!       nnk = nonlinear term in numerical method
!       zzknew = updated vorticity from Crank-Nicholson plus nonlinear term

PROGRAM MAIN

! -----------------------------------------------------------------
! LIBRARIES AND STUFF
! -----------------------------------------------------------------

  use, intrinsic :: iso_c_binding
  use param_John_MPI
  use param_John_fftw_MPI
  implicit none

! -----------------------------------------------------------------
! ADDITIONAL PARAMETERS
! -----------------------------------------------------------------

! DECLARE TIME VARIABLES
  real, save :: tint = 0.5, tstart, tstop

  ! DISSIPATION
  real, parameter :: VISC = 1e-5

! -----------------------------------------------------------------
! VARIABLE DECLARATIONS
! -----------------------------------------------------------------

  ! WORKING VARIABLES
  integer :: ikx, iky, ikya, jj, i, j, k, switch = 0
  integer :: iter = 0
  real    :: kx, ky, kh, k2, Lino, xx, yy

  ! THE FIELDS WHOSE FOURIER TRANSFORM IS CALCULATED
  real(C_DOUBLE), dimension(:, :), pointer            :: ur, vr, zzr, nur, nvr
  complex(C_DOUBLE_COMPLEX), dimension(:, :), pointer :: uk, vk, zzk, nuk, nvk

  ! THE FIELDS WHOSE FOURIER TRANSFORM IS NOT CALCULATED
  complex(C_DOUBLE_COMPLEX), allocatable, dimension(:, :) :: nnk, nnk_nm1, nnk_nm2, zzknew, zzk_nm1, zzk_nm2
  real(C_DOUBLE), dimension(NX, NY) :: masterbuf

  ! TEMP
  complex(C_DOUBLE_COMPLEX), allocatable, dimension(:, :) :: fk
  real(C_DOUBLE), allocatable, dimension(:, :)            :: fr

  ! SUBROUTINES
  external :: initR, WriteFieldk, WriteFieldr

! -----------------------------------------------------------------
! INITIALIZATION OF MPI
! -----------------------------------------------------------------

  ! START FLOW OF INFORMATION

  ! INITIALIZE MPI
  call MPI_INIT(ierr)

  ! START TIMING
  start = MPI_WTIME()


  ! CHECK TO SEE IF MPI_INIT WAS CALLED
  call MPI_INITIALIZED(flag, ierr)

  if (flag.ne.1) then
    write(*, *) 'Error: MPI_INIT not called'
    call MPI_FINALIZE(ierr)
    stop
  end if

  ! PARALLELIZATION TOOLS
  comm = MPI_COMM_WORLD
  call MPI_COMM_SIZE(comm, size_d, ierr)
  call MPI_COMM_RANK(comm, rank_proc, ierr)
  call FFTW_MPI_INIT()

  if (rank_proc.eq.0) then
    write(*, *) 'We are using ', P, 'processes'
  end if

  ! CHECK WHETHER NUMBER OF PROCESSES IS CORRECT
  if (rank_proc.eq.0) then
    if (size_d.ne.P) then
        write(*, *) 'Error: Number of processes is ', size_d, ' while it should be ', P
        call MPI_FINALIZE(ierr)
        stop
    end if
  end if

  ! SET DOMAIN SIZE FOR EACH PROCESS
  if (rank_proc.eq.0) then
    if (mod(NY, P).ne.0) then
        write(*, *) 'Error: Size of domain of each process must be integer valued'
        call MPI_FINALIZE(ierr)
        stop
    end if
  end if

  Nxp   = NX
  Nyp   = NY/P
  iktyp = IKTY/P

! -----------------------------------------------------------------
! ALLOCATE ARRAYS NOT USED IN FFTW
! -----------------------------------------------------------------

  allocate(kxa(IKTX))
  allocate(kya(IKTY))
  allocate(L(IKTX, iktyp))
  allocate(nnk(IKTX, iktyp))
  allocate(nnk_nm1(IKTX, iktyp))
  allocate(nnk_nm2(IKTX, iktyp))
  allocate(zzknew(IKTX, iktyp))
  allocate(zzk_nm1(IKTX, iktyp))
  allocate(zzk_nm2(IKTX, iktyp))
  allocate(fr(Nxp, Nyp))

! -----------------------------------------------------------------
! ALLOCATE FFTW ARRAYS
! -----------------------------------------------------------------

! Use in-place transforms.
! Refer to Section 6.13 of FFTW v. 3.3.3. for details on how inputs are used.

  if(rank_proc.eq.0) then
    write(*, *) 'Initializing FFTW memory allocation'
  end if

  alloc_local = fftw_mpi_local_size_2d(IKTY, IKTX, MPI_COMM_WORLD, locy, locystart)
  u_p         = fftw_alloc_complex(alloc_local)
  call c_f_pointer(u_p, uk, [IKTX, iktyp])
  call c_f_pointer(u_p, ur, [NXD, Nyp])

  alloc_local = fftw_mpi_local_size_2d(IKTY, IKTX, MPI_COMM_WORLD, locy, locystart)
  v_p         = fftw_alloc_complex(alloc_local)
  call c_f_pointer(v_p, vk, [IKTX, iktyp])
  call c_f_pointer(v_p, vr, [NXD, Nyp])

  alloc_local = fftw_mpi_local_size_2d(IKTY, IKTX, MPI_COMM_WORLD, locy, locystart)
  zz_p        = fftw_alloc_complex(alloc_local)
  call c_f_pointer(zz_p, zzk, [IKTX, iktyp])
  call c_f_pointer(zz_p, zzr, [NXD, Nyp])

  alloc_local = fftw_mpi_local_size_2d(IKTY, IKTX, MPI_COMM_WORLD, locy, locystart)
  nu_p         = fftw_alloc_complex(alloc_local)
  call c_f_pointer(nu_p, nuk, [IKTX, iktyp])
  call c_f_pointer(nu_p, nur, [NXD, Nyp])

  alloc_local = fftw_mpi_local_size_2d(IKTY, IKTX, MPI_COMM_WORLD, locy, locystart)
  nv_p         = fftw_alloc_complex(alloc_local)
  call c_f_pointer(nv_p, nvk, [IKTX, iktyp])
  call c_f_pointer(nv_p, nvr, [NXD, Nyp])

  ! -----------------------------------------------------------------
  ! INITIALIZE FFTW PLANS (consider dimension are passed reversely)
  ! -----------------------------------------------------------------


  if (rank_proc.eq.0) then
    write(*, *) 'Initializing FFTW'
  end if

  ! FFTW_MEASURE tells machine to find an optimized fft by taking multiple ffts of the problem and measuring the computational time required by
  plan2_ur_uk   = fftw_mpi_plan_dft_r2c_2d(NY, NX, ur, uk, MPI_COMM_WORLD, FFTW_MEASURE)
  plan2_vr_vk   = fftw_mpi_plan_dft_r2c_2d(NY, NX, vr, vk, MPI_COMM_WORLD, FFTW_MEASURE)
  plan2_zzr_zzk = fftw_mpi_plan_dft_r2c_2d(NY, NX, zzr, zzk, MPI_COMM_WORLD, FFTW_MEASURE)
  plan2_nur_nuk = fftw_mpi_plan_dft_r2c_2d(NY, NX, nur, nuk, MPI_COMM_WORLD, FFTW_MEASURE)
  plan2_nvr_nvk = fftw_mpi_plan_dft_r2c_2d(NY, NX, nvr, nvk, MPI_COMM_WORLD, FFTW_MEASURE)

  plan2_uk_ur   = fftw_mpi_plan_dft_c2r_2d(NY, NX, uk, ur, MPI_COMM_WORLD, FFTW_MEASURE)
  plan2_vk_vr   = fftw_mpi_plan_dft_c2r_2d(NY, NX, vk, vr, MPI_COMM_WORLD, FFTW_MEASURE)
  plan2_zzk_zzr = fftw_mpi_plan_dft_c2r_2d(NY, NX, zzk, zzr, MPI_COMM_WORLD, FFTW_MEASURE)

! -----------------------------------------------------------------
! PRODUCE INITIAL CONDITION
! -----------------------------------------------------------------

  if(rank_proc.eq.0) then
    write(*, *) 'Finding initial condition'
  end if

  if (ICRK.eq.0) then
    ! I.C. described in real space
    call initR(fr, tstart)
  else
    ! I.C. described in spectral space
    ! call initK(zzk, tstart)
  end if

  zzr = fr

  call fftw_mpi_execute_dft_r2c(plan2_zzr_zzk, zzr, zzk)
  zzk=zzk/FFTNORM

 ! -----------------------------------------------------------------
 ! INITIALIZE WAVENUMBERS AND TRUNCATION MASK
 ! -----------------------------------------------------------------

! INITIALIZE THE 2D WAVENUMBER VECTOR K_HAT
! NOTE: 2*PI/LX = 1 and 2*PI/LY = 1

  if (rank_proc.eq.0) then
    write(*, *) 'Initializing wavenumbers'
  end if

  do  ikx = 1, IKTX
     kxa(ikx) = float(ikx-1)*2*PI/LX
  end do

  do iky = 1, iktyp
     ikya = rank_proc*iktyp + iky
     jj = ikya - 1
     if (ikya.gt.KTY)   jj = jj - 2*KTY
     if (ikya.eq.KTY+1) jj = 0
     if (ikya.gt.2*KTY) jj = 0
     kya(ikya) = float(jj)*2*PI/LY
  end do

  L  = 1
  do iky = 1, iktyp
     ikya = rank_proc*iktyp + iky
     ky = kya(ikya)

     do  ikx = 1, IKTX
        kx = kxa(ikx)
        kh = sqrt(kx*kx + ky*ky)
        if (ikya.eq.KTY+1)                        L(ikx, iky) = 0
        if (kx.eq.0 .and. ky.le.0)                L(ikx, iky) = 0
        ! truncation with 2/3 instead of 8/9 truncation
        if (abs(kx*LX/2.0/PI).gt.int(float(NX)*1./3. + 0.5)-0.5)    L(ikx, iky) = 0
        if (abs(ky*LY/2.0/PI).gt.int(float(NY)*1./3. + 0.5)-0.5)    L(ikx, iky) = 0
     end do

  end do

! -----------------------------------------------------------------
! 3-STEP ADAMS BASHFORTH STARTING PROCEDURE
! -----------------------------------------------------------------

  if (METHOD.eq."AB") then

    if (rank_proc.eq.0) then
        write(*, *) 'Using Adams-Bashforth 3-Step method'
    end if

    call wave_numbers(uk, vk, zzk)

    call fftw_mpi_execute_dft_c2r(plan2_zzk_zzr, zzk, zzr)
    call fftw_mpi_execute_dft_c2r(plan2_uk_ur, uk, ur)
    call fftw_mpi_execute_dft_c2r(plan2_vk_vr, vk, vr)

    nur = ur*zzr
    nvr = vr*zzr

    call fftw_mpi_execute_dft_r2c(plan2_nur_nuk, nur, nuk)
    nuk=nuk/FFTNORM
    call fftw_mpi_execute_dft_r2c(plan2_nvr_nvk, nvr, nvk)
    nvk=nvk/FFTNORM
    call fftw_mpi_execute_dft_r2c(plan2_zzr_zzk, zzr, zzk)
    zzk=zzk/FFTNORM

    do iky = 1, iktyp
        ikya = rank_proc*iktyp + iky
        ky = kya(ikya)

        do ikx = 1, IKTX
            kx = kxa(ikx)
            k2 = kx*kx + ky*ky
            nnk(ikx, iky) = - ZI * (kx * nuk(ikx, iky) + ky * nvk(ikx, iky))* L(ikx, iky)
            Lino =  0.5 * VISC * k2
            zzknew(ikx, iky) = zzk(ikx, iky) + delt*(-nnk(ikx, iky) - Lino*zzk(ikx, iky))
        end do

    end do

    time = tstart + delt

    zzk_nm1 = zzk
    zzk = zzknew
    nnk_nm1 = nnk

    call wave_numbers(uk, vk, zzk)

    call fftw_mpi_execute_dft_c2r(plan2_zzk_zzr, zzk, zzr)
    call fftw_mpi_execute_dft_c2r(plan2_uk_ur, uk, ur)
    call fftw_mpi_execute_dft_c2r(plan2_vk_vr, vk, vr)

    nur = ur*zzr
    nvr = vr*zzr

    call fftw_mpi_execute_dft_r2c(plan2_nur_nuk, nur, nuk)
    nuk=nuk/FFTNORM
    call fftw_mpi_execute_dft_r2c(plan2_nvr_nvk, nvr, nvk)
    nvk=nvk/FFTNORM
    call fftw_mpi_execute_dft_r2c(plan2_zzr_zzk, zzr, zzk)
    zzk=zzk/FFTNORM

    do iky = 1, iktyp
        ikya = rank_proc*iktyp + iky
        ky = kya(ikya)

        do ikx = 1, IKTX
            kx = kxa(ikx)
            k2 = kx*kx + ky*ky
            nnk(ikx, iky) = - ZI * (kx * nuk(ikx, iky) + ky * nvk(ikx, iky))* L(ikx, iky)
            Lino =  0.5 * VISC * k2
            zzknew(ikx, iky) = zzk(ikx, iky) + delt*(3./2.)*(-Lino*zzk(ikx, iky) - nnk(ikx, iky)) - delt*(1./2.)*(-Lino*zzk_nm1(ikx, iky) - nnk_nm1(ikx, iky))
        end do

    end do

    time = time + delt
    zzk_nm2 = zzk_nm1
    zzk_nm1 = zzk
    zzk = zzknew
    nnk_nm2 = nnk_nm1
    nnk_nm1 = nnk

  else if (method.eq."Forward") then
    if (rank_proc.eq.0) then
        write(*, *) 'Using forward Euler method'
    end if

    time = tstart           ! Set first time to starting time

  else
    if (rank_proc.eq.0) then
        write(*, *) 'Using Crank-Nicholson (trapezium) method'
    end if

    time = tstart

  end if

!-----------------------------------------------------------------
! EXECUTIONS
!-----------------------------------------------------------------

  tstop = tstart + tint

  do while (time.le.tstop)
    time = time + delt
    iter = iter + 1

    if (mod(iter, PRINTFREQ).eq.0) then
        if (rank_proc.eq.0) write(*, *) 'Time = ', time
    end if

  ! DOING ONE STEP OF TIME-STEPPING
  call wave_numbers(uk, vk, zzk)

  call fftw_mpi_execute_dft_c2r(plan2_zzk_zzr, zzk, zzr)
  call fftw_mpi_execute_dft_c2r(plan2_uk_ur, uk, ur)
  call fftw_mpi_execute_dft_c2r(plan2_vk_vr, vk, vr)

  nur = ur*zzr
  nvr = vr*zzr

  call fftw_mpi_execute_dft_r2c(plan2_nur_nuk, nur, nuk)
  nuk=nuk/FFTNORM
  call fftw_mpi_execute_dft_r2c(plan2_nvr_nvk, nvr, nvk)
  nvk=nvk/FFTNORM
  call fftw_mpi_execute_dft_r2c(plan2_zzr_zzk, zzr, zzk)
  zzk=zzk/FFTNORM

  do iky = 1, iktyp
     ikya = rank_proc*iktyp + iky
     ky = kya(ikya)

     do  ikx = 1, IKTX
        kx = kxa(ikx)
        k2 = kx*kx + ky*ky
        ! nnk is in the RHS of equation. Hence F^-1(nnk) = - ( u z_x + v z_y)
        nnk(ikx, iky) = - ZI * (kx * nuk(ikx, iky) + ky * nvk(ikx, iky))* L(ikx, iky)
        Lino =  0.5 * VISC * k2

        if (METHOD.eq."Forward") then
            zzknew(ikx, iky) = zzk(ikx, iky) + delt*(-nnk(ikx, iky) - Lino*zzk(ikx, iky))
        else if (METHOD.eq."AB") then
            zzknew(ikx, iky) = zzk(ikx, iky) + delt*23./12.*(-Lino*zzk(ikx, iky) - nnk(ikx, iky)) - &
            delt*4./3.*(-Lino*zzk_nm1(ikx, iky) - nnk_nm1(ikx, iky)) + delt*5./12.*(-Lino*zzk_nm2(ikx, iky) - nnk_nm2(ikx, iky))
        else
            zzknew(ikx, iky)=((1/delt - Lino)*zzk(ikx, iky)-nnk(ikx, iky))/(1/delt + Lino) &
                            * L(ikx, iky)
        end if

     end do

  end do

    if ((METHOD.eq."Forward").or.(METHOD.eq."Crank")) then
        zzk = zzknew
    else
        zzk_nm2 = zzk_nm1
        zzk_nm1 = zzk
        nnk_nm2 = nnk_nm1
        nnk_nm1 = nnk
        zzk = zzknew
    end if

  end do

  call fftw_mpi_execute_dft_c2r(plan2_zzk_zzr, zzk, zzr)

  ! TRANSFER ALL DATA BACK TO MASTER PROCESS
  call MPI_GATHER(zzr(1:Nxp, 1:Nyp), Nxp*Nyp, MPI_DOUBLE_PRECISION, masterbuf, Nxp*Nyp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  if (rank_proc.eq.0) then
    print*, 'zzr(1:4,1:4) = ', masterbuf(1:4,1:4)
    print*,  'final time = ', time
  end if

  ! FILE OUTPUT OF REAL VORTICITIES AT FINAL TIMESTEP
  if (rank_proc.eq.0) then

    filename = 'real_vorticities.csv'
    open (10, file = filename, status = 'unknown', action = 'write')
    write(*, *) 'Writing REAL vorticities for physical grid of size ', NX, ' by ', NY, 'to file: ', filename
    do k = 1, Nxp
        write (10, *) masterbuf(k, :)
    end do
    close(10)

  end if

  !-----------------------------------------------------------------
  ! CLEAN-UP AND DEALLOCATIONS
  !-----------------------------------------------------------------

  ! FREE MEMORY OF DYNAMICALLY ALLOCATED VARIABLES
  deallocate(kxa)
  deallocate(kya)
  deallocate(L)
  deallocate(nnk)
  deallocate(nnk_nm1)
  deallocate(nnk_nm2)
  deallocate(zzknew)
  deallocate(zzk_nm1)
  deallocate(zzk_nm2)
  deallocate(fr)

  ! DESTROY FFTW PLANS AND FREE MEMORY LOCATIONS
  call fftw_destroy_plan(plan2_ur_uk)
  call fftw_destroy_plan(plan2_uk_ur)
  call fftw_destroy_plan(plan2_vr_vk)
  call fftw_destroy_plan(plan2_vk_vr)
  call fftw_destroy_plan(plan2_zzr_zzk)
  call fftw_destroy_plan(plan2_zzk_zzr)
  call fftw_destroy_plan(plan2_nur_nuk)
  call fftw_destroy_plan(plan2_nvr_nvk)
  call fftw_free(u_p)
  call fftw_free(v_p)
  call fftw_free(zz_p)
  call fftw_free(nu_p)
  call fftw_free(nv_p)

  ! FIND FINAL WALLTIME
  finish = MPI_WTIME()

  ! End MPI
  call MPI_FINALIZE(ierr)

  ! FIND FINAL TIME
  difference = start - finish

  ! PRINT RESULT
  difference = finish - start

  filename = 'times.csv'
  inquire(file = filename, exist = exist)

  if (exist) then
    open(20, file = filename, status = 'old', position = 'append', action = 'write')
    write(20, *)
  else
    open(20, file = filename, status = 'new', action = 'write')
  end if

  write(20, *) P, difference
  close(20)

  ! LAST OUTPUT

  if (rank_proc.eq.0) then
    write(*, *)
    write(*, *) 'Numerical scheme complete'
  end if

  !-----------------------------------------------------------------
  ! CREATE INTERNAL SUBROUTINE FOR CALCULATING WAVE NUMBERS
  !-----------------------------------------------------------------

  contains

    ! CALCULATES VELOCITY COMPONENTS IN FOURIER SPACE
    subroutine wave_numbers(uk, vk, zzk)

        use param_John_MPI
        implicit none

        integer :: ikx, iky, ikya
        real :: kx, ky, k2
        complex(C_DOUBLE_COMPLEX), intent(in), dimension(:,:), pointer  :: zzk
        complex(C_DOUBLE_COMPLEX), intent(out), dimension(:,:), pointer :: uk, vk


        do iky = 1, iktyp
            ikya =  rank_proc*iktyp + iky
            ky = kya(ikya)

            do ikx = 1, IKTX
                kx = kxa(ikx)
                k2 = kx*kx + ky*ky
                if (L(ikx, iky).eq.1) then
                    uk(ikx, iky) =  ZI * ky * zzk(ikx, iky) / k2
                    vk(ikx, iky) = -ZI * kx * zzk(ikx, iky) / k2
                else
                    uk(ikx, iky) = 0
                    vk(ikx, iky) = 0
                end if
            end do

        end do

    end subroutine wave_numbers

END PROGRAM MAIN

! -----------------------------------------------------------------
! OTHER SUBROUTINES
! -----------------------------------------------------------------

! THIS SUBROUTINE DERIVES THE INITIAL CONDITION
subroutine initR(zzr, tstart)

! DERIVE THE INITIAL CONDITION IN REAL SPACE
  use param_John_MPI
  implicit none

  real(8), intent(out), dimension(Nxp, Nyp)   :: zzr
  real, intent(out)                           :: tstart
  real                                        :: xx, yy, uu, vv, u_y, v_x
  integer                                     :: i, j, ikya

  do j = 1, iktyp
    ikya = rank_proc*iktyp + j

    do i = 1, Nxp

        xx = (float(i-1)/NX)*LX
        yy = (float(ikya-1)/NY)*LY
        uu = sin(xx)*cos(yy)
        vv = -cos(xx)*sin(yy)
        u_y = -sin(xx)*sin(yy)
        v_x =  sin(xx)*sin(yy)
        zzr(i, j) = v_x-u_y

    end do

  end do

  tstart = 0.0

end subroutine initR

! -----------------------------------------------------------------
! DIAGNOSTICS
! -----------------------------------------------------------------

! THIS SUBROUTINE IS NOT USED IN THE CODE
! subroutine WriteFieldk(fk)
!
!  use param_John_MPI
!  implicit none
!
! complex(8), intent(in), dimension(IKTX,IKTY) :: fk
! integer              :: i,j
!  write(6,'(A12)',advance="no") ' ky|  kx -->'
!  write(6,5004) (kxa(i),i=1,IKTX)
!  write(6,'(A12)') ' v        '
!  do j = 1, IKTY
!     write(6,'(F6.2,"   ")',advance="no") kya(j)
!     write(6,'(*(F8.2,"+",F6.2,"i"))') (fk(i,j),i=1,IKTX)
!  end do
!  print*
!5004 format(*(F8.2,8x))
!
!end subroutine WriteFieldk

! THIS SUBROURTINE IS NOT USED IN THE CODE
!subroutine WriteFieldr(fr)
!
!  use param_John_MPI
!  implicit none
!
!  real(8), intent(in), dimension(NX,NY) :: fr
!  real                 :: xx(NX), yy(NY)
!  integer              :: i,j

  ! grid points
!  do i = 1, NX
!     xx(i) = (float(i-1)/NX)*LX
!  end do
!  do j = 1, NY
!     yy(j) = (float(j-1)/NY)*LY
!  end do
!
!  write(6,'(A12)',advance="no") ' yy|  xx -->'
!  write(6,5004) (xx(i),i=1,NX)
!  write(6,'(A12)') ' v        '
!  do j = 1, NY
!     write(6,'(F6.2,"   ")',advance="no") yy(j)
!     write(6,'(*(F9.3))') (fr(i,j),i=1,NX)
!  end do
!  print*
!5004 format(*(F5.2,4x))
!
!end subroutine WriteFieldr
