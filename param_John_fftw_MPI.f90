module param_John_fftw_MPI
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3-mpi.f03'
  include 'mpif.h'

! -----------------------------------------------------------------
! FFTW 2 PLANS
  type(C_PTR), save :: plan2_ur_uk, plan2_vr_vk                           ! Initialize fft plans for "u" and "v" velocity components from real to complex
  type(C_PTR), save :: plan2_zzr_zzk, plan2_nur_nuk, plan2_nvr_nvk        ! Initialize fft plans for vorticity and nonlinear term from real to complex
  type(C_PTR), save :: u_p, v_p, zz_p, nu_p, nv_p                         ! Initialize pointers used as memory locations for fft arrays

  type(C_PTR), save :: plan2_uk_ur, plan2_vk_vr                           ! Initialize fft plans for "u" and "v" velocity components from complex to real
  type(C_PTR), save :: plan2_zzk_zzr                                      ! Initialize fft plan for vorticity from complex to real

! -----------------------------------------------------------------

end module param_John_fftw_MPI
