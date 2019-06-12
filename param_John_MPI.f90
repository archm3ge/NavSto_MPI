! ==========================================================
! ----------------------------------------------------------
! ++++++++++  GLOBAL PARAMETERS AND VARIABLES ++++++++++++++
! ----------------------------------------------------------

! Note:
! "ind" dimensionless wavenum. "dim" dimensional. dim wavenum = 2*pi/L ind wavenum
! --- make sure Nx & Ny are even preferably powers of small primes (like 2,3)

module param_John_MPI
  use, intrinsic :: iso_c_binding
  implicit none

! -----------------------------------------------------------------
! SPECIFY NUMERICAL METHOD
  character(len = 32), parameter :: METHOD = "Crank"
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! SET THE DIMENSION OF SPACE
  integer(C_INTPTR_T), parameter :: NX = 256, NY = 256                 ! Physical grid resolution in x and y
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! NUMBER OF PROCESSORS (must divide NY and IKTY)
  integer, parameter :: P = 1
! -----------------------------------------------------------------

! -----------------------------------------------------------------
  real, parameter :: FFTNORM = float(NX*NY)                            ! Fft-size to normalize the transform
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! MATH CONSTANTS
  real, parameter    :: PI = 2.*asin(1.)
  complex, parameter :: ZI = cmplx(0.,1.)
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! TIME VARIABLES
  real, save :: time
  real, save :: delt = 0.01/256
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! OUTPUTTING DIAGNOSTICS
  integer, parameter       :: PRINTFREQ = 1000
  logical, save            :: exist
  integer, parameter       :: MAXLEN = 32
  character*(MAXLEN), save :: filename
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! DOMAIN SIZE
  real, parameter :: LX = 2*PI, LY = 2*PI                              ! Domain size in physical space
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! PHYSICAL AND FOURIER ARRAY DIMENSIONS
  integer(C_INTPTR_T), parameter :: NXD = NX+2, NYD = NY               ! Size of physical arrays
  integer(C_INTPTR_T), parameter :: KTX = NX/2, KTY = NY/2             ! Max ind wavenumber
  integer(C_INTPTR_T), parameter :: IKTX = KTX+1, IKTY = NY            ! Total number of wavenumbers
  integer(C_INTPTR_T), save      :: Nxp, Nyp                           ! New domain size variables (Nxp = NX and Nyp = NY/P)
  integer(C_INTPTR_T), save      :: iktyp                              ! Number of wavenumbers in y-direction with MPI (IKTY/P)
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! WAVENUMBERS
  real, allocatable, dimension(:), save       :: kxa, kya              ! Dimension wavenumber
  integer, allocatable, dimension(:, :), save :: L                     ! Mask for wavenumber truncation
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! RELATED TO INITIAL CONDITION
  real, parameter :: ICREAD = 0                                        ! Reads the initial condition from an ext. file (UNUSED)
  real, parameter :: ICRK = 0                                          ! 0: I.C. set in real space (UNUSED)
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! OTHER MPI VARIABLE DECLARATIONS
  integer, save :: ierr, comm, size_d, rank_proc, uprank, dnrank       ! Variables for storing error code, communications link, size of domain, process ranks, and up one and down one process ranks
  integer(C_INTPTR_T), save :: locy, locystart                         ! Size and starting position of blocks in y-direction
  integer, save :: status                                              ! Variable for storing information about source of information and information tag
  integer, save :: flag                                                ! Returns TRUE if MPI_INIT called and FALSE if not
  integer(C_INTPTR_T), save :: alloc_local                             ! Local data size for malloc
! -----------------------------------------------------------------

! -----------------------------------------------------------------
! TIMING
  double precision, save :: start, finish, difference
! -----------------------------------------------------------------

end module param_John_MPI
