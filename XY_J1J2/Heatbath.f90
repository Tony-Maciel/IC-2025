!------------------------------------------------------------------------------------------!
! Simulates 2D XY model with nearest neighbors (NN) and next nearest neighbors (NNN).      !
!                                                                                          !
! Its hamiltonian is:                                                                      !
!               H = \sum_{<i,j>} \vec{S}_i \cdot \vec{S}_j +                               !
!                   + x \sum_{<<i,j>>} \vec{S}_i \cdot \vec{S}_j                           !
! (Uses Heat Bath algorithm with Binary Search)                                            !
!                                                                                          !
! To modify parameters, change "XY-J1J2_input.in".                                         !
!------------------------------------------------------------------------------------------!

! COMPILE WITH:
! gfortran -O3 -flto -c rkiss05.f90 -o fastkiss.o
! gfortran -O3 -march=native -flto -c Heatbath.f90 -o Heatbath.o 
! gfortran -O3 -flto Heatbath.o fastkiss.o -o Heatbath

! RUN WITH: 
! ./Heatbath 1 

!NOTE: There is NO LAG in between each burn in MC step.
!NOTE: (make sure to give an integer in CLI for program to initialize)

!NOTE: To speed up computations, this program doesn't calculate the errors of E, M, C or X, because it is meant 
! to be trivially parallelized on many nodes. So the data from these nodes constitute a set to collect averages 
! and errors for the final plots. 

module Tools 
  use , intrinsic :: iso_fortran_env, only: i64=>int64, dp=>real64, compiler_version, compiler_options
  implicit none

  real(dp), parameter :: twopi = 2.0_dp*acos(-1.0_dp)

  ! INFORMATION ABOUT MAGNETIC MATERIAL (Cannot write Total = L**dim or status(Total), ...) 
  type :: System
    integer                :: dim                  ! dimensions of sample. 
    integer                :: q                    ! 0, 1, ..., q-1 possible spin orientations
    integer                :: L                    ! Number of sites per dimension. 
    integer                :: Total                ! L^dim. Serves as one Monte Carlo step (MCS). 
    integer, allocatable   :: Neighbors(:,:)       ! map of neighbors with PBC. 
    integer, allocatable   :: Status(:)            ! integer state (0, q-1)
    real(dp), allocatable  :: sin_table(:)         ! sin lookup table 
    real(dp), allocatable  :: cos_table(:)         ! cos lookup table 
    real(dp), allocatable  :: P(:)                 ! boltzmann weights (for HB)
    real(dp)               :: x                    ! J2 / J1 

    contains 

      procedure, pass :: init_sys ! pass just tells the compiler that the first arg will be "self" 
      procedure, pass :: set_neighbors
      procedure, pass :: MCstep_HB

  end type 

  ! INFORMATION ABOUT MEASUREMENTS TAKEN FROM System 
  type :: Measurements 
    integer               :: N                      ! number of measurements taken (MCS == "Monte Carlo Steps") 
    integer               :: TNC                    ! Total Number of Configurations (for temporal correlation) 
    integer               :: burnin                 ! Burn in value before equilibrium (in MC steps) 
    integer               :: lag                    ! time between measurements (in MCS) 
    real(dp)              :: Ti, xi                 ! initial temperature and J2/J1
    real(dp)              :: Tf, xf                 ! final temperature and J2/J1
    real(dp)              :: dT, dx                 ! temperature step and J2/J1
    integer               :: length                 ! total number of points to plot (here zero) 
    integer               :: nboot                  ! number of resamples for bootstrap method. 
    real(dp), allocatable :: E(:), M(:), C(:), X(:) ! Thermodynamic variables vs temperatures 
    real(dp), allocatable :: m_phi(:)               ! another order parameter for finding BKT Tc 
    real(dp), allocatable :: mc_E(:), mc_M(:)       ! Variables for each temperature to get C and X (auxiliary) 
    real(dp), allocatable :: E_err(:), M_err(:)     ! errors of E and M at each temperature 
    real(dp), allocatable :: C_err(:), X_err(:)     ! errors of C and X at each temperature (bootstrap method) 
    real(dp), allocatable :: corr(:)                ! spatial correlation and temporal correlation 
    real(dp), allocatable :: bincu(:)               ! Binder cumulant for a given L over range of temperatures.
    real(dp), allocatable :: U_phi(:)               ! Another cumulant for a given L over range of temperatures.
    real(dp), allocatable :: U_Msigma(:)            ! Yet another cumulant for a given L over range of temperatures.
    real(dp), allocatable :: Upsilon(:)             ! Helicity modulus or spin stiffness
    real(dp), allocatable :: M_sigma_avg(:)         ! Average Collinear Order Parameter
    real(dp), allocatable :: M_sigma_susc(:)        ! "magnetic susceptibility" of nematic order parameter
    real(dp), allocatable :: ms_time_series(:)      ! M_sigma values at each MC step
    real(dp), allocatable :: ms_autocorr(:)         ! Calculated autocorrelation
    real(dp)              :: E1, E2, M1, M2, M4     ! auxiliary variables 
    real(dp)              :: Mp, Mp2, Mp4           ! more auxiliary variables 
    real(dp)              :: Msigma, Msigma2,Msigma4! even more auxiliary variables 
    real(dp)              :: tot_E, tot_M           ! for multiple histogram method 
    real(dp)              :: sx, jx                 ! Instantaneous Stiffness Energy and Current

    contains 

      procedure, pass :: init_meas

  end type

contains 
  
  ! INITIALIZES SYSTEM VARIABLES 
  subroutine init_sys(self, q, L, dim, random_start, x)
    class(System), intent(in out) :: self 
    real(dp), external            :: rkiss05 ! assumed to already be initialized. 
    integer, intent(in)           :: q, L, dim
    real(dp), intent(in)          :: x
    logical, intent(in)           :: random_start
    integer                       :: i
    
    self%q     = q
    self%L     = L 
    self%dim   = dim
    self%Total = int(self%L ** self%dim)
    self%x     = x

    if (self%dim /= 2) then
      print*, "ERROR: This program is hardcoded for dim=2."
      stop
    end if

    allocate(self%Neighbors(self%Total, 8), &
      self%Status(self%Total)             , &
      self%sin_table(0:self%q-1)          , &
      self%cos_table(0:self%q-1)          , &
      self%P(0:self%q-1))
    
    do i = 0, self%q-1
      self%cos_table(i) = cos(real(i, dp) * twopi / real(self%q, dp))
      self%sin_table(i) = sin(real(i, dp) * twopi / real(self%q, dp))
    end do
  
    if (random_start .eqv. .false.) then ! ordered start
      self%Status(:) = 0 
    else  ! disordered start
    do i = 1, self%Total
      self%Status(i) = int(rkiss05() * self%q)
    end do
    end if
  end subroutine init_sys

  ! INITIALIZES MEASUREMENTS VARIABLES 
  subroutine init_meas(self, L, N, burnin, lag, Ti, Tf, dT, TNC, nboot, xi, xf, dx)
    class(Measurements), intent(in out) :: self 
    integer, intent(in)                 :: L, N, burnin, Lag, TNC, nboot
    real(dp), intent(in)                :: Ti, Tf, dT, xi, xf, dx

    self%N      = N
    self%burnin = burnin 
    self%lag    = lag 
    self%Ti     = Ti
    self%Tf     = Tf
    self%dT     = dT
    self%xi     = xi 
    self%xf     = xf
    self%dx     = dx 
    self%TNC    = TNC
    self%nboot  = nboot
    self%length = int((self%Tf - self%Ti) / self%dT + 1)
    
    allocate(self%E(self%length)       , &
             self%M(self%length)       , &
             self%C(self%length)       , & 
             self%X(self%length)       , &
             self%mc_E(N)              , & 
             self%mc_M(N)              , &
             self%E_err(self%length)   , &
             self%M_err(self%length)   , & 
             self%C_err(self%length)   , &
             self%X_err(self%length)   , & 
             self%corr(0:L/2-1)        , &
             self%bincu(self%length)   , &
             self%m_phi(self%length)   , &
             self%U_phi(self%length)   , & 
             self%U_Msigma(self%length), &
             self%Upsilon(self%length) , &
             self%M_sigma_avg(self%length), &
             self%M_sigma_susc(self%length), &
             self%ms_time_series(N)    , &
             self%ms_autocorr(0:N/2))
  end subroutine init_meas 

  ! CREATES DATA STRUCTURE TO REPRESENT NEIGHBORS OF SPINS (2D ONLY!) 
  subroutine set_neighbors(self) 
    class(System), intent(in out) :: self 
    integer :: i, j, site
    integer :: i_up, i_dn, j_up, j_dn

    if (self%dim /= 2) then
        print*, "ERROR: set_neighbors is hardcoded for dim=2."
        stop
    end if

    do j = 0, self%L - 1  ! y-coordinate
      do i = 0, self%L - 1  ! x-coordinate
        site = i + j*self%L + 1 ! 1D index
        
        ! Periodic Boundary Conditions
        i_up = modulo(i+1, self%L)
        i_dn = modulo(i-1, self%L)
        j_up = modulo(j+1, self%L)
        j_dn = modulo(j-1, self%L)

        ! 4 Nearest Neighbors (NN)
        self%Neighbors(site, 1) = (i_up) + j*self%L + 1  ! right
        self%Neighbors(site, 2) = i + (j_dn)*self%L + 1  ! up 
        self%Neighbors(site, 3) = (i_dn) + j*self%L + 1  ! left
        self%Neighbors(site, 4) = i + (j_up)*self%L + 1  ! down

        ! 4 Next-Nearest Neighbors (NNN)
        self%Neighbors(site, 5) = (i_up) + (j_dn)*self%L + 1 ! up-right
        self%Neighbors(site, 6) = (i_dn) + (j_dn)*self%L + 1 ! up-left
        self%Neighbors(site, 7) = (i_dn) + (j_up)*self%L + 1 ! down-left
        self%Neighbors(site, 8) = (i_up) + (j_up)*self%L + 1 ! down-right
      end do
    end do
  end subroutine set_neighbors

  ! Frees all used memeory. 
  subroutine freeall(sys, meas)
    class(System), intent(in out)       :: sys 
    class(Measurements), intent(in out) :: meas

    deallocate(sys%Neighbors , &
               sys%Status    , &   
               sys%sin_table , &
               sys%cos_table , &
               sys%P         , &
               meas%E        , & 
               meas%M        , &
               meas%C        , & 
               meas%X        , &
               meas%mc_E     , & 
               meas%mc_M     , & 
               meas%E_err    , & 
               meas%M_err    , & 
               meas%C_err    , & 
               meas%X_err    , &
               meas%corr     , &
               meas%bincu    , &
               meas%m_phi    , & 
               meas%U_phi    , &
               meas%U_Msigma , &
               meas%Upsilon  , &
               meas%M_sigma_avg, &
               meas%M_sigma_susc, &
               meas%ms_time_series, &
               meas%ms_autocorr)
              
  end subroutine freeall

  ! CALCULATES ERROR OF CORRELATED SAMPLE (LANDAU) 
  pure real(dp) function error(x, mean, len) result(res)
    real(dp), intent(in) :: x(:) ! set 
    real(dp), intent(in) :: mean ! mean of set 
    integer, intent(in)  :: len  ! length of set 
    integer              :: i, j
    real(dp)             :: minus_mean, term1, term2, temp1, temp2
    
    minus_mean = - mean
    term1 = 0.0_dp 
    term2 = 0.0_dp 
    do i = 1, len 
      temp1 = (x(i) + minus_mean)*(x(i) + minus_mean)
      temp2 = 0.0_dp 
      do j = 1, (len - i)
        temp2 = temp2 + (x(j) + minus_mean)*(x(i+j) + minus_mean)
      end do
      term1 = term1 + temp1
      term2 = term2 + (len-i)*temp2
    end do

    res = sqrt((term1 + 2.0_dp / len * term2) / (len*(len - 1)))
  end function error

  ! ONE MC STEP (MCS) USING HEAT BATH ALGORITHM WITH BINARY SEARCH
  subroutine MCstep_HB(self, beta, exp_table)
    class(System), intent(in out) :: self
    real(dp), intent(in)          :: beta
    real(dp), intent(in)          :: exp_table(0:,0:)
    real(dp), external            :: rkiss05
    real(dp)                      :: weight, total_weight, r
    real(dp)                      :: cdf(0:self%q-1) ! CDF array
    integer                       :: site, i, j, k
    integer                       :: state1, state2, neighbstate1, neighbstate2
    integer                       :: lower, upper, mid, new_state
    
    do i = 1, self%Total 
        site = 1 + int(rkiss05() * self%Total)

        total_weight = 0.0_dp
        cdf = 0.0_dp

        ! Construct CDF for all possible states q
        do k = 0, self%q - 1
            weight = 1.0_dp
            
            ! Loop over 4 pairs of neighbors (NN and NNN)
            ! Neighbors 1-4 are NN (J1), Neighbors 5-8 are NNN (J2)
            ! The exp_table(d1, d2) contains exp(-beta*(cos(d1) + x*cos(d2)))
            do j = 1, 4
                neighbstate1 = self%Status(self%Neighbors(site, j))   ! NN
                neighbstate2 = self%Status(self%Neighbors(site, j+4)) ! NNN
                
                ! Calculate differences modulo q to stay in range [0, q-1]
                state1 = modulo(k - neighbstate1, self%q)
                state2 = modulo(k - neighbstate2, self%q)
                
                weight = weight * exp_table(state1, state2)
            end do

            if (k == 0) then
                cdf(k) = weight
            else
                cdf(k) = cdf(k-1) + weight
            end if
        end do
        
        ! Normalize CDF
        total_weight = cdf(self%q - 1)
        cdf(:) = cdf(:) / total_weight

        ! Sample using Binary Search
        r = rkiss05()
        
        lower = 0
        upper = self%q - 1
        new_state = upper 
        do while (lower < upper)
            mid = (lower + upper) / 2
            if (cdf(mid) > r) then
                new_state = mid
                upper = mid - 1
            else
                lower = mid + 1
            end if
        end do

        ! Update Spin
        self%Status(site) = new_state
    end do
  end subroutine MCstep_HB

  ! GETS THERMODYNAMIC VARIABLES AND WRITES IT TO MEASUREMENTS
  subroutine Get(sys, meas, mcs_idx) 
    class(System), intent(in)          :: sys 
    class(Measurements), intent(in out):: meas
    integer, intent(in)                :: mcs_idx
    real(dp)                           :: temp, Mx, My, E_temp
    integer                            :: site, r, neighb
    integer                            :: neighb_r, neighb_u, neighb_ur, neighb_ul
    integer                            :: diff_idx
    real(dp)                           :: s_i_x, s_i_y, s_tr_x, s_tr_y, s_r_x, s_r_y, s_tl_x, s_tl_y   ! variables for M_sigma
    real(dp)                           :: v1x, v1y, v2x, v2y, dot_prod, M_sigma_accum
    real(dp)                           :: Sx_temp, Jx_temp ! variables for spin stiffness

    Mx            = 0.0_dp 
    My            = 0.0_dp
    E_temp        = 0.0_dp
    M_sigma_accum = 0.0_dp
    Sx_temp       = 0.0_dp
    Jx_temp       = 0.0_dp

    do site = 1, sys%Total 
      ! MAGNETIZATION 
      Mx = Mx + sys%cos_table(sys%Status(site))
      My = My + sys%sin_table(sys%Status(site))
      
      ! ENERGY & STIFFNESS
      ! Use forward/half neighbors to avoid double counting for extensive quantities
      ! Neighbors map: 1:R, 2:U, 3:L, 4:D, 5:UR, 6:UL, 7:DL, 8:DR
      
      neighb_r = sys%Neighbors(site, 1) ! Right
      neighb_u = sys%Neighbors(site, 2) ! Up
      neighb_ur = sys%Neighbors(site, 5) ! Up-Right
      neighb_ul = sys%Neighbors(site, 6) ! Up-Left

      ! NN Interactions (J1)
      ! Right
      diff_idx = modulo(sys%Status(site) - sys%Status(neighb_r), sys%q) ! couldn't be abs()...
      E_temp   = E_temp + sys%cos_table(diff_idx)
      ! Stiffness X: dx = +1
      Sx_temp  = Sx_temp + sys%cos_table(diff_idx) ! dx^2 * cos = 1 * cos
      Jx_temp  = Jx_temp + sys%sin_table(diff_idx) ! dx * sin   = 1 * sin

      ! Up doesn't contribute to stiffness X
      diff_idx = modulo(sys%Status(site) - sys%Status(neighb_u), sys%q)
      E_temp   = E_temp + sys%cos_table(diff_idx)
      
      ! NNN Interactions (J2)
      ! Up-Right 
      diff_idx = modulo(sys%Status(site) - sys%Status(neighb_ur), sys%q)
      E_temp   = E_temp + sys%x * sys%cos_table(diff_idx)
      ! Stiffness X: dx = +1
      Sx_temp  = Sx_temp + sys%x * sys%cos_table(diff_idx) ! dx^2 = 1
      Jx_temp  = Jx_temp + sys%x * sys%sin_table(diff_idx) ! dx = 1
      
      ! Up-Left 
      diff_idx = modulo(sys%Status(site) - sys%Status(neighb_ul), sys%q)
      E_temp   = E_temp + sys%x * sys%cos_table(diff_idx)
      ! Stiffness X: dx = -1
      Sx_temp  = Sx_temp + sys%x * sys%cos_table(diff_idx) ! dx^2 = 1
      Jx_temp  = Jx_temp - sys%x * sys%sin_table(diff_idx) ! dx = -1 (NOTE MINUS SIGN)

      ! M_SIGMA (Collinear Order Parameter)
      ! Plaquette vertices: BL(site), BR(neighb_r), TL(neighb_u), TR(neighb_ur)

      s_i_x  = sys%cos_table(sys%Status(site))
      s_i_y  = sys%sin_table(sys%Status(site))
      s_r_x  = sys%cos_table(sys%Status(neighb_r))
      s_r_y  = sys%sin_table(sys%Status(neighb_r))
      s_tl_x = sys%cos_table(sys%Status(neighb_u))
      s_tl_y = sys%sin_table(sys%Status(neighb_u))
      s_tr_x = sys%cos_table(sys%Status(neighb_ur))
      s_tr_y = sys%sin_table(sys%Status(neighb_ur))

      ! V1: Main Diagonal (BL to TR) -> S_i - S_tr
      v1x = s_i_x - s_tr_x
      v1y = s_i_y - s_tr_y
      
      ! V2: Anti Diagonal (BR to TL) -> S_r - S_tl
      v2x = s_r_x - s_tl_x
      v2y = s_r_y - s_tl_y
      
      dot_prod = v1x*v2x + v1y*v2y
      
      ! M_sigma_accum = M_sigma_accum + dot_prod / abs(dot_prod) ! +/- 1
      if (dot_prod > 1.0e-10_dp) then
        M_sigma_accum = M_sigma_accum + 1.0_dp
      else if (dot_prod < -1.0e-10_dp) then
        M_sigma_accum = M_sigma_accum - 1.0_dp
      end if
    end do
    
    meas%tot_E = E_temp    
    meas%tot_M = Mx*Mx + My*My

    meas%E1  = E_temp / real(sys%Total,dp) 
    meas%E2  = meas%E1 * meas%E1  
    meas%M1  = sqrt(Mx*Mx + My*My) / real(sys%Total,dp)
    meas%M2  = meas%M1 * meas%M1
    meas%M4  = meas%M2 * meas%M2
    meas%Mp  = cos(sys%q * atan2(My, Mx)) 
    meas%Mp2 = meas%Mp * meas%Mp 
    meas%Mp4 = meas%Mp2 * meas%Mp2

    ! Normalization of new observables
    meas%Msigma  = M_sigma_accum / real(sys%Total, dp)
    meas%Msigma2 = meas%Msigma * meas%Msigma
    meas%Msigma4 = meas%Msigma2 * meas%Msigma2
    meas%sx      = Sx_temp / real(sys%Total, dp)
    meas%jx      = Jx_temp / real(sys%Total, dp)

    meas%mc_E(mcs_idx) = meas%E1 
  end subroutine Get

  ! PERFORMS SIMULATION FOR A GIVEN TEMPERATURE 
  subroutine Run(sys, meas, beta, idx, int_param)
    class(System), intent(in out)      :: sys 
    class(Measurements), intent(in out):: meas
    real(dp), intent(in)               :: beta 
    integer, intent(in)                :: idx, int_param 
    real(dp)                           :: Et, Mt, Ct, Xt
    real(dp)                           :: bc, mphi, mphi2, mphi4, energy_diff
    real(dp)                           :: Msig_t, Stiff_T_t, Stiff_I2_t, Msig2, Msig4
    integer                            :: i, j, k, io
    character(len=100)                 :: mcresults, config_results
    real(dp), allocatable              :: exp_table(:,:)
    real(dp)                           :: mean_ms, var_ms, sum_ms ! for autocorrelation
    integer                            :: tau, t, ac_io
    character(len=100)                 :: ac_filename
   
    allocate(exp_table(0:sys%q-1, 0:sys%q-1))
    do i = 0, sys%q-1
      do j = 0, sys%q-1
        energy_diff = sys%cos_table(i) + sys%x * sys%cos_table(j)
        exp_table(i,j) = exp(-beta * energy_diff)
      end do
    end do

    Et               = 0.0_dp 
    Mt               = 0.0_dp 
    Ct               = 0.0_dp 
    Xt               = 0.0_dp
    bc               = 0.0_dp
    mphi             = 0.0_dp 
    mphi2            = 0.0_dp 
    mphi4            = 0.0_dp
    meas%corr(:)     = 0.0_dp
    Msig_t           = 0.0_dp
    Msig2            = 0.0_dp 
    Msig4            = 0.0_dp
    Stiff_T_t        = 0.0_dp
    Stiff_I2_t       = 0.0_dp
    
    ! Burn in time 
    do i = 1,meas%burnin
      call sys%MCstep_HB(beta, exp_table)
    end do

    ! Measurement time 
    do i = 1, meas%N
      do j = 1,meas%lag
        call sys%MCstep_HB(beta, exp_table)
      end do
      call Get(sys, meas, i)

      Et         = Et + meas%E1 
      Mt         = Mt + meas%M1 
      Ct         = Ct + meas%E2
      Xt         = Xt + meas%M2
      bc         = bc + meas%M4 
      mphi       = mphi + meas%Mp
      mphi2      = mphi2 + meas%Mp2 
      mphi4      = mphi4 + meas%Mp4 
      Msig_t     = Msig_t + abs(meas%Msigma)
      Msig2      = Msig2 + meas%Msigma2
      Msig4      = Msig4 + meas%Msigma4
      Stiff_T_t  = Stiff_T_t + meas%sx
      Stiff_I2_t = Stiff_I2_t + (meas%jx * meas%jx)

      ! Store M_sigma for autocorrelation
      meas%ms_time_series(i) = abs(meas%Msigma)

    end do
    deallocate(exp_table)

    ! Autocorrelation
    mean_ms = sum(meas%ms_time_series) / real(meas%N, dp)
    var_ms = 0.0_dp
    do i = 1, meas%N
       var_ms = var_ms + (meas%ms_time_series(i) - mean_ms)**2
    end do
    var_ms = var_ms / real(meas%N, dp)

    do tau = 0, meas%N/2
       sum_ms = 0.0_dp
       do t = 1, meas%N - tau
          sum_ms = sum_ms + (meas%ms_time_series(t) - mean_ms) * &
                            (meas%ms_time_series(t + tau) - mean_ms)
       end do
       if (var_ms > 1.0e-12_dp) then
          meas%ms_autocorr(tau) = (sum_ms / real(meas%N - tau, dp)) / var_ms
       else
          meas%ms_autocorr(tau) = 0.0_dp
       end if
    end do

    ! Save autocorrelation to file
    write(ac_filename,'(A,F0.4,A)') "MS_Autocorr_HB_T_", 1.0_dp/beta, ".dat"
    open(newunit=ac_io, file=ac_filename, action="write")
    write(ac_io,*) "# Lag  Autocorrelation"
    do tau = 0, meas%N/2
       write(ac_io,*) tau, meas%ms_autocorr(tau)
    end do
    close(ac_io)
    
    ! get thermodynamic variables
    Et              = Et / real(meas%N,dp)  ! <E>
    Mt              = Mt / real(meas%N,dp)  ! <M>
    Ct              = Ct / real(meas%N,dp)  ! <E^2>
    Xt              = Xt / real(meas%N,dp)  ! <M^2>
    bc              = bc / real(meas%N,dp)  ! <M^4>
    mphi            = mphi / real(meas%N,dp)! <mphi> 
    mphi2           = mphi2 / real(meas%N,dp) ! <mphi^2>
    mphi4           = mphi4 / real(meas%N,dp) ! <mphi^4>
    Msig_t          = Msig_t / real(meas%N, dp) ! <Msigma>
    Msig2           = Msig2 / real(meas%N, dp)  ! <Msigma^2>
    Msig4           = Msig4 / real(meas%N, dp)  ! <Msigma^4>
    Stiff_T_t       = Stiff_T_t / real(meas%N, dp)  ! <e_x>
    Stiff_I2_t      = Stiff_I2_t / real(meas%N, dp) ! <j_x^2>

    meas%m_phi(idx) = mphi
    meas%corr(0)    = 1.0_dp - Mt
    meas%corr(1:)   = meas%corr(1:) / (real(sys%dim * sys%Total,dp) * real(meas%N,dp)) - Mt
    meas%E(idx)     = Et 
    meas%M(idx)     = Mt 
    meas%C(idx)     = (Ct - Et*Et) * beta * beta * sys%Total  
    meas%X(idx)     = (Xt - Mt*Mt) * beta * sys%Total         
    
    meas%E_err(idx) = 0.0_dp
    meas%M_err(idx) = 0.0_dp
    meas%bincu(idx) = 1.0_dp - bc / (3.0_dp * Xt*Xt) 
    meas%U_phi(idx) = 1.0_dp - mphi4 / (2.0_dp * mphi2*mphi2)
    meas%U_Msigma(idx) = 1.0_dp - Msig4 / (3.0_dp * Msig2 * Msig2)
    meas%M_sigma_avg(idx) = Msig_t
    meas%M_sigma_susc(idx) = (Msig2 - Msig_t*Msig_t) * beta * sys%Total
    ! Upsilon = -<e_x> - beta * N * <j_x^2>
    ! Note: Stiff_T_t (sx) calculates sum(cos).
    ! The stiffness term (second derivative) is sum(-cos) = -sx.
    ! So <T> = -Stiff_T_t.
    ! Result: -Stiff_T_t - beta * N * Stiff_I2_t
    meas%Upsilon(idx) = -Stiff_T_t - (beta * real(sys%Total, dp) * Stiff_I2_t)

    meas%C_err(idx) = 0.0_dp 
    meas%X_err(idx) = 0.0_dp
 end subroutine Run 
 
end module Tools

!######################################################################################################!
program main 
  use Tools   
  implicit none 
  
  real(dp), external  :: rkiss05 ! has to be double precision... 
  type(System)        :: sample
  type(Measurements)  :: vals
  character(len=100)  :: output, corr_output, param
  real(dp)            :: Ti, Tf, dT, current_temp, xi, xf, dx, current_x
  integer(i64)        :: io
  integer             :: i, j, k, L, dim, burnin, lag, mcs, TNC, int_param
  integer             :: x_length, nboot, n_sizes, seed, q, init_L
  logical             :: random_start, exists 

  ! Initialization
  call get_command_argument(1, param)
  read(param,*) int_param
  open(newunit=io, file="input_XYJ1J2.in", action="read")
    read(io,*) dim    ! dimensions of model 
    read(io,*) q      ! spin directions
    read(io,*) init_L ! initial number of sites per dimension 
    read(io,*) n_sizes! total number of system sizes to use (1->32, 2->32,64, ...) 
    read(io,*) random_start ! boolean.
    read(io,*) Ti     ! initial temperature 
    read(io,*) Tf     ! final temperature 
    read(io,*) dT     ! temperature increment 
    read(io,*) mcs    ! total Monte Carlo steps 
    read(io,*) burnin ! mcs before registering variables 
    read(io,*) TNC    ! Total number of configurations saved for tcorr 
    read(io,*) nboot  ! number of resamples ofr bootstrapping.
    read(io,*) lag    ! mcs in between each measurement 
    read(io,*) seed   ! prng seed 
    read(io,*) xi     ! J2 / J1 initial 
    read(io,*) xf     ! J2 / J1 final
    read(io,*) dx     ! J2 / J1 step
  close(io)

  call kissinit(seed + int_param) 
  L         = init_L
  x_length  = int((xf - xi) / dx + 1)
  do k = 1, n_sizes
    print*, "L = ", L
    call vals%init_meas(L=L,N=mcs,burnin=burnin,lag=lag,Ti=Ti,Tf=Tf,dT=dT,TNC=TNC,nboot=nboot,xi=xi,xf=xf,dx=dx)
    call sample%init_sys(q=q,L=L,dim=dim,random_start=random_start, x=xi) 
    call sample%set_neighbors()         

    current_x = xi
    do j = 1, x_length
      print*, "x = ", current_x
      write(output,'(A,I0,A,I0,A,I0,A)') "output_", L, "-J1J2-", int_param, "-", j,".dat" 
      inquire(file=output, exist=exists)
      if (exists) then 
        print*, output, "exists already"
        current_x = current_x + dx
        cycle 
      else

      ! FILE FOR THERMODYNAMIC VARIABLES AND BINDER CUMULANT 
      open(newunit=io, file=output, action="write")  
        write(io,*) "# Compiler version = ", compiler_version()
        write(io,*) "# Compiler flags = ", compiler_options()
        write(io,*) "# dimension = ", sample%dim 
        write(io,*) "# q = ", sample%q
        write(io,*) "# XY J1-J2 Model"
        write(io,*) "# J2 / J1 = ", current_x
        write(io,*) "# L = ", sample%L
        write(io,*) "# Seed = ", seed
        write(io,*) "# temperature interval (Ti, Tf, dT) = ", Ti, Tf, dT 
        write(io,*) "# random start =", random_start 
        write(io,*) "# Burn in time (in MC steps units) = ", vals%burnin
        write(io,*) "# Lag time (in MC steps units) = ", vals%lag
        write(io,*) "# Number of MC steps = ", vals%N
        write(io,*) "# T / E / Eerr / M / Merr / C / Cerr / X / Xerr / BC / Uphi / mphi / Upsilon / Msigma / U_Msigma / Msigma_susc"
      close(io)

      ! FILE FOR SPATIAL CORRELATION 
      ! write(corr_output,'(A,I0,A,I0,A)') "CORR_", L, "_XY-", int_param, ".dat" 
      ! open(newunit=io, file=corr_output, action="write")  
      ! write(io,*) "# Compiler version = ", compiler_version()
      ! write(io,*) "# Compiler flags = ", compiler_options()
      ! write(io,*) "# dimension = ", sample%dim 
      ! write(io,*) "# XY J1-J2 Model"
      ! write(io,*) "# L = ", sample%L
      ! write(io,*) "# Seed = ", seed
      ! write(io,*) "# temperature interval (Ti, Tf, dT) = ", Ti, Tf, dT 
      ! write(io,*) "# random start =", random_start
      ! write(io,*) "# Burn in time (in MC steps units) = ", vals%burnin
      ! write(io,*) "# Lag time (in MC steps units) = ", vals%lag
      ! write(io,*) "# Number of MC steps = ", vals%N
      ! write(io,*) "# Temperature / correlation at distance 0, 1, ..."
      ! close(io)
      
      current_temp = vals%Ti
      sample%x     = current_x
      do i = 1, vals%length 
        print*, "T = ", current_temp, " (", i, "/", vals%length, ")" 
        write(*,"(A, F5.1, A)", advance="no") "Completion: ", real(i-1,dp)*100_dp/vals%length, "%" 
        write(*, "(A)", advance="no") CHAR(13) 
        
        call run(sample, vals, 1.0_dp/current_temp, i, int_param)

        open(newunit=io, file=output, action="write", position="append")
          write(io,*) current_temp, vals%E(i), vals%E_err(i), vals%M(i), vals%M_err(i),&
            vals%C(i), vals%C_err(i), vals%X(i), vals%X_err(i), vals%bincu(i), vals%U_phi(i), &
            vals%m_phi(i), vals%Upsilon(i), vals%M_sigma_avg(i), vals%U_Msigma(i), vals%M_sigma_susc(i)
        close(io)

        ! open(newunit=io, file=corr_output, action="write", position="append")
        ! write(io,*) current_temp, vals%corr(:)
        ! close(io)
        current_temp = current_temp + vals%dT 
      end do ! temperatures
      current_x = current_x + dx
      end if
    end do ! x
    call freeall(sample,vals)
    L = L + L ! 8, 8 + 8 = 16 , 16 + 16 = 32, 32 + 32 = 64, ...
  end do ! L
end program main
