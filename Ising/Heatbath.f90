!--------------------------------------------------------------------------------------------!
! N dimensional Ising model on cubic lattice, w/o external field, only NN, using the         !
! Heat bath algorithm.                                                                       !
!                                                                                            ! 
! To modify parameters, change "Ising_input.in".                                             !
!--------------------------------------------------------------------------------------------!

!NOTE: compile with:
! gfortran -O3 -march=native -flto rng_wrapper.o Heatbath.f90 -L. -lRANLUX++ -lranluxpp_c -Wl,-rpath='$ORIGIN' -o Heatbath
!NOTE: run and call python3 show.py
!NOTE: make sure to give an integer in CLI for program to initialize
!NOTE: gets <|m|> and NOT <m>
!NOTE: To speed up computations, this program doesn't calculate the errors of E, M, C or X, because it is meant 
! to be trivially parallelized on many nodes. So the data from these nodes constitute a set to collect averages 
! and errors for the final plots.

!----------------------------------------------------------------------------------------------------------------------

module Tools ! Ideally would have this in a separate file...
  use , intrinsic :: iso_fortran_env, only: i64=>int64, dp=>real64, compiler_version, compiler_options
  use rng_wrapper
  implicit none

  real(dp), parameter:: Tc = 2.0_dp/log(1.0_dp+sqrt(2.0_dp))

  ! INFORMATION ABOUT MAGNETIC MATERIAL (Cannot write Total = L**dim or status(Total), ...)
  type :: System
    integer                :: L                    ! Number of sites per dimension.
    integer                :: dim                  ! Dimensions of sample. 
    integer                :: Total                ! L^dim. Serves as Monte Carlo Step (MCS).
    integer, allocatable   :: Neighbors(:,:)       ! Map of neighbors with periodic boundary conditions (PBC).
    real(dp), allocatable  :: Status(:)            ! +1 or -1 for each site. In sum(), result is same type as elements...
    integer, allocatable   :: positions(:)         ! Auxiliary vector just to set up neighbors.
    real(dp), allocatable  :: expdus(:)            ! To avoid calculating exp(). effective size is 2*dim + 1.

    contains 

      procedure, pass :: init_sys ! pass just tells the compiler that the first arg will be "self"
      procedure, pass :: f, set_neighbors, MCstep

  end type 

  ! INFORMATION ABOUT MEASUREMENTS TAKEN FROM System
  type :: Measurements 
    integer               :: N                      ! Number of measurements taken
    integer               :: TNC                    ! Total Number of Configurations (for temporal correlation)
    integer               :: burnin                 ! Burn in value before equilibrium (in MCS)
    integer               :: lag                    ! Time between measurements (in MCS)
    real(dp)              :: Ti                     ! Initial temperature 
    real(dp)              :: Tf                     ! Final temperature 
    real(dp)              :: dT                     ! Temperature step 
    integer               :: length                 ! Total number of temperatures simulated
    integer               :: nboot                  ! Number of resamples for bootstrap method.
    real(dp), allocatable :: E(:), M(:), C(:), X(:) ! Thermodynamic variables vs temperatures
    real(dp), allocatable :: mc_E(:), mc_M(:)       ! Variables for each temperature to get C and X (auxiliary)
    real(dp), allocatable :: E_err(:), M_err(:)     ! Errors of E and M at each temperature
    real(dp), allocatable :: C_err(:), X_err(:)     ! Errors of C and X at each temperature (bootstrap method)
    real(dp), allocatable :: corr(:), tcorr(:)      ! Spatial correlation and temporal correlation
    real(dp), allocatable :: bincu(:)               ! Binder cumulant for a given L over range of temperatures.
    real(dp)              :: E1, E2, M1, M2, M4     ! Auxiliary variables for multiple histogram method
    real(dp)              :: tot_E, tot_M           ! ...Same 
    contains 

      procedure, pass :: init_meas

  end type

contains 
  
  ! INITIALIZES SYSTEM VARIABLES 
  subroutine init_sys(self, L, dim, p, rng)
    class(System), intent(in out) :: self 
    integer, intent(in)           :: L, dim
    real(dp), intent(in)          :: p   ! probability of having spin up initially.
    type(rng_t), intent(in out)   :: rng ! asssumes prng is already initialized.
    integer                       :: i
    
    self%L     = L
    self%dim   = dim
    self%Total = int(self%L ** self%dim)
    allocate(self%Neighbors(self%Total, 2*self%dim), &
      self%Status(self%Total)            , &
      self%positions(self%dim)           , &
      self%expdus(-4*self%dim:4*self%dim))
    
    do i = 1, self%Total
      if (next_rng(rng) .lt. p) then 
        self%Status(i) = 1.0_dp 
      else 
        self%Status(i) = -1.0_dp
      end if
    end do
  end subroutine init_sys

  ! INITIALIZES MEASUREMENTS VARIABLES 
  subroutine init_meas(self, L, N, burnin, lag, Ti, Tf, dT, TNC, nboot)
    class(Measurements), intent(in out) :: self 
    integer, intent(in)                 :: L, N, burnin, Lag, TNC, nboot
    real(dp), intent(in)                :: Ti, Tf, dT

    self%N      = N
    self%burnin = burnin 
    self%lag    = lag 
    self%Ti     = Ti*Tc
    self%Tf     = Tf*Tc
    self%dT     = dT
    self%TNC    = TNC
    self%nboot  = nboot
    self%length = int((self%Tf - self%Ti) / self%dT + 1)
    allocate(self%E(self%length),     &
             self%M(self%length),     &
             self%C(self%length),     &
             self%X(self%length),     &
             self%mc_E(N),            & 
             self%mc_M(N),            &
             self%E_err(self%length), &
             self%M_err(self%length), & 
             self%C_err(self%length), &
             self%X_err(self%length), & 
             self%corr(0:L/2-1),      &
             self%tcorr(0:TNC-1),     & 
             self%bincu(self%length))
  end subroutine init_meas 

  ! AUXILIARY TO SETTING UP NEIGHBORS
  pure integer function f(self) result(res)
    class(System), intent(in) :: self
    integer                   :: i
    
    res = 1
    do i = 1, self%dim 
      res = res + (self%positions(i) - 1)*self%L**(self%dim - i)
    end do
  end function f

  ! CREATES DATA STRUCTURE TO REPRESENT NEIGHBORS OF SPINS (N-DIMENSIONAL)
  recursive subroutine set_neighbors(self, depth) 
    class(System), intent(in out) :: self 
    integer, intent(in)           :: depth   ! recursion parameter
    integer                       :: i, j, temp, temp2

    if (depth .gt. self%dim) then ! found a combination
      temp = self%f()
      do j = 1, self%dim          ! for first N
        temp2 = self%positions(j) ! save value
        self%positions(j) = modulo(temp2,self%L) + 1
        self%Neighbors(temp,j) = self%f()  ! coordinates of neighbors
        self%positions(j) = temp2 ! put back value
      end do

      do j = 1, self%dim          ! for second N
        temp2 = self%positions(j) ! save value
        self%positions(j) = modulo(temp2 - 2,self%L) + 1
        self%Neighbors(temp,j + self%dim) = self%f()
        self%positions(j) = temp2 ! put back value
      end do
    else 
      do i = 1, self%L 
        self%positions(depth) = i
        call self%set_neighbors(depth + 1)
      end do
    end if    
  end subroutine set_neighbors

  ! FREES ALL USED MEMORY.
  subroutine freeall(sys, meas, rng)
    class(System), intent(in out)       :: sys 
    class(Measurements), intent(in out) :: meas
    type(rng_t), intent(in out)         :: rng 

    deallocate(sys%Neighbors , &
               sys%Status    , &   
               sys%positions , &
               sys%expdus    , &
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
               meas%tcorr    , &
               meas%bincu)
            
    call destroy_rng(rng)
  end subroutine freeall

  ! CALCULATES ERROR OF CORRELATED SAMPLE (FROM LANDAU AND BINDER)
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

  ! ONE MONTE CARLO STEP USING THE HEAT BATH ALGORITHM
  subroutine MCstep(self, rng)  
    class(System), intent(in out) :: self 
    type(rng_t), intent(in out)   :: rng 
    integer                       :: site, U, i 

    ! For one MC step (Heat bath)
    do i = 1, self%Total
      site = 1 + int(next_rng(rng) * self%Total)
      U    = sum(self%Status(self%Neighbors(site,:))) * self%status(site)
      if (next_rng(rng) .lt. self%expdus(U)) self%Status(site) = -self%Status(site)
    end do
  end subroutine MCstep

  ! just in case want to see configurations...
  ! subroutine save_config(self, output) 
  !   class(System), intent(in)    :: self 
  !   character(len=*), intent(in) :: output
  !   integer                      :: io, i, j
  !
  !   open(newunit=io,file=output,action="write")
  !     do i = 0, (self%L - 1)
  !       write(io,*) (self%Status(i*self%L + j), j=1,self%L)
  !     end do
  !   close(io)
  !
  ! end subroutine save_config

  ! GETS THERMODYNAMIC VARIABLES / SPATIAL CORRELATION AND WRITES IT TO MEASUREMENTS
  subroutine Get(sys, meas, mcs_idx) 
  class(System), intent(in)            :: sys  ! doesn't change system properties, only measures them.
    class(Measurements), intent(in out):: meas
    integer, intent(in)                :: mcs_idx
    real(dp)                           :: E_temp, M_temp, temporary
    integer                            :: site, direction, r, neighb
    
    E_temp = 0.0_dp 
    M_temp = 0.0_dp 
    do site = 1, sys%Total 
      temporary = 0.0_dp
      do direction = 1, sys%dim ! sums across only half the neighbors to avoid repeated counting
        ! ENERGY
        temporary = temporary + sys%Status(sys%Neighbors(site,direction))

        ! CORRELATION
        neighb       = sys%Neighbors(site,direction)
        meas%corr(1) = meas%corr(1) + sys%Status(site) * sys%Status(neighb)
        do r = 2, (sys%L/2 - 1) ! for all distances r between sites... -> optimization: L/2-1 = 15 (exponential decay)
          neighb       = sys%Neighbors(neighb,direction)
          meas%corr(r) = meas%corr(r) + sys%Status(site) * sys%Status(neighb)
        end do
      end do
      temporary = temporary * (-sys%Status(site))
      E_temp    = E_temp + temporary
      M_temp    = M_temp + sys%Status(site) 
    end do
    
    meas%corr(0) = meas%corr(0) + sys%dim * sys%Total ! r = 0: (+/- 1)^2 = 1

    meas%tot_E = E_temp 
    meas%tot_M = abs(M_temp)

    meas%E1 = E_temp / sys%Total 
    meas%E2 = meas%E1 * meas%E1  
    meas%M1 = abs(M_temp) / sys%Total
    meas%M2 = meas%M1 * meas%M1
    meas%M4 = meas%M1**4

    meas%mc_E(mcs_idx) = meas%E1
    meas%mc_M(mcs_idx) = meas%M1
  end subroutine Get

  ! PERFORMS SIMULATION FOR A GIVEN TEMPERATURE 
  subroutine Run(sys, meas, beta, rng, idx, int_param)
    class(System), intent(in out)      :: sys 
    class(Measurements), intent(in out):: meas
    real(dp), intent(in)               :: beta 
    type(rng_t), intent(in out)        :: rng
    integer, intent(in)                :: idx, int_param  ! idx is for indexing thermodynamic variables
    real(dp)                           :: samples, mean_mt, mean_mtp, Et, Mt, Ct, Xt, bc, temp ! temporary variables
    real(dp)                           :: re_C1, re_C2, re_X1, re_X2, C1, C2, X1, X2, rE, rm   ! bootstrapping
    integer                            :: i, j, k, dutemp, io
    character(len=100)                 :: mcresults

    ! set up exponential vector for Heat Bath algorithm
    do i = 1 , (2*sys%dim + 1)
      dutemp = 2*(sys%dim - i + 1)
      temp   = real(dutemp,dp)*beta
      sys%expdus(dutemp) = exp(-temp) / (exp(temp) + exp(-temp))
    end do

    Et             = 0.0_dp 
    Mt             = 0.0_dp 
    Ct             = 0.0_dp 
    Xt             = 0.0_dp
    bc             = 0.0_dp
    meas%corr(:)   = 0.0_dp
    meas%tcorr(:)  = 0.0_dp
    ! Burn in time 
    do i = 1,meas%burnin
      call sys%MCstep(rng)
    end do

    ! Measurement time
    do i = 1, meas%N
      do j = 1,meas%lag
        call sys%MCstep(rng)
      end do
      call Get(sys, meas, i)

      Et = Et + meas%E1 
      Mt = Mt + meas%M1 
      Ct = Ct + meas%E2
      Xt = Xt + meas%M2
      bc = bc + meas%M4
    end do

    ! get thermodynamic variables and spatial correlation
    Et              = Et / real(meas%N,dp) ! <E>
    Mt              = Mt / real(meas%N,dp) ! <m> 
    Ct              = Ct / real(meas%N,dp) ! <E^2> 
    Xt              = Xt / real(meas%N,dp) ! <m^2>
    bc              = bc / real(meas%N,dp) ! <m^4>
    meas%tcorr(0)   = Xt - Mt*Mt ! C(t=0) = <m^2> - <m>^2, temporal correlation at t=0
    do i = 0, (sys%L/2-1)
      meas%corr(i)    = meas%corr(i) / (real(sys%dim * sys%Total,dp) * real(meas%N,dp)) - Mt**2
   end do
    meas%E(idx)     = Et 
    meas%M(idx)     = Mt 
    meas%C(idx)     = (Ct - Et*Et) * beta * beta * sys%Total ! error needs E / Total.
    meas%X(idx)     = (Xt - Mt*Mt) * beta * sys%Total
    ! meas%E_err(idx) = error(meas%mc_E, Et, meas%N) 
    ! meas%M_err(idx) = error(meas%mc_M, Mt, meas%N)
    meas%E_err(idx) = 0.0_dp
    meas%M_err(idx) = 0.0_dp
    meas%bincu(idx) = 1.0_dp - bc / (3.0_dp * Xt**2) ! binder cumulant

    ! bootstrapping to get errors for C and X (unnecessary if trivially parallelized)
    ! C1 = 0.0_dp 
    ! C2 = 0.0_dp 
    ! X1 = 0.0_dp 
    ! X2 = 0.0_dp
    ! do i = 1, meas%nboot 
    !   re_C1 = 0.0_dp 
    !   re_C2 = 0.0_dp 
    !   re_X1 = 0.0_dp 
    !   re_X2 = 0.0_dp 
    !   do j = 1, meas%N 
    !     rE    = meas%mc_E(int(next_rng(rng) * meas%N) + 1) ! resample energy randomly 
    !     rm    = meas%mc_M(int(next_rng(rng) * meas%N) + 1) ! resample magnetization randomly 
    !     re_C1 = re_C1 + rE 
    !     re_C2 = re_C2 + rE * rE 
    !     re_X1 = re_X1 + rm 
    !     re_X2 = re_X2 + rm * rm 
    !   end do 
    !   re_C1 = re_C1 / real(meas%N,dp)
    !   re_C2 = re_C2 / real(meas%N,dp)
    !   re_X1 = re_X1 / real(meas%N,dp)
    !   re_X2 = re_X2 / real(meas%N,dp)
    !   Ct    = sys%Total * beta * beta * (re_C2 - re_C1 * re_C1)
    !   Xt    = sys%Total * beta * (re_X2 - re_X1 * re_X1)
    !   C1    = C1 + Ct
    !   C2    = C2 + Ct * Ct
    !   X1    = X1 + Xt
    !   X2    = X2 + Xt * Xt
    ! end do 
    ! meas%C_err(idx) = sqrt((C2 / real(meas%nboot,dp)) - (C1 / real(meas%nboot,dp))**2)
    ! meas%X_err(idx) = sqrt((X2 / real(meas%nboot,dp)) - (X1 / real(meas%nboot,dp))**2)
    meas%C_err(idx) = 0.0_dp 
    meas%X_err(idx) = 0.0_dp

    ! get temporal correlation (p. 63 Newman and Barkema)
    samples       = real(meas%TNC - 1,dp)
    for_every_dif: do j = 1, (meas%TNC  - 1) ! maximum possible time difference is TNC-1 (first and last)
      mean_mt  = 0.0_dp 
      mean_mtp = 0.0_dp
      for_every_pair: do k = 1, (meas%TNC - j) 
        meas%tcorr(j) = meas%tcorr(j) + meas%mc_M(k) * meas%mc_M(k+j)
        mean_mt       = mean_mt  + meas%mc_M(k) 
        mean_mtp      = mean_mtp + meas%mc_M(k+j) 
      end do for_every_pair
      meas%tcorr(j) = meas%tcorr(j) / samples - mean_mt*mean_mtp/samples**2
      samples       = samples - 1.0_dp
    end do for_every_dif

    !call sys%save_config(output)
 end subroutine Run 
 
end module Tools


program main
  use Tools   
  implicit none 
  
  type(rng_t)         :: rng ! Ranlux++, C++ implementation.
  type(System)        :: sample
  type(Measurements)  :: vals
  character(len=100)  :: output, corr_output, tcorr_output, param
  real(dp)            :: Ti, Tf, dT, p, current_temp
  integer(i64)        :: io, seed
  integer             :: i, k, L, dim, burnin, lag, mcs, TNC, int_param, nboot, n_sizes

  ! Initialization
  call get_command_argument(1, param)
  read(param,*) int_param
  open(newunit=io, file="Ising_input.in", action="read")
    read(io,*) dim    ! dimensions of model 
    read(io,*) L      ! initial number of sites per dimension 
    read(io,*) n_sizes! total number of system sizes to use (1->32, 2->32,64, ...)
    read(io,*) p      ! probability of a spin pointing up for initial configuration
    read(io,*) Ti     ! initial temperature (in percentage of Tc)
    read(io,*) Tf     ! final temperature (in percentage of Tc)
    read(io,*) dT     ! temperature increment
    read(io,*) mcs    ! total Monte Carlo steps
    read(io,*) burnin ! mcs before registering variables 
    read(io,*) TNC    ! Total number of configurations saved for tcorr
    read(io,*) nboot  ! number of resamples ofr bootstrapping.
    read(io,*) lag    ! mcs in between each measurement 
    read(io,*) seed   ! prng seed
  close(io)

  do k = 1, n_sizes
    print*, L
    rng = init_rng(seed*k + int_param) 
    write(output,'(A,I0,A,I0,A,I0,A)') "HB_EMCX_", L, "_", dim, "-", int_param, ".dat"  ! prefix M_ for Metropolis
    call vals%init_meas(L=L,N=mcs,burnin=burnin,lag=lag,Ti=Ti,Tf=Tf,dT=dT,TNC=TNC,nboot=nboot)
    call sample%init_sys(L=L,dim=dim,p=p,rng=rng)              ! ordered => p=1 (all 1) or p=0 (all -1) , random => p=0.5
    call sample%set_neighbors(1)                               ! needs agrument to be 1 to start recursion

    ! FILE FOR THERMODYNAMIC VARIABLES AND BINDER CUMULANT
    open(newunit=io, file=output, action="write")  
      write(io,*) "# Compiler version = ", compiler_version()
      write(io,*) "# Compiler flags = ", compiler_options()
      write(io,*) "# dimension = ", sample%dim 
      write(io,*) "# L = ", sample%L
      write(io,*) "# Seed = ", seed
      write(io,*) "# temperature interval (Ti, Tf, dT) = ", Ti, Tf, dT 
      write(io,*) "# initial spin up probability =", p
      write(io,*) "# Burn in time (in MC steps units) = ", vals%burnin
      write(io,*) "# Lag time (in MC steps units) = ", vals%lag
      write(io,*) "# Number of MC steps = ", vals%N
      write(io,*) "# Temperature / E / E error / M / M error / C / C error / X / X error / BC"
    close(io)

    ! FILE FOR SPATIAL CORRELATION
    write(corr_output,'(A,I0,A,I0,A,I0,A)') "HB_CORR_", L, "_", dim, "-", int_param, ".dat" 
    open(newunit=io, file=corr_output, action="write")  
      write(io,*) "# Compiler version = ", compiler_version()
      write(io,*) "# Compiler flags = ", compiler_options()
      write(io,*) "# dimension = ", sample%dim 
      write(io,*) "# L = ", sample%L
      write(io,*) "# Seed = ", seed
      write(io,*) "# temperature interval (Ti, Tf, dT) = ", Ti, Tf, dT 
      write(io,*) "# initial spin up probability =", p
      write(io,*) "# Burn in time (in MC steps units) = ", vals%burnin
      write(io,*) "# Lag time (in MC steps units) = ", vals%lag
      write(io,*) "# Number of MC steps = ", vals%N
      write(io,*) "# Temperature / correlation at distance 0, 1, ..."
    close(io)

    ! FILE FOR TEMPORAL CORRELATION
    write(tcorr_output,'(A,I0,A,I0,A,I0,A)') "HB_TCORR_", L, "_", dim, "-", int_param, ".dat" 
    open(newunit=io, file=tcorr_output, action="write")  
      write(io,*) "# Compiler version = ", compiler_version()
      write(io,*) "# Compiler flags = ", compiler_options()
      write(io,*) "# dimension = ", sample%dim 
      write(io,*) "# L = ", sample%L
      write(io,*) "# Seed = ", seed
      write(io,*) "# temperature interval (Ti, Tf, dT) = ", Ti, Tf, dT 
      write(io,*) "# initial spin up probability =", p
      write(io,*) "# Burn in time (in MC steps units) = ", vals%burnin
      write(io,*) "# Lag time (in MC steps units) = ", vals%lag
      write(io,*) "# Number of MC steps = ", vals%N
      write(io,*) "# TNC = ", vals%TNC
      write(io,*) "# Temperature / correlation at time 0, 1, ..."
    close(io)

    current_temp = vals%Ti
    do i = 1, vals%length 
      print*, i, vals%length
      ! prints progress on a single line
      ! write(*,"(A, F5.1, A)", advance="no") "Completion: ", real(i-1,dp)*100_dp/vals%length, "%"
      ! write(*, "(A)", advance="no") CHAR(13) ! Move the cursor to the start of the line. (carriage return)
      call run(sample, vals, 1.0_dp/current_temp, rng, i, int_param)

      open(newunit=io, file=output, action="write", position="append")
        write(io,*) current_temp, vals%E(i), vals%E_err(i), vals%M(i), vals%M_err(i),&
          vals%C(i), vals%C_err(i), vals%X(i), vals%X_err(i), vals%bincu(i)
      close(io)

      open(newunit=io, file=corr_output, action="write", position="append")
        write(io,*) current_temp, vals%corr(:)
      close(io)

      open(newunit=io, file=tcorr_output, action="write", position="append")
        ! sample%cluster_size == 0.0_dp to avoid any more modifications to show.py (no clusters in Heat bath)
        write(io,*) current_temp, 0.0_dp, vals%tcorr(:)  ! Not divided by C(t=0)
      close(io)
      current_temp = current_temp + vals%dT
    end do
    call freeall(sample,vals,rng)
    L = 2*L
  end do
end program main
