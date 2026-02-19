!------------------------------------------------------------------------------------------!
! Simulates 2D clock model and XY model with hamiltonian:                                  !
!               H = -J \sum_{<i,j>} \vec{S}_i \cdot \vec{S}_j                              !
!                                                                                          !
! To modify parameters, change "clock_input.in".                                           !
!------------------------------------------------------------------------------------------!

! COMPILE WITH:
! gfortran -O3 -flto -c rkiss05.f90 -o fastkiss.o
! gfortran -O3 -march=native -flto -c Wolff.f90 -o Wolff.o 
! gfortran -O3 -flto Wolff.o fastkiss.o -o exec          

! RUN WITH: 
! ./exec 1 

module Tools
  use , intrinsic :: iso_fortran_env, only: i64=>int64, dp=>real64, compiler_version, compiler_options
  implicit none

  real(dp), parameter :: twopi = 2.0_dp*acos(-1.0_dp)

  ! INFORMATION ABOUT MAGNETIC MATERIAL
  type :: System
    integer                :: dim                  ! dimensions of sample. 
    integer                :: q                    ! number of possible spin orientations (0, 1, ..., q-1). 
    integer                :: L                    ! Number of sites per dimension.
    integer                :: Total                ! L^dim. Serves as one Monte Carlo step (MCS).
    integer, allocatable   :: Neighbors(:,:)       ! map of neighbors with PBC.
    integer, allocatable   :: Status(:)            ! Stores integer state of each spin (0 to q-1)
    real(dp), allocatable  :: cos_table(:)         ! Lookup table for cos(k*2pi/q)
    real(dp), allocatable  :: sin_table(:)         ! Lookup table for sin(k*2pi/q)
    integer, allocatable   :: positions(:)         ! just to set up neighbors.
    integer, allocatable   :: stack(:), cluster(:) ! data structures for Wolff.
    logical, allocatable   :: in_cluster(:)        ! Boolean mask to see if a spin is in cluster.
    real(dp)               :: cluster_size         ! total number of sites in an average cluster

    contains 
      procedure, pass :: init_sys
      procedure, pass :: f, set_neighbors, MCstep, save_config
  end type 

  ! INFORMATION ABOUT MEASUREMENTS
  type :: Measurements 
    integer               :: N                      ! number of measurements taken (MCS == "Monte Carlo Steps")
    integer               :: TNC                    ! Total Number of Configurations (for temporal correlation)
    integer               :: burnin                 ! Burn in value before equilibrium (in MC steps)
    integer               :: lag                    ! time between measurements (in MCS)
    real(dp)              :: Ti                     ! initial temperature 
    real(dp)              :: Tf                     ! final temperature 
    real(dp)              :: dT                     ! temperature step 
    integer               :: length                 ! total number of points to plot
    integer               :: nboot                  ! number of resamples for bootstrap method.
    real(dp), allocatable :: E(:), M(:), C(:), X(:) ! Thermodynamic variables vs temperatures
    real(dp), allocatable :: m_phi(:)               ! another order parameter for finding BKT Tc
    real(dp), allocatable :: Upsilon(:)             ! helicity modulus
    real(dp), allocatable :: mc_E(:), mc_M(:)       ! Variables for each temperature to get C and X (auxiliary)
    real(dp), allocatable :: E_err(:), M_err(:)     ! errors of E and M at each temperature
    real(dp), allocatable :: C_err(:), X_err(:)     ! errors of C and X at each temperature
    real(dp), allocatable :: corr(:)                ! spatial correlation and temporal correlation
    real(dp), allocatable :: bincu(:)               ! Binder cumulant for a given L over range of temperatures.
    real(dp), allocatable :: U_phi(:)               ! Another cumulant for a given L over range of temperatures.
    real(dp)              :: E1, E2, M2, M4, M8     ! auxiliary variables
    real(dp)              :: Mp, Mp2, Mp4, Ex, Sx2  ! more auxiliary variables
    real(dp)              :: tot_E, tot_M           ! for multiple histogram method

    contains 
      procedure, pass :: init_meas
  end type

contains 
  
  ! INITIALIZES SYSTEM VARIABLES 
  subroutine init_sys(self, L, dim, q, random_start)
    class(System), intent(in out) :: self 
    real(dp), external            :: rkiss05 
    integer, intent(in)           :: L, dim, q
    logical, intent(in)           :: random_start
    integer                       :: i
    
    self%L     = L
    self%dim   = dim
    self%q     = q
    self%Total = int(self%L ** self%dim)

    allocate(self%Neighbors(self%Total, 2*self%dim), &
      self%Status(self%Total)                  , &
      self%positions(self%dim)                 , &
      self%stack(2*self%dim*self%Total)        , &
      self%cluster(self%Total)                 , &
      self%in_cluster(self%Total)              , &
      self%cos_table(0:self%q-1)               , &
      self%sin_table(0:self%q-1))
    
    ! Populate the lookup tables
    do i = 0, self%q-1
      self%cos_table(i) = cos(real(i, dp) * twopi / real(self%q, dp))
      self%sin_table(i) = sin(real(i, dp) * twopi / real(self%q, dp))
    end do
    
    self%in_cluster(:) = .false.
    if (random_start .eqv. .false.) then ! ordered start
      self%Status(:) = 0 
    else  ! disordered start
      do i = 1, self%Total
        self%Status(i) = int(rkiss05() * self%q) 
      end do
    end if
  end subroutine init_sys

  ! INITIALIZES MEASUREMENTS VARIABLES 
  subroutine init_meas(self, L, N, burnin, lag, Ti, Tf, dT, TNC, nboot)
    class(Measurements), intent(in out) :: self 
    integer, intent(in)                 :: L, N, burnin, Lag, TNC, nboot
    real(dp), intent(in)                :: Ti, Tf, dT

    self%N      = N
    self%burnin = burnin 
    self%lag    = lag 
    self%Ti     = Ti
    self%Tf     = Tf
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
             self%bincu(self%length), &
             self%m_phi(self%length), &
             self%U_phi(self%length), & 
             self%Upsilon(self%length))
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
    integer, intent(in)           :: depth
    integer                       :: i, j, temp, temp2

    if (depth .gt. self%dim) then
      temp = self%f()
      do j = 1, self%dim          
        temp2 = self%positions(j)
        self%positions(j) = modulo(temp2,self%L) + 1
        self%Neighbors(temp,j) = self%f()
        self%positions(j) = temp2
      end do
      do j = 1, self%dim         
        temp2 = self%positions(j)
        self%positions(j) = modulo(temp2 - 2,self%L) + 1
        self%Neighbors(temp,j + self%dim) = self%f()
        self%positions(j) = temp2
      end do
    else 
      do i = 1, self%L 
        self%positions(depth) = i
        call self%set_neighbors(depth + 1)
      end do
    end if    
  end subroutine set_neighbors

  ! Frees all used memeory.
  subroutine freeall(sys, meas)
    class(System), intent(in out)       :: sys 
    class(Measurements), intent(in out) :: meas
    deallocate(sys%Neighbors , &
               sys%Status    , &
               sys%cos_table , &
               sys%sin_table , &
               sys%positions , &
               sys%stack     , &
               sys%cluster   , &
               sys%in_cluster, &
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
               meas%Upsilon)
  end subroutine freeall

  ! ONE ITERATION OF WOLFF FOR EVERY SPIN
  subroutine MCstep(self, P_add_table)
    class(System), intent(in out) :: self
    real(dp), intent(in)          :: P_add_table(0:, 0:)
    real(dp), external            :: rkiss05
    real(dp)                      :: P_add
    integer                       :: site, i, idx_cl, idx_st, neighb
    integer                       :: state_s, state_n, state_s_prime, rand_rot, k_old, k_new

    idx_cl = 1 ! index for cluster
    idx_st = 1 ! index for stack
    site   = int(rkiss05() * self%Total) + 1 ! randomly choose spin

    ! Define a random rotation for the entire cluster step (1 to q-1)
    rand_rot = int(rkiss05() * (self%q - 1)) + 1
    
    self%cluster(idx_cl)  = site  ! add it to cluster
    idx_cl                = idx_cl + 1
    self%in_cluster(site) = .true.
    self%stack(idx_st)    = site
    idx_st                = idx_st + 1
    
    do while(idx_st .gt. 1) ! while still have spins in stack...
      ! pop spin from stack
      idx_st = idx_st - 1
      site   = self%stack(idx_st)
      state_s = self%Status(site) ! Get integer state of current spin

      do i = 1, 4 ! for each of its neighbors...
        neighb = self%Neighbors(site,i)
        if (self%in_cluster(neighb)) cycle

        ! Get integer state of neighbor
        state_n = self%Status(neighb)

        ! Proposed new state for the current spin 'site' after rotation
        state_s_prime = modulo(state_s + rand_rot, self%q)
        
        ! Indices for the probability lookup table
        k_old = modulo(state_s - state_n, self%q)
        k_new = modulo(state_s_prime - state_n, self%q)
        
        ! Look up the probability directly
        P_add = P_add_table(k_old, k_new)

        if (rkiss05() .lt. P_add) then     ! add neighbor to stack and cluster
          self%cluster(idx_cl)    = neighb
          idx_cl                  = idx_cl + 1
          self%in_cluster(neighb) = .true.
          self%stack(idx_st)      = neighb
          idx_st                  = idx_st + 1
        end if
      end do
    end do 

    ! rotate all spins in cluster and reset it
    idx_cl = idx_cl - 1
    do i = 1, idx_cl
      site = self%cluster(i)
      ! Apply the same integer rotation to all spins in the cluster
      self%Status(site) = modulo(self%Status(site) + rand_rot, self%q)
      self%in_cluster(site) = .false.
    end do
    self%cluster_size = self%cluster_size + idx_cl
  end subroutine MCstep

  ! just in case want to see configurations...
  subroutine save_config(self, output) 
    class(System), intent(in)    :: self 
    character(len=*), intent(in) :: output
    integer                      :: io, i, j
 
    open(newunit=io,file=output,action="write")
      do i = 0, (self%L - 1)
        write(io,*) (self%Status(i*self%L + j), j=1,self%L)
      end do
    close(io)
  end subroutine save_config

  ! GETS THERMODYNAMIC VARIABLES / SPATIAL CORRELATION AND WRITES IT TO MEASUREMENTS
  subroutine Get(sys, meas, mcs_idx) 
    class(System), intent(in)           :: sys
    class(Measurements), intent(in out) :: meas
    integer, intent(in)                 :: mcs_idx
    real(dp)                            :: temp, Mx, My, E_temp, Exprime, Sxprime
    integer                             :: site, direction, r, neighb, state_s, state_n
    
    Mx         = 0.0_dp 
    My         = 0.0_dp
    E_temp     = 0.0_dp
    Exprime    = 0.0_dp 
    Sxprime    = 0.0_dp
    do site = 1, sys%Total
      state_s = sys%Status(site) ! Get integer state of the site
      
      Mx = Mx + sys%cos_table(state_s)
      My = My + sys%sin_table(state_s)

      ! right neighbor
      neighb = sys%Neighbors(site,1)
      state_n = sys%Status(neighb)
      temp = sys%cos_table(modulo(state_s - state_n, sys%q))
      E_temp = E_temp + temp 
      Exprime = Exprime + temp 
      temp = sys%sin_table(modulo(state_s - state_n, sys%q))
      Sxprime = Sxprime + temp

      ! up neighbor
      neighb = sys%Neighbors(site,2)
      state_n = sys%Status(neighb)
      temp = sys%cos_table(modulo(state_s - state_n, sys%q))
      E_temp = E_temp + temp 

      ! do direction = 1, 2 ! sys%dim (avoids repeated counting)
        ! neighb = sys%Neighbors(site,direction)
        ! state_n = sys%Status(neighb) ! Get integer state of neighbor

        ! temp   = sys%cos_table(modulo(state_s - state_n, sys%q))
        ! E_temp = E_temp + temp 

        ! CORRELATION
        ! meas%corr(1) = meas%corr(1) + temp
        ! do r = 2, (sys%L/2 - 1)
        !   neighb       = sys%Neighbors(neighb,direction)
        !   state_n      = sys%Status(neighb)
        !   meas%corr(r) = meas%corr(r) + sys%cos_table(modulo(state_s - state_n, sys%q))
        ! end do
      ! end do
    end do
    
    ! divide by Total for calculation of error
    meas%E1  = -E_temp / real(sys%Total,dp)
    meas%E2  = meas%E1 * meas%E1  
    meas%M2  = sqrt(Mx*Mx + My*My) / (real(sys%Total,dp))  ! actually M, mot M2
    meas%M4  = meas%M2 * meas%M2                           ! M2 
    meas%M8  = meas%M4 * meas%M4                           ! M4
    meas%Mp  = cos(sys%q * atan(My/Mx))
    meas%Mp2 = meas%Mp * meas%Mp 
    meas%Mp4 = meas%Mp2 * meas%Mp2

    meas%Ex  = Exprime / real(sys%Total,dp)
    meas%Sx2  = Sxprime * Sxprime / real(sys%Total,dp)

    meas%mc_E(mcs_idx) = meas%E1
    meas%mc_M(mcs_idx) = meas%M2
  end subroutine Get

  ! PERFORMS SIMULATION FOR A GIVEN TEMPERATURE 
  subroutine Run(sys, meas, beta, idx, int_param)
    class(System), intent(in out)      :: sys 
    class(Measurements), intent(in out):: meas
    real(dp), intent(in)               :: beta 
    integer, intent(in)                :: idx, int_param
    real(dp)                           :: Et, Mt, Ct, Xt, bc, mphi, mphi2, mphi4, energy_diff, mExp, mSx2
    integer                            :: i, j
    character(len=100)                 :: config_results
    real(dp), allocatable              :: P_add_table(:,:) 

    ! Allocate and populate the probability table for the current beta
    allocate(P_add_table(0:sys%q-1, 0:sys%q-1))
    do i = 0, sys%q-1
      do j = 0, sys%q-1
        energy_diff = sys%cos_table(i) - sys%cos_table(j)
        if (energy_diff > 0.0_dp) then
          P_add_table(i,j) = 1.0_dp - exp(-beta * energy_diff)
        else
          P_add_table(i,j) = 0.0_dp
        endif
      end do
    end do
    
    Et    = 0.0_dp
    Mt    = 0.0_dp
    Ct    = 0.0_dp
    Xt    = 0.0_dp
    bc    = 0.0_dp
    mphi  = 0.0_dp
    mphi2 = 0.0_dp
    mphi4 = 0.0_dp
    mExp  = 0.0_dp 
    mSx2  = 0.0_dp
    meas%corr(:)     = 0.0_dp
    sys%cluster_size = 0.0_dp

    ! Burn in time 
    do i = 1,meas%burnin
      do j = 1, meas%lag
        call sys%MCstep(P_add_table) 
      end do
    end do

    sys%cluster_size = 0.0_dp
    ! Measurement time
    do i = 1, meas%N
      do j = 1,meas%lag
        call sys%MCstep(P_add_table) 
      end do
      call Get(sys, meas, i)

      Et    = Et + meas%E1
      Mt    = Mt + meas%M2
      Ct    = Ct + meas%E2
      Xt    = Xt + meas%M4
      bc    = bc + meas%M8
      mphi  = mphi + meas%Mp
      mphi2 = mphi2 + meas%Mp2
      mphi4 = mphi4 + meas%Mp4
      mExp  = mExp + meas%Ex
      mSx2  = mSx2 + meas%Sx2
    end do
    deallocate(P_add_table) 

    
    ! get thermodynamic variables and spatial correlation
    Et              = Et / real(meas%N,dp); Mt = Mt / real(meas%N,dp)
    Ct              = Ct / real(meas%N,dp); Xt = Xt / real(meas%N,dp)
    bc              = bc / real(meas%N,dp); mphi = mphi / real(meas%N,dp)
    mphi2           = mphi2 / real(meas%N,dp); mphi4 = mphi4 / real(meas%N,dp)
    mSx2 = mSx2 / real(meas%N, dp)
    mExp = mExp / real(meas%N, dp)
    meas%m_phi(idx) = mphi
    ! meas%corr(0)    = 1.0_dp - Mt
    ! meas%corr(1:)   = meas%corr(1:) / (real(sys%dim * sys%Total,dp) * real(meas%N,dp)) - Mt
    meas%E(idx)     = Et 
    meas%M(idx)     = Mt 
    meas%C(idx)     = (Ct - Et*Et) * beta * beta * sys%Total
    meas%X(idx)     = (Xt - Mt*Mt) * beta * sys%Total
    meas%E_err(idx) = 0.0_dp
    meas%M_err(idx) = 0.0_dp
    meas%bincu(idx) = 1.0_dp - bc / (3.0_dp * Xt*Xt)
    meas%U_phi(idx) = 1.0_dp - mphi4 / (2.0_dp * mphi2*mphi2)
    meas%C_err(idx) = 0.0_dp 
    meas%X_err(idx) = 0.0_dp
    meas%Upsilon(idx) = mExp - mSx2 * beta

    sys%cluster_size = sys%cluster_size / (real(meas%N, dp) * real(meas%lag, dp))
 end subroutine Run 
 
end module Tools

program main
  use Tools   
  implicit none 
  
  real(dp), external  :: rkiss05
  type(System)        :: sample
  type(Measurements)  :: vals
  character(len=100)  :: output, corr_output, param
  real(dp)            :: Ti, Tf, dT, current_temp
  integer(i64)        :: io
  integer             :: i, k, L, dim, burnin, lag, mcs, TNC, int_param, nboot, n_sizes, q, seed, init_L
  logical             :: random_start, exists

  ! Initialization
  call get_command_argument(1, param)
  read(param,*) int_param
  open(newunit=io, file="clock_input.in", action="read")
    read(io,*) dim
    read(io,*) q
    read(io,*) init_L
    read(io,*) n_sizes
    read(io,*) random_start
    read(io,*) Ti
    read(io,*) Tf
    read(io,*) dT
    read(io,*) mcs
    read(io,*) burnin
    read(io,*) TNC
    read(io,*) nboot
    read(io,*) lag
    read(io,*) seed
  close(io)

  call kissinit(seed + int_param) 
  L = init_L
  do k = 1, n_sizes
    print*, "Simulating L = ", L
    write(output,'(A,I0,A,I0,A,I0,A)') "EMCX_", L, "_", q, "-", int_param, ".dat" 
    inquire(file=output, exist=exists)
    if (exists) then 
      print*, output , " exists already. Skipping."
      L = L + init_L  ! init_L, 2*init_L , ...
      cycle 
    else
      lag = L * L ! One MC step is L^2 Wolff cluster flips
      call vals%init_meas(L=L,N=mcs,burnin=burnin,lag=lag,Ti=Ti,Tf=Tf,dT=dT,TNC=TNC,nboot=nboot)
      call sample%init_sys(L=L,dim=dim,q=q,random_start=random_start) 
      call sample%set_neighbors(1)

      ! FILE FOR THERMODYNAMIC VARIABLES AND BINDER CUMULANT
      open(newunit=io, file=output, action="write")  
        write(io,*) "# Compiler version = ", compiler_version()
        write(io,*) "# Compiler flags = ", compiler_options()
        write(io,*) "# dimension = ", sample%dim 
        write(io,*) "# q = ", sample%q
        write(io,*) "# L = ", sample%L
        write(io,*) "# Seed = ", seed + int_param
        write(io,*) "# temperature interval (Ti, Tf, dT) = ", Ti, Tf, dT 
        write(io,*) "# random start =", random_start
        write(io,*) "# Burn in time (in MC steps units) = ", vals%burnin
        write(io,*) "# Lag time (in MC steps units) = ", vals%lag
        write(io,*) "# Number of MC steps = ", vals%N
        write(io,*) "# Temperature / E / E error / M / M error / C / C error / X / X error / BC / Uphi / mphi / Upsilon"
      close(io)

      ! FILE FOR SPATIAL CORRELATION
      ! write(corr_output,'(A,I0,A,I0,A,I0,A)') "CORR_", L, "_", q, "-", int_param, ".dat" 
      ! open(newunit=io, file=corr_output, action="write")  
      !   write(io,*) "# Compiler version = ", compiler_version()
      !   write(io,*) "# Compiler flags = ", compiler_options()
      !   write(io,*) "# dimension = ", sample%dim 
      !   write(io,*) "# q = ", sample%q
      !   write(io,*) "# L = ", sample%L
      !   write(io,*) "# Seed = ", seed + int_param
      !   write(io,*) "# temperature interval (Ti, Tf, dT) = ", Ti, Tf, dT 
      !   write(io,*) "# random start =", random_start
      !   write(io,*) "# Burn in time (in MC steps units) = ", vals%burnin
      !   write(io,*) "# Lag time (in MC steps units) = ", vals%lag
      !   write(io,*) "# Number of MC steps = ", vals%N
      !   write(io,*) "# Temperature / correlation at distance 0, 1, ..."
      ! close(io)
      !
      current_temp = vals%Ti
      do i = 1, vals%length 
        print*, "Progress: ", i, "/", vals%length, " (T = ", current_temp, ")"
        call run(sample, vals, 1.0_dp/current_temp, i, int_param)

        open(newunit=io, file=output, action="write", position="append")
          write(io,*) current_temp, vals%E(i), vals%E_err(i), vals%M(i), vals%M_err(i),&
            vals%C(i), vals%C_err(i), vals%X(i), vals%X_err(i), vals%bincu(i), vals%U_phi(i), vals%m_phi(i), vals%Upsilon(i)
        close(io)

        ! open(newunit=io, file=corr_output, action="write", position="append")
        !   write(io,*) current_temp, vals%corr(:)
        ! close(io)
        
        current_temp = current_temp + vals%dT
      end do
      call freeall(sample,vals)
      L = L + init_L
    end if
  end do
end program main
