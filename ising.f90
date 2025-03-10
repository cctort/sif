! Metropolis algorithm to calculate <E>, <M>, in the canonical ensemble
! (fix T,N,V) with a 2D Ising model

module module

  implicit none
  public :: metro_step, metro_step_ordered, DeltaE
  integer, public, parameter :: double = 8
  integer, public :: N, L, accept
  character(len=3), public :: boundary
  real(kind=double), public, dimension(-8:8) :: w
  real(kind=double), public :: E, M

contains
    
  subroutine put_seed(seed_number)

     ! Sets the seed for the rng

     integer, intent(in) :: seed_number
     integer, dimension(:), allocatable :: seed
     integer :: i, size

     call random_seed(size=size)
     allocate(seed(size))

     do i=1, size
       seed(i) = seed_number
     end do
     call random_seed(put=seed)
     
  end subroutine put_seed

  
  subroutine metro_step(spin)

    ! N random spin flips

    integer, dimension(:,:), intent(inout) :: spin
    integer :: ispin, x, y, dE
    real :: rnd

    do ispin = 1,N
       
       ! Random coordinates x and y
       call random_number(rnd)
       x = int(L*rnd) + 1
       call random_number(rnd)
       y = int(L*rnd) + 1

       ! Energy difference due to spin flip in x,y
       dE = DeltaE(x,y,spin)

       ! Flip the spin with a probability of w(dE)
       call random_number(rnd)
       if (rnd <= w(dE)) then
          spin(x,y) = -spin(x,y)
          accept = accept + 1
          M = M + 2*spin(x,y)
          E = E + dE
       end if
    end do

  end subroutine metro_step


  subroutine metro_step_ordered(spin)
   
    ! N ordered spin flips

    integer, dimension(:,:), intent(inout) :: spin
    integer :: x, y, dE
    real :: rnd

    do x = 1,L
      do y = 1,L

        ! Energy difference due to spin flip in x,y
        dE = DeltaE(x,y,spin)

        ! Flip the spin with a probability of w(dE)
        call random_number(rnd)
        if (rnd <= w(dE)) then
            spin(x,y) = -spin(x,y)
            accept = accept + 1
            M = M + 2*spin(x,y)
            E = E + dE
        end if
      end do
    end do

  end subroutine metro_step_ordered


  function DeltaE(x,y,spin) result (DeltaE_result)
    
    ! Calculate energy difference due to spin flip

    integer, dimension(:,:), intent(inout) :: spin
    integer, intent (in) :: x,y
    integer :: DeltaE_result
    integer :: left
    integer :: right
    integer :: up
    integer :: down
 
    ! nearest neighbors' sites
    left = x - 1
    right = x + 1
    down = y - 1
    up = y + 1
    
    ! Periodic boundary conditions
    if (boundary == 'pbc') then
      left = left + merge(L, 0, x == 1)
      right = right - merge(L, 0, x == L)
      down = down + merge(L, 0, y == 1)
      up = up - merge(L, 0, y == L)
      DeltaE_result = 2*spin(x,y) * (spin(x,up) + spin(x,down) + spin(left,y) + spin(right,y))
   
    ! Open boundary conditions
    else  
      DeltaE_result = 2*spin(x,y) * &
      (merge(spin(x,up), 0, up <= L) &
      + merge(spin(x,down), 0, down >= 1) &
      + merge(spin(left,y), 0, left >= 1) &
      + merge(spin(right,y), 0, right <= L))  
   end if

 end function DeltaE

end module module