!   ======= Lambda_Prepare =========
! Purpose:
!          Create a new format of covariance matrix
!     for lambda ambiguity fixing.
! Inputs:
!      Q0                   the old format covariance matrix
!      amb0               the old format ambiguities
!      Q                     the new format covariance matrix
!      amb                 the new format ambiguities

subroutine LAMBDA_Prepare(Q0, amb0, N0, Q, amb, N, PRN)
use MOD_VAR
implicit none
    integer :: N0, N
    real(8) :: Q0(N0, N0), amb0(N0)
    real(8) :: Q(N0, N0), amb(N0)
    integer :: PRN(N0)
!    real(8), allocatable :: Q(:,:)), amb(:)
    ! Local variables
    integer :: i, PRN_S

    Q=Q0; amb=amb0

!    N=0
!    do i=1,N0
!        if (amb0(i)==0.d0) then
!            N=N+1
!        end if
!    end do

!    if (allocated(Q)) deallocate(Q)
!    if (allocated(Q)) deallocate(amb)
!    allocate(Q(N,N))
!    allocate(amb(N))
    
    N=N0
    do i=N0,  1, -1
        PRN_S=PRN(i)
        if ((amb(i)==0.d0) .or. (GloParaNum>0 .and. PRN_S>GNum .and. PRN_S<=GNum+RNum)) then  ! exclude GLONASS satellite in Lambda
            PRN(i)=0           
            Q(1:N0, i:N0-1)=Q(1:N0, i+1:N0)
            Q(1:N0, N0) = 0.d0
            Q(i:N0-1, 1:N0)=Q(i+1:N0, 1:N0)
            Q(N0, 1:N0) = 0.d0
            amb(i:N0-1)=amb(i+1:N0)
            amb(N0)=0.d0
            PRN(i:N0-1)=PRN(i+1:N0)
            PRN(N0)=0
            N=N-1
        end if
    end do
    

    return
    
end subroutine