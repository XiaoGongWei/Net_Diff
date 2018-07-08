! ==============  Cross  ============
!
! PURPOSE:
!      Vector product of two vectors of 3-dim.
!
! INPUTS:
!      a       3*1 vector
!      b       3*1 vector
!
! OUTPUTS:
!      c        a cross b, 3*1 vector
!
! Note: a cross b not equals to b cross a
!
! WRITTEN BY: Yize Zhang, zhyize@163.com, Tongji & SHAO
! ==========  End of Header  ========

subroutine Cross(a,b,c)
implicit none
    ! Intents
    real(8) :: a(3), b(3)
    real(8) :: c(3)
    
    c(1)=a(2)*b(3)-a(3)*b(2)
    c(2)=a(3)*b(1)-a(1)*b(3)
    c(3)=a(1)*b(2)-a(2)*b(1)
    
    return
end subroutine