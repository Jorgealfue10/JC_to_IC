program get_fit
    implicit none
    integer :: i,j,npoints
    real(kind=8) :: rg,drg,rp,theta,rgini
    real(kind=8) :: r12,r13,r23,e,der

    interface
        subroutine fit3d(r1,r2,r3,e,der)
            implicit none
            real(8) :: r1,r2,r3,e,der
        end subroutine
    end interface

    open(15,file="jac-param.dat",action="read")
    open(16,file="pot-jc.dat",action="write")

    read(15,*) npoints,drg,rp,theta

    rgini=1.0
    do i=1,npoints
        rg=rgini+drg*(i-1)
        call get_ic(rg,rp,theta,r12,r13,r23)
        call fit3d(r12,r13,r23,e,der)
        write(16,*) rg,rp,theta,e
    enddo

end program get_fit

