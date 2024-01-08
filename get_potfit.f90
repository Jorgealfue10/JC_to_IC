program get_fit
    implicit none
    integer :: i,j,npoints
    real(kind=8) :: rg,drg,rp,theta,rgini
    real(kind=8) :: r12,r13,r23,e,de
    real(kind=8) :: der(3)
    character(len=1) :: diattype

    interface
        subroutine fit3d(r1,r2,r3,e,der)
            implicit none
            real(8) :: r1,r2,r3,e
            real(8) :: der(3)
        end subroutine
        subroutine triabb(r1,r2,r3,e123,der)
            implicit none
            real(8) :: r1,r2,r3,e123
            real(8) :: der(3)
        end subroutine
        subroutine diat12(r,e,der)
            implicit none
            real(8) :: r,e,der
        end subroutine
        subroutine diat23(r,e,der)
            implicit none
            real(8) :: r,e,der
        end subroutine
    end interface

    open(15,file="jac-param.dat",action="read")
    open(16,file="pot-jc.dat",action="write")
!   open(20,file="pot-ic.dat",action="write")
    open(17,file="pot-12.dat",action="write")
    open(18,file="pot-13.dat",action="write")
    open(19,file="pot-23.dat",action="write")
    open(20,file="pot-3body.dat",action="write")

    write(*,*) "-----------------------------------------------"
    write(*,*) "For the diatomic molecule."
    write(*,fmt='(A19)',advance='no') "Homonuclear (y/n): "
    !read(*,*) diattype
    write(*,*) "-----------------------------------------------"
    write(*,*)

    diattype='y'
    read(15,*) npoints,drg,rp,theta

    rgini=1.5
    do i=1,npoints
        rg=rgini+drg*(i-1)
        call get_ic(rg,rp,theta,r12,r13,r23,diattype)
        call diat12(r12,e,de)
        write(17,*) r12,e
        call diat12(r13,e,de)
        write(18,*) r13,e
        call diat23(r23,e,de)
        write(19,*) r23,e
        call triabb(r12,r13,r23,e,der)
        write(20,*) r12,r13,r23,rg,rp,theta,e
        call fit3d(r12,r13,r23,e,der)
        write(16,*) rg,rp,theta,e
    enddo

end program get_fit

