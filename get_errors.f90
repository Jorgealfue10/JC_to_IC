program errors 
    implicit none
    integer :: i,j,npoints
    real(kind=8) :: r12,r13,r23,e,de,en,aux
    real(kind=8) :: der(3)
    character(len=1) :: diattype
    integer :: io

    interface
        subroutine fit3d(r1,r2,r3,e,der)
            implicit none
            real(8) :: r1,r2,r3,e
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

    open(15,file="plot_ABC.dat",action="read")
    open(16,file="new_inp.dat",action="write")
    open(17,file="discard.dat",action="write")

    write(*,*) "-----------------------------------------------"
    write(*,*) "For the diatomic molecule."
    write(*,fmt='(A19)',advance='no') "Homonuclear (y/n): "
    read(*,*) diattype
    write(*,*) "-----------------------------------------------"
    write(*,*)

    do
        read(15,*,iostat=io) r12,r13,r23,e,aux,aux,aux
        if (io.ne.0) exit
        call get_ic(rg,rp,theta,r12,r13,r23,diattype)
        call fit3d(r12,r13,r23,en,der)
        if (abs(en-e)*219474>20) then
            write(17,*) r12,r13,r23,e
        else
            write(16,*) r12,r13,r23,e
        endif
    enddo

end program errors