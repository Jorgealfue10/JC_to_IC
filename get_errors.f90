program errors 
    implicit none
    integer :: i,j,npoints
    real(kind=8) :: r12,r13,r23,e,de,en,aux,eh
    real(kind=8) :: der(3)
    character(len=1) :: diattype
    integer :: io

    interface
        subroutine fit3d(r1,r2,r3,e,eh,der)
            implicit none
            real(8) :: r1,r2,r3,eh,e
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
    write(*,*)
    write(*,*) "-----------------------------------------------"

    i=0
    do
        read(15,*,iostat=io) r12,r13,r23,e!,aux,aux,aux
        if (io.ne.0) exit
        ! call get_ic(rg,rp,theta,r12,r13,r23,diattype)
        call fit3d(r12,r13,r23,en,eh,der)
        !e=eh+e
        if (abs(en-e)*219474.ge.2000) then
            write(17,*) r12,r23,r13,e,en,abs(en-e)*219474
            print*,i
        else
            write(16,*) r12,r23,r13,e
            print*,"A",i
        endif
        i=i+1
    enddo
    close(15)
    close(16)
    close(17)

end program errors