program errors 
    implicit none
    integer :: i,j,npoints
    real(kind=8) :: r12,r13,r23,e,de,en,aux,xp
    real(kind=8) :: emin,eminfit
    real(kind=8) :: der(3),emax
    character(len=1) :: diattype
    integer :: io

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

    open(15,file="suf.dat",action="read")
    open(16,file="new_inp.dat",action="write")
    open(17,file="discard.dat",action="write")

    ! write(*,*) "-----------------------------------------------"
    ! write(*,*) "For the diatomic molecule."
    ! write(*,fmt='(A19)',advance='no') "Homonuclear (y/n): "
    ! read(*,*) diattype
    ! write(*,*)
    ! write(*,*) "-----------------------------------------------"

    rewind(15)

    emax=20000
    i=0
    do 
        read(15,*,iostat=io) r12,r23,r13,e!,aux,aux,aux
        if (i.eq.0)then
            emin=e
            call fit3d(r12,r13,r23,en,der)
            eminfit=en
        endif
        e=e-emin
        if (io.ne.0) exit
        ! call get_ic(rg,rp,theta,r12,r13,r23,diattype)
        call fit3d(r12,r13,r23,en,der)
        en=en-eminfit
        !e=eh+e
        if (abs(en-e)*219474.ge.emax) then
            write(17,*) r12,r23,r13,e,en,abs(en-e)*219474
            xp=0.01d0
            write(16,*) r12,r23,r13,e+emin,xp
            print*,i
        else
            if (r12.gt.10.0d0.or.r13.gt.10.0d0.or.r23.gt.10.0d0) xp=0.5d0
            if ((r12.gt.10.0d0.or.r13.gt.10.0d0).and.r23.lt.0.80d0) xp=0.03d0
            if ((r23.gt.10.0d0.or.r13.gt.10.0d0).and.r12.lt.1.70d0) xp=0.03d0
            if ((r12.lt.5.0d0.and.r12.gt.1.8d0).and.r23.lt.3.30d0) xp=1.0d0
            if ((r13.lt.5.0d0.and.r13.gt.1.8d0).and.r12.lt.3.30d0) xp=1.0d0
            if (abs(en-e)*219474.le.500.and.((r12.lt.5.0d0.and.r23.lt.5.0d0).or.(r12.lt.5.0d0.and.r23.lt.5.0d0))) xp=4.0d0
            write(16,*) r12,r23,r13,e+emin,xp
            print*,"A",i,emin,eminfit
        endif
        i=i+1
    enddo
    close(15)
    close(16)
    close(17)

end program errors