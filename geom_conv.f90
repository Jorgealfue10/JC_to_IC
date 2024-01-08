program jcconv 
    implicit none
    integer :: i,j,npoints
    real(kind=8) :: rg,rp,theta
    real(kind=8) :: r12,r13,r23
    real(kind=8) :: e(6)
    character(len=1) :: diattype
    integer :: io
    
    open(15,file="ci.dat",action="read" )
    open(16,file="icci.dat",action="write")

    write(*,*) "-----------------------------------------------"
    write(*,*) "For the diatomic molecule."
    write(*,fmt='(A19)',advance='no') "Homonuclear (y/n): "
    read(*,*) diattype
    write(*,*) "-----------------------------------------------"
    write(*,*)

    do
        read(15,*,iostat=io) rg,rp,theta,e(:)
        if (io.ne.0) exit
        call get_ic(rg,rp,theta,r12,r13,r23,diattype)
        write(16,*) rg,rp,theta,r12,r13,r23,e(:)
    enddo
    close(15)
    close(16)

end program jcconv