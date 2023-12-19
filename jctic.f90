program jctransic
        implicit none
        real(kind=8),parameter :: pi=3.141592653589793
        integer :: i,j,nl,nst
        integer :: nargs,io
        real(kind=8) :: rg,rp,theta,th2
        real(kind=8) :: r1,r2,r3,beta
        real(kind=8),allocatable,dimension(:) :: est
        character(len=3),allocatable,dimension(:) :: args
        
        nargs=command_argument_count()
        allocate(args(nargs))

        if (nargs.lt.1) then
          write(*,fmt='(A33)') "Please provide theta (000 format)"
          stop
        endif
        
        do i=1, nargs
         call get_command_argument(i,args(i))
        end do      
        
        open(15,file="pes-"//args(1)//".dat",action="read")
        open(16,file="ic-"//args(1)//".dat",action="write")
        open(17,file="jc-"//args(1)//".dat",action="write")

        print*,"pes-"//args(1)//".dat","ic-"//args(1)//".dat"

        do
          read(15,*,iostat=io)
          if (io/=0) exit
          nl=nl + 1
        end do
        rewind(15)

        print*,nl

        write(*,fmt='(A21)',advance='no') "# electronic states: "
        read(*,*) nst

        allocate(est(nst))

        write(16,*) "r1  r2  r3  e(:)"
        do i=1,nl
          read(15,*) rg,rp,theta,est(:)
          th2=theta
          theta=theta*pi/180
          r1=((rp*0.5)*(rp*0.5)+rg*rg-2*0.5*rg*rp*cos(theta))**0.5
          r2=((rp*0.5)*(rp*0.5)+rg*rg-2*0.5*rg*rp*cos(pi-theta))**0.5
          r3=rp
          write(16,*) r1,r3,r2,est(:)
!         write(*,*) r1,r3,r2,est(:)
          rp=r3
          beta=acos(-(r2*r2-r1*r1-r3*r3)/(2*r1*r3))
          rg=(r1*r1+(rp*0.5)*(rp*0.5)-2*0.5*r1*rp*cos(beta))**0.5
          theta=acos(-(r1*r1-rg*rg-(rp*0.5)*(rp*0.5))/(2*0.5*rg*rp))
          theta=theta*180/pi
          write(17,*) rg,rp,th2,est(:)
        enddo
        
        deallocate(est,args)
        close(15)
        close(16)
        close(17)
 
end program jctransic
