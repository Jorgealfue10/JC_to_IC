subroutine get_ic(rg,rp,theta,r12,r13,r23)
    implicit none
    real(kind=8),parameter :: pi=3.141592653589793
    real(kind=8),intent(in) :: rg,rp,theta
    real(kind=8) :: thrad
    real(kind=8),intent(out) :: r12,r13,r23

    thrad=theta*pi/180
    r12=((rp*0.5)*(rp*0.5)+rg*rg-2*0.5*rg*rp*cos(thrad))**0.5
    r13=((rp*0.5)*(rp*0.5)+rg*rg-2*0.5*rg*rp*cos(pi-thrad))**0.5
    r23=rp
end