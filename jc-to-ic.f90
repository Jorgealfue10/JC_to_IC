subroutine get_ic(rg,rp,theta,r12,r13,r23,ditype)
    implicit none
    real(kind=8),parameter :: pi=3.141592653589793
    real(kind=8),intent(in) :: rg,rp,theta
    real(kind=8) :: thrad
    real(kind=8),intent(out) :: r12,r13,r23
    real(kind=8) :: mp,mh,rpx,rhx
    character(len=1) :: ditype

    if (ditype.eq.'y') then
        thrad=theta*pi/180
        r12=((rp*0.5)*(rp*0.5)+rg*rg-2*0.5*rg*rp*cos(thrad))**0.5
        r13=((rp*0.5)*(rp*0.5)+rg*rg-2*0.5*rg*rp*cos(pi-thrad))**0.5
        r23=rp
    else
        mp=30.973762
        mh=1.007
        write(*,*) 'Mass atom 1 (P): ',mp
        write(*,*) 'Mass atom 2 (H): ',mh
        rpx=(rp*mh)/(mh+mp)
        rhx=(rp*mp)/(mh+mp)
        write(*,*) 'R P->X: ', rpx
        write(*,*) 'R H->X: ', rhx

        thrad=theta*pi/180

        r12=rp
        r13=(rpx**2+rg**2-2*rpx*rg*cos(thrad))**0.5
        r23=(rhx**2+rg**2-2*rhx*rg*cos(180-thrad))**0.5
    endif

end
