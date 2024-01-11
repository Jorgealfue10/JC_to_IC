subroutine get_ic(rg,rp,theta,r12,r13,r23,ditype)
    implicit none
    real(kind=8),intent(in) :: rg,rp,theta
    real(kind=8) :: thrad,pi
    real(kind=8),intent(out) :: r12,r13,r23
    real(kind=8) :: mp,mh,rpx,rhx
    character(len=1) :: ditype

    pi=dacos(-1.d0)

    if (ditype.eq.'y') then
        thrad=theta*pi/180.d0
        r12=((rp*0.5d0)*(rp*0.5d0)+rg*rg-2.d0*0.5d0*rg*rp*dcos(thrad))**0.5d0
        r13=((rp*0.5d0)*(rp*0.5d0)+rg*rg-2.d0*0.5d0*rg*rp*dcos(pi-thrad))**0.5d0
        r23=rp
    else
        mp=30.973762d0
        mh=1.007d0
        ! write(*,*) 'Mass atom 1 (P): ',mp
        ! write(*,*) 'Mass atom 2 (H): ',mh
        rpx=(rp*mh)/(mh+mp)
        rhx=(rp*mp)/(mh+mp)
        write(*,*) 'R P->X: ',rg,rp,rpx
        write(*,*) 'R H->X: ', rhx

        thrad=theta*pi/180.d0

        r12=rp
        r13=(rpx**2+rg**2-2.d0*rpx*rg*dcos(pi-thrad))**0.5d0
        r23=(rhx**2+rg**2-2.d0*rhx*rg*dcos(thrad))**0.5d0
    endif

end
