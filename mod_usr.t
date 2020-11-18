module mod_usr
  use mod_hd
  implicit none
  double precision :: rhoj, eta, vj, Re,Ma,RR,rho_0,beta

contains

  subroutine usr_init()

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_internal_bc   => no_vel
   !usr_special_bc    => specialbound_usr
    usr_refine_grid   => specialrefine_grid

    call hd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr

    !Ma = 0.1d0
    !Re = 100.0d0
    !vc_mu=1.0/Re
    hd_gamma=5.0d0/3.0d0
    RR = 0.1d0
    beta=0.5d0
    rho_0=1.0d0
    rhoj=0.01d0
    if(iprob==1)then
        eta=10.d0
    endif
    if(iprob==2)then
        eta=2.d0
    endif
    if(iprob==3)then
        eta=0.1d0
    endif
    vj=50.d0

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)

    ! initialize one grid

    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: rad(ixI^S),costhe(ixI^S),sinthe(ixI^S),cos2theta(ixI^S),rad2(ixI^S),RR2
    integer :: idims

    
    RR2=RR**2
    rad2(ixO^S)=x(ixO^S,1)**2+x(ixO^S,2)**2
    rad(ixO^S)=dsqrt(x(ixO^S,1)**2+x(ixO^S,2)**2)
    costhe(ixO^S)=x(ixO^S,1)/rad(ixO^S)
    sinthe(ixO^S)=x(ixO^S,2)/rad(ixO^S)
    cos2theta(ixO^S)=2.0d0*costhe(ixO^S)**2-1.0d0
    
    !w(ix^S,rho_)  =rhoj
    w(ixO^S,rho_) = rho_0*(1.0d0+((rad(ixO^S)/RR)**2.0d0)**(-3.0d0*beta/2.0d0))
    w(ixO^S,mom(1))=0.0d0
    w(ixO^S,mom(2))=0.0d0
    !w(ix^S,p_)=p0
    !w(ix^S,rho_) = rhoj/eta
    w(ixO^S,e_) = 1.0d0/(hd_gamma-1.0d0)

    where ((rad(ixO^S)<RR) .and. (dabs(costhe(ixO^S))<0.259d0))
      w(ixO^S,mom(1))=rhoj*vj*costhe(ixO^S)
      w(ixO^S,mom(2))=rhoj*vj*sinthe(ixO^S)
      w(ixO^S,rho_)  =rhoj
      !w(ix^S,p_)    =p0
      w(ixO^S,e_)=one/(hd_gamma-one)+0.5d0*rhoj*vj**2.0d0
    end where
  end subroutine initonegrid_usr


  subroutine specialrefine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)

    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

    ! you must set consistent values for integers refine/coarsen:

    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement

    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen

    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    double precision:: R(ixI^S)

    ! always refine the jet inlet zone
    R(ixO^S)=dsqrt(x(ixO^S,1)**2+x(ixO^S,2)**2)

    if (any(R(ixO^S) <= 2.0d0*RR)) refine=1

  end subroutine specialrefine_grid

  subroutine no_vel(level,qt,ixI^L,ixO^L,w,x)
    integer, intent(in) :: ixI^L,ixO^L,level
    double precision, intent(in) :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision :: rad2(ixI^S),costhe(ixI^S),sinthe(ixI^S),rad(ixI^S)
    rad2(ixO^S)=x(ixO^S,1)**2+x(ixO^S,2)**2
    rad(ixO^S)=dsqrt(x(ixO^S,1)**2+x(ixO^S,2)**2)
    costhe(ixO^S)=x(ixO^S,1)/rad(ixO^S)
    sinthe(ixO^S)=x(ixO^S,2)/rad(ixO^S)

    where (rad(ixO^S)<RR)
        w(ixO^S,rho_) = rho_0*(1.0d0+((rad(ixO^S)/RR)**2.0d0)**(-3.0d0*beta/2.0d0))
        w(ixO^S,mom(1)) = 0.0d0
        w(ixO^S,mom(2)) = 0.0d0
        !w(ixO^S,p_)    =p0
        w(ixO^S,e_)=1.0d0/(hd_gamma-1)
    endwhere 
    where ((rad(ixO^S)<RR) .and. (dabs(costhe(ixO^S))<0.259d0))
        w(ixO^S,mom(1)) = rhoj*vj*costhe(ixO^S)
        w(ixO^S,mom(2)) = rhoj*vj*sinthe(ixO^S)
        w(ixO^S,rho_)  =rhoj
        !w(ixO^S,p_)    =p0
        w(ixO^S,e_)=one/(hd_gamma-one)+0.5d0*rhoj*vj**2.0d0
    end where
  end subroutine no_vel
end module mod_usr
