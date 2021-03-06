module mod_usr
  use mod_hd
  implicit none
  double precision :: rhoj, eta, vj, Re,Ma,p0,RR

contains

  subroutine usr_init()

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_internal_bc   => no_vel
   !usr_special_bc    => specialbound_usr
    usr_refine_grid   => specialrefine_grid
    usr_modify_output  => set_internal_cylinder

    call hd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr

    Ma = 0.1d0
    Re = 100.0d0
    !vc_mu=1.0/Re
    p0=1.0d0/(hd_gamma*Ma**2)
    RR = 0.1d0
    hd_gamma=1.4d0
    rhoj=hd_gamma
    if(iprob==1)then
        eta=10.d0
    endif
    if(iprob==2)then
        eta=2.d0
    endif
    if(iprob==3)then
        eta=0.1d0
    endif
    vj=100.d0

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)

    ! initialize one grid

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision :: rad(ixG^S),costhe(ixG^S),cos2theta(ixG^S),rad2(ixG^S),RR2
    integer :: idims

    w(ix^S,rho_)  =rhoj
    RR2=RR**2
    rad2(ix^S)=x(ix^S,1)**2+x(ix^S,2)**2
    rad(ix^S)=dsqrt(x(ix^S,1)**2+x(ix^S,2)**2)
    costhe(ix^S)=x(ix^S,1)/rad(ix^S)
    cos2theta(ix^S)=2.0d0*costhe(ix^S)**2-1.0d0
    w(ix^S,mom(1))=0.0d0
    w(ix^S,mom(2))=0.0d0
    w(ix^S,p_)=p0
    w(ix^S,rho_) = rhoj/eta
    w(ix^S,e_) = one/(hd_gamma-one)

    where(rad(ix^S)<RR .and. dabs(x(ix^S,1))<0.025d0 .and. x(ix^S,2)<-0.05d0)
      w(ix^S,mom(1))=0.0d0
      w(ix^S,mom(2))=-rhoj*vj
      w(ix^S,rho_)  =rhoj
      w(ix^S,p_)    =p0
      w(ix^S,e_)=one/(hd_gamma-one)+0.5d0*rhoj*vj**2.0d0
    else where(rad(ix^S)<RR .and. dabs(x(ix^S,1))<0.025d0 .and. x(ix^S,2)>0.05d0)
      w(ix^S,mom(1))=0.0d0
      w(ix^S,mom(2))=rhoj*vj
      w(ix^S,rho_)  =rhoj
      w(ix^S,p_)    =p0
      w(ix^S,e_)=one/(hd_gamma-one)+0.5d0*rhoj*vj**2.0d0
    endwhere

    call hd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine initonegrid_usr


  subroutine specialrefine_grid(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)

    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

    ! you must set consistent values for integers refine/coarsen:

    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement

    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen

    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    double precision:: R(ixG^S)

    ! always refine the jet inlet zone
    R(ix^S)=dsqrt(x(ix^S,1)**2+x(ix^S,2)**2)

    if (any(R(ix^S) <= 2.0d0*RR)) refine=1

  end subroutine specialrefine_grid

  subroutine no_vel(level,qt,ixI^L,ixO^L,w,x)
    integer, intent(in) :: ixI^L,ixO^L,level
    double precision, intent(in) :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision :: rad2(ixI^S)
    rad2(ixO^S)=x(ixO^S,1)**2+x(ixO^S,2)**2
    where (rad2(ixO^S)<RR**2)
        w(ixO^S,mom(1)) = 0.0d0
        w(ixO^S,mom(2)) = 0.0d0
        w(ixO^S,rho_)  =rhoj
        w(ixO^S,p_)    =p0
        w(ixO^S,e_)=one/(hd_gamma-one)+0.5d0*rhoj*vj**2.0d0
    else where (rad2(ixO^S)<RR**2 .and. dabs(x(ixI^S,1))<0.025d0 .and. x(ixI^S,2)<-0.05d0)
        w(ixO^S,mom(1)) = 0.0d0
        w(ixO^S,mom(2)) = -rhoj*vj
        w(ixO^S,rho_)  =rhoj
        w(ixO^S,p_)    =p0
        w(ixO^S,e_)=one/(hd_gamma-one)+0.5d0*rhoj*vj**2.0d0
    else where (rad2(ixO^S)<RR**2 .and. dabs(x(ixI^S,1))<0.025d0 .and. x(ixI^S,2)>0.05d0)
        w(ixO^S,mom(1)) = 0.0d0
        w(ixO^S,mom(2)) = rhoj*vj
        w(ixO^S,rho_)  =rhoj
        w(ixO^S,p_)    =p0
        w(ixO^S,e_)=one/(hd_gamma-one)+0.5d0*rhoj*vj**2.0d0
    end where
  end subroutine no_vel


  subroutine set_internal_cylinder(ixI^L,ixO^L,qt,w,x)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision :: rad(ixI^S)

    rad(ixO^S)=dsqrt(x(ixO^S,1)**2+(x(ixO^S,2))**2)
    where (rad(ixO^S)<RR)
        w(ixO^S,mom(1)) = 0.d0
        w(ixO^S,mom(2)) = 0.d0
        w(ixO^S,p_)     = p0
        w(ixO^S,rho_)   = rhoj
        w(ixO^S,e_)=one/(hd_gamma-one)+0.5d0*rhoj*vj**2.0d0
    else where (rad(ixO^S)<RR .and. dabs(x(ixI^S,1))<0.025d0 .and. x(ixI^S,2)<-0.05d0)
        w(ixO^S, mom(1)) = 0.0d0
        w(ixO^S,mom(2)) = -rhoj*vj
        w(ixO^S,rho_)  =rhoj
        w(ixO^S,p_)    =p0
        w(ixO^S,e_)=one/(hd_gamma-one)+0.5d0*rhoj*vj**2.0d0

    else where (rad(ixO^S)<RR .and. dabs(x(ixI^S,1))<0.025d0 .and. x(ixI^S,2)>0.05d0)
        w(ixO^S,mom(1)) = 0.0d0
        w(ixO^S,mom(2)) = rhoj*vj
        w(ixO^S,rho_)  =rhoj
        w(ixO^S,p_)    =p0
        w(ixO^S,e_)=one/(hd_gamma-one)+0.5d0*rhoj*vj**2.0d0

    end where
  end subroutine set_internal_cylinder

end module mod_usr
