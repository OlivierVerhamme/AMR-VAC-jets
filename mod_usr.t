module mod_usr
  use mod_hd
  implicit none
  double precision :: rhoj, eta, vj, Re,Ma,RR,rho_0,beta,g_cor

contains


  subroutine usr_init()

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_internal_bc   => no_vel
    !usr_special_bc    => specialbound_usr
    usr_refine_grid   => specialrefine_grid
    usr_gravity       => gravity

    ! Physical scaling factors

    unit_length        = 50*1.0d3*3.0857d18 ! 50 kilo parsec in cm
    unit_numberdensity = 1.0d-2 ! cm^-3
    unit_velocity = unit_length/(50d6)

    call hd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr

    !Ma = 0.1d0
    !Re = 100.0d0
    !vc_mu=1.0/Re

    ! Physical parameters
    hd_gamma=5.0d0/3.0d0
    RR = 0.05d0
    beta=0.5d0
    rho_0=1.0d0
    rhoj=0.01d0
    !g_cor = 2d5!6.67408d-11
    g_cor = 2.74d4*unit_length/(unit_velocity**2)
    if(iprob==1)then
        eta=10.d0
    endif
    if(iprob==2)then
        eta=2.d0
    endif
    if(iprob==3)then
        eta=0.1d0
    endif

    ! Jet velocity
    vj=70.d0


    print *,'unit_time: ', unit_time





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

    where ((rad(ixO^S)<RR) .and. (dabs(costhe(ixO^S))<0.131d0))
      w(ixO^S,mom(1))=rhoj*vj*costhe(ixO^S)
      w(ixO^S,mom(2))=rhoj*vj*sinthe(ixO^S)
      w(ixO^S,rho_)  =rhoj
      !w(ix^S,p_)    =p0
      w(ixO^S,e_)=one/(hd_gamma-one)+0.5d0*rhoj*vj**2.0d0
    end where

    !where((rad(ixO^S)<RR))
    !  w(ixO^S,rho_)  =rhoj*1000
    !  w(ixO^S,mom(1))=0.0d0
    !  w(ixO^S,mom(2))=0.0d0
    !  w(ixO^S,e_) = 1.0d0/(hd_gamma-1.0d0)

      !w(ixO^S,mom(1))=rhoj*vj*costhe(ixO^S)
      !w(ixO^S,e_)=one/(hd_gamma-one)+0.5d0*rhoj*vj**2.0d0

    !end where

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

    where ((rad(ixO^S)<RR) .and. (dabs(costhe(ixO^S))<0.131d0))
        w(ixO^S,mom(1)) = 0.0d0
        w(ixO^S,mom(2)) = 0.0d0
        w(ixO^S,rho_)  =rhoj
        !w(ixO^S,p_)    =p0
        w(ixO^S,e_)=one/(hd_gamma-one)+0.5d0*rhoj*vj**2.0d0
    end where

    where ((rad(ixO^S)<RR) .and. (dabs(costhe(ixO^S))<0.131d0).and.(qt<0.002))
        w(ixO^S,mom(1)) = rhoj*vj*costhe(ixO^S)
        w(ixO^S,mom(2)) = rhoj*vj*sinthe(ixO^S)
        w(ixO^S,rho_)  =rhoj
        !w(ixO^S,p_)    =p0
        w(ixO^S,e_)=one/(hd_gamma-one)+0.5d0*rhoj*vj**2.0d0
    end where


    !where((rad(ixO^S)<RR))
    !  w(ixO^S,rho_)  =rhoj*1000
    !  w(ixO^S,mom(1))=rho*vj
    !  w(ixO^S,mom(2))=0.0d0
    !  w(ixO^S,e_) = 1.0d0/(hd_gamma-1.0d0)

      !w(ixO^S,mom(1))=rhoj*vj*costhe(ixO^S)
      !w(ixO^S,e_)=one/(hd_gamma-one)+0.5d0*rhoj*vj**2.0d0

    !end where
  end subroutine no_vel

  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
     integer, intent(in)             :: ixI^L, ixO^L
     double precision, intent(in)    :: x(ixI^S,1:ndim)
     double precision, intent(in)    :: wCT(ixI^S,1:nw)
     double precision, intent(out)   :: gravity_field(ixI^S,ndim)
     double precision                :: ggrid(ixI^S),costhe(ixI^S),sinthe(ixI^S),rad(ixI^S)
     rad(ixO^S)=dsqrt(x(ixO^S,1)**2+x(ixO^S,2)**2)
     gravity_field=0.d0
     costhe(ixO^S)=x(ixO^S,1)/rad(ixO^S)
     sinthe(ixO^S)=x(ixO^S,2)/rad(ixO^S)
     call getggrav(ggrid,ixI^L,ixO^L,x)
     gravity_field(ixO^S,1)=-ggrid(ixO^S)*costhe(ixO^S)
     gravity_field(ixO^S,2)=-ggrid(ixO^S)*sinthe(ixO^S)
     !gravity_field(ixO^S,2)=g_cor
   end subroutine gravity
   subroutine getggrav(ggrid,ixI^L,ixO^L,x)
     integer, intent(in)             :: ixI^L, ixO^L
     double precision, intent(in)    :: x(ixI^S,1:ndim)
     double precision, intent(out)   :: ggrid(ixI^S)
     ggrid(ixO^S)=g_cor*(1/(x(ixO^S,1)**2+x(ixO^S,2)**2))
   end subroutine

end module mod_usr
