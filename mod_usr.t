module mod_usr
  use mod_hd
  implicit none
  double precision :: rhoj, eta, vj, Re,Ma,RR,rho_0,beta,k,mu,massproton
contains
  subroutine usr_init()
    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_internal_bc   => no_vel
    !usr_special_bc    => specialbound_usr
    usr_refine_grid   => specialrefine_grid
    usr_var_for_errest => laplacian_for_errest
    usr_gravity      => gravity
    usr_aux_output     => specialvar_output
    usr_add_aux_names  => specialvarnames_output

    ! Physical scaling factors

    unit_length = 50*1.0d3*3.0857d18 !50 kilo parsec in cm
    unit_numberdensity = 1.0d-2 ! cm^-3
    unit_velocity = unit_length/(50d6*31556926) ! now time in seconds!

    call hd_activate()
  end subroutine usr_init
  subroutine initglobaldata_usr
    !Ma = 0.1d0
    !Re = 100.0d0
    !vc_mu=1.0/Re
    ! Physical parameters
    hd_gamma=5.0d0/3.0d0
    RR = 8*3.26d-3
    beta=0.5d0
    rho_0=1.0d0
    rhoj=0.01d0
    k = 1.38d-23
    mu= 0.6
    massproton = 1.67*10d-24
    !g_cor = 2d5!6.67408d-11
    !g_cor = 2.74d4*unit_length/(unit_velocity**2)

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

    !Print *,'unit_time: ', unit_time

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


  subroutine laplacian_for_errest(ixI^L,ixO^L,iflag,w,x,var)
      integer, intent(in)           :: ixI^L,ixO^L,iflag
      double precision, intent(in)  :: x(ixI^S,1:ndim)
      double precision, intent(in)  :: w(ixI^S,1:nw)
      double precision, intent(out) :: var(ixI^S)
      double precision :: wlocal(ixI^S,1:nw)
      integer :: idir
      double precision :: laplacianpx(ixI^S), laplacianrhox(ixI^S), laplacianpy(ixI^S), laplacianrhoy(ixI^S)
      double precision :: laplacianp(ixI^S), laplacianrho(ixI^S)
      double precision :: lapl(ixI^S, 1:2)

      wlocal(ixO^S,1:nw) = w(ixO^S,1:nw)
      call hd_to_primitive(ixI^L,ixO^L,wlocal,x)

      if (iflag == nw+1) then
          idir = 1 ! in the x-direction
          call gradient(wlocal(ixO^S, p_), ixI^L, ixO^L, idir, laplacianpx(ixO^S))
          call gradient(laplacianpx(ixO^S), ixI^L, ixO^L, idir, laplacianpx(ixO^S))
          idir = 2 ! in the y-direction
          call gradient(wlocal(ixO^S, p_), ixI^L, ixO^L, idir, laplacianpy(ixO^S))
          call gradient(laplacianpy(ixO^S), ixI^L, ixO^L, idir, laplacianpy(ixO^S))
          var(ixI^S)=laplacianpx(ixO^S)+laplacianpy(ixO^S)
      endif
      if (iflag == nw+2) then
          idir = 1 ! in the x-direction
          call gradient(wlocal(ixO^S, rho_), ixI^L, ixO^L, idir, laplacianrhox(ixO^S))
          call gradient(laplacianrhox(ixO^S), ixI^L, ixO^L, idir, laplacianrhox(ixO^S))
          idir = 2 ! in the y-direction
          call gradient(wlocal(ixO^S, rho_), ixI^L, ixO^L, idir, laplacianrhoy(ixO^S))
          call gradient(laplacianrhoy(ixO^S), ixI^L, ixO^L, idir, laplacianrhoy(ixO^S))
          var(ixI^S)=laplacianrhox(ixO^S)+laplacianrhoy(ixO^S)
      endif
      call hd_to_conserved(ixI^L,ixO^L,wlocal,x)
    end subroutine laplacian_for_errest


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
        w(ixO^S,mom(1)) = 0.0d0
        w(ixO^S,mom(2)) = 0.0d0
        w(ixO^S,rho_)  =rhoj
        !w(ixO^S,p_)    =p0
        w(ixO^S,e_)=one/(hd_gamma-one)+0.5d0*rhoj*vj**2.0d0
    end where
    where ((rad(ixO^S)<RR) .and. (dabs(costhe(ixO^S))<0.259d0).and.(qt<1.0d0))
        w(ixO^S,mom(1)) = rhoj*vj*costhe(ixO^S)
        w(ixO^S,mom(2)) = rhoj*vj*sinthe(ixO^S)
        w(ixO^S,rho_)  =rhoj
        !w(ixO^S,p_)    =p0
        w(ixO^S,e_)=one/(hd_gamma-one)+0.5d0*rhoj*vj**2.0d0
    end where
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
    call getggrav(ggrid,ixI^L,ixO^L,x,wCT)
    gravity_field(ixO^S,1)=-ggrid(ixO^S)*costhe(ixO^S)
    gravity_field(ixO^S,2)=-ggrid(ixO^S)*sinthe(ixO^S)
    !gravity_field(ixO^S,2)=g_cor
  end subroutine gravity

  subroutine getggrav(ggrid,ixI^L,ixO^L,x,w)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision :: wlocal(1:nw)
    double precision :: pth(ixI^S)
    double precision :: t
    double precision :: g_cor
    wlocal(1:nw)=w(ixImax1 /2, ixImax2 /2,1:nw)

    !call hd_to_primitive(ixI^L,ixIInt^L,w,x)
    call hd_get_pthermal(wlocal,x,ixI^L,ixO^L,pth)
    !t = pth(ixImax1/2, ixImax2/2)/w(ixImax1/2,ixImax2/2,rho_)
    t = 1/rho_0


    g_cor = -3.0d0*beta*t/2.0d0
    ggrid(ixO^S)=g_cor*((dsqrt(x(ixO^S,1)**2+x(ixO^S,2)**2)/2.0d0)/(1+(x(ixO^S,1)**2+x(ixO^S,2)**2)/4.0d0))
    !g_cor = -3.0d0*beta*t*unit_temperature/(mu*massproton*2.0d0*unit_length)
    !ggrid(ixO^S)=g_cor*((dsqrt((x(ixO^S,1)*unit_length)**2+(x(ixO^S,2)*unit_length)**2)/(2.0d0*unit_length))/(1+((x(ixO^S,1)*unit_length)**2+(x(ixO^S,2)*unit_length)**2)/((2.0d0*unit_length)**2)))*(unit_time**2/unit_length)
  end subroutine getggrav

  subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)

    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: pth(ixI^S),gradrho(ixI^S),drho(ixI^S),vrot(ixI^S),tmp(ixI^S)
    double precision                   :: kk,grhomax,kk1
    double precision :: wlocal(ixI^S,1:nw)
    integer                            :: idims

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    ! output temperature
    call hd_get_pthermal(wlocal,x,ixI^L,ixO^L,pth)
    w(ixO^S,nw+1)=pth(ixO^S)/w(ixO^S,rho_)

    ! output Mach number V/c_s
    w(ixO^S,nw+2)=dsqrt(wlocal(ixO^S,mom(1))**2+wlocal(ixO^S,mom(2))**2) /dsqrt(hd_gamma*pth(ixO^S)*w(ixO^S,rho_))

    ! output vorticity
    vrot(ixO^S)=zero
    idims=1
    tmp(ixI^S)=wlocal(ixI^S,mom(2))/wlocal(ixI^S,rho_)
    call gradient(tmp,ixI^L,ixO^L,idims,drho)
    vrot(ixO^S)=vrot(ixO^S)+drho(ixO^S)
    idims=2
    tmp(ixI^S)=wlocal(ixI^S,mom(1))/wlocal(ixI^S,rho_)
    call gradient(tmp,ixI^L,ixO^L,idims,drho)
    vrot(ixO^S)=vrot(ixO^S)-drho(ixO^S)
    w(ixO^S,nw+3)=vrot(ixO^S)

    ! output schlieren plot
    gradrho(ixO^S)=zero
    do idims=1,ndim
       select case(typegrad)
          case("central")
           call gradient(wlocal(ixI^S,rho_),ixI^L,ixO^L,idims,drho)
          case("limited")
           call gradientS(wlocal(ixI^S,rho_),ixI^L,ixO^L,idims,drho)
       end select
       gradrho(ixO^S)=gradrho(ixO^S)+drho(ixO^S)**2.0d0
    enddo
    gradrho(ixO^S)=dsqrt(gradrho(ixO^S))
    kk=5.0d0
    kk1=0.0001d0
    ! need the global maximum here, otherwise see the block structure reflected...
    !grhomax=max(1000.0d0,maxval(gradrho(ixO^S)))
    grhomax=100.0d0
    w(ixO^S,nw+4)=dexp(-kk*(gradrho(ixO^S)/grhomax-kk1))

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
    ! newly added variables need to be concatenated with the w_names/primnames string
    character(len=*) :: varnames

    varnames='Te Mach omega Schlier'

  end subroutine specialvarnames_output

end module mod_usr
