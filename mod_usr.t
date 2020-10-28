module mod_usr
  use mod_hd
  implicit none
  double precision :: rhoj, eta, vj, p0, RR, Re, Ma

contains

  subroutine usr_init()

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_special_bc    => specialbound_usr
    usr_refine_grid   => specialrefine_grid
    usr_internal_bc    => no_vel
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

    RR2=RR**2
    rad2(ix^S)=(x(ix^S,1))**2+(x(ix^S,2)-1)**2
    rad(ix^S)=dsqrt((x(ix^S,1))**2+(x(ix^S,2)-1)**2)
    costhe(ix^S)=x(ix^S,1)/rad(ix^S)
    cos2theta(ix^S)=2.0d0*costhe(ix^S)**2-1.0d0
    !where( dabs(x(ix^S,1))<0.05d0.and.x(ix^S,2)<0.05d0)
    !   w(ix^S,rho_)=rhoj
    !   w(ix^S,mom(1))=zero
    !   w(ix^S,mom(2))=rhoj*vj
    !   w(ix^S,e_)=one/(hd_gamma-one)+0.5d0*rhoj*vj**2.0d0

    where(rad(ix^S)<RR)
      w(ix^S,mom(1))=0.0d0
      w(ix^S,mom(2))=rhoj*vj
      w(ix^S,rho_)  =rhoj
      w(ix^S,p_)    =p0
      w(ix^S,e_)=one/(hd_gamma-one)+0.5d0*rhoj*vj**2.0d0
    else where
       w(ix^S,rho_) = rhoj/eta
       w(ix^S,e_) = one/(hd_gamma-one)
       w(ix^S,mom(1)) = zero
       w(ix^S,mom(2)) = zero
    end where


    !w(ix^S,mom(1))=1.0d0+RR2/rad2(ix^S)-2.0d0*x(ix^S,1)**2*RR2/rad2(ix^S)**2
    !w(ix^S,mom(2))=-2.0d0*x(ix^S,1)*x(ix^S,2)*RR2/rad2(ix^S)**2
    !w(ix^S,p_)=0.5d0*(2.0d0*cos2theta(ix^S)*RR2/rad2(ix^S)-RR2**2/rad2(ix^S)**2)+p0


    !endwhere
    !where(  dabs(x(ix^S,1))<0.05d0.and. 0.0d0>x(ix^S,2)>-0.05d0)
    !   w(ix^S,rho_)=rhoj
    !   w(ix^S,mom(1))=zero
    !   w(ix^S,mom(2))=-rhoj*vj
    !   w(ix^S,e_)=one/(hd_gamma-one)+0.5d0*rhoj*vj**2.0d0
    !end where

  end subroutine initonegrid_usr

  subroutine specialbound_usr(qt,ixG^L,ixO^L,iB,w,x)

    ! special boundary types, user defined

    integer, intent(in) :: ixG^L, ixO^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    integer :: ixOInt^L, ix2
    double precision :: rad(ixG^S),costhe(ixG^S),cos2theta(ixG^S),rad2(ixG^S),RR2

   select case(iB)
     ! implementation of special bottom boundary
     case(3)
      ! extrapolate all primitives continuously, and in jet region: fix jet
      !
      ! first switch internal zone above boundary zone to primitive variables
      ixOInt^L=ixO^L;
      ixOIntmin2=ixOmax2+1
      ixOIntmax2=ixOmax2+1
      RR2=RR**2
      rad2(ixO^S)=(x(ixO^S,1))**2+(x(ixO^S,2))**2
      rad(ixO^S)=dsqrt((x(ixO^S,1))**2+(x(ixO^S,2))**2)
      costhe(ixO^S)=x(ixO^S,1)/rad(ixO^S)
      cos2theta(ixO^S)=2.0d0*costhe(ixO^S)**2-1.0d0
      call hd_to_primitive(ixG^L,ixOInt^L,w,x)
      ! extrapolate primitives, first everywhere on boundary
      do ix2 = ixOmin2,ixOmax2
         w(ixOmin1:ixOmax1,ix2,rho_)  = w(ixOmin1:ixOmax1,ixOmax2+1,rho_)
         w(ixOmin1:ixOmax1,ix2,mom(1))= w(ixOmin1:ixOmax1,ixOmax2+1,mom(1))
         w(ixOmin1:ixOmax1,ix2,mom(2))= w(ixOmin1:ixOmax1,ixOmax2+1,mom(2))
         w(ixOmin1:ixOmax1,ix2,e_)    = w(ixOmin1:ixOmax1,ixOmax2+1,e_)
      enddo
      ! in jet zone: fix all primitives to the jet values
      where(dabs(x(ixO^S,1))<0.05d0)
         w(ixO^S,rho_)=rhoj
         w(ixO^S,mom(1))=zero
         w(ixO^S,mom(2))=vj
         w(ixO^S,e_)=one
      endwhere
      ! switch to conservative variables in internal zone
      call hd_to_conserved(ixG^L,ixOInt^L,w,x)
      ! switch to conservative variables in ghost cells
      call hd_to_conserved(ixG^L,ixO^L,w,x)
    case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine specialbound_usr

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
    R(ix^S)=dsqrt((x(ix^S,1))**2+(x(ix^S,2)-1)**2)
    if (any(R(ix^S) <= 2.0d0*RR)) refine=1

    ! always refine the jet inlet zone
    if (minval(dabs(x(ix^S,1))) < 0.1.and.minval(dabs(x(ix^S,2))) < 0.1) refine=1

  end subroutine specialrefine_grid

  subroutine no_vel(level,qt,ixI^L,ixO^L,w,x)
    integer, intent(in) :: ixI^L,ixO^L,level
    double precision, intent(in) :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision :: rad2(ixI^S)
    rad2(ixO^S)=x(ixO^S,1)**2+(x(ixO^S,2)-1)**2
    where (rad2(ixO^S)<RR**2)
        w(ixO^S,mom(1)) = 0.d0
        w(ixO^S,mom(2)) = 0.d0
    end where
  end subroutine no_vel

  subroutine set_internal_cylinder(ixI^L,ixO^L,qt,w,x)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: qt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision :: rad(ixI^S)

    rad(ixO^S)=dsqrt(x(ixO^S,1)**2+(x(ixO^S,2)-1)**2)
    where (rad(ixO^S)<RR)
        w(ixO^S,mom(1)) = 0.d0
        w(ixO^S,mom(2)) = 0.d0
        w(ixO^S,p_)     = p0/(hd_gamma-1.0d0)
        w(ixO^S,rho_)   = 1.d0
    end where
  end subroutine set_internal_cylinder
end module mod_usr
