! **********************************************************************
! *********** Abaqus/ STANDARD USER ELEMENT SUBROUTINE (UEL) ***********
! **********************************************************************
!           fully coupled chemo-mechanics of non-ionic hydrogels
!   the formulation uses PK-II stress based total Lagrangian framework
! with F-bar modification for fully-integrated HEX8 and QUAD4-PE element
! currently supports first-order coupled elements for 3D and 2D-PE cases
!    FUTURE TODO: add first-order elements for 2D-PS and 2D-AX cases
! **********************************************************************
!                     BIBEKANANDA DATTA (C) MAY 2024
!                 JOHNS HOPKINS UNIVERSITY, BALTIMORE, MD
! **********************************************************************
!
!                       JTYPE DEFINITION
!
!     U1                THREE-DIMENSIONAL TET4 ELEMENT
!     U2                THREE-DIMENSIONAL HEX8 ELEMENT
!     U3                PLANE STRAIN TRI3 ELEMENT
!     U4                PLANE STRAIN QUAD4 ELEMENT
!
! **********************************************************************
!
!          VOIGT NOTATION CONVENTION FOR STRESS/ STRAIN TENSORS
!
!       In this subroutine we adopted the following convention for
!       symmetric stress and strain tensor following Voigt notation
!       This is different than what is followed by Abaqus/ Standard
!
!          sigma11, sigma22, sigma33, sigma23, sigma13, sigma12
!       strain11, strain22, strain33, strain23, strain13, strain12
!
! **********************************************************************
!
!                       LIST OF MATERIAL PROPERTIES
!
!         props(1)-props(12): universal parameters and gel properties
!
!     Rgas        = props(1)        Universal gas constant
!     theta       = props(2)        Absolute temperature (K)
!     phi0        = props(3)        Initial polymer volume fraction
!     rho         = props(4)        Density of the gel
!     Gshear      = props(5)        Shear modulus
!     Kappa       = props(6)        Bulk modulus
!     lam_L       = props(7)        Locking stretch (for AB model)
!     Vp          = props(8)        Molar volume of polymer
!     mu0         = props(9)        Chemical potential of pure solvent
!     Vw          = props(10)       Molar volume of the solvent
!     chi         = props(11)       Flory-Huggins parameter
!     Dw          = props(12)       Diffusion coefficient of the solvent
!
! **********************************************************************
!
!                        LIST OF ELEMENT PROPERTIES
!
!     jprops(1)   = nInt            no of integration points in element
!     jprops(2)   = matID           constitutive model for elastomeric network
!     jprops(3)   = nPostVars       no of local (int pt) post-processing variables
!
! **********************************************************************
!
!                        POST-PROCESSED VARIABLES
!                     (follows the convention above)
!
!     uvar(1:nStress)               Cauchy stress tensor components
!     uvar(nStress+1:2*nStress)     Euler-Almansi strain tensor components
!     uvar(2*nStress+1)             Polymer volume fraction (phi)
!
! **********************************************************************
!
!               VARIABLES TO BE UPDATED WITHIN THE SUBROUTINE
!
!     RHS(i,NRHS)                   Right hand side vector
!     AMATRX(i,j)                   Stiffness matrix (NDOFEL x NDOFEL)
!     SVARS(1:NSVARS)               Element state variables.  Must be updated in this routine
!     ENERGY(1:8)                   Energy(1) Kinetic Energy
!                                   Energy(2) Elastic Strain Energy
!                                   Energy(3) Creep Dissipation
!                                   Energy(4) Plastic Dissipation
!                                   Energy(5) Viscous Dissipation
!                                   Energy(6) Artificial strain energy
!                                   Energy(7) Electrostatic energy
!                                   Energy(8) Incremental work done by loads applied to the element
!     PNEWDT                        Allows user to control ABAQUS time increments.
!                                   If PNEWDT<1 then time step is abandoned and computation is restarted with
!                                   a time increment equal to PNEWDT*DTIME
!                                   If PNEWDT>1 ABAQUS may increase the time increment by a factor PNEWDT
!
!                       VARIABLES PROVIDED FOR INFORMATION
!
!     NDOFEL                        Total # DOF for the element
!     NRHS                          Dimension variable
!     NSVARS                        Total # element state variables
!     PROPS(1:NPROPS)               User-specified properties of the element
!     NPROPS                        No. properties
!     JPROPS(1:NJPROPS)             Integer valued user specified properties for the element
!     NJPROPS                       No. integer valued properties
!     COORDS(i,N)                   ith coordinate of Nth node on element
!     MCRD                          Maximum of (# coords,minimum of (3,#DOF)) on any node
!     Uall                          Vector of DOF at the end of the increment
!     DUall                         Vector of DOF increments
!     Vel                           Vector of velocities (defined only for implicit dynamics)
!     Accn                          Vector of accelerations (defined only for implicit dynamics)
!     JTYPE                         Integer identifying element type (the number n in the Un specification in the input file)
!     TIME(1:2)                     TIME(1)   Current value of step time
!                                   TIME(2)   Total time
!     DTIME                         Time increment
!     KSTEP                         Current step number
!     KINC                          Increment number
!     JELEM                         User assigned element number in ABAQUS
!     PARAMS(1:3)                   Time increment parameters alpha, beta, gamma for implicit dynamics
!     NDLOAD                        Number of user-defined distributed loads defined for this element
!     JDLTYP(1:NDLOAD)              Integers n defining distributed load types defined as Un or (if negative) UnNU in input file
!     ADLMAG(1:NDLOAD)              Distributed load magnitudes
!     DDLMAG(1:NDLOAD)              Increment in distributed load magnitudes
!     PREDEF(1:2,1:NPREDF,1:NNODE)  Predefined fields.
!     PREDEF(1,...)                 Value of predefined field
!     PREDEF(2,...)                 Increment in predefined field
!     PREDEF(1:2,1,k)               Value of temperature/temperature increment at kth node
!     PREDEF(1:2,2:NPREDF,k)        Value of user defined field/field increment at kth node
!     NPREDF                        Number of predefined fields
!     LFLAGS                        Load type control variable
!     MLVARX                        Dimension variable
!     PERIOD                        Time period of the current step
!
! **********************************************************************
! **********************************************************************

      ! make sure the relative directory is correct
      include 'global_parameters.for'     ! global parameters module
      include 'error_logging.for'         ! error/ debugging module
      include 'linear_algebra.for'        ! linear algebra module
      include 'nonlinear_solver.for'      ! nonlinear solver module
      include 'lagrange_element.for'      ! Lagrange element module
      include 'gauss_quadrature.for'      ! Guassian quadrature module
      include 'solid_mechanics.for'       ! solid mechanics module
      include 'post_processing.for'       ! post-processing module

! **********************************************************************
! **********************************************************************

      module user_element

      ! This module contains subroutines related to element formulation
      ! and constitutive calculation. Abaqus user subroutines can not
      ! be included in a module. Instead we extended the list of arguments
      ! of the Abaqus UEL subroutine and wrote another subroutine of
      ! similar kind which is included in the user_element module.
      ! Compilers can perform additional checks on the arguments when
      ! any modularized subroutines are called. The first subroutine is
      ! called by UEL subroutine of Abaqus with an extended set of
      ! input arguments. The first subroutine calls other subroutines.

      contains

! **********************************************************************
! **********************************************************************

      SUBROUTINE uel_hydrogel(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,
     & NSVARS,PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     & TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD,
     & NDIM,ANALYSIS,NSTRESS,NINT,UDOF,UDOFEL,MDOF,MDOFEL)

      use global_parameters
      use lagrange_element
      use gauss_quadrature
      use linear_algebra
      use solid_mechanics
      use error_logging

      implicit none

      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     & SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),UAll(NDOFEL),
     & DUAll(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     & JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     & PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      ! input arguments to the subroutine
      integer, intent(in)   :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD
      integer, intent(in)   :: NNODE, JTYPE, KSTEP, KINC, JELEM
      integer, intent(in)   :: NDLOAD, JDLTYP, NPREDF, LFLAGS
      integer, intent(in)   :: MLVARX, MDLOAD, JPROPS, NJPROPS

      real(wp), intent(in)  :: PROPS, COORDS, DUall, Uall, Vel, Accn
      real(wp), intent(in)  :: TIME, DTIME, PARAMS, ADLMAG, PREDEF
      real(wp), intent(in)  :: DDLMAG, PERIOD

      character(len=2), intent(in)    :: analysis
      integer, intent(in)             :: nDim, nStress
      integer, intent(in)             :: uDOF, uDOFEL, mDOF, mDOFEL
      integer, intent(in)             :: nInt

      ! output of the suboutine
      real(wp), intent(out)           :: RHS, AMATRX
      real(wp), intent(out), optional :: SVARS, ENERGY, PNEWDT

      ! variables local to the subroutine
      real(wp)          :: ID(nDim,nDim)

      ! nodal degrees of freedom in the matrix form
      real(wp)          :: UallMat(uDOF+mDOF,nNode)
      real(wp)          :: DUallMat(uDOF+mDOF,nNODE)
      real(wp)          :: uNode(nDim,nNode), muNode(mDOFEL,1)
      real(wp)          :: duNode(nDim,nNode), dmuNode(mDOFEL,1)

      ! finite element matrices and parameters
      real(wp)          :: w(nInt), xi(nInt,nDim)
      real(wp)          :: Nxi(nNode), dNdxi(nNode,nDim)
      real(wp)          :: dXdxi(nDim,nDim), dxidX(nDim,nDim)
      real(wp)          :: detJ, dNdX(nNode,nDim)

      ! finite element operator matrices
      real(wp)          :: Na(nDim,nDim), Nmat(nDim,uDOFEl)
      real(wp)          :: Ga(nDim*nDim,nDim)
      real(wp)          :: Gmat(nDim*nDim,nDim*nNode)
      real(wp)          :: Ba(nStress,nDim), Bmat(nStress,uDOFEl)
      real(wp)          :: BNLmat(nStress,nDim*nNode)
      real(wp)          :: stressTensorPK2(nDim,nDim)
      real(wp)          :: SIGMA_S(nDim*nDim,nDim*nDim)
      real(wp)          :: SIGMA_F(nDim*nNode,nDim*nNode)

      ! integration point quantities (variables)
      real(wp)          :: coord_ip(nDim,1)
      real(wp)          :: F(3,3), detF, Fbar(3,3)
      real(wp)          :: FInv(3,3), FInvT(3,3)
      real(wp)          :: mu, dMUdX(nDim,1)

      ! additional field variables at the nodes and integration point
      real(wp)          :: fieldNode(npredf,nNode)
      real(wp)          :: dfieldNode(npredf,nNode)
      real(wp)          :: fieldVar(npredf), dfieldVar(npredf)

      ! element and material properties used in this subroutine
      real(wp), allocatable :: statev(:)  ! state variables per int pt
      integer               :: matID      ! material constitutive law

      ! output from the material point subroutine (UMAT)
      real(wp)          :: stressPK2(nStress,1)
      real(wp)          :: Dmat(nStress,nStress), Cmat(3,3,3,3)
      real(wp)          :: dCwdt
      real(wp)          :: Jw(nDim,1), dCwdotdmu, dJwdMu(nDim,1)
      real(wp)          :: dVectUM(nDim*nDim,1), DmatMU(nDim,nDim*nDim)
      real(wp)          :: MmatW(nDim,nDim), dCwdotdF(1,nDim*nDim)

      ! additional variables for F-bar method (element and material)
      real(wp)          :: centroid(nDim)
      real(wp)          :: Nxi0(nNode), dNdxi0(nNode,nDim)
      real(wp)          :: dXdxi0(nDim,nDim), dxidX0(nDim,nDim)
      real(wp)          :: dNdX0(nNode,nDim), detJ0
      real(wp)          :: Ga0(nDim**2,nDim), Gmat0(nDim**2,uDOFEl)
      real(wp)          :: F0(3,3), detF0
      real(wp)          :: F0Inv(3,3), F0InvT(3,3)
      real(wp)          :: QR0Tensor(nDim,nDim,nDim,nDim)
      real(wp)          :: QRTensor(nDim,nDim,nDim,nDim)
      real(wp)          :: QR0mat(nDim*nDim,nDim*nDim)
      real(wp)          :: QRmat(nDim*nDim,nDim*nDim)
      real(wp)          :: tanFac1, tanFac2, resFac

      ! tangent matrix components and residual vectors
      real(wp)          :: Kuu(uDOFEl,uDOFEl), Kum(uDOFEl,mDOFEL)
      real(wp)          :: Kmu(mDOFEL,uDOFEl), Kmm(mDOFEL,mDOFEL)
      real(wp)          :: Ru(uDOFEL,1), Rm(mDOFEL,1)
      real(wp)          :: Kelem(nDOFEL,nDOFEL), Relem(nDOFEL,1)

      integer           :: i, j, k, l, m, n, p, q, intPt
      integer           :: nstatev
      type(element)     :: hydrogel
      type(logger)      :: msg

      ! initialize hydrogel element
      hydrogel  = element(nDim=nDim, analysis=analysis,
     &                    nNode=nNode, nInt=nInt)


      ! N_mu = transpose(Nxi), B_mu = transpose(dNdX) in the theory
      F0      = zero
      Fbar    = zero
      Ga0     = zero
      Gmat0   = zero
      F       = zero
      Na      = zero
      Ba      = zero
      Ga      = zero
      Nmat    = zero
      Bmat    = zero
      Gmat    = zero
      SIGMA_F = zero
      SIGMA_S = zero
      Kuu     = zero
      Kum     = zero
      Kmu     = zero
      Kmm     = zero
      Ru      = zero
      Rm      = zero
      Kelem   = zero
      Relem   = zero

      matID     = jprops(2)
      nstatev   = nsvars/nint       ! state variables per integration point

      if ( .not. allocated(statev) ) allocate(statev(nstatev))

      ! reshape the degrees of freedom and their change into matrices
      uAllMat   = reshape( uAll, shape=[uDOF+mDOF,nNode] )
      duAllMat  = reshape( duAll(:,1), shape=[uDOF+mDOF,nNode] )

      ! seperate the displacement and chemical potential field
      uNode(1:uDOF,1:nNode)   = uAllMat(1:uDOF,1:nNode)
      muNode(1:mDOFEL,1)      = uAllMat(uDOF+1,1:mDOFEL)

      ! seperate the change of displacement and chemical potential field
      duNode(1:uDOF,1:nNode)  = duAllMat(1:uDOF,1:nNode)
      dmuNode(1:mDOFEL,1)     = duAllMat(uDOF+1,1:mDOFEL)

      call eyeMat(ID)                       ! create an identity matrix w/ size of nDim

      ! For fully-integrated linear quad and hex elements, calculate Gmat0.
      ! These calculations are done to evaluate volumetric deformation
      ! gradient at centroid which will be used in to define F-bar later.
      if ( ((jtype .eq. 2) .and. (nInt .eq. 8))
     &    .or. ((jtype .eq. 4) .and. (nInt .eq. 4)) ) then

        centroid = zero

        ! evaluate the interpolation functions and derivates at centroid
        call calcInterpFunc(hydrogel, centroid, Nxi0, dNdxi0)

        ! calculate element jacobian and global shape func gradient at centroid
        dXdxi0  = matmul(coords,dNdxi0)       ! calculate the jacobian (dXdxi) at centroid
        detJ0   = det(dXdxi0)                 ! calculate jacobian determinant at centroid

        if (detJ0 .le. zero) then
          call msg%ferror( flag=warn, src='uelHydrogel',
     &      msg='Negative element jacobian at centroid.', ia=jelem)
        end if

        dxidX0 = inv(dXdxi0)                  ! calculate jacobian inverse
        dNdX0  = matmul(dNdxi0,dxidX0)        ! calculate dNdX at centroid

        do i = 1, nNode

          ! form the nodal-level matrix: [Ga0] at the centroid
          do j = 1, nDim
            Ga0(nDim*(j-1)+1:nDim*j,1:nDim) = dNdX0(i,j)*ID
          end do

          ! form the [G0] matrix at the centroid
          Gmat0(1:nDim**2,nDim*(i-1)+1:nDim*i) = Ga0(1:nDim**2,1:nDim)
        end do                             ! end of nodal point loop

        F0(1:nDim,1:nDim) = ID + matmul(uNode,dNdX0)

        if (analysis .eq. 'PE') F0(3,3) = one

        detF0   = det(F0)
        F0Inv   = inv(F0)
        F0InvT  = transpose(F0Inv)

      end if
      ! end of centroid level calculation for F-bar

      ! obtain the standard gauss quadrature points and weights
      call getGaussQuadrtr(hydrogel,w,xi)


      do intPt = 1, nInt

        call calcInterpFunc(hydrogel,xi(intPt,:),Nxi,dNdxi)

        ! calculate element jacobian and global shape function gradient
        dxdxi = matmul(coords,dNdxi)        ! calculate dxdxi
        detJ  = det(dxdxi)                  ! calculate determinant

        if (detJ .lt. zero) then
          call msg%ferror( flag=warn, src='uelHydrogel',
     &          msg='Negative element jacobian.', ivec=[jelem, intpt])
        end if

        dxidX = inv(dxdxi)                  ! calculate inverse
        dNdX  = matmul(dNdxi,dxidX)         ! calculate dNdX


        ! loop over all the nodes (internal loop)
        do i=1,nNode

          ! form the nodal-level matrices: [Na] and [Ga]
          do j = 1, nDim
            Na(j,j) = Nxi(i)
            Ga(nDim*(j-1)+1:nDim*j,1:nDim) = dNdX(i,j)*ID
          end do

          ! form [Ba] matrix: 3D case
          if (analysis .eq. '3D') then
            Ba(1,1)       = dNdx(i,1)
            Ba(2,2)       = dNdx(i,2)
            Ba(3,3)       = dNdx(i,3)
            Ba(4,1:nDim)  = [  zero,      dNdx(i,3),  dNdx(i,2)]
            Ba(5,1:nDim)  = [dNdx(i,3),     zero,     dNdx(i,1)]
            Ba(6,1:nDim)  = [dNdx(i,2),   dNdx(i,1),    zero   ]

          ! form [Ba] matrix: plane stress/ plane strain case
          else if (analysis .eq. 'PE') then
            Ba(1,1)       = dNdx(i,1)
            Ba(2,2)       = dNdx(i,2)
            Ba(3,1:nDim)  = [dNdx(i,2), dNdx(i,1)]
          else 
            call msg%ferror( flag=error, src='uelHydrogel',
     &                msg='Wrong analysis.', ch=analysis )
            call xit
          end if

          ! form the [N], [B], and [G] matrices
          Nmat(1:nDim,nDim*(i-1)+1:nDim*i)    = Na(1:nDim,1:nDim)
          Bmat(1:nStress,nDim*(i-1)+1:nDim*i) = Ba(1:nStress,1:nDim)
          Gmat(1:nDim**2,nDim*(i-1)+1:nDim*i) = Ga(1:nDim**2,1:nDim)
        end do                             ! end of nodal point loop

      !!!!!!!!!!!!!! COMPLETE ELEMENT RELATED OPERATIONS !!!!!!!!!!!!!!!



      !!!!!!!!!!!!!!!!!!!!!!! CONSTITUTIVE MODEL !!!!!!!!!!!!!!!!!!!!!!!

        ! calculate the coordinate of integration point
        coord_ip = matmul(Nmat, reshape(coords, [nDOFEL, 1]))


        ! calculate deformation gradient and deformation tensors
        F(1:nDim,1:nDim) = ID + matmul(uNode,dNdX)

        if (analysis .eq. 'PE')  F(3,3) = one

        ! calculate chemical potential and its gradient
        mu    = dot_product( Nxi, reshape(muNode, [mDOFEl]) )

        dMUdX = matmul(transpose(dNdX), muNode)

        ! F-bar modification
        ! calculate material point jacobian (volume change)
        detF    = det(F)
        FInv    = inv(F)
        FInvT   = transpose(FInv)

        ! definition of modified deformation gradient
        if ((jtype .eq. 2) .and. (nInt .eq. 8)) then
          ! fully-integrated three-dimensional trilinear hex element
          Fbar    = (detF0/detF)**(third) * F
          resFac  = (detF0/detF)**(-two/three)
          tanFac1 = (detF0/detF)**(-third)
          tanFac2 = (detF0/detF)**(-two/three)
        else if ((jtype .eq. 4) .and. (nInt .eq. 4))  then
          ! fully-integrated plane strain bilinear quad element
          Fbar(3,3)           = one
          Fbar(1:nDim,1:nDim) = (detF0/detF)**(half) * F(1:nDim,1:nDim)
          resFac              = (detF0/detF)**(-half)
          tanFac1             = one
          tanFac2             = (detF0/detF)**(-half)
        else
          ! standard F for all other types of implemented elements
          Fbar    = F
          resFac  = one
          tanFac1 = one
        end if


        ! call material point subroutine (UMAT) for the neutral gel
        call umat_hydrogel(kstep,kinc,time,dtime,nDim,analysis,
     &          nStress,nNode,jelem,intpt,coord_ip,props,nprops,
     &          jprops,njprops,matID,Fbar,mu,dMudX,
     &          svars,nsvars,fieldVar,dfieldVar,npredf,
     &          stressPK2,Jw,dCwdt,
     &          Dmat,Cmat,dVectUM,dCwdotdF,DmatMU,dCwdotdMu,dJwdMu,
     &          MmatW)

      !!!!!!!!!!!!!!!!!!!! END CONSTITUTIVE MODEL !!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!! TANGENT MATRIX AND RESIDUAL VECTOR !!!!!!!!!!!!!!!!!

        call voigtVectorScatter(stressPK2,stressTensorPK2)

        ! form [SIGMA_F] matrix for material stiffness
        do i = 1, nNode
          do j = 1, nNode
            if (i .eq. j) then
              SIGMA_F(nDim*(i-1)+1:nDim*i,nDim*(j-1)+1:nDim*j)
     &                           = Fbar(1:nDim,1:nDim)
            end if
          end do
        end do

        ! form the [SIGMA_S] matrix for geometric stiffness
        do i = 1, nDim
          do j = 1, nDim
            SIGMA_S(nDim*(i-1)+1:nDim*i,nDim*(j-1)+1:nDim*j) =
     &                    stressTensorPK2(i,j)*ID
          end do
        end do

        BNLmat = matmul(Bmat,SIGMA_F)


        ! form the tangent matrices
        ! mechanical tangent matrix
        Kuu = Kuu + w(intPt) * detJ * tanFac1 *
     &      (
     &      matmul( transpose(BNLmat), matmul(Dmat,BNLmat) )
     &      + matmul( matmul( transpose(Gmat),SIGMA_S ), Gmat )
     &      )


        ! perform F-bar modification on linear hex element
        if ((jtype .eq. 2) .and. (nInt .eq. 8)) then
          ! form fourth-order QR0 and QR tensor
          QR0Tensor = zero
          QRTensor  = zero

          do i = 1,nDim
            do j = 1,nDim
              do k = 1,nDim
                do l = 1,nDim
                  do m = 1,nDim
                    do n = 1,nDim
                      do p = 1,nDim
                        do q = 1,nDim

                          QR0Tensor(i,j,k,l) = QR0Tensor(i,j,k,l)
     &                        + third * F0InvT(k,l) *
     &                          (
     &                            Fbar(i,p) * Cmat(p,j,m,n)
     &                            * Fbar(q,m) * Fbar(q,n)
     &                            - Fbar(i,q) * stressTensorPK2(q,j)
     &                          )

                          QRTensor(i,j,k,l) = QRTensor(i,j,k,l)
     &                        + third * FInvT(k,l) *
     &                          (
     &                            Fbar(i,p) * Cmat(p,j,m,n)
     &                            * Fbar(q,m) * Fbar(q,n)
     &                            - Fbar(i,q) * stressTensorPK2(q,j)
     &                          )
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do

          ! reshape QR and QR0 tensor into matrix form
          call unsymmMatrix(QR0Tensor,QR0mat)
          call unsymmMatrix(QRTensor,QRmat)

          ! modify the element tangent matrix
          Kuu   = Kuu + w(intPt) * detJ * tanFac2  *
     &              (
     &              matmul(transpose(Gmat), matmul(QR0mat,Gmat0))
     &              - matmul(transpose(Gmat), matmul(QRmat,Gmat))
     &              )

        end if


        ! do F-bar modification on plane strain linear quad element
        if ( (jtype .eq. 4) .and. (nInt .eq. 4) ) then
          ! form fourth-order QR0 and QR tensor
          QR0Tensor = zero
          QRTensor  = zero

          do i = 1,nDim
            do j = 1,nDim
              do k = 1,nDim
                do l = 1,nDim
                  do m = 1,nDim
                    do n = 1,nDim
                      do p = 1,nDim
                        do q = 1,nDim
                          QR0Tensor(i,j,k,l) = QR0Tensor(i,j,k,l)
     &                        + half * Fbar(i,p) * Cmat(p,j,q,n)
     &                        * Fbar(m,q) * Fbar(m,n) * F0InvT(k,l)

                          QRTensor(i,j,k,l) = QRTensor(i,j,k,l)
     &                        + half * Fbar(i,p) * Cmat(p,j,q,n)
     &                        * Fbar(m,q) * Fbar(m,n) * FInvT(k,l)
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do

          ! reshape QR and QR0 tensor into matrix form
          call unsymmMatrix(QR0Tensor,QR0mat)
          call unsymmMatrix(QRTensor,QRmat)

          ! modify the element tangent matrix
          Kuu = Kuu + w(intPt) * detJ * tanFac2  *
     &              (
     &              matmul(transpose(Gmat), matmul(QR0mat,Gmat0))
     &              - matmul(transpose(Gmat), matmul(QRmat,Gmat))
     &              )

        end if


        ! mechano-chemical tangent matrix
        Kum = Kum + w(intpt) * detJ *
     &            matmul( matmul( transpose(Gmat), dVectUM ),
     &                    reshape(Nxi, [1, nNode]) )


        ! chemo-mechanical tangent matrix
        Kmu = Kmu + w(intPt) * detJ *
     &            (
     &            matmul( reshape( Nxi, [nNode, 1] ),
     &                    matmul(dCwdotdF, Gmat) )
     &            - matmul( matmul(dNdX, DmatMU), Gmat)
     &            )

        ! chemical tangent matrix
        Kmm = Kmm + w(intPt) * detJ *
     &            (
     &            matmul( reshape(Nxi, [nNode, 1]),
     &                    reshape(Nxi, [1, nNode]) ) * dCwdotdMU
     &            - matmul( matmul(dNdX, dJwdMu),
     &                      reshape(Nxi, [1, nNode]) )
     &            + matmul( matmul(dNdX, MmatW), transpose(dNdX) )
     &            )

        ! mechanical residual (TODO: add body force)
        Ru = Ru - w(intPt) * detJ * resFac *
     &            matmul( transpose(BNLmat), stressPK2 )

        ! chemical residual
        Rm = Rm + w(intPt) * detJ *
     &          ( - reshape(Nxi, [nNode,1]) * dCwdt + matmul(dNdX,Jw) )

      end do                           ! end of integration point loop


      ! apply traction boundary condition on the overlaying mechanical elements
      ! TODO: add code snippet to calculate the flux for the solvent (modify Rm)


      ! assemble the element tangent sub-matrices and residual sub-vectors
      call assembleElement(nDim,nNode,uDOFEL,mDOFEL,nDOFEl,       ! input to the subroutine
     &     Kuu,Kum,Kmu,Kmm,Ru,Rm,                                 ! input to the subroutine
     &     Kelem,Relem)                                           ! output from the subroutine


      ! assign the element tangent matrices to Abaqus defined variables
      amatrx(1:NDOFEL,1:NDOFEL)   = Kelem(1:NDOFEL,1:NDOFEL)
      rhs(1:NDOFEL,1)             = Relem(1:NDOFEL,1)


      end subroutine uel_hydrogel

! **********************************************************************
! **********************************************************************

      subroutine umat_hydrogel(kstep,kinc,time,dtime,nDim,analysis,
     &          nStress,nNode,jelem,intpt,coord_ip,props,nprops,
     &          jprops,njprops,matID,F,mu,dMudX,
     &          svars,nsvars,fieldVar,dfieldVar,npredf,
     &          stressPK2,Jw,dCwdt,
     &          Dmat,Cmat,dVectUM,dCwdotdF,DmatMU,dCwdotdMu,dJwdMu,
     &          MmatW)

      ! This material point subroutine calculates constitutive response
      ! of an elastomeric hydrogel and returns the flux and tangents.
      ! Optionally, it can return some other strain and stress quantities
      ! in vector form (needs to be programmed).
      ! This material subroutine also stores the user-defined element
      ! output in a global array for post=processing in Abaqus/Viewer.

      use global_parameters
      use solid_mechanics
      use linear_algebra
      use nonlinear_solver
      use post_processing
      use error_logging

      implicit none

      ! minimum and maximum value for polymer volume fraction
      real(wp), parameter :: phiMin = 0.05_wp, phiMax = 0.9999_wp

      ! input arguments to the subroutine
      character(len=2), intent(in)  :: analysis

      integer, intent(in)   :: kstep, kinc, nDim, nstress
      integer, intent(in)   :: nNode, jelem, intpt, nprops
      integer, intent(in)   :: njprops, nsvars, npredf
      integer, intent(in)   :: matID

      real(wp), intent(in)  :: time(2), dtime
      real(wp), intent(in)  :: coord_ip(nDim,1)
      real(wp), intent(in)  :: props(nprops)
      integer,  intent(in)  :: jprops(njprops)

      real(wp), intent(in)  :: F(3,3), mu, dMUdX(nDim,1)
      real(wp), intent(in)  :: fieldVar(npredf)
      real(wp), intent(in)  :: dfieldVar(npredf)

      ! output from the subroutine
      real(wp), intent(out) :: stressPK2(nStress,1)
      real(wp), intent(out) :: Jw(nDim,1)
      real(wp), intent(out) :: dCwdt
      real(wp), intent(out) :: Dmat(nstress,nstress)
      real(wp), intent(out) :: Cmat(3,3,3,3)
      real(wp), intent(out) :: dVectUM(nDim*nDim,1)
      real(wp), intent(out) :: dCwdotdF(1,nDim*nDim)
      real(wp), intent(out) :: DmatMU(nDim,nDim*nDim)
      real(wp), intent(out) :: dCwdotdMu, dJwdMu(nDim,1)
      real(wp), intent(out) :: MmatW(nDim,nDim)

      real(wp), intent(inout), optional   :: svars(nsvars)


      ! local variables (kinematic quantities)
      real(wp)          :: detF, FInv(3,3), FInvT(3,3)
      real(wp)          :: C(3,3), CInv(3,3), detC, trC
      real(wp)          :: B(3,3), Binv(3,3), detB
      real(wp)          :: strainTensorLagrange(3,3)
      real(wp)          :: strainTensorEuler(3,3)

      ! local variables (internal variables)
      real(wp)          :: phi_old, phi_new, dPhidt, Cw_old, Cw_new
      real(wp)          :: vars(nprops+2)
      real(wp)          :: detFe, detFs

      ! local variables (stress tensors)
      real(wp)          :: stressTensorPK1(3,3)
      real(wp)          :: stressTensorPK2(3,3)
      real(wp)          :: stressTensorCauchy(3,3)

      ! tangent matrices and related quantities
      real(wp)          :: dCwdPhi, dPhidG, dGdPhi, dCwdMU
      real(wp)          :: dGdFTensor(3,3)
      real(wp)          :: dCwdFTensor(3,3)
      real(wp)          :: dCwdCTensor(3,3)
      real(wp)          :: dFdCwTensor(3,3)
      real(wp)          :: dSdCwTensor(3,3)
      real(wp)          :: FSTensorUM(3,3)
      real(wp)          :: JwTensor(nDim,nDim,nDim)
      real(wp)          :: dCwdotdFTensor(3,3)

      ! intermeidate variables for post-processing and output
      real(wp)          :: strainVectLagrange(nSymm,1)
      real(wp)          :: strainVectEuler(nSymm,1)
      real(wp)          :: stressVectPK1(nUnsymmm,1)
      real(wp)          :: stressVectPK2(nSymm,1)
      real(wp)          :: stressVectCauchy(nSymm,1)
      real(wp)          :: VoigtMat(nSymm,nSymm)

      ! strain and stress vectors for output purposes
      real(wp)          :: strainLagrange(nStress,1)
      real(wp)          :: strainEuler(nStress,1)
      real(wp)          :: stressPK1(nDim*nDim,1)
      real(wp)          :: stressCauchy(nStress,1)

      ! material properties
      real(wp)          :: Rgas, theta, phi0, rho, Gshear, Kappa
      real(wp)          :: lam_L, Vp, mu0, Vw, chi, Dw, RT

      real(wp)          :: lam_c, lam_r, beta_0, beta_c, dBeta_c

      integer           :: i, j, k, l

      type(logger)      :: msg

      ! initialize matrial stiffness tensors
      Cmat        = zero
      Dmat        = zero
      JwTensor    = zero

      ! assign material properties to variables
      Rgas        = props(1)        ! universal gas constant
      theta       = props(2)        ! absolute temperature (K)
      phi0        = props(3)        ! initial polymer volume fraction
      rho         = props(4)        ! density of the hydrogel
      Gshear      = props(5)        ! shear modulus
      Kappa       = props(6)        ! bulk modulus
      lam_L       = props(7)        ! locking stretch (only for AB model)
      Vp          = props(8)        ! molar volume of polymer
      mu0         = props(9)        ! chemical potential of pure solvent
      Vw          = props(10)       ! molar volume of the solvent
      chi         = props(11)       ! flory-huggins mixing parameter
      Dw          = props(12)       ! diffusion coefficient of the solvent

      RT          = Rgas*theta

      !!!!!!!!!!!!!!!!!!!!!!!!! KINEMATIC PART !!!!!!!!!!!!!!!!!!!!!!!!!

      detF        = det(F)

      if (detF .le. zero) then
        call msg%ferror(flag=error, src='umatHydrogel',
     &        msg='Issue with volume change (detF, jelem, intPt)',
     &        ra= detF, ivec=[jelem, intpt])
        call xit
      end if

      FInv  = inv(F)
      FInvT = transpose(Finv)

      C     = matmul(transpose(F),F)
      B     = matmul(F,transpose(F))
      CInv  = inv(C)
      Binv  = inv(B)
      trC   = trace(C)

      ! calculate Euler-Almansi strain tensor
      strainTensorLagrange  = half*(C-ID3)
      strainTensorEuler     = half*(ID3-Binv)

      !!!!!!!!!!!!!!!!!!!!!!!!! KINEMATIC PART !!!!!!!!!!!!!!!!!!!!!!!!!



      !!!!!!!!!!!!!!!!!!!!! SOLVE INT PT VARIABLES !!!!!!!!!!!!!!!!!!!!!

      ! (1.1) get the internal variable, phi (polymer volume fraction)
      if( (kinc .le. 1) .and. (kstep .eq. 1) ) then
        phi_old  = phi0               ! read the initial value at the first step
      else
        phi_old  = svars(intPt)       ! read stored state variables
      end if


      ! (1.2) local iteration to solve integration pt variable
      vars(1:nprops)  = props
      vars(nprops+1)  = detF
      vars(nprops+2)  = mu

      call fzero( chemicalState, phi_old, phi_new,
     &            phiMin, phiMax, jac=.true., vars=vars)

      ! update the state variable for next step (or iteration)
      svars(intPt)    = phi_new

      ! (1.3) calculate the solution-dependent variables
      detFs     = one/phi_new
      detFe     = detF/(phi0*detFs)

      ! (1.4) calculate solvent concentration
      Cw_old    = phi0*(one/phi_old-one)/Vw
      Cw_new    = phi0*(one/phi_new-one)/Vw

      ! (1.5) calculate time derivatives: dPhi/dt and dCw/dt
      dPhidt    = (phi_new-phi_old)/dtime

      dCwdt     = - phi0*dPhidt/( Vw*(phi_new)**two )
      ! dCwdt   = (Cw_new - Cw_old)/dtimes

      !!!!!!!!!!!!!!!!!!!!! SOLVE INT PT VARIABLES !!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!! TANGENT CALCULATION !!!!!!!!!!!!!!!!!!!!!!!

      ! (2.1) calculate dCw/dPhi
      dCwdPhi = - (phi0/ (Vw*(phi_new)**two) )

      ! (2.2) calculate dG/dPhi (same as dMu/dPhi)
      dGdPhi = RT * (one - one/(one-phi_new) + two*chi*phi_new)
     &            - Kappa*Vw/phi_new
     &            + (Kappa*Vw/phi_new)*log(detF*phi_new/phi0)

      ! (2.3) calculate dPhi/dG (same as dPhi/dMu)
      dPhidG = one/dGdPhi

      ! (2.4) calculate dCw/dMu
      dCwdMU = dCwdPhi*dPhidG


      ! (3.1) calculate dG/dF = dMu/dF
      dGdFTensor  = Kappa * Vw * (log(detFe)-one) * FInvT

      ! (3.2) calculate dCw/dF (second order tensor)
      dCwdFTensor = - dCwdPhi * dPhidG * dGdFTensor

      ! (3.3) calculate dCw/dC
      dCwdCTensor = half * dCwdPhi * dPhidG
     &              * Kappa * Vw * (log(detFe)-one) * CInv


      ! (4) dS/dCw (a symmetric second order tensor)
      dSdCwTensor = Kappa * Vw * (log(detFe)-one) * CInv


      ! (5) stress tensor and material tangent
      if (matID .eq. 1) then                  ! neo-hookean elstomeric gel
        ! (5.1) stress tensors
        stressTensorPK2 = Gshear * ( ID3 - (phi0)**(two/three) * CInv )
     &                  + Kappa * phi0 * detFs * log(detFe) * CInv

        stressTensorCauchy = (one/detF)
     &                      * ( Gshear * (B - (phi0)**(two/three) * ID3)
     &                      + Kappa * phi0 * detFs * log(detFe) * ID3 )


        ! (5.2) calculate material tangent
        do i = 1,3
          do j = 1,3
            do k = 1,3
              do l = 1,3
                Cmat(i,j,k,l) = Cmat(i,j,k,l)
     &              + Kappa * phi0 * detFs * CInv(i,j) * CInv(k,l)
     &              + ( (phi0)**(two/three) * Gshear
     &              - Kappa * phi0 * detFs * log(detFe) )
     &              * ( CInv(i,k) * CInv(j,l) + CInv(i,l) * CInv(j,k) )
     &              + two * dSdCwTensor(i,j) * dCwdCTensor(k,l)
              end do
            end do
          end do
        end do

      else if (matID .eq. 2) then          ! Arruda-Boyce elastomeric gel
        lam_c     = sqrt(trC/three)
        lam_r     = lam_c/lam_L
        beta_0    = InvLangevin( (phi0**third)/ lam_L) 
        beta_c    = InvLangevin(lam_r)
        dBeta_c   = DInvLangevin(lam_r)

        ! (5.1) stress tensors
        stressTensorCauchy  = (1/detF) * ( (Gshear/three)*lam_r*beta_c*B
     &        - ( (Gshear*lam_L)/three * beta_0 * (phi0)**(two/three) 
     &        - Kappa*phi0*detFs*log(detFe) ) * ID3 )

        stressTensorPK2     = (Gshear/three) * lam_r * beta_c * ID3
     &        - ( ( (phi0)**(two/three) *Gshear*lam_L )/three 
     &        - Kappa*phi0*detFs*log(detFe) ) * Cinv

        ! (5.2) calculate material tangent
        do i = 1,3
          do j = 1,3
            do k = 1,3
              do l = 1,3
                Cmat(i,j,k,l) = Cmat(i,j,k,l)
     &              + Gshear/(nine*lam_c**two)
     &              * ( dBeta_c- lam_r*beta_c ) * ID3(i,j)*ID3(k,l)
     &              + Kappa * phi0 * detFs * Cinv(i,j)*Cinv(k,l)
     &              + ( (Gshear/three) * lam_L 
     &              - Kappa * phi0 * detFs * log(detFe) )
     &              * ( Cinv(i,k)*Cinv(j,l) + Cinv(i,l)*Cinv(j,k) )
     &              + two * dSdCwTensor(i,j) * dCwdCTensor(k,l)
              end do
            end do
          end do
        end do

      else
        call msg%ferror(flag=error, src='umatHydrogel',
     &        msg='Wrong material ID.', ia=matID)
        call xit
      end if

      ! (6.1) mechano-chemical tangent: F*dSdCwTensor*dCwdMU
      FSTensorUM  = matmul(F,dSdCwTensor) * dCwdMU

      ! (6.2) reshape into vector form
      dVectUM     = reshape( FSTensorUM(1:nDim,1:nDim), shape(dVectUM) )



      ! (7.1) mobility matrix : Mw = Dw*Cw/RT*Inv(C) (dimension-dependent)
      MmatW = (Dw*Cw_new/RT)*CInv(1:nDim,1:nDim)


      ! (7.2) calculate the solvent molar flux: Jw = - Mw*Grad(mu)
      Jw    = - matmul(MmatW,dMUdX)


      ! (7.3) calculate Jw_tensor
      JwTensor = zero
      do i = 1,nDim
        do k = 1,nDim
          do l = 1,nDim
            do j = 1,nDim               ! summation over dummy index j
              JwTensor(i,k,l) = JwTensor(i,k,l)
     &            + (Dw*Cw_new)/RT 
     &            * ( FInv(i,k)*CInv(l,j) ) * dMUdX(j,1)
     &            - (Dw/RT) * CInv(i,j) * dMUdX(j,1) * dCwdFTensor(k,l)
            end do
          end do
        end do
      end do


      ! (7.4) map the third-order tensor, JwTensor, to a rank-2 matrix
      do i = 1, nDim
        do l = 1, nDim
          do k = 1, nDim
            DmatMU(i,(l-1)*nDim+k) = JwTensor(i,k,l)
          end do
        end do
      end do


      ! (8) Calculate dJw/dMU = (dJw/dCw)*(dCw/dMU) (dimension-dependent)
      dJwdMU  = - (Dw/RT) * matmul( CInv(1:nDim,1:nDim),dMUdX ) * dCwdMU


      ! (9.1) calculate dCw_dot/dMu (scalar)
      dCwdotdMU       = dCwdMu/dtime


      ! (9.2) calculate of dCw_dot/dF and reshape it (dimension-dependent)
      dCwdotdFTensor  = dCwdFTensor/dtime
      dCwdotdF        = reshape(dCwdotdFTensor(1:nDim,1:nDim),
     &                          shape(dCwdotdF))


      !!!!!!!!!!!!!!!!!!!!!! TANGENT CALCULATION !!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!! ANALYSIS-BASED RESHAPING OF TENSORS !!!!!!!!!!!!!!!

      ! transforms the stiffness tensor (3x3x3x3) to Voigt matrix (6x6)
      ! transform the strain/ stress tensor (3x3) to Voigt vector form (6x1)
      call voigtVector(stressTensorPK2,stressVectPK2)
      call voigtMatrix(Cmat,VoigtMat)

      ! truncate the Voigt vector and matrix for FEM formulation (if needed)
      ! for 2D PE case, it performs truncation: Dmat (3x3), stressPK2 (3x1)
      ! for 3D case, it returns the same output as input argument
      call voigtVectorTruncate(stressVectPK2,stressPK2)
      call voigtMatrixTruncate(VoigtMat,Dmat)


      ! perform reshape and truncation (if needed) for post-processing
      call voigtVector(strainTensorLagrange,strainVectLagrange)
      call voigtVector(strainTensorEuler,strainVectEuler)
      call voigtVector(stressTensorCauchy,stressVectCauchy)

      call voigtVectorTruncate(strainVectLagrange,strainLagrange)
      call voigtVectorTruncate(strainVectEuler,strainEuler)
      call voigtVectorTruncate(stressVectCauchy,stressCauchy)

      !!!!!!!!!!!!!! ANALYSIS-BASED RESHAPING OF TENSORS !!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!! POST-PROCESSING SECTION !!!!!!!!!!!!!!!!!!!!!

      ! save the variables to be post-processed in globalPostVars
      ! more variables can be added for element level output
      globalPostVars(jelem,intPt,1:nStress) = stressCauchy(1:nStress,1)
      globalPostVars(jelem,intPt,nStress+1:2*nStress)
     &                                      = strainEuler(1:nStress,1)
      globalPostVars(jElem,intPt,2*nStress+1) = phi_new

      !!!!!!!!!!!!!!!!!!!! POST-PROCESSING SECTION !!!!!!!!!!!!!!!!!!!!!

! **********************************************************************

      contains

      function InvLangevin(x)

      ! calculates an approximation of the inverse Langevin function
      ! reference: Bergstorm (PhD thesis, MIT, 1999)

      implicit none

      real(wp),intent(in)  :: x
      real(wp)             :: InvLangevin

      if (abs(x) .lt. 0.84136_wp) then
        InvLangevin = 1.31446_wp*tan(1.58986_wp*x) + 0.91209_wp*x

      else if ((abs(x) .ge. 0.84136_wp) .and. (abs(x) .lt. one)) then
        InvLangevin = one/(sign(one,x)-x)
        
      else
        call msg%ferror(flag=error, src='umatArrudaBoyce:InvLangevin',
     &                  msg='Unbound argument.', ra = x)
        call xit
      end if

      end function InvLangevin

! **********************************************************************

      function DInvLangevin(x)

      ! calculates an approximation of derivative of
      ! the inverse Langevin function
      ! reference: Bergstorm (PhD thesis, MIT, 1999)

      implicit none

      real(wp), intent(in)   :: x
      real(wp)               :: DInvLangevin, sec

      if (abs(x) .lt. 0.84136_wp) then
        DInvLangevin = 2.0898073756_wp*(tan(1.58986_wp*x))**two
     &                + 3.0018973756_wp
      else if ((abs(x) .ge. 0.84136_wp) .and. (abs(x) .lt. one)) then
        DInvLangevin = one/( (sign(one,x)-x)**two )
      else
       call msg%ferror(flag=error, src='umatArrudaBoyce:DInvLangevin',
     &                  msg='Unbound argument.', ra = x)
        call xit
      end if

      end function DInvLangevin

      end subroutine umat_hydrogel

! **********************************************************************

      subroutine chemicalState(phi, f, df, vars)

      use global_parameters, only: wp, zero, one, two

      implicit none

      real(wp), intent(in)              :: phi
      real(wp), intent(out)             :: f
      real(wp), intent(out), optional   :: df
      real(wp), intent(in), optional    :: vars(:)

      ! material properties
      real(wp)    :: Rgas, theta, phi0, rho, Gshear, Kappa
      real(wp)    :: lam_L, Vp, mu0, Vw, chi, Dw, RT

      ! state variables
      real(wp)    :: detF, mu

      ! assign material properties to variables
      Rgas        = vars(1)         ! universal gas constant
      theta       = vars(2)         ! absolute temperature (K)
      phi0        = vars(3)         ! initial polymer volume fraction
      rho         = vars(4)         ! density of the hydrogel
      Gshear      = vars(5)         ! shear modulus
      Kappa       = vars(6)         ! bulk modulus
      lam_L       = vars(7)         ! locking stretch (only for AB model)
      Vp          = vars(8)         ! molar volume of polymer
      mu0         = vars(9)         ! chemical potential of pure solvent
      Vw          = vars(10)        ! molar volume of the solvent
      chi         = vars(11)        ! flory-huggins parameter
      Dw          = vars(12)        ! diffusion coefficient of the solvent
      detF        = vars(13)        ! volume change
      mu          = vars(14)        ! chemical potential

      RT          = Rgas*theta

      ! calculate the residual
      f = (mu0 - mu)
     &      + RT*(phi + log(one-phi) + chi*phi**two)
     &      - Kappa*Vw * log(detF*phi/phi0)
     &      +  ( (Kappa*Vw)/two ) * ( log(detF*phi/phi0) )**two

        ! analytical gradient of the residual
        if ( present(df) ) then
          df = RT * ( one - ( one/(one - phi) ) + two*chi*phi )
     &            - (Kappa*Vw)/phi
     &            + ( (Kappa*Vw)/phi ) * log(detF*phi/phi0)
        endif

      end subroutine chemicalState

! **********************************************************************
! **********************************************************************

      subroutine assembleElement(nDim, nNode, uDOFEL, mDOFEL, nDOFEL,
     &     Kuu, Kum, Kmu, Kmm, Ru, Rm, Kelem, Relem)

      ! This subroutine performs assembly of the element residual
      ! vectors and the element tangent matrix of a linear element.

      use global_parameters, only: wp, zero

      implicit none

      integer,  intent(in)  :: nDim, nNode, uDOFEL, mDOFEL, nDOFEL

      real(wp), intent(in)  :: Kuu(uDOFEL,uDOFEL), Kum(uDOFEl,mDOFEL)
      real(wp), intent(in)  :: Kmu(mDOFEL,uDOFEL), Kmm(mDOFEL,mDOFEL)
      real(wp), intent(in)  :: Ru(uDOFEL,1), Rm(mDOFEL,1)

      real(wp), intent(out) :: Kelem(ndofel,ndofel), Relem(ndofel,1)

      ! loop counters and indexes
      integer               :: i, j, nDOF, R11, R12, C11, C12

      ! total number of degrees of freedom per node
      nDOF = nDOFEl/nNode

      ! assemble the element level residuals and tangent matrices
      do i = 1, nNode

        ! residual vector assembly
        R11 = nDOF*(i-1)+1
        R12 = nDim*(i-1)+1

        Relem(R11:R11+(nDim-1),1) = Ru(R12:R12+(nDim-1),1)
        Relem(R11+nDim,1)         = Rm(i,1)

        do j = 1, nNode

          ! tangent matrix assembly
          C11 = nDOF*(j-1)+1
          C12 = nDim*(j-1)+1

          ! mechanical tangent
          Kelem(R11:R11+(nDim-1), C11:C11+(nDim-1)) =
     &            Kuu( R12:R12+(nDim-1),C12:C12+(nDim-1) )

          ! mechano-chemical tangent
          Kelem(R11:R11+(nDim-1),C11+nDim) = Kum(R12:R12+nDim-1,j)

          ! chemo-mechanical tangent
          Kelem(R11+nDim,C11:C11+(nDim-1)) = Kmu(i,C12:C12+(nDim-1))

          ! chemical tangent
          Kelem(R11+nDim,C11+nDim)         = Kmm(i,j)

        end do
      end do

      end subroutine assembleElement

      end module user_element

! **********************************************************************
! ****************** ABAQUS USER ELEMENT SUBROUTINE ********************
! **********************************************************************

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD)

      ! This subroutine is called by Abaqus with following arguments
      ! for each user elements defined in an Abaqus model. Users are
      ! responsible for programming the element tangent/ stiffness
      ! matrix and residual vectors which will be then assembled and
      ! solved by Abaqus after applying the boundary conditions.

      use global_parameters
      use error_logging
      use user_element
      use post_processing

      INCLUDE 'ABA_PARAM.INC'

      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     & SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),UAll(NDOFEL),
     & DUAll(MLVARX,*),Vel(NDOFEL),Accn(NDOFEL),TIME(2),PARAMS(*),
     & JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     & PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      ! user coding to define RHS, AMATRX, SVARS, ENERGY, and PNEWDT
      integer, intent(in)   :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD
      integer, intent(in)   :: NNODE, JTYPE, KSTEP, KINC, JELEM
      integer, intent(in)   :: NDLOAD, JDLTYP, NPREDF, LFLAGS
      integer, intent(in)   :: MLVARX, MDLOAD, JPROPS, NJPROPS

      real(wp), intent(in)  :: PROPS, COORDS, DUall, Uall, Vel, Accn
      real(wp), intent(in)  :: TIME, DTIME, PARAMS, ADLMAG, PREDEF
      real(wp), intent(in)  :: DDLMAG, PERIOD

      real(wp), intent(out)           :: RHS, AMATRX
      real(wp), intent(out), optional :: SVARS, ENERGY, PNEWDT


      character(len=2)    :: analysis
      character(len=8)    :: abqProcedure
      logical             :: nlgeom
      integer             :: nInt, nPostVars
      integer             :: nDim, nStress
      integer             :: uDOF, uDOFEL, mDOF, mDOFEL
      
      
      integer             :: lenJobName,lenOutDir
      character(len=256)  :: outDir
      character(len=256)  :: jobName
      character(len=512)  :: errFile, dbgFile
      type(logger)        :: msg


      ! initialize primary output variable to be zero
      amatrx        = zero
      rhs(:,nrhs)   = zero
      energy        = zero


      ! open a debug file for the current job
      call getJobName(jobName,lenJobName)
      call getOutDir(outDir,lenOutDir)
      errFile = trim(outDir)//'\aaERR_'//trim(jobName)//'.dat'
      dbgFile = trim(outDir)//'\aaDBG_'//trim(jobName)//'.dat'
      call msg%fopen( errfile=errFile, dbgfile=dbgFile )


      ! change the LFLAGS criteria as needed (check abaqus UEL manual)
      if((lflags(1) .eq. 72) .or. (lflags(1) .eq. 73)) then
        abqProcedure = 'COUPLED'
      else
        call msg%ferror(flag=error, src='UEL',
     &            msg='Incorrect Abaqus procedure.', ia=lflags(1))
        call xit
      end if

      ! check if the procedure is linear or nonlinear
      ! it does not matter except for post processing
      if (lflags(2) .eq. 0) then
        NLGEOM = .false.
      else if (lflags(2) .eq. 1) then
        NLGEOM = .true.
      end if


      ! check to see if it's a general step or a linear purturbation step
      if(lflags(4) .eq. 1) then
        call msg%ferror(flag=error, src='UEL',
     &        msg='The step should be a GENERAL step.', ia=lflags(4))
        call xit
      end if

      ! assign parameter specific to analysis and element types
      ! more anlysis type can be added here
      if ( (JTYPE .eq. 1) .or. (JTYPE .eq. 2) ) then
        analysis  = '3D'         ! three-dimensional analysis
        nDim      = 3
        nStress   = 6
      else if ( (JTYPE .eq. 3) .or. (JTYPE .eq. 4) ) then
        analysis  = 'PE'         ! plane strain analysis
        nDim      = 2
        nStress   = 3
      else
        call msg%ferror(error,'UEL','Element is unavailable.')
        call xit
      end if

      uDOF   = nDim
      uDOFEl = uDOF*nNode
      mDOF   = 1
      mDOFEL = mDOF*nNODE

      nInt      = jprops(1)
      matID     = jprops(2)
      nPostVars = jprops(3)

      ! array containing variables for post-processing
      if (.not. allocated(globalPostVars)) then

        allocate(globalPostVars(numElem,nInt,nPostVars))

        call msg%finfo('---------------------------------------')
        call msg%finfo('--------- ABAQUS HYDROGEL UEL ---------')
        call msg%finfo('---------------------------------------')
        call msg%finfo('--- Abaqus Job: ', ch=trim(jobName))
        call msg%finfo('---------------------------------------')
        call msg%finfo('------- PROCEDURE       = ', ch=abqProcedure)
        call msg%finfo('------- ANALYSIS TYPE   = ', ch=analysis)
        call msg%finfo('---------- NLGEOM       = ', la=nlgeom)
        call msg%finfo('------- MODEL DIMENSION = ', ia=nDim)
        call msg%finfo('------- ELEMENT NODES   = ', ia=nNode)
        call msg%finfo('---------------------------------------')
        call msg%finfo('-------- INTEGRATION SCHEME -----------')
        call msg%finfo('------------ NINT  = ', ia=nInt)
        call msg%finfo('---------------------------------------')
        call msg%finfo('---------- POST-PROCESSING ------------')
        call msg%finfo('--- NO OF ELEMENTS            = ',ia=numElem)
        call msg%finfo('--- OVERLAY ELEMENT OFFSET    = ',ia=elemOffset)
        call msg%finfo('--- NO OF VARIABLES AT INT PT = ',ia=nPostVars)
        call msg%finfo('---------------------------------------')

      end if

      ! return when Abaqus performs dummy step calculation with dt = 0
      if( dtime .eq. zero) return

      
      call uel_hydrogel(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     &    PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     &    DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     &    NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,
     &    PERIOD,NDIM,ANALYSIS,NSTRESS,NINT,UDOF,UDOFEL,MDOF,MDOFEL)

      END SUBROUTINE UEL

! **********************************************************************
! ************** ABAQUS USER OUTPUT VARIABLES SUBROUTINE ***************
! **********************************************************************

       SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     & NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     & JMAC,JMATYP,MATLAYO,LACCFLA)

      ! This subroutine is called by Abaqus at each material point (int pt)
      ! to obtain the user defined output variables for standard Abaqus
      ! elements. We used an additional layer of standard Abaqus elements
      ! with same topology (same number of nodes and int pts) on top of
      ! the user elements with an offset in the numbering between the user
      ! elements and standard elements. This number is defined in the
      ! post_processing module and should match with Abaqus input file.

      use global_parameters
      use post_processing

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

      ! the dimensions of the variables FLGRAY, ARRAY and JARRAY
      ! must be set equal to or greater than 15.

      ! explicityly define the type for uvar to avoid issues
      real(wp)        :: uvar

      ! assign the stored global variables to the UVAR for Abaqus to process
      uvar(1:nuvarm)  = globalPostVars(noel-elemOffset,npt,1:nuvarm)

      END SUBROUTINE UVARM

! **********************************************************************
! **********************************************************************