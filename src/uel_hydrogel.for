! **********************************************************************
! *********** Abaqus/ STANDARD USER ELEMENT SUBROUTINE (UEL) ***********
! **********************************************************************
!           fully coupled chemo-mechanics of non-ionic hydrogels
!   the formulation uses PK-II stress based total Lagrangian framework
! with F-bar modification for fully-integrated HEX8 and QUAD4-PE element
! currently supports first-order elements for 3D, 2D-AX, and 2D-PE cases
! **********************************************************************
!                     BIBEKANANDA DATTA (C) MAY 2024
!                 JOHNS HOPKINS UNIVERSITY, BALTIMORE, MD
! **********************************************************************
!
!                       JTYPE DEFINITION
!
!     U1                THREE-DIMENSIONAL TET4 ELEMENT
!     U2                THREE-DIMENSIONAL HEX8 ELEMENT
!     U3                AXISYMMETRIC TRI3 ELEMENT
!     U4                AXISYMMETRIC QUAD4 ELEMENT
!     U5                PLANE STRAIN TRI3 ELEMENT
!     U6                PLANE STRAIN QUAD4 ELEMENT
!
!       FUTURE TODO: first-order elements for 2D plane stress
!     U7                PLANE STRESS TRI3 ELEMENT
!     U8                PLANE STRESS QUAD4 ELEMENT
!         REMARK: U7 and U8 are not currently available
! **********************************************************************
!          VOIGT NOTATION CONVENTION FOR STRESS/ STRAIN TENSORS
!
!       In this subroutine we adopted the following convention for
!       symmetric stress and strain tensor following Voigt notation
!       This is different than what is followed by Abaqus/ Standard
!
!          sigma11, sigma22, sigma33, sigma23, sigma13, sigma12
!       strain11, strain22, strain33, strain23, strain13, strain12
! **********************************************************************
!                        LIST OF ELEMENT PROPERTIES
!
!     jprops(1)   = nInt            no of integration points in element
!     jprops(2)   = fbarFlag        flag to use F-bar modification or not
!     jprops(3)   = matID           constitutive model for elastomeric network
!     jprops(4)   = nPostVars       no of local (int pt) post-processing variables
! **********************************************************************
!                        POST-PROCESSED VARIABLES
!                     (follows the convention above)
!
!     uvar(1:nStress)               Cauchy stress tensor components
!     uvar(nStress+1:2*nStress)     Euler-Almansi strain tensor components
!     uvar(2*nStress+1)             Polymer volume fraction (phi)
! **********************************************************************
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
!     MDLOAD                        Total number of distributed loads and/or fluxes defined on this element
!     PERIOD                        Time period of the current step
! **********************************************************************
! **********************************************************************

      ! make sure the relative directory is correct
      include 'global_parameters.for'     ! global parameters module
      include 'error_logging.for'         ! error/ debugging module
      include 'linear_algebra.for'        ! linear algebra module
      include 'nonlinear_solver.for'      ! nonlinear solver module
      include 'lagrange_element.for'      ! Lagrange element module
      include 'gauss_quadrature.for'      ! Guassian quadrature module
      include 'surface_integration.for'   ! surface integration module
      include 'solid_mechanics.for'       ! solid mechanics module
      include 'post_processing.for'       ! post-processing module

! **********************************************************************
! **********************************************************************

      module hydrogel_material

      contains

      subroutine neohookean_flory(kstep,kinc,time,dtime,nDim,analysis,
     &          nStress,nNode,jelem,intpt,coord_ip,props,nprops,
     &          jprops,njprops,matID,F,mu,dMudX,
     &          svars,nsvars,fieldVar,dfieldVar,npredf,
     &          stressTensorPK2,dCwdt,Jw,
     &          CTensor,stressTensorCauchy,
     &          FSTensorUM,dCwdFTensor,dJwdFTensor,dCwdMu,dJwdMu,MmatW)

! **********************************************************************
!
!     This material point subroutine calculates constitutive response
!     of an elastomeric hydrogel and returns the flux and tangents.
!     This material subroutine also stores the user-defined element
!     output in a global array for post-processing in Abaqus/Viewer.
!
! **********************************************************************
!                       LIST OF MATERIAL PROPERTIES
!
!     props(1)-props(12): universal parameters and gel properties
!
!     Rgas        = props(1)        Universal gas constant
!     theta       = props(2)        Absolute temperature (K)
!     phi0        = props(3)        Initial polymer volume fraction
!     rho         = props(4)        Density of the gel
!     Gshear      = props(5)        Shear modulus
!     Kappa       = props(6)        Bulk modulus
!     lam_L       = props(7)        Locking stretch (only for AB model)
!     Vp          = props(8)        Molar volume of polymer
!     mu0         = props(9)        Chemical potential of pure solvent
!     Vw          = props(10)       Molar volume of the solvent
!     chi         = props(11)       Flory-Huggins parameter
!     Dw          = props(12)       Diffusion coefficient of the solvent
! **********************************************************************

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
      real(wp), intent(out) :: stressTensorPK2(3,3)
      real(wp), intent(out) :: Jw(nDim,1)
      real(wp), intent(out) :: dCwdt
      real(wp), intent(out) :: CTensor(3,3,3,3)
      real(wp), intent(out) :: stressTensorCauchy(3,3)
      real(wp), intent(out) :: FSTensorUM(3,3)
      real(wp), intent(out) :: dCwdFTensor(3,3)
      real(wp), intent(out) :: dCwdMu, dJwdMu(nDim,1)
      real(wp), intent(out) :: dJwdFTensor(nDim,3,3)
      real(wp), intent(out) :: MmatW(nDim,nDim)

      real(wp), intent(inout), optional   :: svars(nsvars)


      ! local variables (kinematic quantities)
      real(wp)          :: detF, FInv(3,3), FInvT(3,3)
      real(wp)          :: C(3,3), CInv(3,3), detC, trC
      real(wp)          :: B(3,3), Binv(3,3), detB
      real(wp)          :: strainTensorLagrange(3,3)
      real(wp)          :: strainTensorEuler(3,3)

      ! local variables (internal variables)
      real(wp)          :: phi_old, phi, dPhidt, Cw_old, Cw
      real(wp)          :: vars(nprops+2)
      real(wp)          :: detFe, detFs
      logical           :: ivarFlag

      ! local variables (stress tensors)
      real(wp)          :: stressTensorPK1(3,3)
      real(wp)          :: pressure

      ! tangent matrices and related quantities
      real(wp)          :: dCwdPhi, dPhidG, G, dGdPhi
      real(wp)          :: dGdFTensor(3,3)
      real(wp)          :: dGdCTensor(3,3)
      real(wp)          :: dCwdCTensor(3,3)
      real(wp)          :: dSdCwTensor(3,3)
      real(wp)          :: dPdCwTensor(3,3)
      real(wp)          :: dPdFTensor(3,3,3,3)

      ! intermeidate variables for post-processing and output
      real(wp)          :: strainVectLagrange(nSymm,1)
      real(wp)          :: strainVectEuler(nSymm,1)
      real(wp)          :: stressVectPK1(nUnSymm,1)
      real(wp)          :: stressVectCauchy(nSymm,1)

      ! strain and stress vectors for output purposes
      real(wp)          :: strainLagrange(nStress,1)
      real(wp)          :: strainEuler(nStress,1)
      real(wp)          :: stressPK1(nDim*nDim,1)
      real(wp)          :: stressCauchy(nStress,1)

      ! material properties
      real(wp)          :: Rgas, theta, phi0, rho, Gshear, Kappa
      real(wp)          :: lam_L, Vp, mu0, Vw, chi, Dw, RT

      real(wp)          :: lam_c, lam_r, beta_0, beta_c, dBeta_c

      integer           :: i, j, k, l, m, n
      type(options)     :: solverOpts


      ! initialize matrial stiffness tensors
      CTensor     = zero
      dJwdFTensor = zero

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

      solverOpts%maxIter  = 500
      solverOpts%tolfx    = 1.0e-9_wp
      solverOpts%tolx     = 1.0e-9_wp
      solverOpts%algo     = 'Newton'

      call fzero( chemicalState, phi_old, phi, phiMin, phiMax,
     &            jac=.true., vars=vars, opts=solverOpts,
     &            sflag=ivarFlag)


      ! update the state variable for next step (or iteration)
      svars(intPt)    = phi

      ! (1.3) calculate the solution-dependent variables
      detFs     = one/phi
      detFe     = detF/(phi0*detFs)

      ! (1.4) calculate solvent concentration
      Cw_old    = phi0*(one/phi_old-one)/Vw
      Cw        = phi0*(one/phi-one)/Vw


      !!!!!!!!!!!!!!!!!!!!! SOLVE INT PT VARIABLES !!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!! ELEMENT RESIDUAL QUANTITITES !!!!!!!!!!!!!!!!!!
      ! (2.1) stress tensors
      stressTensorPK2   = Gshear * ( ID3 - (phi0)**(two/three) * CInv )
     &                    + Kappa * phi0 * detFs * log(detFe) * CInv

      stressTensorCauchy = (one/detF)
     &                      * ( Gshear * (B - (phi0)**(two/three) * ID3)
     &                      + Kappa * phi0 * detFs * log(detFe) * ID3 )

      ! (2.2) calculate the mean pressure, p = -(Je/3)*trace(sigma)
      pressure  = -(detFe/three) * trace(stressTensorCauchy)

      ! (3) calculate time derivative: dCw/dt
      dCwdt     = (Cw - Cw_old)/dtime


      ! (4.1) mobility matrix : Mw = Dw*Cw/RT*Inv(C) (dimension-dependent)
      MmatW = (Dw*Cw/RT)*CInv(1:nDim,1:nDim)


      ! (4.2) calculate the solvent molar flux: Jw = - Mw*Grad(mu)
      Jw    = - matmul(MmatW,dMUdX)

      !!!!!!!!!!!!!!!! END ELEMENT RESIDUAL QUANTITITES !!!!!!!!!!!!!!!!



      !!!!!!!!!!!!!!!!!!! ELEMENT TANGENT QUANTITITES !!!!!!!!!!!!!!!!!!

      ! (5.1) calculate dCw/dPhi
      dCwdPhi = - (phi0/ (Vw*(phi)**two) )


      ! (5.2) calculate dG/dPhi (same as dMu/dPhi)
      call chemicalState(phi, G, dGdPhi, vars)

      ! (5.3) calculate dPhi/dG (same as dPhi/dMu)
      dPhidG = one/dGdPhi

      ! (5.4) calculate dCw/dMu
      dCwdMU = dCwdPhi*dPhidG


      ! (6.1) calculate dG/dF = dMu/dF
      dGdFTensor  = Kappa * Vw * (log(detFe)-one) * FInvT

      ! (6.2) calculate dG/dC = dMu/dC
      dGdCTensor  = Kappa * Vw * (log(detFe)-one) * CInv/two

      ! (6.3) calculate dCw/dF (using implicit function theorem)
      dCwdFTensor = - dCwdPhi * dPhidG * dGdFTensor

      ! (6.4) calculate dCw/dC (using implicit function theorem)
      dCwdCTensor = - dCwdPhi * dPhidG * dGdCTensor


      ! (7.1) dP/dCw
      dPdCwTensor = Kappa * Vw * (log(detFe)-one) * FInvT

      ! (7.2) dS/dCw (a symmetric second order tensor)
      dSdCwTensor = Kappa * Vw * (log(detFe)-one) * CInv

      ! (8) calculate material tangent (CTensor = 2*dS/dC)
      if (analysis .eq. 'AX') then

        ! for axisymmetry we are defining dP/dF and then calculating
        ! spatial tangent tensor by performing a transformation
        dPdFTensor = zero

        do i=1,3
          do j = 1,3
            do k = 1,3
              do l = 1,3
                dPdFTensor(i,j,k,l) = dPdFTensor(i,j,k,l)
     &          + Gshear * ID3(i,k) * ID3(j,l)
     &          + Gshear * (phi0)**(two/three) * Finv(l,i) * Finv(j,k)
     &          + Kappa * phi0 * detFs * Finv(j,i) * Finv(l,k)
     &          - Kappa*phi0*detFs*log(detFe) * Finv(l,i) * Finv(j,k)
     &          + dPdCwTensor(i,j) * dCwdFTensor(k,l)
              end do
            end do
          end do
        end do

        ! now transforming dp/dF to a_ijkl = 1/J * F_jm * (dP/dF)_imkn * F_ln
        ! follow the appendix of Shawn Chester (IJSS 2015)
        CTensor = zero

        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                do m=1,3
                  do n=1,3
                    CTensor(i,j,k,l) = CTensor(i,j,k,l) +
     &                 ( dPdFTensor(i,m,k,n) * F(j,m) * F(l,n) ) /detF
                  end do
                end do
              end do
            end do
          end do
        end do

      ! for all other analysis we are defining C_IJKL = 2*dS_IJ/dC_KL
      else if ( (analysis .eq. '3D') .or. (analysis .eq. 'PE') ) then
        CTensor   =   zero

        do i = 1,3
          do j = 1,3
            do k = 1,3
              do l = 1,3
                CTensor(i,j,k,l) = CTensor(i,j,k,l)
     &              + Kappa * phi0 * detFs * CInv(i,j) * CInv(k,l)
     &              + ( (phi0)**(two/three) * Gshear
     &              - Kappa * phi0 * detFs * log(detFe) )
     &              * ( CInv(i,k) * CInv(j,l) + CInv(j,k) * CInv(i,l) )
     &              + two * dSdCwTensor(i,j) * dCwdCTensor(k,l)
              end do
            end do
          end do
        end do

      end if


      ! (9) mechano-chemical tangent: F*dSdCwTensor*dCwdMU
      FSTensorUM  = matmul(F,dSdCwTensor) * dCwdMU


      ! (10.1) calculate Jw_tensor
      dJwdFTensor = zero
      do i = 1,nDim
        do k = 1,3
          do l = 1,3
            do j = 1,nDim                 ! summation over dummy index j
              dJwdFTensor(i,k,l) = dJwdFTensor(i,k,l)
     &            + (Dw*Cw)/RT
     &            * ( FInv(i,k)*CInv(l,j) ) * dMUdX(j,1)
     &            - (Dw/RT) * CInv(i,j) * dMUdX(j,1) * dCwdFTensor(k,l)
            end do
          end do
        end do
      end do


      ! (10.2) Calculate dJw/dMU = (dJw/dCw)*(dCw/dMU) (dimension-dependent)
      dJwdMU  = - (Dw/RT) * matmul( CInv(1:nDim,1:nDim),dMUdX ) * dCwdMU


      !!!!!!!!!!!!!!!!!!!!!! TANGENT CALCULATION !!!!!!!!!!!!!!!!!!!!!!!



      !!!!!!!!!!!!!!!!!!!! POST-PROCESSING SECTION !!!!!!!!!!!!!!!!!!!!!

      ! perform reshape and truncation (if needed) for post-processing
      call voigtVector(strainTensorLagrange,strainVectLagrange)
      call voigtVector(strainTensorEuler,strainVectEuler)
      call voigtVector(stressTensorCauchy,stressVectCauchy)

      call voigtVectorTruncate(strainVectLagrange,strainLagrange)
      call voigtVectorTruncate(strainVectEuler,strainEuler)
      call voigtVectorTruncate(stressVectCauchy,stressCauchy)

      ! save the variables to be post-processed in globalPostVars
      ! more variables can be added for element level output
      globalPostVars(jelem,intPt,1:nStress) = stressCauchy(1:nStress,1)
      globalPostVars(jelem,intPt,nStress+1:2*nStress)
     &                                      = strainEuler(1:nStress,1)
      globalPostVars(jElem,intPt,2*nStress+1) = phi

      !!!!!!!!!!!!!!!!!!!! POST-PROCESSING SECTION !!!!!!!!!!!!!!!!!!!!!

! **********************************************************************
! **********************************************************************

      contains

      subroutine chemicalState(phi, G, dG, vars)

! **********************************************************************
!     this subroutine solves the algorithmic internal variable (phi)
!     for Neo-Hookean and Flory-Huggins based elastomeric gel.
!     This subroutine is called at each integration point before
!     constitutive relations and their tangents are computed.
! **********************************************************************

      use global_parameters, only: wp, zero, one, two

      implicit none

      real(wp), intent(in)              :: phi
      real(wp), intent(out)             :: G
      real(wp), intent(out), optional   :: dG
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
      detFe       = detF*phi/phi0

      ! calculate the residual
      G = (mu0 - mu)
     &      + RT*(phi + log(one-phi) + chi*phi**two)
     &      - Kappa*Vw * log(detF*phi/phi0)
     &      +  ( (Kappa*Vw)/two ) * ( log(detF*phi/phi0) )**two

        ! analytical gradient of the residual
        if ( present(dG) ) then
          dG = RT * ( one - ( one/(one - phi) ) + two*chi*phi )
     &            + ( (Kappa*Vw)/phi ) * log(detFe) - (Kappa*Vw)/phi
        endif

      end subroutine chemicalState

      end subroutine neohookean_flory


      end module hydrogel_material

! **********************************************************************
! **********************************************************************
! **********************************************************************

      module hydrogel_element

! **********************************************************************
!     This module contains subroutines related to element formulation
!     Abaqus user subroutines can not be included in a module.
!     Instead we extended the list of arguments of the Abaqus UEL subroutine
!     and wrote two other subroutines contained in this module.
!     Compilers can perform additional checks on the arguments when
!     any modularized subroutines are called. The first subroutine is
!     called by UEL subroutine of Abaqus with an extended set of
!     input arguments. The first subroutine calls other subroutines.
! **********************************************************************

      contains

      SUBROUTINE gel_general(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,
     & NSVARS,PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     & TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD,
     & NDIM,ANALYSIS,NSTRESS,NINT,UDOF,UDOFEL,MDOF,MDOFEL)

      ! this subroutine computes the AMATRX and RHS for 3D elements and
      ! 2D plane strain elements. axisymmetric elements are placed in a
      ! seperate subroutine as that requires different matrix operators 
      ! and tensor conversion

      use global_parameters
      use error_logging
      use lagrange_element
      use gauss_quadrature
      use surface_integration
      use solid_mechanics
      use linear_algebra
      use hydrogel_material

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

      ! nodal degrees of freedom
      real(wp)          :: UallMat(uDOF+mDOF,nNode)
      real(wp)          :: DUallMat(uDOF+mDOF,nNODE)
      real(wp)          :: uNode(nDim,nNode), muNode(mDOFEL,1)
      real(wp)          :: duNode(nDim,nNode), dmuNode(mDOFEL,1)

      ! additional field variables at the nodes and integration point
      real(wp)          :: fieldNode(npredf,nNode)
      real(wp)          :: dfieldNode(npredf,nNode)


      ! finite element parameters (integration and shape functions)
      real(wp)          :: w(nInt), xi(nInt,nDim)
      real(wp)          :: Nxi(nNode), dNdxi(nNode,nDim)
      real(wp)          :: dXdxi(nDim,nDim), dxidX(nDim,nDim)
      real(wp)          :: detJ, dNdX(nNode,nDim)

      ! finite element matrix operators
      real(wp)          :: Na(nDim,nDim)
      real(wp)          :: Ba(nStress,nDim)
      real(wp)          :: Ga(nDim*nDim,nDim)
      real(wp)          :: Nmat(nDim,uDOFEl)
      real(wp)          :: NmatT(uDOFEL,nDim)
      real(wp)          :: Bmat(nStress,uDOFEl)
      real(wp)          :: BmatT(uDOFEL,nStress)
      real(wp)          :: Gmat(nDim*nDim,nDim*nNode)
      real(wp)          :: GmatT(nDim*nNode,nDim*nDim)
      real(wp)          :: NmatScalar(1,nNode)
      real(wp)          :: NmatScalarT(nNode,1)
      real(wp)          :: BmatScalar(nDim,nNode)
      real(wp)          :: BmatScalarT(nNode,nDim)


      ! additional variables for F-bar method (element and material)
      logical           :: fbarFlag
      real(wp)          :: centroid(nDim)
      real(wp)          :: Nxi0(nNode), dNdxi0(nNode,nDim)
      real(wp)          :: dXdxi0(nDim,nDim), dxidX0(nDim,nDim)
      real(wp)          :: dNdX0(nNode,nDim), detJ0
      real(wp)          :: Ga0(nDim*nDim,nDim)
      real(wp)          :: Gmat0(nDim*nDim,uDOFEl)
      real(wp)          :: F0(3,3), detF0
      real(wp)          :: F0Inv(3,3), F0InvT(3,3)
      real(wp)          :: QR0Tensor(nDim,nDim,nDim,nDim)
      real(wp)          :: QRTensor(nDim,nDim,nDim,nDim)
      real(wp)          :: QR0mat(nDim*nDim,nDim*nDim)
      real(wp)          :: QRmat(nDim*nDim,nDim*nDim)
      real(wp)          :: tanFac1, tanFac2, resFac


      ! integration point quantities (variables)
      real(wp)          :: coord_ip(nDim,1)
      real(wp)          :: F(3,3), detF, Fbar(3,3)
      real(wp)          :: FInv(3,3), FInvT(3,3)
      real(wp)          :: mu, dMUdX(nDim,1)
      real(wp)          :: fieldVar(npredf), dfieldVar(npredf)


      ! material model identifier and state variables
      integer               :: matID      ! material constitutive law
      real(wp), allocatable :: statev(:)  ! state variables per int pt


      ! constitutive output from the material point subroutine (UMAT)
      real(wp)          :: stressTensorPK2(3,3)
      real(wp)          :: dCwdt, Jw(nDim,1)
      real(wp)          :: CTensor(3,3,3,3)
      real(wp)          :: stressTensorCauchy(3,3)
      real(wp)          :: FSTensorUM(3,3)
      real(wp)          :: dJwdFTensor(nDim,3,3)
      real(wp)          :: dCwdFTensor(3,3)
      real(wp)          :: dCwdMu, dJwdMu(nDim,1)
      real(wp)          :: MmatW(nDim,nDim)

      ! additional matrix operator for element formulation
      real(wp)          :: VoigtMat(nSymm,nSymm)
      real(wp)          :: Dmat(nStress,nStress)
      real(wp)          :: stressVectPK2(nSymm,1)
      real(wp)          :: stressPK2(nStress,1)
      real(wp)          :: SIGMA_S(nDim*nDim,nDim*nDim)
      real(wp)          :: SIGMA_F(nDim*nNode,nDim*nNode)
      real(wp)          :: BNLmat(nStress,nDim*nNode)
      real(wp)          :: BNLmatT(nDim*nNode,nStress)
      real(wp)          :: aVectUM(nDim*nDim,1)
      real(wp)          :: DmatMU(nDim,nDim*nDim)
      real(wp)          :: dCwdotdFTensor(3,3)
      real(wp)          :: dCwdotdF(1,nDim*nDim)
      real(wp)          :: dCwdotdMu


      ! element residual vectors and tangent matrix components
      real(wp)          :: Ru(uDOFEL,1)
      real(wp)          :: Rm(mDOFEL,1)

      real(wp)          :: Kuu(uDOFEl,uDOFEl)
      real(wp)          :: Kum(uDOFEl,mDOFEL)
      real(wp)          :: Kmu(mDOFEL,uDOFEl)
      real(wp)          :: Kmm(mDOFEL,mDOFEL)

      real(wp)          :: Kelem(nDOFEL,nDOFEL)
      real(wp)          :: Relem(nDOFEL,1)

      integer           :: i, j, k, l, m, n, p, q, intPt
      integer           :: nstatev
      type(element)     :: hydrogel

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

      matID     = jprops(3)
      nstatev   = nsvars/nint       ! state variables per integration point

      ! set the F-bar flag based on the input
      if (jprops(2) .eq. 0) then
        fbarFlag = .false.
      else
        fbarFlag = .true.
      end if


      if ( .not. allocated(statev) ) allocate(statev(nstatev) )

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


      ! For fully-integrated QUAD4 and HEX8 element, calculate Gmat0.
      ! These calculations are done to evaluate volumetric deformation
      ! gradient at the element centroid to calculate F-bar.
      if (fbarFlag .eq. .true.) then

        if ( ((jtype .eq. 2) .and. (nInt .eq. 8))
     &      .or. ((jtype .eq. 6) .and. (nInt .eq. 4)) ) then

          centroid = zero

          ! evaluate the interpolation functions and derivates at centroid
          call calcInterpFunc(hydrogel, centroid, Nxi0, dNdxi0)

          ! calculate element jacobian and global shape func gradient at centroid
          dXdxi0  = matmul(coords,dNdxi0)       ! calculate the jacobian (dXdxi) at centroid
          detJ0   = det(dXdxi0)                 ! calculate jacobian determinant at centroid

          if (detJ0 .le. zero) then
            call msg%ferror( flag=warn, src='gel_general',
     &      msg='Negative element jacobian at centroid: ', ia=jelem)
          call xit
          end if

          dxidX0 = inv(dXdxi0)                  ! calculate jacobian inverse
          dNdX0  = matmul(dNdxi0,dxidX0)        ! calculate dNdX at centroid

          do i=1,nNode

            ! form the nodal-level matrix: [Ga0] at the centroid
            do j = 1, nDim
              Ga0(nDim*(j-1)+1:nDim*j, 1:nDim) = dNdX0(i,j)*ID
            end do

            ! form the [G0] matrix at the centroid
            Gmat0(1:nDim**2,nDim*(i-1)+1:nDim*i) = Ga0(1:nDim**2,1:nDim)
          end do                             ! end of nodal point loop


          F0(1:nDim,1:nDim) = ID + matmul(uNode,dNdX0)

          if (analysis .eq. 'PE') F0(3,3) = one

          detF0   = det(F0)
          F0Inv   = inv(F0)
          F0InvT  = transpose(F0Inv)

        else
          call msg%ferror( flag=warn, src='gel_general',
     &      msg='F-bar is not available: ', ivec=[jtype, nInt])
          call xit
        end if

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
          call msg%ferror( flag=warn, src='gel_general',
     &          msg='Negative element jacobian: ', ivec=[jelem, intpt])
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
            call msg%ferror( flag=error, src='gel_general',
     &                msg='Wrong analysis: ', ch=analysis )
            call xit
          end if

          ! form the [N], [B], and [G] matrices
          Nmat(1:nDim,nDim*(i-1)+1:nDim*i)    = Na(1:nDim,1:nDim)
          Bmat(1:nStress,nDim*(i-1)+1:nDim*i) = Ba(1:nStress,1:nDim)
          Gmat(1:nDim**2,nDim*(i-1)+1:nDim*i) = Ga(1:nDim**2,1:nDim)
        end do                             ! end of nodal point loop


        ! transpose the vector field matrix operators
        NmatT       = transpose(Nmat)
        BmatT       = transpose(Bmat)
        GmatT       = transpose(Gmat)

        ! all the scalar matrix operators for the element
        NmatScalar  = reshape(Nxi, [1, nNode])
        NmatScalarT = transpose(NmatScalar)
        BmatScalar  = transpose(dNdX)
        BmatScalarT = dNdX


        !!!!!!!!!!!!! END CALCULATING ELEMENT OPERATORS !!!!!!!!!!!!!!!


        !!!!!!!!!!!!!!!!!! CONSTITUTIVE CALCULATION !!!!!!!!!!!!!!!!!!!

        ! calculate the coordinate of integration point
        coord_ip = matmul(Nmat, reshape(coords, [nDOFEL, 1]))


        ! calculate deformation gradient and deformation tensors
        F(1:nDim,1:nDim)  = ID + matmul(uNode,dNdX)

        if (analysis .eq. 'PE')  F(3,3) = one

        ! calculate jacobian (volume change) at the current integration pt
        detF    = det(F)
        FInv    = inv(F)
        FInvT   = transpose(FInv)

        ! calculate chemical potential and its gradient
        mu    = dot_product( Nxi, reshape(muNode, [mDOFEl]) )

        dMudX = matmul( BmatScalar, muNode )

        !! definition of modified deformation gradient, F-bar
        if (fbarFlag .eq. .true.) then
          if ( (jtype .eq. 2) .and. (nInt .eq. 8) ) then
            ! fully-integrated HEX8 element
            Fbar    = (detF0/detF)**(third) * F
            resFac  = (detF0/detF)**(-two/three)
            tanFac1 = (detF0/detF)**(-one/three)
            tanFac2 = (detF0/detF)**(-two/three)

          else if ( (jtype .eq. 6) .and. (nInt .eq. 4) )  then
            ! fully-integrated QUAD4-PE element
            Fbar(3,3)           = one
            Fbar(1:nDim,1:nDim) = (detF0/detF)**(half)*F(1:nDim,1:nDim)
            resFac              = (detF0/detF)**(-half)
            tanFac1             = one
            tanFac2             = (detF0/detF)**(-half)
          else
            ! standard F for all other available elements
            Fbar    = F
            resFac  = one
            tanFac1 = one
            tanFac2 = one
            call msg%ferror( flag=warn, src='gel_general',
     &      msg='F-bar is not available: ', ivec=[jtype, nInt])
          call xit
          end if
        else
          ! set F-bar = F if fbarFlag is .false. for all element
          Fbar    = F
          resFac  = one
          tanFac1 = one
          tanFac2 = one
        end if


        ! call material point subroutine (UMAT) for the neutral gel
        select case(matID)
        case (1)
          call neohookean_flory(kstep,kinc,time,dtime,nDim,analysis,
     &          nStress,nNode,jelem,intpt,coord_ip,props,nprops,
     &          jprops,njprops,matID,Fbar,mu,dMudX,
     &          svars,nsvars,fieldVar,dfieldVar,npredf,
     &          stressTensorPK2,dCwdt,Jw,
     &          CTensor,stressTensorCauchy,
     &          FSTensorUM,dCwdFTensor,dJwdFTensor,dCwdMu,dJwdMu,MmatW)

        case default
          call msg%ferror(flag=error, src='gel_general',
     &              msg='Material model is not available', ia=matID)
          call xit
      end select

        !!!!!!!!!!!!!!! END CONSTITUTIVE CALCULATION !!!!!!!!!!!!!!!!!!


        !!!!!!!!!!!!! FORM ADDITIONAL ELEMENT OPERATORS !!!!!!!!!!!!!!!

        call voigtMatrix(CTensor,VoigtMat)
        call voigtMatrixTruncate(VoigtMat,Dmat)

        call voigtVector(stressTensorPK2,stressVectPK2)
        call voigtVectorTruncate(stressVectPK2,stressPK2)


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

        BNLmat  = matmul(Bmat,transpose(SIGMA_F))
        BNLmatT = transpose(BNLmat)


        ! reshape FSTensorUM into vector form
        aVectUM = reshape( FSTensorUM(1:nDim,1:nDim),  [nDim*nDim,1] )

        ! map the third-order tensor, dJw/dF, to a rank-2 matrix
        do i = 1, nDim
          do l = 1, nDim
            do k = 1, nDim
              DmatMU(i,(l-1)*nDim+k) = dJwdFTensor(i,k,l)
            end do
          end do
        end do

        ! calculate dCw_dot/dMu (scalar)
        dCwdotdMU       = dCwdMu/dtime

        ! calculate of dCw_dot/dF and reshape it (dimension-dependent)
        dCwdotdFTensor  = dCwdFTensor/dtime
        dCwdotdF        = reshape( dCwdotdFTensor(1:nDim,1:nDim),
     &                          [1,nDim*nDim] )

        !!!!!!!!!! END FORMING ADDITIONAL ELEMENT OPERATORS !!!!!!!!!!!



        !!!!!!!!!!!!!!!! RESIDUAL VECTOR CALCULATION !!!!!!!!!!!!!!!!!!


        ! mechanical residual
        Ru = Ru - w(intPt) * detJ * resFac *
     &            matmul( transpose(BNLmat), stressPK2 )

        ! chemical residual
        Rm = Rm + w(intPt) * detJ *
     &           (
     &            - NmatScalarT * dCwdt + matmul( BmatScalarT, Jw )
     &           )

        !!!!!!!!!!!!!! END RESIDUAL VECTOR CALCULATION !!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!! TANGENT MATRIX CALCULATION !!!!!!!!!!!!!!!!!!!

        ! mechanical tangent matrix
        Kuu = Kuu + w(intPt) * detJ * tanFac1 *
     &      (
     &        matmul( matmul( GmatT, SIGMA_S ), Gmat )
     &        + matmul( BNLmatT, matmul( Dmat,BNLmat ) )
     &      )


        ! mechanical-solvent tangent matrix
        Kum = Kum + w(intpt) * detJ * tanFac2 *
     &        matmul( matmul( GmatT, aVectUM ), NmatScalar )


        ! solvent-mechanical tangent matrix
        Kmu = Kmu + w(intPt) * detJ *
     &        (
     &        matmul( matmul( NmatScalarT, dCwdotdF ), Gmat)
     &        - matmul( matmul( BmatScalarT, DmatMU ), Gmat )
     &        )


        ! solvent tangent matrix
        Kmm = Kmm + w(intPt) * detJ *
     &        (
     &        matmul( NmatScalarT, NmatScalar ) * dCwdotdMu
     &        - matmul( matmul(BmatScalarT, dJwdMu), NmatScalar )
     &        + matmul( matmul(BmatScalarT, MmatW), BmatScalar )
     &        )


        !! F-bar modification block for Kuu
        if (fbarFlag .eq. .true.) then

          ! form fourth-order QR0 and QR tensor
          QR0Tensor = zero
          QRTensor  = zero

          !! fully-integrated HEX8 element
          if ( (jtype .eq. 2) .and. (nInt .eq. 8) ) then

            do i = 1,nDim
              do j = 1,nDim
                do k = 1,nDim
                  do l = 1,nDim
                    do m = 1,nDim
                      do n = 1,nDim
                        do p = 1,nDim
                          do q = 1,nDim
                            QR0Tensor(i,j,k,l) = QR0Tensor(i,j,k,l)
     &                          + third * F0InvT(k,l) *
     &                            (
     &                              Fbar(i,p) * CTensor(p,j,m,n)
     &                              * Fbar(q,m) * Fbar(q,n)
     &                              - Fbar(i,q) * stressTensorPK2(q,j)
     &                            )

                            QRTensor(i,j,k,l) = QRTensor(i,j,k,l)
     &                          + third * FInvT(k,l) *
     &                            (
     &                              Fbar(i,p) * CTensor(p,j,m,n)
     &                              * Fbar(q,m) * Fbar(q,n)
     &                              - Fbar(i,q) * stressTensorPK2(q,j)
     &                            )
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do

          ! fully-integrated QUAD4-PE element
          else if ( (jtype .eq. 6) .and. (nInt .eq. 4) ) then

            do i = 1,nDim
              do j = 1,nDim
                do k = 1,nDim
                  do l = 1,nDim
                    do m = 1,nDim
                      do n = 1,nDim
                        do p = 1,nDim
                          do q = 1,nDim
                            QR0Tensor(i,j,k,l) = QR0Tensor(i,j,k,l)
     &                          + half * Fbar(i,p) * CTensor(p,j,m,n)
     &                          * Fbar(q,m) * Fbar(q,n) * F0InvT(k,l)

                            QRTensor(i,j,k,l) = QRTensor(i,j,k,l)
     &                          + half * Fbar(i,p) * CTensor(p,j,m,n)
     &                          * Fbar(q,m) * Fbar(q,n) * FInvT(k,l)
                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do

          end if

          ! reshape QR and QR0 tensor into matrix form
          call unsymmMatrix(QR0Tensor,QR0mat)
          call unsymmMatrix(QRTensor,QRmat)

          ! modify the element tangent matrix
          Kuu = Kuu + w(intPt) * detJ * tanFac2  *
     &              (
     &                matmul( GmatT, matmul(QR0mat,Gmat0) )
     &                - matmul( GmatT, matmul(QRmat,Gmat) )
     &              )
        end if

        !!!!!!!!!!!!!! END TANGENT MATRIX CALCULATION !!!!!!!!!!!!!!!!!

      end do                           ! end of integration point loop



      ! assemble the element tangent sub-matrices and residual sub-vectors
      call assembleElement(nDim,nNode,uDOFEL,mDOFEL,nDOFEl,       ! input to the subroutine
     &     Kuu,Kum,Kmu,Kmm,Ru,Rm,                                 ! input to the subroutine
     &     Kelem,Relem)                                           ! output from the subroutine


      ! assign the element tangent matrices to Abaqus defined variables
      amatrx(1:NDOFEL,1:NDOFEL)   = Kelem(1:NDOFEL,1:NDOFEL)
      rhs(1:NDOFEL,1)             = Relem(1:NDOFEL,1)


      end subroutine gel_general

! **********************************************************************
! **********************************************************************

      SUBROUTINE gel_axisymmetric(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,
     & NSVARS,PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,
     & TIME,DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD,
     & NDIM,ANALYSIS,NSTRESS,NINT,UDOF,UDOFEL,MDOF,MDOFEL)

      ! this subroutine computes the AMATRX and RHS for axisymmetric elements

      use global_parameters
      use error_logging
      use lagrange_element
      use gauss_quadrature
      use surface_integration
      use solid_mechanics
      use linear_algebra
      use hydrogel_material

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

      ! nodal degrees of freedom
      real(wp)          :: UallMat(uDOF+mDOF,nNode)
      real(wp)          :: DUallMat(uDOF+mDOF,nNODE)
      real(wp)          :: uNode(nDim,nNode), muNode(mDOFEL,1)
      real(wp)          :: duNode(nDim,nNode), dmuNode(mDOFEL,1)

      ! additional field variables at the nodes and integration point
      real(wp)          :: fieldNode(npredf,nNode)
      real(wp)          :: dfieldNode(npredf,nNode)


      ! finite element parameters (integration and shape functions)
      real(wp)          :: w(nInt), xi(nInt,nDim)
      real(wp)          :: Nxi(nNode), dNdxi(nNode,nDim)
      real(wp)          :: dXdxi(nDim,nDim), dxidX(nDim,nDim)
      real(wp)          :: detJ, dNdX(nNode,nDim)

      ! finite element matrix operators
      real(wp)          :: Na(nDim,nDim)
      real(wp)          :: Ba(nStress,nDim)
      real(wp)          :: Ga(nDim*nDim+1,nDim)
      real(wp)          :: Nmat(nDim,uDOFEl)
      real(wp)          :: NmatT(uDOFEL,nDim)
      real(wp)          :: Bmat(nStress,uDOFEl)
      real(wp)          :: BmatT(uDOFEL,nStress)
      real(wp)          :: Gmat(nDim*nDim+1,nDim*nNode)
      real(wp)          :: GmatT(nDim*nNode,nDim*nDim+1)
      real(wp)          :: NmatScalar(1,nNode)
      real(wp)          :: NmatScalarT(nNode,1)
      real(wp)          :: BmatScalar(nDim,nNode)
      real(wp)          :: BmatScalarT(nNode,nDim)

      ! integration point quantities (variables)
      real(wp)          :: coord_ip(nDim,1)
      real(wp)          :: F(3,3), detF, Fbar(3,3)
      real(wp)          :: FInv(3,3), FInvT(3,3)
      real(wp)          :: mu, dMUdX(nDim,1)
      real(wp)          :: fieldVar(npredf), dfieldVar(npredf)


      ! current coordinate variables
      real(wp)          :: coords_t(ndim,nnode)
      real(wp)          :: R_t, AR_t, R, AR
      real(wp)          :: dxdxi_t(nDim,nDim), dxidx_t(nDim,nDim)
      real(wp)          :: detJ_t, dNdX_t(nNode,nDim)
      real(wp)          :: Nmat_t(nDim,uDOFEl)
      real(wp)          :: NmatT_t(uDOFEL,nDim)
      real(wp)          :: Bmat_t(nStress,uDOFEl)
      real(wp)          :: BmatT_t(uDOFEL,nStress)
      real(wp)          :: Gmat_t(nDim*nDim+1,nDim*nNode)
      real(wp)          :: GmatT_t(nDim*nNode,nDim*nDim+1)
      real(wp)          :: NmatScalar_t(1,nNode)
      real(wp)          :: NmatScalarT_t(nNode,1)
      real(wp)          :: BmatScalar_t(nDim,nNode)
      real(wp)          :: BmatScalarT_t(nNode,nDim)
      real(wp)          :: dMudx_t(nDim,1)


      ! material model identifier and state variables
      integer               :: matID      ! material constitutive law
      real(wp), allocatable :: statev(:)  ! state variables per int pt


      ! constitutive output from the material point subroutine (UMAT)
      real(wp)          :: stressTensorPK2(3,3)
      real(wp)          :: dCwdt, Jw(nDim,1)
      real(wp)          :: CTensor(3,3,3,3)
      real(wp)          :: stressTensorCauchy(3,3)
      real(wp)          :: FSTensorUM(3,3)
      real(wp)          :: dJwdFTensor(nDim,3,3)
      real(wp)          :: dCwdFTensor(3,3)
      real(wp)          :: dCwdMu, dJwdMu(nDim,1)
      real(wp)          :: MmatW(nDim,nDim)


      ! additional matrix operator for element formulation
      real(wp)          :: stressCauchy(4,1)
      real(wp)          :: stressPK1(5,1)
      real(wp)          :: ATensor(3,3,3,3)
      real(wp)          :: Amat(5,5)
      real(wp)          :: aVectUM(nDim*nDim+1,1)
      real(wp)          :: DmatMU(nDim,nDim*nDim+1)
      real(wp)          :: dCwdotdFTensor(3,3)
      real(wp)          :: dCwdotdF(1,nDim*nDim+1)
      real(wp)          :: dCwdotdMu

      ! additional variables for F-bar method (element and material)
      logical           :: fbarFlag
      real(wp)          :: centroid(nDim)
      real(wp)          :: Nxi0(nNode), dNdxi0(nNode,nDim)
      real(wp)          :: dXdxi0(nDim,nDim), dxidX0(nDim,nDim)
      real(wp)          :: dNdX0(nNode,nDim), detJ0
      real(wp)          :: dXdxi0_t(nDim,nDim), dxidX0_t(nDim,nDim)
      real(wp)          :: dNdX0_t(nNode,nDim), detJ0_t
      real(wp)          :: Ga0(nDim*nDim+1,nDim)
      real(wp)          :: Gmat0(nDim*nDim+1,uDOFEl)
      real(wp)          :: Gmat0_t(nDim*nDim+1,uDOFEl)
      real(wp)          :: Gmat0T(uDOFEl,nDim*nDim+1)
      real(wp)          :: Gmat0T_t(uDOFEl,nDim*nDim+1)
      real(wp)          :: R0, r0_t
      real(wp)          :: F0(3,3), detF0
      real(wp)          :: F0Inv(3,3), F0InvT(3,3)
      real(wp)          :: Qmat(nDim*nDim+1,nDim*nDim+1)
      real(wp)          :: tanFac1, tanFac2, resFac


      ! element residual vectors and tangent matrix components
      real(wp)          :: Ru(uDOFEL,1)
      real(wp)          :: Rm(mDOFEL,1)

      real(wp)          :: Kuu(uDOFEl,uDOFEl)
      real(wp)          :: Kum(uDOFEl,mDOFEL)
      real(wp)          :: Kmu(mDOFEL,uDOFEl)
      real(wp)          :: Kmm(mDOFEL,mDOFEL)

      real(wp)          :: Kelem(nDOFEL,nDOFEL)
      real(wp)          :: Relem(nDOFEL,1)

      integer           :: i, j, k, l, m, n, p, q, intPt
      integer           :: nstatev
      type(element)     :: hydrogel

      ! initialize hydrogel element
      hydrogel  = element(nDim=nDim, analysis=analysis,
     &                    nNode=nNode, nInt=nInt)

      F0      = zero
      Fbar    = zero
      Ga0     = zero
      Gmat0   = zero
      Gmat0_t = zero
      F       = zero
      Na      = zero
      Ba      = zero
      Ga      = zero
      Nmat    = zero
      Bmat    = zero
      Gmat    = zero
      Kuu     = zero
      Kum     = zero
      Kmu     = zero
      Kmm     = zero
      Ru      = zero
      Rm      = zero
      Kelem   = zero
      Relem   = zero

      matID     = jprops(3)
      nstatev   = nsvars/nint       ! state variables per integration point

      ! set the F-bar flag based on the input
      if (jprops(2) .eq. 0) then
        fbarFlag = .false.
      else
        fbarFlag = .true.
      end if


      if ( .not. allocated(statev) ) allocate(statev(nstatev) )

      call eyeMat(ID)                       ! create an identity matrix w/ size of nDim

      ! reshape the degrees of freedom and their change into matrices
      uAllMat   = reshape( uAll, shape=[uDOF+mDOF,nNode] )
      duAllMat  = reshape( duAll(:,1), shape=[uDOF+mDOF,nNode] )

      ! seperate the displacement and chemical potential field
      uNode(1:uDOF,1:nNode)   = uAllMat(1:uDOF,1:nNode)
      muNode(1:mDOFEL,1)      = uAllMat(uDOF+1,1:mDOFEL)

      ! seperate the change of displacement and chemical potential field
      duNode(1:uDOF,1:nNode)  = duAllMat(1:uDOF,1:nNode)
      dmuNode(1:mDOFEL,1)     = duAllMat(uDOF+1,1:mDOFEL)


      ! calculate the current/ deformed coordinate
      coords_t  = coords + uNode



      !!!!!!!!!!!!!!!!!!! CENTROID LEVEL CALCULATION !!!!!!!!!!!!!!!!!!

      ! For fully-integrated QUAD4 and HEX8 element, calculate Gmat0.
      ! These calculations are done to evaluate volumetric deformation
      ! gradient at the element centroid to calculate F-bar.
      if (fbarFlag .eq. .true.) then

        ! fully-integrated QUAD4-AX element only
        if ( (jtype .eq. 4) .and. (nInt .eq. 4) ) then

          centroid = zero

          ! evaluate the interpolation functions and derivates at centroid
          call calcInterpFunc(hydrogel, centroid, Nxi0, dNdxi0)

          ! calculate element jacobian and global shape func gradient at centroid
          dXdxi0  = matmul(coords,dNdxi0)       ! calculate the jacobian (dXdxi) at centroid
          detJ0   = det(dXdxi0)                 ! calculate jacobian determinant at centroid
          dxidX0  = inv(dXdxi0)                 ! calculate jacobian inverse
          dNdX0   = matmul(dNdxi0,dxidX0)       ! calculate dNdX0 at centroid

          if (detJ0 .le. zero) then
            call msg%ferror( flag=warn, src='gel_axisymmetric',
     &      msg='Negative element jacobian at centroid: ', ia=jelem)
            pnewdt = fourth
            return
          end if


          ! shape functions and their gradients in current coordinate
          dxdxi0_t   = matmul(coords_t,dNdxi0)        ! calculate dxdxi
          detJ0_t    = det(dxdxi0_t)                  ! calculate determinant
          dxidx0_t   = inv(dxdxi0_t)                  ! calculate inverse
          dNdx0_t    = matmul(dNdxi0,dxidx0_t)        ! calculate dNdX0_t

          if (detJ0_t .lt. zero) then
            call msg%ferror( flag=warn, src='gel_axisymmetric',
     &          msg='Negative element jacobian: ', ivec=[jelem, intpt])
            pnewdt = fourth
            return
          end if

          ! calculate the centroid radius and circumference (current and old)
          R0    = dot_product( Nxi0, coords(1,:) )
          r0_t  = dot_product( Nxi0, coords_t(1,:) )


          do i=1,nNode

            ! form the nodal-level matrix: [Ga0] at the centroid
            do j = 1, nDim
              Ga0(nDim*(j-1)+1:nDim*j, 1:nDim) = dNdx0_t(i,j)*ID
            end do
            Ga0(nDim**2+1,1)   = Nxi0(i)/r0_t

            ! form the [G0_t] matrix at the centroid
            Gmat0_t(1:nDim**2+1,nDim*(i-1)+1:nDim*i)
     &                        = Ga0(1:nDim**2+1,1:nDim)
          end do                             ! end of nodal point loop


          Gmat0T_t  = transpose(Gmat0_t)

          F0                = zero
          F0(1:nDim,1:nDim) = ID + matmul(uNode,dNdX0)
          F0(3,3)           = r0_t/R0


          detF0   = det(F0)
          F0Inv   = inv(F0)
          F0InvT  = transpose(F0Inv)

        else
          call msg%ferror( flag=error, src='gel_axisymmetric',
     &      msg='F-bar is not available: ', ivec=[jtype, nInt])
          call xit
        end if

      end if

      !!!!!!!!!!!!!!!!! END CENTROID LEVEL CALCULATION !!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!! INTEGRATION POINT LOOP !!!!!!!!!!!!!!!!!!!!!

      ! obtain the standard gauss quadrature points and weights
      call getGaussQuadrtr(hydrogel,w,xi)


      do intPt = 1, nInt

        call calcInterpFunc(hydrogel,xi(intPt,:),Nxi,dNdxi)

        ! shape functions and their gradient in undeformed coordinate
        dXdxi = matmul(coords,dNdxi)        ! calculate dxdxi
        detJ  = det(dXdxi)                  ! calculate determinant

        if (detJ .lt. zero) then
          call msg%ferror( flag=warn, src='gel_axisymmetric',
     &          msg='Negative element jacobian: ', ivec=[jelem, intpt])
        end if

        dxidX = inv(dXdxi)                  ! calculate inverse
        dNdX  = matmul(dNdxi,dxidX)         ! calculate dNdX



        ! shape functions and their gradient in current coordinate
        dxdxi_t   = matmul(coords_t,dNdxi)        ! calculate dxdxi
        detJ_t    = det(dxdxi_t)                  ! calculate determinant

        if (detJ_t .lt. zero) then
          call msg%ferror( flag=warn, src='gel_general',
     &          msg='Negative element jacobian: ', ivec=[jelem, intpt])
        end if

        dxidx_t   = inv(dxdxi_t)                  ! calculate inverse
        dNdx_t    = matmul(dNdxi,dxidx_t)         ! calculate dNdX


        ! calculate the centroid radius and circumference (current and old)
        R     = dot_product( Nxi, coords(1,:) )
        r_t   = dot_product( Nxi, coords_t(1,:) )
        AR    = two * pi * R
        Ar_t  = two * pi * r_t


        ! loop over all the nodes (internal loop)
        do i=1, nNode

          ! form the matrix operator in reference coordinate
          do j = 1, nDim
            Na(j,j) = Nxi(i)
            Ga(nDim*(j-1)+1:nDim*j,1:nDim) = dNdX(i,j)*ID
          end do
            Ga(nDim**2+1,1)   = Nxi(i)/R

            Ba(1,1)           = dNdX(i,1)
            Ba(2,2)           = dNdX(i,2)
            Ba(3,1:nDim)      = [dNdX(i,2), dNdX(i,1)]
            Ba(4,1)           = Nxi(i)/R

          ! form the [N], [B], and [G] matrices
          Nmat(1:nDim,nDim*(i-1)+1:nDim*i)      = Na(1:nDim,1:nDim)
          Bmat(1:nStress,nDim*(i-1)+1:nDim*i)   = Ba(1:nStress,1:nDim)
          Gmat(1:nDim**2+1,nDim*(i-1)+1:nDim*i) = Ga(1:nDim**2+1,1:nDim)
        end do

        ! transpose the vector field matrix operators
        NmatT       = transpose(Nmat)
        BmatT       = transpose(Bmat)
        GmatT       = transpose(Gmat)

        ! all the scalar matrix operators for the element
        NmatScalar  = reshape(Nxi, [1, nNode])
        NmatScalarT = transpose(NmatScalar)
        BmatScalar  = transpose(dNdX)
        BmatScalarT = dNdX


        ! loop over all the nodes (internal loop)
        do i=1, nNode

          ! form the matrix operator in current coordinate
          do j = 1, nDim
            Na(j,j) = Nxi(i)
            Ga(nDim*(j-1)+1:nDim*j,1:nDim) = dNdX_t(i,j)*ID
          end do
            Ga(nDim**2+1,1)   = Nxi(i)/r_t

            Ba(1,1)           = dNdX_t(i,1)
            Ba(2,2)           = dNdX_t(i,2)
            Ba(3,1:nDim)      = [dNdX_t(i,2), dNdX_t(i,1)]
            Ba(4,1)           = Nxi(i)/r_t

          ! form the [N], [B], and [G] matrices
          Nmat_t(1:nDim,nDim*(i-1)+1:nDim*i)      = Na(1:nDim,1:nDim)
          Bmat_t(1:nStress,nDim*(i-1)+1:nDim*i)   = Ba(1:nStress,1:nDim)
          Gmat_t(1:nDim**2+1,nDim*(i-1)+1:nDim*i)
     &                                          = Ga(1:nDim**2+1,1:nDim)
        end do

        ! transpose the vector field matrix operators
        NmatT_t     = transpose(Nmat_t)
        BmatT_t     = transpose(Bmat_t)
        GmatT_t     = transpose(Gmat_t)

        !!!!!!!!!!!!! END CALCULATING ELEMENT OPERATORS !!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!!!! CONSTITUTIVE CALCULATION !!!!!!!!!!!!!!!!!!!

         ! calculate the coordinate of integration point
        coord_ip = matmul(Nmat, reshape(coords, [nDOFEL, 1]))


        ! calculate deformation gradient and deformation tensors
        F                 = zero
        F(1:nDim,1:nDim)  = ID + matmul(uNode,dNdX)
        F(3,3)            = r_t/R

        ! calculate jacobian (volume change) at the current integration pt
        detF    = det(F)
        FInv    = inv(F)
        FInvT   = transpose(FInv)

        ! calculate chemical potential and its gradient
        mu      = dot_product( Nxi, reshape(muNode, [mDOFEl]) )

        dMudX   = matmul( BmatScalar, muNode )


        !! definition of modified deformation gradient, F-bar
        if (fbarFlag .eq. .true.) then
          if ( (jtype .eq. 4) .and. (nInt .eq. 4) ) then
            ! fully-integrated HEX8 element
            Fbar    = (detF0/detF)**(third) * F
            tanFac1 = (detF0/detF)**(-one/three)
            tanFac2 = (detF0/detF)**(-two/three)

          else
            ! standard F for all other available elements
            Fbar    = F
            tanFac1 = one
            tanFac2 = one
            call msg%ferror( flag=error, src='gel_axisymmetric',
     &          msg='F-bar is not available: ', ivec=[jtype, nInt])
          call xit
          end if
        else
          ! set F-bar = F if fbarFlag is .false. for all element
          Fbar    = F
          tanFac1 = one
          tanFac2 = one
        end if


        select case(matID)
        case (1)
          call neohookean_flory(kstep,kinc,time,dtime,nDim,analysis,
     &          nStress,nNode,jelem,intpt,coord_ip,props,nprops,
     &          jprops,njprops,matID,Fbar,mu,dMudX,
     &          svars,nsvars,fieldVar,dfieldVar,npredf,
     &          stressTensorPK2,dCwdt,Jw,
     &          CTensor,stressTensorCauchy,
     &          FSTensorUM,dCwdFTensor,dJwdFTensor,dCwdMu,dJwdMu,MmatW)

        case default
          call msg%ferror(flag=error, src='gel_axisymmetric',
     &              msg='Material model is not available', ia=matID)
          call xit
      end select

        !!!!!!!!!!!!!!! END CONSTITUTIVE CALCULATION !!!!!!!!!!!!!!!!!!



        !!!!!!!!!!!!! FORM ADDITIONAL ELEMENT OPERATORS !!!!!!!!!!!!!!!

        stressCauchy(1,1)   = stressTensorCauchy(1,1)
        stressCauchy(2,1)   = stressTensorCauchy(2,2)
        stressCauchy(3,1)   = stressTensorCauchy(1,2)
        stressCauchy(4,1)   = stressTensorCauchy(3,3)



        ! transform unsymmetric A tensor to A matrix
        call unsymmMatrix( CTensor(1:nDim,1:nDim,1:nDim,1:nDim),
     &        Amat(1:nDim*nDim,1:nDim*nDim) )

        k = 1
        do i = 1, 2
          do j = 1, 2
            Amat(5,k)   = CTensor(3,3,j,i)
            Amat(k,5)   = CTensor(j,i,3,3)
            k           = k + 1
          end do
        end do
        Amat(5,5)       = CTensor(3,3,3,3)

        ! reshape FSTensorUM into vector form
        aVectUM   = reshape( FSTensorUM(1:nDim,1:nDim), [nDim*nDim,1] )
        aVectUM(nDim*nDim+1,1)  = FSTensorUM(3,3)

        ! map the third-order tensor, dJw/dF, to a rank-2 matrix
        DmatMU  = zero
        do i = 1, nDim
          do l = 1, nDim
            do k = 1, nDim
              DmatMU(i,(l-1)*nDim+k) = dJwdFTensor(i,k,l)
            end do
          end do
          DmatMU(i,nDim*nDim+1)      = dJwdFTensor(i,3,3)
        end do

        ! time derivatives
        dCwdotdMu       = dCwdMu/dtime

        dCwdotdFTensor  = dCwdFTensor/dtime

        ! reshape dCwdotdFtensor into a vector
        dCwdotdF(:,:)   =
     &        reshape( dCwdotdFTensor(1:nDim,1:nDim),[1,nDim*nDim] )

        dCwdotdF(1,5)   = dCwdotdFTensor(3,3)

        !!!!!!!!!! END FORMING ADDITIONAL ELEMENT OPERATORS !!!!!!!!!!!



        !!!!!!!!!!!!!!!! RESIDUAL VECTOR CALCULATION !!!!!!!!!!!!!!!!!!

        Ru = Ru - w(intPt) * detJ_t * AR_t *
     &            matmul( BmatT_t, stressCauchy )

        ! chemical residual
        Rm = Rm + w(intPt) * detJ * AR *
     &            (
     &            - NmatScalarT * dCwdt + matmul( BmatScalarT, Jw )
     &            )

        !!!!!!!!!!!!!! END RESIDUAL VECTOR CALCULATION !!!!!!!!!!!!!!!!



        !!!!!!!!!!!!!!!! TANGENT MATRIX CALCULATION !!!!!!!!!!!!!!!!!!!

        ! mechanical tangent matrix
        Kuu = Kuu + w(intPt) * detJ_t * AR_t *
     &        matmul( matmul( GmatT_t, Amat ), Gmat_t )


        ! mechanical-solvent tangent matrix
        Kum = Kum + w(intpt) * detJ * AR *
     &        matmul( matmul( GmatT, aVectUM ), NmatScalar )


        ! solvent-mechanical tangent matrix
        Kmu = Kmu + w(intPt) * detJ * AR *
     &        (
     &        matmul( matmul( NmatScalarT, dCwdotdF ), Gmat )
     &        - matmul( matmul( BmatScalarT, DmatMU ), Gmat )
     &        )


        ! solvent tangent matrix
        Kmm = Kmm + w(intPt) * detJ * AR *
     &        (
     &        matmul( NmatScalarT, NmatScalar ) * dCwdotdMu
     &        - matmul( matmul(BmatScalarT, dJwdMu), NmatScalar )
     &        + matmul( matmul(BmatScalarT, MmatW), BmatScalar )
     &        )

        !!!!!!!!!!!!!! END TANGENT MATRIX CALCULATION !!!!!!!!!!!!!!!!!

        !! F-bar modification block
        if (fbarFlag .eq. .true.) then

          ! adopted from Neto et al., IJSS (1996) and Chester et al. (2015)
          if ( (jtype .eq. 4) .and. (nInt .eq. 4) ) then

            Qmat = zero

            Qmat(1,1) = third*(Amat(1,1)+Amat(1,4)+Amat(1,5))
     &        - (two/three)*stressTensorCauchy(1,1)
            Qmat(2,1) = third*(Amat(2,1)+Amat(2,4)+Amat(2,5))
     &        - (two/three)*stressTensorCauchy(1,2)
            Qmat(3,1) = third*(Amat(3,1)+Amat(3,4)+Amat(3,5))
     &        - (two/three)*stressTensorCauchy(1,2)
            Qmat(4,1) = third*(Amat(4,1)+Amat(4,4)+Amat(4,5))
     &        - (two/three)*stressTensorCauchy(2,2)
            Qmat(5,1) = third*(Amat(5,1)+Amat(5,4)+Amat(5,5))
     &        - (two/three)*stressTensorCauchy(3,3)

            do j = 4, 5
              do i = 1, 5
                Qmat(i,j)   = Qmat(i,1)
              end do
            end do

            Kuu   = Kuu + w(intPt) * detJ_t * Ar_t *
     &              matmul( GmatT_t, matmul(Qmat,Gmat0_t-Gmat_t) )

          end if

        end if

      end do                           ! end of integration point loop


      ! assemble the element tangent sub-matrices and residual sub-vectors
      call assembleElement(nDim,nNode,uDOFEL,mDOFEL,nDOFEl,       ! input to the subroutine
     &     Kuu,Kum,Kmu,Kmm,Ru,Rm,                                 ! input to the subroutine
     &     Kelem,Relem)                                           ! output from the subroutine


      ! assign the element tangent matrices to Abaqus defined variables
      amatrx(1:NDOFEL,1:NDOFEL)   = Kelem(1:NDOFEL,1:NDOFEL)
      rhs(1:NDOFEL,1)             = Relem(1:NDOFEL,1)


      end subroutine gel_axisymmetric

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

      end module hydrogel_element

! **********************************************************************
! **********************************************************************
! ****************** ABAQUS USER ELEMENT SUBROUTINE ********************
! **********************************************************************
! **********************************************************************

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     & PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     & DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     & NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,PERIOD)

! **********************************************************************
!     This subroutine is called by Abaqus with above arguments
!     for each user elements defined in an Abaqus model. Users are
!     responsible for programming the element tangent/ stiffness
!     matrix and residual vectors which will be then assembled and
!     solved by Abaqus after applying the boundary conditions.
! **********************************************************************

      use global_parameters
      use error_logging
      use hydrogel_element
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


      character(len=2)      :: analysis
      character(len=8)      :: abqProcedure
      logical               :: nlgeom
      integer               :: nInt, nPostVars
      integer               :: nDim, nStress
      integer               :: uDOF, uDOFEL, mDOF, mDOFEL

      logical, parameter    :: devMode = .false.
      integer               :: lenJobName,lenOutDir
      character(len=256)    :: outDir
      character(len=256)    :: jobName
      character(len=512)    :: errFile


      ! open a log file for the current job from Abaqus job
      if (devMode .eq. .false.) then
        call getJobName(jobName, lenJobName)
        call getOutDir(outDir, lenOutDir)
        errFile = trim(outDir)//'\aaERR_'//trim(jobName)//'.dat'
        call msg%fopen( errfile=errFile )
      end if


      ! initialize primary output variable to be zero
      amatrx        = zero
      rhs(:,nrhs)   = zero
      energy        = zero



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
      ! define different element parameters
      if ( (jtype .eq. 1) .or. (jtype .eq. 2) ) then
        nDim      = 3
        analysis  = '3D'          ! three-dimensional analysis
        nStress   = 6
      else if ( (jtype .eq. 3) .or. jtype .eq. 4) then
        nDim      = 2
        analysis  = 'AX'          ! plane axisymmetric analysis
        nStress   = 4
      else if ( (jtype .eq. 5) .or. (jtype .eq. 6) ) then
        nDim      = 2
        analysis  = 'PE'          ! 2D plane-strain analysis
        nStress   = 3
      else if ( (jtype .eq. 7) .or. (jtype .eq. 8) ) then
        nDim      = 2
        analysis  = 'PS'          ! 2D plane-stress analysis (currently unavailable)
        nStress   = 3
      else
        call msg%ferror(error,src='uel',
     &            msg='Element type is unavailable: ', ia=jtype)
        call xit
      end if


      uDOF   = nDim
      uDOFEl = uDOF*nNode
      mDOF   = 1
      mDOFEL = mDOF*nNODE

      nInt      = jprops(1)
      matID     = jprops(3)
      nPostVars = jprops(4)

      ! array containing variables for post-processing
      if (.not. allocated(globalPostVars)) then

        allocate(globalPostVars(numElem,nInt,nPostVars))

        call msg%finfo('---------------------------------------')
        call msg%finfo('--------- ABAQUS HYDROGEL UEL ---------')
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
      if( dtime .eq. zero ) return

      ! call the first-order polyelectrolyte gel element subroutine
      select case (jtype)
      case (1, 2, 5, 6)
        ! (1) TET4, (2) HEX8, (5) TRI3-PE, and (6) QUAD4-PE
        call gel_general(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     &    PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     &    DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     &    NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,
     &    PERIOD,NDIM,ANALYSIS,NSTRESS,NINT,UDOF,UDOFEL,MDOF,MDOFEL)

      case (3, 4)
        ! axisymmetric subroutine for (3) TRI3 and (4) QUAD4-AX
        call gel_axisymmetric(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     &    PROPS,NPROPS,COORDS,MCRD,NNODE,Uall,DUall,Vel,Accn,JTYPE,TIME,
     &    DTIME,KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     &    NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROPS,
     &    PERIOD,NDIM,ANALYSIS,NSTRESS,NINT,UDOF,UDOFEL,MDOF,MDOFEL)

      case default
        call msg%ferror(flag=error, src='UEL',
     &                  msg='Wrong element type: ', ia=jtype)
        call xit
      end select

      END SUBROUTINE UEL

! **********************************************************************
! **********************************************************************
! ************** ABAQUS USER OUTPUT VARIABLES SUBROUTINE ***************
! **********************************************************************
! **********************************************************************

       SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     & NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     & JMAC,JMATYP,MATLAYO,LACCFLA)

! **********************************************************************
!     This subroutine is called by Abaqus at each material point (int pt)
!     to obtain the user defined output variables for standard Abaqus
!     elements. We used an additional layer of standard Abaqus elements
!     with same topology (same number of nodes and int pts) on top of
!     the user elements with an offset in the numbering between the user
!     elements and standard elements. This number is defined in the
!     post_processing module and should match with Abaqus input file.
! **********************************************************************

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