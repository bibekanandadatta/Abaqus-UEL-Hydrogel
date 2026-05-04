! **********************************************************************
! *************** POLYELECTROLYTE HYDROGEL ELEMENT MODULE **************
! **********************************************************************
!     Author: Bibekananda Datta (C) May 2024. All Rights Reserved.
!  This module and dependencies are shared under 3-clause BSD license
! **********************************************************************
!   This module contains the material point calculation and returns all
!   the constitutive output to element formulation subroutine at each
!   integration point. It also stores the post-processing global variables.
!   Currently, one one material model is available:
!   Neo-Hookean elastomer + Flory-Huggins potential
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
!                        LIST OF ELEMENT PROPERTIES
!
!     jprops(1)   = nInt            no of integration points in element
!     jprops(2)   = fbarFlag        flag to use F-bar modification or not
!     jprops(3)   = matID           constitutive model for elastomeric network
!     jprops(4)   = nPostVars       no of local (int pt) post-processing variables
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
!                        POST-PROCESSED VARIABLES
!                     (follows the convention above)
!
!     uvar(1:nStress)                       Cauchy stress tensor components
!     uvar(nStress+1:2*nStress)             Euler-Almansi strain tensor components
!     uvar(2*nStress+1)                     Polymer volume fraction (phi)
!
!     nStress                               = 6 for 3D elements
!                                           = 4 for axisymmetric elements
!                                           = 3 for plane strain elements
! **********************************************************************

      module hydrogel_material

      use global_parameters
      use error_logging
      use linear_algebra
      use solid_mechanics
      use nonlinear_solver
      use post_processing

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