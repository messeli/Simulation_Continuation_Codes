!“I don't know what the language of the year 2000 will look like, but I know it will be called Fortran.” —Tony Hoare, winner of the 1980 Turing Award, in 1982.
! From an article by Lee Philips at arstechnica.com .
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!   frc modfd2 duffingBN_forced modfd2 ZilliCubic:      
!   ZILLI CUBIC ROT (AUTONOMOUS - ISOTROPIC STIFFNESS)
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
        DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
        DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
        DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)
        DOUBLE PRECISION ZETA,Q1,Q2,Q3,Q4,U5,S4,C4,R2,GAMMA,GAMMA2,THETA,THETAP,THETAPP,MH,EPSH,JPH,X,Y
      
        GAMMA = PAR(1) 
        THETAP = PAR(2)
        MH = PAR(3) 
        GAMMA2 = PAR(4) 
        EPSH = PAR(5)
        ZETA = PAR(6) 
        JPH  = PAR(7)
        THETAPP = PAR(8)

        Q1=U(1)
        Q2=U(2)
        Q3=U(3)
        Q4=U(4)
 
        R2 = (Q1**2)+(Q3**2)
  
        F(1) = Q2
        F(2) = -THETAP*(JPH-2)*Q4          &
               -2*ZETA*Q2                  &
               +((THETAP**2)*(1-JPH)-1)*Q1 &
               +2*ZETA*THETAP*Q3           &
               +MH*EPSH*(THETAP**2)        &
               -GAMMA*R2*Q1 !!!-GAMMA*Q1**3 - GAMMA2*Q3**2*Q1
               ! -GAMMA/4*( 3*Q1*R2 + C4*Q1*(R2-4*(Q3**2)) + S4*Q3*(R2-4*(Q1**2)) ) 
        F(3) = Q4
        F(4) = +THETAP*(JPH-2)*Q2          &
               -2*ZETA*Q4                  &
               +((THETAP**2)*(1-JPH)-1)*Q3 &
               -2*ZETA*THETAP*Q1           &
               -MH*EPSH*THETAPP            &
               -GAMMA*R2*Q3  !!!-GAMMA*Q3**3 - GAMMA2*Q1**2*Q3
               ! -GAMMA/4*( 3*Q3*R2 + C4*Q3*(R2-4*(Q1**2)) - S4*Q1*(R2-4*(Q3**2)) ) 
        

      ! IF (IJAC.EQ.1) RETURN
      !   DFDU(1,1)=0
      !   DFDU(1,2)=1
      !   DFDU(1,3)=0
      !   DFDU(1,4)=0
      !   DFDU(1,5)=0
      !   DFDU(1,6)=0

      !   DFDU(2,1)=+((THETAP**2)*(1-JPH)-1) &
      !             -GAMMA/4*( 3*R2 + 6*Q1**2 + C4*(R2-4*Q3**2) + 2*C4*Q1**2 - 6*S4*Q3*Q1 )
      !   DFDU(2,2)=-2*ZETA
      !   DFDU(2,3)=+2*ZETA*THETAP + &
      !             -GAMMA/4*( +6*Q1*Q3 - 6*C4*Q1*Q3 + S4*(R2-4*Q1**2) + 2*S4*Q3**2  )
      !   DFDU(2,4)=-THETAP*(JPH-2)
      !   DFDU(2,5)=-GAMMA/4*Q3*(R2-4*Q1**2) 
      !   DFDU(2,6)=-GAMMA/4*Q1*(R2-4*Q3**2)  

      !   DFDU(3,1)=0
      !   DFDU(3,2)=0
      !   DFDU(3,3)=0
      !   DFDU(3,4)=1
      !   DFDU(3,5)=0
      !   DFDU(3,6)=0

      !   DFDU(4,1)=-2*ZETA*THETAP &
      !             -GAMMA/4*( 6*Q3*Q1 - 6*C4*Q3*Q1 - S4*(R2-4*Q3**2) - 2*S4*Q1**2 )
      !   DFDU(4,2)=+THETAP*(JPH-2)
      !   DFDU(4,3)=+((THETAP**2)*(1-JPH)-1) &
      !             -GAMMA/4*( 3*R2 + 6*Q3**2 + C4*(R2-4*Q1**2) + 2*C4*Q3**2 +6*S4*Q1*Q3 )
      !   DFDU(4,4)=-2*ZETA
      !   DFDU(4,5)=+GAMMA/4*Q1*(R2-4*Q3**2)  
      !   DFDU(4,6)=-GAMMA/4*Q3*(R2-4*Q1**2) 

      !   DFDU(5,1)=0
      !   DFDU(5,2)=0
      !   DFDU(5,3)=0
      !   DFDU(5,4)=0
      !   DFDU(5,5)=1- ( 1*(X**2 + Y**2) + X*(2*X) )    
      !   DFDU(5,6)=4*THETAP-X*(2*Y)    

      !   DFDU(6,1)=0
      !   DFDU(6,2)=0
      !   DFDU(6,3)=0
      !   DFDU(6,4)=0
      !   DFDU(6,5)=-4*THETAP-Y*(2*X)
      !   DFDU(6,6)=1- ( 1*(X**2 + Y**2) + Y*(2*Y))
      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)  
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM
        DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
        DOUBLE PRECISION, INTENT(IN) :: T
        DOUBLE PRECISION ZETA,Q1,Q2,Q3,Q4,S4,C4,R2,GAMMA,GAMMA2,THETA,THETAP,THETAPP,MH,EPSH,JPH,X,Y,TPI!WN,F_M,K3_M,OMEG,

        ! GIVEN BACKWARDS (MATLAB LOWER AMPLITUDE END DIMENSIONS)
        GAMMA = 0.25 !2nd continue gamma from 0 to 0.25
        THETAP = 7.0 ! Continue THETAP from 0 to 7. 
        MH = 0.9 !1st continue MH from 0 to 0.9
  
        EPSH =0.353
        ZETA = 1e-2
        JPH  = 0.143
        THETAPP = 0.0

        PAR(1)=GAMMA 
        PAR(2)=THETAP 
        PAR(3)=MH 
        PAR(4)=GAMMA2
        PAR(5)=EPSH 
        PAR(6)=ZETA 
        PAR(7)=JPH  
        PAR(8)=THETAPP 

        !| Orig one used for orig givenBackwards, where it was zeta=0.01 
        ! U(1)=-0.380009000620448
        ! U(2)= 0.000022150720618
        ! U(3)=-0.001297860009550
        ! U(4)= 0.000206860716715

        !| New ones from "Google Drive\PhD\Zilli ISO cubic stiffness\Resutls ISO cubic\
        !| ...7p0 end lowAmp datum, differentZeta for AUTO givenBackwards.txt"
        U(1)=-0.380085547357884
        U(2)=0.000000023340008
        U(3)=-0.001299216246013
        U(4)=0.000000244301170

      END SUBROUTINE STPNT


      ! SUBROUTINE PVLS(NDIM,U,PAR)
      ! IMPLICIT NONE
      ! INTEGER, INTENT(IN) :: NDIM
      ! DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      ! DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
      ! DOUBLE PRECISION R2,Q1,Q3
      !   Q1 = U(1)
      !   Q3 = U(3)

      !   R2 = Q1**2 + Q3**2 
      !   ! Which one is this?
      !   ! ?? max( U(1) )^2 + max( U(3) )^2
      !   ! ?? max( U(1)^2 + U(3)^2 )

      !   PAR(9) = R2**(0.5)
      ! END SUBROUTINE PVLS

      SUBROUTINE PVLS(NDIM,U,PAR)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NDIM
        DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
        DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
        DOUBLE PRECISION, ALLOCATABLE :: UU(:,:),R2(:),R(:)
        DOUBLE PRECISION, EXTERNAL :: GETP, GETU
        INTEGER :: NDX,NCOL,NTST,i,j,k
        !|SEE the explanation of demo pvl.f90 : 
        !|:"For algebraic problems the argument U is, as usual, the state vector.
        !|:For differential equations the argument U represents the approximate 
        !|:solution on the entire interval [0,1]. In this case its values can
        !|:be accessed indirectly by calls to GETP, as illustrated below, or
        !|:by obtaining NDIM, NCOL, NTST via GETP and then dimensioning U as
        !|:U(NDIM,0:NCOL*NTST) in a seperate subroutine that is called by PVLS. "

        ! NDX  = NINT(GETP('NDX',0,U))
        ! NTST = NINT(GETP('NTST',0,U))
        ! NCOL = NINT(GETP('NCOL',0,U))

        !|AAA/BBB METHOD A/B
        ! PAR(9) = GETU(i,j,U,NDX,NTST,NCOL)
        !|AAA/BBB END
        
        !|CCC METHOD C
        ! DO i=1,NDX
        !   DO j=0,NCOL*NTST
        !     UU(i,j) = GETU(i,j,U,NDX,NTST,NCOL)
        !   END DO
        ! END DO
        ! R2(:) = UU(1,:)**2 + UU(3,:)**2  !El by El multip
        ! R(:) = R2**0.5 
        ! PAR(9) = MAXVAL(R)
        !|CCC END

        !|DDD DEFAULT METHOD - PHASED AMP
        !|:CORRECT AMPLITUDE FOR IPS=1
        PAR(9) = (U(1)**2 + U(3)**2)**0.5 
        !|:The state vector is directly used to calc amplitude.
        !|:See the qoute above.
        !|DDD END 
      END SUBROUTINE PVLS


      DOUBLE PRECISION FUNCTION GETU(i,j,U,NDX,NTST,NCOL)
        INTEGER, INTENT(IN) :: NDX,NCOL,NTST,i,j
        INTEGER :: p
        DOUBLE PRECISION, INTENT(IN) :: U(NDX,0:NCOL*NTST)
        DOUBLE PRECISION, ALLOCATABLE :: R2(:),R(:)
        DOUBLE PRECISION :: M,N

        !|AAA METHOD A
        M = 0.0
        DO p=0,NCOL*NTST
          N = ( U(1,p)**2+U(3,p)**2 )**0.5
          if (N.GT.M) THEN
            M = N
          END IF 
        END DO 

        GETU = M
        !|AAA END

        !|BBB METHOD B
        ! R2(:) = U(1,:)**2 + U(3,:)**2  !
        ! R(:) = R2**0.5 
        ! GETU = MAXVAL(R)
        !|BBB END

        !|CCC METHOD C
        ! GETU = U(i,j)
        !|CCC END
      END FUNCTION GETU



      SUBROUTINE BCND
      END SUBROUTINE BCND

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT

