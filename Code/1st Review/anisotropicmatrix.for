      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
C
      DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS, NTENS),
     1 DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),
     2 PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3), DROT(3, 3),
     3 DFGRD0(3, 3), DFGRD1(3, 3)
      
      DIMENSION EELAS(6), EPLAS(6), FLOW(6)
C
      PARAMETER(ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0, THREE=3.0D0, SIX=6.0D0,
     1 ENUMAX=.4999D0, NEWTON=10, TOLER=1.0D-6)

C   http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_orthotropic.cfm
Cg11=2*g23
Cg22=2*g31
Cg33=2*g12
      pr12=PROPS(1)
      pr13=PROPS(2)
      pr23=PROPS(3)
      e11=PROPS(4)
      e22=PROPS(5)
      e33=PROPS(6)   
      g11=PROPS(7)
      g22=PROPS(8)
      g33=PROPS(9)
      
      pr21=e22*pr12/e11
      pr31=e33*pr13/e11
      pr32=e33*pr23/e22
      
      tri=(ONE-pr12*pr21-pr23*pr32-pr31*pr13-TWO*pr12*pr23*pr31)/(e11*e22*e33)   

      DDSDDE(1,1)=(ONE-pr23*pr32)/(e22*e33*tri)
      DDSDDE(1,2)=(pr21+pr31*pr23)/(e22*e33*tri)
      DDSDDE(1,3)=(pr31+pr21*pr32)/(e22*e33*tri)
      DDSDDE(1,4)=0
      DDSDDE(1,5)=0
      DDSDDE(1,6)=0
      
      DDSDDE(2,1)=(pr12+pr13*pr32)/(e33*e11*tri)
      DDSDDE(2,2)=(ONE-pr31*pr13)/(e33*e11*tri)
      DDSDDE(2,3)=(pr32+pr31*pr12)/(e33*e11*tri)
      DDSDDE(2,4)=0
      DDSDDE(2,5)=0
      DDSDDE(2,6)=0
      
      DDSDDE(3,1)=(pr13+pr12*pr23)/(e11*e22*tri)
      DDSDDE(3,2)=(pr23+pr13*pr21)/(e11*e22*tri)
      DDSDDE(3,3)=(ONE-pr12*pr21)/(e11*e22*tri)
      DDSDDE(3,4)=0
      DDSDDE(3,5)=0
      DDSDDE(3,6)=0
      
      DDSDDE(4,1)=0
      DDSDDE(4,2)=0
      DDSDDE(4,3)=0
      DDSDDE(4,4)=ONE*g11
      DDSDDE(4,5)=0
      DDSDDE(4,6)=0
      
      DDSDDE(5,1)=0
      DDSDDE(5,2)=0
      DDSDDE(5,3)=0
      DDSDDE(5,4)=0
      DDSDDE(5,5)=ONE*g22
      DDSDDE(5,6)=0
      
      DDSDDE(6,1)=0
      DDSDDE(6,2)=0
      DDSDDE(6,3)=0
      DDSDDE(6,4)=0
      DDSDDE(6,5)=0
      DDSDDE(6,6)=ONE*g33
      
      DO K1=1, NTENS
          DO K2=1, NTENS
              STRESS(K2)=STRESS(K2)+DDSDDE(K2, K1)*DSTRAN(K1)
          END DO
      END DO
      
C
      print *,"***********STRESS***********" 
      print *, STRESS
      RETURN
      END