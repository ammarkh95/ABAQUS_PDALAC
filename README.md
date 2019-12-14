# Development of the Failure Criteria for Composites 

<img src="resources/Picture1.png" width="400">
<img src="resources/Picture2.png" width="220">

- [Development of the Failure Criteria for Composites](#development-of-the-failure-criteria-for-composites)
  - [Overview &amp; Motivation](#overview-amp-motivation)
  - [User Material Subroutines](#user-material-subroutines)


## Overview & Motivation
<p align="left">
<img src="resources/Motivation5.jpg" width="150">
<img src="resources/2019-07-11_15-37-18.jpg" width="135">
</p>

<p style="text-align: justify;">
Composite materials are a challenge to analyze. This challenge is characterized by complex
anisotropic behavior of the material and its heterogenous microstructure. This in turn leads to
complex failure and damage modes which are of primary concern for the safety and
performance of composite structures design. In context of computational mechanics, the
available material models in commercial nonlinear finite element analysis tools for modeling
these new material forms are often lagging behind the material science developments. In this
regard, this project aims to develop custom material subroutines for composite materials, the
developed models predict the onset of damage and damage progression in composite structures
according to wide range of failure theories and damage models that can be customized by users
accordingly to fit their requirements. The developed material models are verified using number
of benchmark problems from literature. Finally, realizing the challenges associated with
writing user defined material models, a standalone code “PDALAC” is developed to help
researchers to visualize the progressive damage analysis process without the need to conduct a
FEA, the developed code is written in vectorized form which can be easily translated to work
with multiple standard FEM packages (e.g. ABAQUS, ANSYS, LS-DYNA)</p>

<img src="resources/introduction.png" width="500">

More text here

- 1
- 2
- 3
  
Then we have numbered list

1. Num1
2. Num2
3. Num3

Here is a link to [google](https://www.google.com)

This is **bold**

This is _italic_

This is underline

## User Material Subroutines

![Anything here](resources/Progressive_Damage.png)


**UMAT Header**

```fortran
SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)


      user coding to define DDSDDE, STRESS, STATEV, SSE, SPD, SCD
      and, if necessary, RPL, DDSDDT, DRPLDE, DRPLDT, PNEWDT


      RETURN
      END
```


**VUMAT Header**

```fortran
      subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock),
C
      character*80 cmname
C

      do 100 km = 1,nblock
        user coding
  100 continue

      return
      end
```

**User-defined property data for the UMAT subroutine**

|**PROPS Array Entry**|**Variable Name**|**Description**|
|---|---|---|
| 1,3|E11,E22,E33|Laminate Orthotropic Young's Moduli values (E11 fiber direction)|
| 4,6|ANU12,ANU13,ANU23|Poisson ratios v12,v13,v23|
| 7,9|G12,G13,G23|Laminate Shear Moduli values|
| 10,11|Xt,Xc|Laminate fiber allowable strengths (Tension/Compression)|
| 12,13|Yt,Yc|Laminate matrix allowable strengths (Tension/Compression)|
| 14,15|Zt,Zc|(Laminate interlaminar allowable strengths (Tension/Compression)|
| 16,17,18|S12,S13,S23|(Laminate allowable shear strengths (1-2,1-3,2-3 planes)|
| 19,20|EpsXt, EpsXc|Laminate fiber allowable strains (Tension/Compression)|
| 21,22|EpsYt, EpsYc|Laminate matrix allowable strains (Tension/Compression)|
| 23,34|EpsZt, EpsZc|(Laminate Interlaminar allowable strains (Tension/Compression)|
|25,26,27|GamS12,Gams13,GamS23|Laminate allowable engineering shear strains (1-2,1-3,2-3 planes)|
|28|failure_id|Failure criteria selection ID: (1) Max Stress (2)Max Strain (3)Tsai-Wu (4)Hoffman (5)Hashin (6)Hashin-Rotem (7)Puck|
|29|damage_id|Damage Model selection ID: (1) Instant (2)Recursive (3)Exponential (4)Constant Stress (5)(CDM)Crack-Band Theory|
|30,31|beta_ft, beta_fc|Instant/Recursive: fiber degradation factors (Tension/Compression)|
|32,33|beta_mt, beta_mc|Instant/Recursive: matrix degradation factors (Tension/Compression)|
|34|beta_s|Instant/Recursive: shear degradation factor (Tension/Compression)|
|35,36|a_ft, n_ft|Exponential: fiber degradation factors (Tension)|
|37,38|a_fc, n_fc|Exponential: fiber degradation factors (Compression)|
|39,40|a_mt, n_mt|Exponential: matrix degradation factors (Tension)|
|41,42|a_mc, n_mc|Exponential: matrix degradation factors (Compression)|
|43,44|a_s, n_s|Exponential: shear degradation factors|
|45,46|G_ft, G_fc|CDM: fiber fracture energies (Tension/Compression)|
|47,48|G_mt, G_mc|CDM: matrix fracture energies (Tension/Compression)|
|49,50|G_IC,G_IIC|CDM: shear fracture energies (Tension/Compression)|
|51,52|le,let|CDM: characteristic element length and thickness|
|53,54|alpha, beta|CDM: nonlinear shear degradation factors|
|55|THETAF|Puck: Maximum fracture angle in radians|
|56|MGF|Puck: Magnification factor|
|57,58|E11F, ANU12F|Puck: Fiber elastic modulus and Poisson ratio|



**UMAT-defined solution-dependent variables**

|**STATEV Array Entry**|**Variable Name**|**Description**|
|---|---|---|
|1| dmg(1)|Degradation factor for the σ11 stress component|
|2| dmg(2)|Degradation factor for the σ22 stress component|
|3| dmg(3)|Degradation factor for the σ33 stress component|
|4| dmg(4)| Degradation factor for the σ12 stress component|
|5| dmg(5)| Degradation factor for the σ13 stress component|
|6| dmg(6)| Degradation factor for the σ23 stress component|
|7|fflags(1)| Failure flag for first failure mode|
|8| fflags(2)| Failure flag for second failure mode|
|9| fflags(3)| Failure flag for third failure mode|
|10| fflags(4)| Failure flag for fourth failure mode|
|11| fflags(5)| Failure flag for fifth failure mode|
|12| fflags(6)| Failure flag for sixth failure mode|
|13| DelEl| Element deletion variable|


This is inline code `To Be Continued...`