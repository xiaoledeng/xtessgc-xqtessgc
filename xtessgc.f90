!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program xtessgc
!
!    Test driver of "tessgc", the double precision Fortran subroutine to compute
!    the gravitational potential/acceleration vector/gradient tensor/curvatures of a tesseroid.
!
!    Reference: 
!           Fukushima T (2018) Accurate computation of gravitational field of a tesseroid.
!                Journal of Geodesy 92(12):1371–1386, doi:10.1007/s00190-018-1126-2
!           Deng XL (2023) Corrections to: “Accurate computation of gravitational field 
!                of a tesseroid” by Fukushima (2018) in J. Geod. 92(12):1371–1386, Journal of Geodesy,
!                doi: 10.1007/s00190-022-01673-2
!
implicit integer (i-n)
implicit real*8 (a-h,o-z)
!
PI=atan(1.d0)*4.d0
raddeg=PI/180.d0
delta=1.d-16
write (*,"(a15,1pe15.7)") "# delta=", delta
!
!    A sample tesseroid covering an eastern Himalaya: size 1 deg x 1 deg x 50 km
!      (Note) Mt. Everest: h = 8.848 km, phi = 27 deg 59 min, lat = 86 deg 55 min
!
R0=6380.d0
HB=-40.d0
HT=10.d0
PhiS=27.d0*raddeg
PhiN=28.d0*raddeg
FlamW=86.d0*raddeg
FlamE=87.d0*raddeg
!
!   Test evaluation points are along a radial straight line
!    passing through the geometrical center of the tesseroid
!
phi=27.5d0*raddeg
flam=86.5d0*raddeg
!
write (*,"(a15,1p3e15.7)") "# R0,HT,HB=",R0,HT,HB
write (*,"(a15,1p2e15.7)") "# PhiN,PhiS=",PhiN/raddeg,PhiS/raddeg
write (*,"(a15,1p2e15.7)") "# LamdaE,LamdaW=",FlamE/raddeg,FlamW/raddeg
write (*,"(a15,1p2e15.7)") "# phi,lambda=",phi/raddeg,flam/raddeg
!
write (*,"(a15,a25)") "# H","V"
write (*,"(a15,3a25)") "#","gPhi","gLambda","gH"
write (*,"(a15,3a25)") "#","GammaPhiPhi","GammaPhiLambda","GammaPhiH"
write (*,"(a15,3a25)") "#","GammaLambdaLambda","GammaLambdaH","GammaHH"
write (*,"(a15,3a25)") "#","gggPhiPhiPhi","gggPhiPhiLam","gggPhiPhiH"
write (*,"(a15,1a25)") "#","gggPhiLamH"
write (*,"(a15,3a25)") "#","gggLamLamPhi","gggLamLamLam","gggLamLamH"
write (*,"(a15,3a25)") "#","gggHHPhi","gggHHLam","gggHHH"
!
!    Compute the field along a radial line passing through the tesseroid
!
call cpu_time(time_begin)
dH=10.d0
do m=-10,29
    h=(dble(m)+0.5d0)*dH
!
    call tess(PhiN,PhiS,FlamE,FlamW,HT,HB,R0,delta, &
        phi,flam,h, &
        v,gPhi,gLam,gH,ggPhiPhi,ggPhiLam,ggPhiH,ggLamLam,ggLamH,ggHH, &
        gggPhiPhiPhi,gggPhiPhiLam,gggPhiPhiH,gggPhiLamH,gggLamLamPhi, &
        gggLamLamLam,gggLamLamH,gggHHPhi,gggHHLam,gggHHH)
    write (*,"(1pe15.7,1pe25.15)") h,v
    write (*,"(15x,1p3e25.15)") gPhi,gLam,gH
    write (*,"(15x,1p3e25.15)") ggPhiPhi,ggPhiLam,ggPhiH
    write (*,"(15x,1p3e25.15)") ggLamLam,ggLamH,ggHH
    write (*,"(15x,1p3e25.15)") gggPhiPhiPhi,gggPhiPhiLam,gggPhiPhiH
    write (*,"(15x,1p1e25.15)") gggPhiLamH
    write (*,"(15x,1p3e25.15)") gggLamLamPhi,gggLamLamLam,gggLamLamH
    write (*,"(15x,1p3e25.15)") gggHHPhi,gggHHLam,gggHHH
enddo
call cpu_time(time_end)
write(*,*) "The total time is",time_end-time_begin,"seconds."
!
stop
end program xtessgc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tess(PhiN0,PhiS0,FlamE0,FlamW0,HT0,HB0,R00,delta0,phi,flam,h, &
    v,gPhi,gLam,gH,ggPhiPhi,ggPhiLam,ggPhiH,ggLamLam,ggLamH,ggHH, &
        gggPhiPhiPhi,gggPhiPhiLam,gggPhiPhiH,gggPhiLamH,gggLamLamPhi, &
        gggLamLamLam,gggLamLamH,gggHHPhi,gggHHLam,gggHHH)
!
!    The double precision Fortran subroutine to compute
!    the gravitational potential, V, the gravitational acceleration vector, g,
!    and the gravity gradient tensor, Gamma,
!    and the gravitational curvatures, ggg, of a tesseroid.
!
!    Input parameters:
!      PhiN0: Latitude of the north end point of the tesseroid
!      PhiS0: Latitude of the south end point of the tesseroid, PhiS0 < PhiN0
!      FlamE0: Longitude of the east end point of the tesseroid
!      FlamW0: Longitude of the west end point of the tesseroid, FlamW0 < FlamE0
!      HT0: Height of the top end point of the tesseroid
!      HB0: Height of the bottom end point of the tesseroid, HB0 < HT0
!      R00: Radius of the reference spherical surface from which H is counted
!      delta0: the relative error tolerance (typical value is 1d-9 or 1d-12)
!    Input variables:
!      phi: Latitude of the evaluation point
!      flam: Longitude of the evaluation point
!      h: Height from the reference sphere of the evaluation point
!        (Note) We adopted "km" as the unit of h as well as HT0 and HB0.
!               But any unit is OK if R00 is expressed in the same unit.
!    Output variables:
!      v: Normalized gravitational potential at the evaluation point
!        (Note) normalization constant: V0 = G rho R00^2
!      gPhi: Latitude component of normalized gravitational acceleration vector
!      gLam: Longitude component of normalized gravitational acceleration vector
!      gH: Height component of normalized gravitational acceleration vector
!        (Note) normalization constant: g0 = G rho R00^2
!      ggPhiPhi: Latitude-latitude component of normalized Gamma
!      ggPhiLam: Latitude-longitude component of normalized Gamma
!      ggPhiH: Latitude-height component of normalized Gamma
!      ggLamLam: Longitude-longitude component of normalized Gamma
!      ggLamH: Longitude-height component of normalized Gamma
!      ggHH: Height-height component of normalized Gamma
!        (Note) normalization constant: Gamma0 = G rho R00^2
!      gggPhiPhiPhi:    Latitude-latitude-latitude component of normalized GGG
!      gggPhiPhiLam:    Latitude-latitude-longitude component of normalized GGG
!      gggPhiPhiH  :    Latitude-latitude-height component of normalized GGG
!      gggPhiLamH  :    Latitude-longitude-height component of normalized GGG
!      gggLamLamPhi:    Longitude-longitude-latitude component of normalized GGG
!      gggLamLamLam:    Longitude-longitude-longitude component of normalized GGG
!      gggLamLamH  :    Longitude-longitude-height component of normalized GGG
!      gggHHPhi    :    Height-height-latitude component of normalized GGG
!      gggHHLam    :    Height-height-longitude component of normalized GGG
!      gggHHH      :    Height-height-height component of normalized GGG
!                       (Note) normalization constant: GGG0 = G rho R00^2
!
implicit integer (i-n)
implicit real*8 (a-h,o-z)
!
!    One should declare "vtess" as an external function name
!
external vtess
!
!    Common blocks for passing the parameters and indirect variables
!
real*8 eps
common /epsV/eps
real*8 R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
common /paramV/R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
real*8 delta,delta1,delta2,d1phi,d1lam,d1h,d2phi,d2lam,d2h,delta3,d3phi,d3lam,d3h
common /deltaV/delta,delta1,delta2,d1phi,d1lam,d1h,d2phi,d2lam,d2h,delta3,d3phi,d3lam,d3h
!
!    Store the input variables into common blocks
!    so that the subprograms can refer them indirectly
!
PhiN=PhiN0
PhiS=PhiS0
FlamE=FlamE0
FlamW=FlamW0
HT=HT0
HB=HB0
R0=R00
!
!    Calculate parameters of numerical integration and differentiation
!
eps=delta0
delta=delta0
delta1=delta**(1.d0/3.d0)
delta2=sqrt(sqrt(delta))
delta3=delta**(1.d0/5.d0)
!
!    Compute tesseroid-dependent constants beforehand
!
deltaH=HT-HB
deltaPhi=PhiN-PhiS
deltaFlam=FlamE-FlamW
beta=deltaH/R0
RC=R0+(HT+HB)*0.5d0
PhiC=(PhiN+PhiS)*0.5d0
FlamC=(FlamE+FlamW)*0.5d0
sinPhiC=sin(PhiC)
cosPhiC=cos(PhiC)
!
!    Compute evaluation-point-dependent constants beforehand
!
r=R0+h
sinphi=sin(phi)
cosphi=cos(phi)
cosdlam=cos(flam-FlamC)
scale0=0.5d0*sqrt(deltaH*deltaH+RC*RC*(deltaPhi*deltaPhi &
    +cosPhiC*cosPhiC*deltaFlam*deltaFlam))
scale1=sqrt(r*r+RC*RC-2.d0*RC*r*(sinphi*sinPhiC+cosphi*cosPhiC*cosdlam))
scale=max(scale0,scale1)
d1h=scale*delta1
! add d1phi, d1lam
d1phi=d1h/r
d1lam=d1h/(r*cosphi)
d2h=scale*delta2
! add d2phi, d2lam
d2phi=d2h/r
d2lam=d2h/(r*cosphi)
! add the d3h, d3phi, d3lam
d3h=scale*delta3
d3phi=d3h/r
d3lam=d3h/(r*cosphi)
!
!    Compute the gravitational potential by numerical integration
!
v=vtess(phi,flam,h)
!
!    Compute the gravitational acceleration vector by numerical differentiation
!
call gtess(vtess,phi,flam,h,v,gPhi,gLam,gH)
!
!    Compute the gravity gradient tensor by numerical differentiation
!
call ggtess(vtess,phi,flam,h,v,gPhi,gLam,gH, &
    ggPhiPhi,ggPhiLam,ggPhiH,ggLamLam,ggLamH,ggHH)
!
!    Compute the gravitational curvatures by numerical differentiation
!
call gggtess(vtess,phi,flam,h,v,gPhi,gLam,gH, &
    ggPhiPhi,ggPhiLam,ggPhiH,ggLamLam,ggLamH,ggHH, &
    gggPhiPhiPhi,gggPhiPhiLam,gggPhiPhiH,gggPhiLamH,gggLamLamPhi, &
    gggLamLamLam,gggLamLamH,gggHHPhi,gggHHLam,gggHHH)
!
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8 function vtess(phi,flam,h)
!
!    The double precision Fortran function to compute
!    the gravitational potential by the conditional split quadrature method
!
implicit integer (i-n)
implicit real*8 (a-h,o-z)
!
!    One should declare "stess" as an external function name
!
external stess
!
!    Common blocks for passing the parameters and indirect variables
!
real*8 eps
common /epsV/eps
real*8 R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
common /paramV/R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
real*8 phi0,flam0,h0,xiN,xiS,etaE,etaW
common /argV/phi0,flam0,h0,xiN,xiS,etaE,etaW
real*8 alpha,gamma,zetaT,zetaB,sigma,fmu
common /paramV2/alpha,gamma,zetaT,zetaB,sigma,fmu
!
!    Store the input variables into common blocks
!    so that the subprograms can refer them indirectly
!
phi0=phi
flam0=flam
h0=h
!
!    Compute tesseroid-dependent constants beforehand
!
xiN=PhiN-phi
xiS=PhiS-phi
etaE=FlamE-flam
etaW=FlamW-flam
!
!    Compute evaluation-point-dependent constants beforehand
!
alpha=1.d0+h/R0
gamma=cos(phi)
zetaB=(HB-h)/R0
zetaT=zetaB+beta
!
!    Conditional split quadrature of latitude line integration 
!
if(xiS.lt.0.d0.and.xiN.gt.0.d0) then
    call dqde(stess,xiS,0.d0,eps,v1)
    call dqde(stess,0.d0,xiN,eps,v2)
    v=v1+v2
else
    call dqde(stess,xiS,xiN,eps,v)
endif
vtess=v
!
return; end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8 function stess(xi)
!
!    The double precision Fortran function to compute
!    the latitude line integral requried in "vtess"
!
implicit integer (i-n)
implicit real*8 (a-h,o-z)
!
!    One should declare "vkern" as an external function name
!
external vkern
!
!    Common blocks for passing the parameters and indirect variables
!
real*8 eps
common /epsV/eps
real*8 phi0,flam0,h0,xiN,xiS,etaE,etaW
common /argV/phi0,flam0,h0,xiN,xiS,etaE,etaW
real*8 alpha,gamma,zetaT,zetaB,sigma,fmu
common /paramV2/alpha,gamma,zetaT,zetaB,sigma,fmu
!
!    Compute latitude-dependent constants beforehand
!
sigma=cos(phi0+xi)
fmu=sin(xi*0.5d0)
!
!    Conditional split quadrature of longitude line integration 
!
if(etaW.lt.0.d0.and.etaE.gt.0.d0) then
    call dqde1(vkern,etaW,0.d0,eps,s1)
    call dqde1(vkern,0.d0,etaE,eps,s2)
    s=s1+s2
else
    call dqde1(vkern,etaW,etaE,eps,s)
endif
stess=s
!
return; end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8 function vkern(eta)
!
!    The double precision Fortran function to compute
!    the kernel function expressed in a cancellation-free form
!
implicit integer (i-n)
implicit real*8 (a-h,o-z)
!
!    Common blocks for passing the parameters and indirect variables
!
real*8 R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
common /paramV/R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
real*8 alpha,gamma,zetaT,zetaB,sigma,fmu
common /paramV2/alpha,gamma,zetaT,zetaB,sigma,fmu
!
!    Cancellation-free computation of the kernel function
!
fnu=sin(eta*0.5d0)
B=2.d0*alpha*(fmu*fmu+gamma*sigma*fnu*fnu)
A=B*(2.d0*alpha-B)
C=(alpha-B)*(alpha-B)-A*0.5d0
D=0.5d0*(zetaT+4.d0*alpha-3.d0*B)
ST=sqrt(A+(B+zetaT)*(B+zetaT))
SB=sqrt(A+(B+zetaB)*(B+zetaB))
T=(zetaT+zetaB+2.d0*B)/(ST+SB)
if(zetaB+B.lt.0.d0) then
    fl=dlog1p((1.d0+T)*beta*(-zetaB-B+SB)/A)
else
    fl=dlog1p((1.d0+T)*beta/(zetaB+B+SB))
endif
vkern=sigma*(C*fl+beta*(D*T+SB*0.5d0))
!
return; end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gtess(vtess,phi,flam,h,v,gPhi,gLam,gH)
!
!    Compute the gravitational acceleration vector by the conditional switch
!    of the central and single-sided second order finite difference formulas
!
implicit integer (i-n)
implicit real*8 (a-h,o-z)
!
!    Common blocks for passing the parameters and indirect variables
!
real*8 R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
common /paramV/R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
real*8 delta,delta1,delta2,d1phi,d1lam,d1h,d2phi,d2lam,d2h,delta3,d3phi,d3lam,d3h
common /deltaV/delta,delta1,delta2,d1phi,d1lam,d1h,d2phi,d2lam,d2h,delta3,d3phi,d3lam,d3h
!
!    Compute conmponent-independent variables beorehand
!
r=R0+h
c=cos(phi)
!
!   Compute gPhi
!
if((flam.ge.FlamW.and.flam.le.FlamE).and.(h.ge.HB.and.h.le.HT)) then
    if((phi.gt.PhiS.and.phi-d1phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d1phi.lt.PhiN)) then
        gPhi=(-vtess(phi+2.d0*d1phi,flam,h)+4.d0*vtess(phi+d1phi,flam,h)-3.d0*v)/(2.d0*d1phi)
    elseif((phi.lt.PhiS.and.phi+d1phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d1phi.gt.PhiN)) then
        gPhi=(vtess(phi-2.d0*d1phi,flam,h)-4.d0*vtess(phi-d1phi,flam,h)+3.d0*v)/(2.d0*d1phi)
    else
        gPhi=(vtess(phi+d1phi,flam,h)-vtess(phi-d1phi,flam,h))/(2.d0*d1phi)
    endif
else
    gPhi=(vtess(phi+d1phi,flam,h)-vtess(phi-d1phi,flam,h))/(2.d0*d1phi)
endif
gPhi=gPhi/r
!
!   Compute gLam
!
if((phi.ge.PhiS.and.phi.le.PhiN).and.(h.ge.HB.and.h.le.HT)) then
    if((flam.gt.FlamW.and.flam-d1lam.lt.FlamW).or.(flam.gt.FlamE.and.flam-d1lam.lt.FlamE)) then
        gLam=(-vtess(phi,flam+2.d0*d1lam,h)+4.d0*vtess(phi,flam+d1lam,h)-3.d0*v)/(2.d0*d1lam)
    elseif((flam.lt.FlamW.and.flam+d1lam.gt.FlamW).or.(flam.lt.FlamE.and.flam+d1lam.gt.FlamE)) then
        gLam=(vtess(phi,flam-2.d0*d1lam,h)-4.d0*vtess(phi,flam-d1lam,h)+3.d0*v)/(2.d0*d1lam)
    else
        gLam=(vtess(phi,flam+d1lam,h)-vtess(phi,flam-d1lam,h))/(2.d0*d1lam)
    endif
else
    gLam=(vtess(phi,flam+d1lam,h)-vtess(phi,flam-d1lam,h))/(2.d0*d1lam)
endif
gLam=gLam/(r*c)
!
!   Compute gH
!
if((phi.ge.PhiS.and.phi.le.PhiN).and.(flam.ge.FlamW.and.flam.le.FlamE)) then
    if((h.gt.HB.and.h-d1h.lt.HB).or.(h.gt.HT.and.h-d1h.lt.HT)) then
        gH=(-vtess(phi,flam,h+2.d0*d1h)+4.d0*vtess(phi,flam,h+d1h)-3.d0*v)/(2.d0*d1h)
    elseif((h.lt.HB.and.h+d1h.gt.HB).or.(h.lt.HT.and.h+d1h.gt.HT)) then
        gH=(vtess(phi,flam,h-2.d0*d1h)-4.d0*vtess(phi,flam,h-d1h)+3.d0*v)/(2.d0*d1h)
    else
        gH=(vtess(phi,flam,h+d1h)-vtess(phi,flam,h-d1h))/(2.d0*d1h)
    endif
else
    gH=(vtess(phi,flam,h+d1h)-vtess(phi,flam,h-d1h))/(2.d0*d1h)
endif
return; end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ggtess(vtess,phi,flam,h,v,gPhi,gLam,gH, &
    ggPhiPhi,ggPhiLam,ggPhiH,ggLamLam,ggLamH,ggHH)
!
!    Compute the gravity gradient tensor by the conditional switch
!    of the central and single-sided second order finite difference formulas
!
implicit integer (i-n)
implicit real*8 (a-h,o-z)
!
!    Common blocks for passing the parameters and indirect variables
!
real*8 R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
common /paramV/R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
real*8 delta,delta1,delta2,d1phi,d1lam,d1h,d2phi,d2lam,d2h,delta3,d3phi,d3lam,d3h
common /deltaV/delta,delta1,delta2,d1phi,d1lam,d1h,d2phi,d2lam,d2h,delta3,d3phi,d3lam,d3h
!
!    Compute conmponent-independent variables beorehand
!
r=R0+h
c=cos(phi)
t=tan(phi)
!
!   Compute ggPhiPhi
!
if((flam.ge.FlamW.and.flam.le.FlamE).and.(h.ge.HB.and.h.le.HT)) then
    if((phi.gt.PhiS.and.phi-d2phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d2phi.lt.PhiN)) then
        ggPhiPhi=(-vtess(phi+3.d0*d2phi,flam,h)+4.d0*vtess(phi+2.d0*d2phi,flam,h) &
            -5.d0*vtess(phi+d2phi,flam,h)+2.d0*v)/(d2phi*d2phi)
    elseif((phi.lt.PhiS.and.phi+d2phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d2phi.gt.PhiN)) then
        ggPhiPhi=(-vtess(phi-3.d0*d2phi,flam,h)+4.d0*vtess(phi-2.d0*d2phi,flam,h) &
            -5.d0*vtess(phi-d2phi,flam,h)+2.d0*v)/(d2phi*d2phi)
    else
        ggPhiPhi=(vtess(phi+d2phi,flam,h)-2.d0*v+vtess(phi-d2phi,flam,h))/(d2phi*d2phi)
    endif
else
    ggPhiPhi=(vtess(phi+d2phi,flam,h)-2.d0*v+vtess(phi-d2phi,flam,h))/(d2phi*d2phi)
endif
ggPhiPhi=ggPhiPhi/(r*r)+gH/r
!
!   Compute ggLamLam
!
if((phi.ge.PhiS.and.phi.le.PhiN).and.(h.ge.HB.and.h.le.HT)) then
    if((flam.gt.FlamW.and.flam-d2lam.lt.FlamW).or.(flam.gt.FlamE.and.flam-d2lam.lt.FlamE)) then
        ggLamLam=(-vtess(phi,flam+3.d0*d2lam,h)+4.d0*vtess(phi,flam+2.d0*d2lam,h) &
            -5.d0*vtess(phi,flam+d2lam,h)+2.d0*v)/(d2lam*d2lam)
    elseif((flam.lt.FlamW.and.flam+d2lam.gt.FlamW).or.(flam.lt.FlamE.and.flam+d2lam.gt.FlamE)) then
        ggLamLam=(-vtess(phi,flam-3.d0*d2lam,h)+4.d0*vtess(phi,flam-2.d0*d2lam,h) &
            -5.d0*vtess(phi,flam-d2lam,h)+2.d0*v)/(d2lam*d2lam)
    else
        ggLamLam=(vtess(phi,flam+d2lam,h)-2.d0*v+vtess(phi,flam-d2lam,h))/(d2lam*d2lam)
    endif
else
    ggLamLam=(vtess(phi,flam+d2lam,h)-2.d0*v+vtess(phi,flam-d2lam,h))/(d2lam*d2lam)
endif
ggLamLam=ggLamLam/(r*r*c*c)+(gH-gPhi*t)/r
!
!   Compute ggHH
!
if((phi.ge.PhiS.and.phi.le.PhiN).and.(flam.ge.FlamW.and.flam.le.FlamE)) then
    if((h.gt.HB.and.h-d2h.lt.HB).or.(h.gt.HT.and.h-d2h.lt.HT)) then
        ggHH=(-vtess(phi,flam,h+3.d0*d2h)+4.d0*vtess(phi,flam,h+2.d0*d2h) &
            -5.d0*vtess(phi,flam,h+d2h)+2.d0*v)/(d2h*d2h)
    elseif((h.lt.HB.and.h+d2h.gt.HB).or.(h.lt.HT.and.h+d2h.gt.HT)) then
        ggHH=(-vtess(phi,flam,h-3.d0*d2h)+4.d0*vtess(phi,flam,h-2.d0*d2h) &
            -5.d0*vtess(phi,flam,h-d2h)+2.d0*v)/(d2h*d2h)
    else
        ggHH=(vtess(phi,flam,h+d2h)-2.d0*v+vtess(phi,flam,h-d2h))/(d2h*d2h)
    endif
else
    ggHH=(vtess(phi,flam,h+d2h)-2.d0*v+vtess(phi,flam,h-d2h))/(d2h*d2h)
endif
!
!   Compute ggPhiLam
!
if(h.ge.HB.and.h.le.HT) then
    if(((phi.gt.PhiS.and.phi-d2phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d2phi.lt.PhiN)).and. &
        ((flam.gt.FlamW.and.flam-d2lam.lt.FlamW).or.(flam.gt.FlamE.and.flam-d2lam.lt.FlamE))) then
        ggPhiLam=(vtess(phi+2.d0*d2phi,flam+2.d0*d2lam,h) &
            -4.d0*vtess(phi+2.d0*d2phi,flam+d2lam,h)+3.d0*vtess(phi+2.d0*d2phi,flam,h) &
            -4.d0*vtess(phi+d2phi,flam+2.d0*d2lam,h)+16.d0*vtess(phi+d2phi,flam+d2lam,h) &
            -12.d0*vtess(phi+d2phi,flam,h)+3.d0*vtess(phi,flam+2.d0*d2lam,h) &
            -12.d0*vtess(phi,flam+d2lam,h)+9.d0*v)/(4.d0*d2phi*d2lam)
    elseif(((phi.gt.PhiS.and.phi-d2phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d2phi.lt.PhiN)).and. &
        ((flam.lt.FlamW.and.flam+d2lam.gt.FlamW).or.(flam.lt.FlamE.and.flam+d2lam.gt.FlamE))) then
        ggPhiLam=-(vtess(phi+2.d0*d2phi,flam-2.d0*d2lam,h) &
            -4.d0*vtess(phi+2.d0*d2phi,flam-d2lam,h)+3.d0*vtess(phi+2.d0*d2phi,flam,h) &
            -4.d0*vtess(phi+d2phi,flam-2.d0*d2lam,h)+16.d0*vtess(phi+d2phi,flam-d2lam,h) &
            -12.d0*vtess(phi+d2phi,flam,h)+3.d0*vtess(phi,flam-2.d0*d2lam,h) &
            -12.d0*vtess(phi,flam-d2lam,h)+9.d0*v)/(4.d0*d2phi*d2lam)
    elseif(((phi.lt.PhiS.and.phi+d2phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d2phi.gt.PhiN)).and. &
        ((flam.gt.FlamW.and.flam-d2lam.lt.FlamW).or.(flam.gt.FlamE.and.flam-d2lam.lt.FlamE))) then
        ggPhiLam=-(vtess(phi-2.d0*d2phi,flam+2.d0*d2lam,h) &
            -4.d0*vtess(phi-2.d0*d2phi,flam+d2lam,h)+3.d0*vtess(phi-2.d0*d2phi,flam,h) &
            -4.d0*vtess(phi-d2phi,flam+2.d0*d2lam,h)+16.d0*vtess(phi-d2phi,flam+d2lam,h) &
            -12.d0*vtess(phi-d2phi,flam,h)+3.d0*vtess(phi,flam+2.d0*d2lam,h) &
            -12.d0*vtess(phi,flam+d2lam,h)+9.d0*v)/(4.d0*d2phi*d2lam)
    elseif(((phi.lt.PhiS.and.phi+d2phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d2phi.gt.PhiN)).and. &
        ((flam.lt.FlamW.and.flam+d2lam.gt.FlamW).or.(flam.lt.FlamE.and.flam+d2lam.gt.FlamE))) then
        ggPhiLam=(vtess(phi-2.d0*d2phi,flam-2.d0*d2lam,h) &
            -4.d0*vtess(phi-2.d0*d2phi,flam-d2lam,h)+3.d0*vtess(phi-2.d0*d2phi,flam,h) &
            -4.d0*vtess(phi-d2phi,flam-2.d0*d2lam,h)+16.d0*vtess(phi-d2phi,flam-d2lam,h) &
            -12.d0*vtess(phi-d2phi,flam,h)+3.d0*vtess(phi,flam-2.d0*d2lam,h) &
            -12.d0*vtess(phi,flam-d2lam,h)+9.d0*v)/(4.d0*d2phi*d2lam)
    elseif(((phi.gt.PhiS.and.phi-d2phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d2phi.lt.PhiN)).and. &
        ((flam+d2lam.lt.FlamW).or.(flam-d2lam.gt.FlamE).or. &
        (flam-d2lam.gt.FlamW.and.flam+d2lam.lt.FlamE))) then
        ggPhiLam=(-vtess(phi+2.d0*d2phi,flam+d2lam,h)+vtess(phi+2.d0*d2phi,flam-d2lam,h) &
            +4.d0*vtess(phi+d2phi,flam+d2lam,h)-4.d0*vtess(phi+d2phi,flam-d2lam,h) &
            -3.d0*vtess(phi,flam+d2lam,h)+3.d0*vtess(phi,flam-d2lam,h))/(4.d0*d2phi*d2lam)
    elseif(((phi.lt.PhiS.and.phi+d2phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d2phi.gt.PhiN)).and. &
        ((flam+d2lam.lt.FlamW).or.(flam-d2lam.gt.FlamE).or. &
        (flam-d2lam.gt.FlamW.and.flam+d2lam.lt.FlamE))) then
        ggPhiLam=(vtess(phi-2.d0*d2phi,flam+d2lam,h)-vtess(phi-2.d0*d2phi,flam-d2lam,h) &
            -4.d0*vtess(phi-d2phi,flam+d2lam,h)+4.d0*vtess(phi-d2phi,flam-d2lam,h) &
            +3.d0*vtess(phi,flam+d2lam,h)-3.d0*vtess(phi,flam-d2lam,h))/(4.d0*d2phi*d2lam)
    elseif(((flam.gt.FlamW.and.flam-d2lam.lt.FlamW).or. &
        (flam.gt.FlamE.and.flam-d2lam.lt.FlamE)).and. &
        ((phi+d2phi.lt.PhiS).or.(phi-d2phi.gt.PhiN).or.(phi-d2phi.gt.PhiS.and.phi+d2phi.lt.PhiN))) then
        ggPhiLam=(-vtess(phi+d2phi,flam+2.d0*d2lam,h)+vtess(phi-d2phi,flam+2.d0*d2lam,h) &
            +4.d0*vtess(phi+d2phi,flam+d2lam,h)-4.d0*vtess(phi-d2phi,flam+d2lam,h) &
            -3.d0*vtess(phi+d2phi,flam,h)+3.d0*vtess(phi-d2phi,flam,h))/(4.d0*d2phi*d2lam)
    elseif(((flam.lt.FlamW.and.flam+d2lam.gt.FlamW).or. &
        (flam.lt.FlamE.and.flam+d2lam.gt.FlamE)).and. &
        ((phi+d2phi.lt.PhiS).or.(phi-d2phi.gt.PhiN).or.(phi-d2phi.gt.PhiS.and.phi+d2phi.lt.PhiN))) then
        ggPhiLam=-(-vtess(phi+d2phi,flam-2.d0*d2lam,h)+vtess(phi-d2phi,flam-2.d0*d2lam,h) &
            +4.d0*vtess(phi+d2phi,flam-d2lam,h)-4.d0*vtess(phi-d2phi,flam-d2lam,h) &
            -3.d0*vtess(phi+d2phi,flam,h)+3.d0*vtess(phi-d2phi,flam,h))/(4.d0*d2phi*d2lam)
    else
        ggPhiLam=(vtess(phi+d2phi,flam+d2lam,h)-vtess(phi-d2phi,flam+d2lam,h) &
            -vtess(phi+d2phi,flam-d2lam,h)+vtess(phi-d2phi,flam-d2lam,h))/(4.d0*d2phi*d2lam)
    endif
else
    ggPhiLam=(vtess(phi+d2phi,flam+d2lam,h)-vtess(phi-d2phi,flam+d2lam,h) &
        -vtess(phi+d2phi,flam-d2lam,h)+vtess(phi-d2phi,flam-d2lam,h))/(4.d0*d2phi*d2lam)
endif
ggPhiLam=ggPhiLam/(r*r*c)+gLam*t/r
!
!   Compute ggPhiH
!
if(flam.ge.FlamW.and.flam.le.FlamE) then
    if(((phi.gt.PhiS.and.phi-d2phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d2phi.lt.PhiN)).and. &
        ((h.gt.HB.and.h-d2h.lt.HB).or.(h.gt.HT.and.h-d2h.lt.HT))) then
        ggPhiH=(vtess(phi+2.d0*d2phi,flam,h+2.d0*d2h)-4.d0*vtess(phi+2.d0*d2phi,flam,h+d2h) &
            +3.d0*vtess(phi+2.d0*d2phi,flam,h)-4.d0*vtess(phi+d2phi,flam,h+2.d0*d2h) &
            +16.d0*vtess(phi+d2phi,flam,h+d2h)-12.d0*vtess(phi+d2phi,flam,h) &
            +3.d0*vtess(phi,flam,h+2.d0*d2h)-12.d0*vtess(phi,flam,h+d2h)+9.d0*v)/(4.d0*d2phi*d2h)
    elseif(((phi.gt.PhiS.and.phi-d2phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d2phi.lt.PhiN)).and. &
        ((h.lt.HB.and.h+d2h.gt.HB).or.(h.lt.HT.and.h+d2h.gt.HT))) then
        ggPhiH=-(vtess(phi+2.d0*d2phi,flam,h-2.d0*d2h)-4.d0*vtess(phi+2.d0*d2phi,flam,h-d2h) &
            +3.d0*vtess(phi+2.d0*d2phi,flam,h)-4.d0*vtess(phi+d2phi,flam,h-2.d0*d2h) &
            +16.d0*vtess(phi+d2phi,flam,h-d2h)-12.d0*vtess(phi+d2phi,flam,h) &
            +3.d0*vtess(phi,flam,h-2.d0*d2h)-12.d0*vtess(phi,flam,h-d2h)+9.d0*v)/(4.d0*d2phi*d2h)
    elseif(((phi.lt.PhiS.and.phi+d2phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d2phi.gt.PhiN)).and. &
        ((h.gt.HB.and.h-d2h.lt.HB).or.(h.gt.HT.and.h-d2h.lt.HT))) then
        ggPhiH=-(vtess(phi-2.d0*d2phi,flam,h+2.d0*d2h)-4.d0*vtess(phi-2.d0*d2phi,flam,h+d2h) &
            +3.d0*vtess(phi-2.d0*d2phi,flam,h)-4.d0*vtess(phi-d2phi,flam,h+2.d0*d2h) &
            +16.d0*vtess(phi-d2phi,flam,h+d2h)-12.d0*vtess(phi-d2phi,flam,h) &
            +3.d0*vtess(phi,flam,h+2.d0*d2h)-12.d0*vtess(phi,flam,h+d2h)+9.d0*v)/(4.d0*d2phi*d2h)
    elseif(((phi.lt.PhiS.and.phi+d2phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d2phi.gt.PhiN)).and. &
        ((h.lt.HB.and.h+d2h.gt.HB).or.(h.lt.HT.and.h+d2h.gt.HT))) then
        ggPhiH=(vtess(phi-2.d0*d2phi,flam,h-2.d0*d2h)-4.d0*vtess(phi-2.d0*d2phi,flam,h-d2h) &
            +3.d0*vtess(phi-2.d0*d2phi,flam,h)-4.d0*vtess(phi-d2phi,flam,h-2.d0*d2h) &
            +16.d0*vtess(phi-d2phi,flam,h-d2h)-12.d0*vtess(phi-d2phi,flam,h) &
            +3.d0*vtess(phi,flam,h-2.d0*d2h)-12.d0*vtess(phi,flam,h-d2h)+9.d0*v)/(4.d0*d2phi*d2h)
    elseif(((phi.gt.PhiS.and.phi-d2phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d2phi.lt.PhiN)).and. &
        ((h+d2h.lt.HB).or.(h-d2h.gt.HT).or.(h-d2h.gt.HB.and.h+d2h.lt.HT))) then
        ggPhiH=(-vtess(phi+2.d0*d2phi,flam,h+d2h)+vtess(phi+2.d0*d2phi,flam,h-d2h) &
            +4.d0*vtess(phi+d2phi,flam,h+d2h)-4.d0*vtess(phi+d2phi,flam,h-d2h) &
            -3.d0*vtess(phi,flam,h+d2h)+3.d0*vtess(phi,flam,h-d2h))/(4.d0*d2phi*d2h)
    elseif(((phi.lt.PhiS.and.phi+d2phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d2phi.gt.PhiN)).and. &
        ((h+d2h.lt.HB).or.(h-d2h.gt.HT).or.(h-d2h.gt.HB.and.h+d2h.lt.HT))) then
        ggPhiH=-(-vtess(phi-2.d0*d2phi,flam,h+d2h)+vtess(phi-2.d0*d2phi,flam,h-d2h) &
            +4.d0*vtess(phi-d2phi,flam,h+d2h)-4.d0*vtess(phi-d2phi,flam,h-d2h) &
            -3.d0*vtess(phi,flam,h+d2h)+3.d0*vtess(phi,flam,h-d2h))/(4.d0*d2phi*d2h)
    elseif(((h.gt.HB.and.h-d2h.lt.HB).or.(h.gt.HT.and.h-d2h.lt.HT)).and.((phi+d2phi.lt.PhiS) &
        .or.(phi-d2phi.gt.PhiN).or.(phi-d2phi.gt.PhiS.and.phi+d2phi.lt.PhiN))) then
        ggPhiH=(-vtess(phi+d2phi,flam,h+2.d0*d2h)+vtess(phi-d2phi,flam,h+2.d0*d2h) &
            +4.d0*vtess(phi+d2phi,flam,h+d2h)-4.d0*vtess(phi-d2phi,flam,h+d2h) &
            -3.d0*vtess(phi+d2phi,flam,h)+3.d0*vtess(phi-d2phi,flam,h))/(4.d0*d2phi*d2h)
    elseif(((h.lt.HB.and.h+d2h.gt.HB).or.(h.lt.HT.and.h+d2h.gt.HT)).and.((phi+d2phi.lt.PhiS) &
        .or.(phi-d2phi.gt.PhiN).or.(phi-d2phi.gt.PhiS.and.phi+d2phi.lt.PhiN))) then
        ggPhiH=-(-vtess(phi+d2phi,flam,h-2.d0*d2h)+vtess(phi-d2phi,flam,h-2.d0*d2h) &
            +4.d0*vtess(phi+d2phi,flam,h-d2h)-4.d0*vtess(phi-d2phi,flam,h-d2h) &
            -3.d0*vtess(phi+d2phi,flam,h)+3.d0*vtess(phi-d2phi,flam,h))/(4.d0*d2phi*d2h)
    else
        ggPhiH=(vtess(phi+d2phi,flam,h+d2h)-vtess(phi-d2phi,flam,h+d2h) &
            -vtess(phi+d2phi,flam,h-d2h)+vtess(phi-d2phi,flam,h-d2h))/(4.d0*d2phi*d2h)
    endif
else
    ggPhiH=(vtess(phi+d2phi,flam,h+d2h)-vtess(phi-d2phi,flam,h+d2h) &
        -vtess(phi+d2phi,flam,h-d2h)+vtess(phi-d2phi,flam,h-d2h))/(4.d0*d2phi*d2h)
endif
ggPhiH=ggPhiH/r-gPhi/r
!
!   Compute ggLamH
!
if(phi.ge.PhiS.and.phi.le.PhiN) then
    if(((flam.gt.FlamW.and.flam-d2flam.lt.FlamW).or.(flam.gt.FlamE.and.flam-d2flam.lt.FlamE)).and. &
        ((h.gt.HB.and.h-d2h.lt.HB).or.(h.gt.HT.and.h-d2h.lt.HT))) then
        ggLamH=(vtess(phi,flam+2.d0*d2lam,h+2.d0*d2h)-4.d0*vtess(phi,flam+2.d0*d2lam,h+d2h) &
            +3.d0*vtess(phi,flam+2.d0*d2lam,h)-4.d0*vtess(phi,flam+d2lam,h+2.d0*d2h) &
            +16.d0*vtess(phi,flam+d2lam,h+d2h)-12.d0*vtess(phi,flam+d2lam,h) &
            +3.d0*vtess(phi,flam,h+2.d0*d2h)-12.d0*vtess(phi,flam,h+d2h)+9.d0*v)/(4.d0*d2lam*d2h)
    elseif(((flam.gt.FlamW.and.flam-d2lam.lt.FlamW).or.(flam.gt.FlamE.and.flam-d2lam.lt.FlamE)).and. &
        ((h.lt.HB.and.h+d2h.gt.HB).or.(h.lt.HT.and.h+d2h.gt.HT))) then
        ggLamH=-(vtess(phi,flam+2.d0*d2lam,h-2.d0*d2h)-4.d0*vtess(phi,flam+2.d0*d2lam,h-d2h) &
            +3.d0*vtess(phi,flam+2.d0*d2lam,h)-4.d0*vtess(phi,flam+d2lam,h-2.d0*d2h) &
            +16.d0*vtess(phi,flam+d2lam,h-d2h)-12.d0*vtess(phi,flam+d2lam,h) &
            +3.d0*vtess(phi,flam,h-2.d0*d2h)-12.d0*vtess(phi,flam,h-d2h)+9.d0*v)/(4.d0*d2lam*d2h)
    elseif(((flam.lt.FlamW.and.flam+d2lam.gt.FlamW).or.(flam.lt.FlamE.and.flam+d2lam.gt.FlamE)).and. &
        ((h.gt.HB.and.h-d2h.lt.HB).or.(h.gt.HT.and.h-d2h.lt.HT))) then
        ggLamH=-(vtess(phi,flam-2.d0*d2lam,h+2.d0*d2h)-4.d0*vtess(phi,flam-2.d0*d2lam,h+d2h) &
            +3.d0*vtess(phi,flam-2.d0*d2lam,h)-4.d0*vtess(phi,flam-d2lam,h+2.d0*d2h) &
            +16.d0*vtess(phi,flam-d2lam,h+d2h)-12.d0*vtess(phi,flam-d2lam,h) &
            +3.d0*vtess(phi,flam,h+2.d0*d2h)-12.d0*vtess(phi,flam,h+d2h)+9.d0*v)/(4.d0*d2lam*d2h)
    elseif(((flam.lt.FlamW.and.flam+d2lam.gt.FlamW).or.(flam.lt.FlamE.and.flam+d2lam.gt.FlamE)).and. &
        ((h.lt.HB.and.h+d2h.gt.HB).or.(h.lt.HT.and.h+d2h.gt.HT))) then
        ggLamH=(vtess(phi,flam-2.d0*d2lam,h-2.d0*d2h)-4.d0*vtess(phi,flam-2.d0*d2lam,h-d2h) &
            +3.d0*vtess(phi,flam-2.d0*d2lam,h)-4.d0*vtess(phi,flam-d2lam,h-2.d0*d2h) &
            +16.d0*vtess(phi,flam-d2lam,h-d2h)-12.d0*vtess(phi,flam-d2lam,h) &
            +3.d0*vtess(phi,flam,h-2.d0*d2h)-12.d0*vtess(phi,flam,h-d2h)+9.d0*v)/(4.d0*d2lam*d2h)
    elseif(((flam.gt.FlamW.and.flam-d2lam.lt.FlamW).or.(flam.gt.FlamE.and.flam-d2lam.lt.FlamE)).and. &
        ((h+d2h.lt.HB).or.(h-d2h.gt.HT).or.(h-d2h.gt.HB.and.h+d2h.lt.HT))) then
        ggLamH=(-vtess(phi,flam+2.d0*d2lam,h+d2h)+vtess(phi,flam+2.d0*d2lam,h-d2h) &
            +4.d0*vtess(phi,flam+d2lam,h+d2h)-4.d0*vtess(phi,flam+d2lam,h-d2h) &
            -3.d0*vtess(phi,flam,h+d2h)+3.d0*vtess(phi,flam,h-d2h))/(4.d0*d2lam*d2h)
    elseif(((flam.lt.FlamW.and.flam+d2lam.gt.FlamW).or.(flam.lt.FlamE.and.flam+d2lam.gt.FlamE)).and. &
        ((h+d2h.lt.HB).or.(h-d2h.gt.HT).or.(h-d2h.gt.HB.and.h+d2h.lt.HT))) then
        ggLamH=-(-vtess(phi,flam-2.d0*d2lam,h+d2h)+vtess(phi,flam-2.d0*d2lam,h-d2h) &
            +4.d0*vtess(phi,flam-d2lam,h+d2h)-4.d0*vtess(phi,flam-d2lam,h-d2h) &
            -3.d0*vtess(phi,flam,h+d2h)+3.d0*vtess(phi,flam,h-d2h))/(4.d0*d2lam*d2h)
    elseif(((h.gt.HB.and.h-d2h.lt.HB).or.(h.gt.HT.and.h-d2h.lt.HT)).and.((flam+d2lam.lt.FlamW).or. &
        (flam-d2lam.gt.FlamE).or.(flam-d2lam.gt.FlamW.and.flam+d2lam.lt.FlamE))) then
        ggLamH=(-vtess(phi,flam+d2lam,h+2.d0*d2h)+vtess(phi,flam-d2lam,h+2.d0*d2h) &
            +4.d0*vtess(phi,flam+d2lam,h+d2h)-4.d0*vtess(phi,flam-d2lam,h+d2h) &
            -3.d0*vtess(phi,flam+d2lam,h)+3.d0*vtess(phi,flam-d2lam,h))/(4.d0*d2lam*d2h)
    elseif(((h.lt.HB.and.h+d2h.gt.HB).or.(h.lt.HT.and.h+d2h.gt.HT)).and.((flam+d2lam.lt.FlamW).or. &
        (flam-d2lam.gt.FlamE).or.(flam-d2lam.gt.FlamW.and.flam+d2lam.lt.FlamE))) then
        ggLamH=-(-vtess(phi,flam+d2lam,h-2.d0*d2h)+vtess(phi,flam-d2lam,h-2.d0*d2h) &
            +4.d0*vtess(phi,flam+d2lam,h-d2h)-4.d0*vtess(phi,flam-d2lam,h-d2h) &
            -3.d0*vtess(phi,flam+d2lam,h)+3.d0*vtess(phi,flam-d2lam,h))/(4.d0*d2lam*d2h)
    else
        ggLamH=(vtess(phi,flam+d2lam,h+d2h)-vtess(phi,flam-d2lam,h+d2h) &
            -vtess(phi,flam+d2lam,h-d2h)+vtess(phi,flam-d2lam,h-d2h))/(4.d0*d2lam*d2h)
    endif
else
    ggLamH=(vtess(phi,flam+d2lam,h+d2h)-vtess(phi,flam-d2lam,h+d2h) &
        -vtess(phi,flam+d2lam,h-d2h)+vtess(phi,flam-d2lam,h-d2h))/(4.d0*d2lam*d2h)
endif
ggLamH=ggLamH/(r*c)-gLam/r
return; end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gggtess(vtess,phi,flam,h,v,gPhi,gLam,gH, &
    ggPhiPhi,ggPhiLam,ggPhiH,ggLamLam,ggLamH,ggHH, &
    gggPhiPhiPhi,gggPhiPhiLam,gggPhiPhiH,gggPhiLamH,gggLamLamPhi, &
    gggLamLamLam,gggLamLamH,gggHHPhi,gggHHLam,gggHHH)
!
!   Compute the gravitational curvatures by the conditional switch
!   of the third-order central and single-sided finite difference formulas
!
implicit integer (i-n)
implicit real*8 (a-h,o-z)
!
!    Common blocks for passing the parameters and indirect variables
!
real*8 R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta,r,s,c,t
common /paramV/R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta,r,s,c,t
real*8 delta,delta1,delta2,d1phi,d1lam,d1h,d2phi,d2lam,d2h,delta3,d3phi,d3lam,d3h
common /deltaV/delta,delta1,delta2,d1phi,d1lam,d1h,d2phi,d2lam,d2h,delta3,d3phi,d3lam,d3h
!
!    Compute conmponent-independent variables beorehand
!
r=R0+h
s=sin(phi)
c=cos(phi)
t=tan(phi)
!
!   Compute gggPhiPhiPhi
!
if((flam.ge.FlamW.and.flam.le.FlamE).and.(h.ge.HB.and.h.le.HT)) then
    if((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) then
        gggPhiPhiPhi=(-3.d0*vtess(phi+4.d0*d3phi,flam,h)+14.d0*vtess(phi+3.d0*d3phi,flam,h) &
            -24.d0*vtess(phi+2.d0*d3phi,flam,h)+18.d0*vtess(phi+d3phi,flam,h) &
            -5.0d0*v)/(2.d0*d3phi*d3phi*d3phi)
    elseif((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) then
        gggPhiPhiPhi=(3.d0*vtess(phi-4.d0*d3phi,flam,h)-14.d0*vtess(phi-3.d0*d3phi,flam,h) &
            +24.d0*vtess(phi-2.d0*d3phi,flam,h)-18.d0*vtess(phi-d3phi,flam,h) &
            +5.0d0*v)/(2.d0*d3phi*d3phi*d3phi)
    else
        gggPhiPhiPhi=(vtess(phi+2.d0*d3phi,flam,h)-2.d0*vtess(phi+d3phi,flam,h) &
            +2.d0*vtess(phi-d3phi,flam,h)-vtess(phi-2.d0*d3phi,flam,h) &
            )/(2.d0*d3phi*d3phi*d3phi)
    endif
else
    gggPhiPhiPhi=(vtess(phi+2.d0*d3phi,flam,h)-2.d0*vtess(phi+d3phi,flam,h) &
            +2.d0*vtess(phi-d3phi,flam,h)-vtess(phi-2.d0*d3phi,flam,h) &
            )/(2.d0*d3phi*d3phi*d3phi)
endif
gggPhiPhiPhi=(gggPhiPhiPhi+r*(gPhi + 3.d0*r*ggPhiH))/(r*r*r)
!
!   Compute gggPhiPhiLam
!
if(h.ge.HB.and.h.le.HT) then
! case 1
    if(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE))) then
        gggPhiPhiLam=(vtess(phi+3.d0*d3phi,flam+2.d0*d3lam,h)-4.d0*vtess(phi+3.d0*d3phi,flam+d3lam,h) &
            +3.d0*vtess(phi+3.d0*d3phi,flam,h)-4.d0*vtess(phi+2.d0*d3phi,flam+2.d0*d3lam,h) &
            +16.d0*vtess(phi+2.d0*d3phi,flam+d3lam,h)-12.d0*vtess(phi+2.d0*d3phi,flam,h) &
            +5.d0*vtess(phi+d3phi,flam+2.d0*d3lam,h)-20.d0*vtess(phi+d3phi,flam+d3lam,h) &
            +15.d0*vtess(phi+d3phi,flam,h)-2.d0*vtess(phi,flam+2.d0*d3lam,h) &
            +8.d0*vtess(phi,flam+d3lam,h)-6.d0*v &
            )/(2.d0*d3phi*d3phi*d3lam)
! case 2
    elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE))) then
        gggPhiPhiLam=-(vtess(phi+3.d0*d3phi,flam-2.d0*d3lam,h)-4.d0*vtess(phi+3.d0*d3phi,flam-d3lam,h) &
            +3.d0*vtess(phi+3.d0*d3phi,flam,h)-4.d0*vtess(phi+2.d0*d3phi,flam-2.d0*d3lam,h) &
            +16.d0*vtess(phi+2.d0*d3phi,flam-d3lam,h)-12.d0*vtess(phi+2.d0*d3phi,flam,h) &
            +5.d0*vtess(phi+d3phi,flam-2.d0*d3lam,h)-20.d0*vtess(phi+d3phi,flam-d3lam,h) &
            +15.d0*vtess(phi+d3phi,flam,h)-2.d0*vtess(phi,flam-2.d0*d3lam,h) &
            +8.d0*vtess(phi,flam-d3lam,h)-6.d0*v &
            )/(2.d0*d3phi*d3phi*d3lam)
! case 3
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE))) then
        gggPhiPhiLam=(vtess(phi-3.d0*d3phi,flam+2.d0*d3lam,h)-4.d0*vtess(phi-3.d0*d3phi,flam+d3lam,h) &
            +3.d0*vtess(phi-3.d0*d3phi,flam,h)-4.d0*vtess(phi-2.d0*d3phi,flam+2.d0*d3lam,h) &
            +16.d0*vtess(phi-2.d0*d3phi,flam+d3lam,h)-12.d0*vtess(phi-2.d0*d3phi,flam,h) &
            +5.d0*vtess(phi-d3phi,flam+2.d0*d3lam,h)-20.d0*vtess(phi-d3phi,flam+d3lam,h) &
            +15.d0*vtess(phi-d3phi,flam,h)-2.d0*vtess(phi,flam+2.d0*d3lam,h) &
            +8.d0*vtess(phi,flam+d3lam,h)-6.d0*v &
            )/(2.d0*d3phi*d3phi*d3lam)
! case 4
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE))) then
        gggPhiPhiLam=-(vtess(phi-3.d0*d3phi,flam-2.d0*d3lam,h)-4.d0*vtess(phi-3.d0*d3phi,flam-d3lam,h) &
            +3.d0*vtess(phi-3.d0*d3phi,flam,h)-4.d0*vtess(phi-2.d0*d3phi,flam-2.d0*d3lam,h) &
            +16.d0*vtess(phi-2.d0*d3phi,flam-d3lam,h)-12.d0*vtess(phi-2.d0*d3phi,flam,h) &
            +5.d0*vtess(phi-d3phi,flam-2.d0*d3lam,h)-20.d0*vtess(phi-d3phi,flam-d3lam,h) &
            +15.d0*vtess(phi-d3phi,flam,h)-2.d0*vtess(phi,flam-2.d0*d3lam,h) &
            +8.d0*vtess(phi,flam-d3lam,h)-6.d0*v &
            )/(2.d0*d3phi*d3phi*d3lam)
! case 5
    elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE))) then
        gggPhiPhiLam=(-vtess(phi+3.d0*d3phi,flam+d3lam,h)+vtess(phi+3.d0*d3phi,flam-d3lam,h) &
            +4.d0*vtess(phi+2.d0*d3phi,flam+d3lam,h)-4.d0*vtess(phi+2.d0*d3phi,flam-d3lam,h) &
            -5.d0*vtess(phi+d3phi,flam+d3lam,h)+5.d0*vtess(phi+d3phi,flam-d3lam,h) &
            +2.d0*vtess(phi,flam+d3lam,h)-2.d0*vtess(phi,flam-d3lam,h) &
            )/(2.d0*d3phi*d3phi*d3lam)
! case 6
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE))) then
        gggPhiPhiLam=(-vtess(phi-3.d0*d3phi,flam+d3lam,h)+vtess(phi-3.d0*d3phi,flam-d3lam,h) &
            +4.d0*vtess(phi-2.d0*d3phi,flam+d3lam,h)-4.d0*vtess(phi-2.d0*d3phi,flam-d3lam,h) &
            -5.d0*vtess(phi-d3phi,flam+d3lam,h)+5.d0*vtess(phi-d3phi,flam-d3lam,h) &
            +2.d0*vtess(phi,flam+d3lam,h)-2.d0*vtess(phi,flam-d3lam,h) &
            )/(2.d0*d3phi*d3phi*d3lam)
! case 7
    elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE))) then
        gggPhiPhiLam=(-vtess(phi+d3phi,flam+2.d0*d3lam,h)+4.d0*vtess(phi+d3phi,flam+d3lam,h) &
            -3.d0*vtess(phi+d3phi,flam,h)+2.d0*vtess(phi,flam+2.d0*d3lam,h) &
            -8.d0*vtess(phi,flam+d3lam,h)+6.d0*v-vtess(phi-d3phi,flam+2.d0*d3lam,h) &
            +4.d0*vtess(phi-d3phi,flam+d3lam,h)-3.d0*vtess(phi-d3phi,flam,h) &
            )/(2.d0*d3phi*d3phi*d3lam)
! case 8
    elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE))) then
        gggPhiPhiLam=(-vtess(phi+d3phi,flam-2.d0*d3lam,h)+4.d0*vtess(phi+d3phi,flam-d3lam,h) &
            -3.d0*vtess(phi+d3phi,flam,h)+2.d0*vtess(phi,flam-2.d0*d3lam,h) &
            -8.d0*vtess(phi,flam-d3lam,h)+6.d0*v-vtess(phi-d3phi,flam-2.d0*d3lam,h) &
            +4.d0*vtess(phi-d3phi,flam-d3lam,h)-3.d0*vtess(phi-d3phi,flam,h) &
            )/(2.d0*d3phi*d3phi*d3lam)
! case 9
    else
        gggPhiPhiLam=(vtess(phi+d3phi,flam+d3lam,h)-vtess(phi+d3phi,flam-d3lam,h) &
            -2.d0*vtess(phi,flam+d3lam,h)+2.d0*vtess(phi,flam-d3lam,h) &
            +vtess(phi-d3phi,flam+d3lam,h)-vtess(phi-d3phi,flam-d3lam,h) &
            )/(2.d0*d3phi*d3phi*d3lam)
    endif
else
    gggPhiPhiLam=(vtess(phi+d3phi,flam+d3lam,h)-vtess(phi+d3phi,flam-d3lam,h) &
            -2.d0*vtess(phi,flam+d3lam,h)+2.d0*vtess(phi,flam-d3lam,h) &
            +vtess(phi-d3phi,flam+d3lam,h)-vtess(phi-d3phi,flam-d3lam,h) &
            )/(2.d0*d3phi*d3phi*d3lam)
endif
gggPhiPhiLam=gggPhiPhiLam/(r*r*r*c)+(gLam+r*ggLamH+2.d0*r*t*ggPhiLam)/(r*r)
!
!   Compute gggPhiPhiH
!
if(flam.ge.FlamW.and.flam.le.FlamE) then
! case 1
    if(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggPhiPhiH=(vtess(phi+3.d0*d3phi,flam,h+2.d0*d3h)-4.d0*vtess(phi+3.d0*d3phi,flam,h+d3h) &
            +3.d0*vtess(phi+3.d0*d3phi,flam,h)-4.d0*vtess(phi+2.d0*d3phi,flam,h+2.d0*d3h) &
            +16.d0*vtess(phi+2.d0*d3phi,flam,h+d3h)-12.d0*vtess(phi+2.d0*d3phi,flam,h) &
            +5.d0*vtess(phi+d3phi,flam,h+2.d0*d3h)-20.d0*vtess(phi+d3phi,flam,h+d3h) &
            +15.d0*vtess(phi+d3phi,flam,h)-2.d0*vtess(phi,flam,h+2.d0*d3h) &
            +8.d0*vtess(phi,flam,h+d3h)-6.d0*v &
            )/(2.d0*d3phi*d3phi*d3h)
! case 2
    elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggPhiPhiH=-(vtess(phi+3.d0*d3phi,flam,h-2.d0*d3h)-4.d0*vtess(phi+3.d0*d3phi,flam,h-d3h) &
            +3.d0*vtess(phi+3.d0*d3phi,flam,h)-4.d0*vtess(phi+2.d0*d3phi,flam,h-2.d0*d3h) &
            +16.d0*vtess(phi+2.d0*d3phi,flam,h-d3h)-12.d0*vtess(phi+2.d0*d3phi,flam,h) &
            +5.d0*vtess(phi+d3phi,flam,h-2.d0*d3h)-20.d0*vtess(phi+d3phi,flam,h-d3h) &
            +15.d0*vtess(phi+d3phi,flam,h)-2.d0*vtess(phi,flam,h-2.d0*d3h) &
            +8.d0*vtess(phi,flam,h-d3h)-6.d0*v &
            )/(2.d0*d3phi*d3phi*d3h)
! case 3
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggPhiPhiH=(vtess(phi-3.d0*d3phi,flam,h+2.d0*d3h)-4.d0*vtess(phi-3.d0*d3phi,flam,h+d3h) &
            +3.d0*vtess(phi-3.d0*d3phi,flam,h)-4.d0*vtess(phi-2.d0*d3phi,flam,h+2.d0*d3h) &
            +16.d0*vtess(phi-2.d0*d3phi,flam,h+d3h)-12.d0*vtess(phi-2.d0*d3phi,flam,h) &
            +5.d0*vtess(phi-d3phi,flam,h+2.d0*d3h)-20.d0*vtess(phi-d3phi,flam,h+d3h) &
            +15.d0*vtess(phi-d3phi,flam,h)-2.d0*vtess(phi,flam,h+2.d0*d3h) &
            +8.d0*vtess(phi,flam,h+d3h)-6.d0*v &
            )/(2.d0*d3phi*d3phi*d3h)
! case 4
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggPhiPhiH=-(vtess(phi-3.d0*d3phi,flam,h-2.d0*d3h)-4.d0*vtess(phi-3.d0*d3phi,flam,h-d3h) &
            +3.d0*vtess(phi-3.d0*d3phi,flam,h)-4.d0*vtess(phi-2.d0*d3phi,flam,h-2.d0*d3h) &
            +16.d0*vtess(phi-2.d0*d3phi,flam,h-d3h)-12.d0*vtess(phi-2.d0*d3phi,flam,h) &
            +5.d0*vtess(phi-d3phi,flam,h-2.d0*d3h)-20.d0*vtess(phi-d3phi,flam,h-d3h) &
            +15.d0*vtess(phi-d3phi,flam,h)-2.d0*vtess(phi,flam,h-2.d0*d3h) &
            +8.d0*vtess(phi,flam,h-d3h)-6.d0*v &
            )/(2.d0*d3phi*d3phi*d3h)
! case 5
    elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
        gggPhiPhiH=(-vtess(phi+3.d0*d3phi,flam,h+d3h)+vtess(phi+3.d0*d3phi,flam,h-d3h) &
            +4.d0*vtess(phi+2.d0*d3phi,flam,h+d3h)-4.d0*vtess(phi+2.d0*d3phi,flam,h-d3h) &
            -5.d0*vtess(phi+d3phi,flam,h+d3h)+5.d0*vtess(phi+d3phi,flam,h-d3h) &
            +2.d0*vtess(phi,flam,h+d3h)-2.d0*vtess(phi,flam,h-d3h) &
            )/(2.d0*d3phi*d3phi*d3h)
! case 6
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
        gggPhiPhiH=(-vtess(phi-3.d0*d3phi,flam,h+d3h)+vtess(phi-3.d0*d3phi,flam,h-d3h) &
            +4.d0*vtess(phi-2.d0*d3phi,flam,h+d3h)-4.d0*vtess(phi-2.d0*d3phi,flam,h-d3h) &
            -5.d0*vtess(phi-d3phi,flam,h+d3h)+5.d0*vtess(phi-d3phi,flam,h-d3h) &
            +2.d0*vtess(phi,flam,h+d3h)-2.d0*vtess(phi,flam,h-d3h) &
            )/(2.d0*d3phi*d3phi*d3h)
! case 7
    elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggPhiPhiH=(-vtess(phi+d3phi,flam,h+2.d0*d3h)+4.d0*vtess(phi+d3phi,flam,h+d3h) &
            -3.d0*vtess(phi+d3phi,flam,h)+2.d0*vtess(phi,flam,h+2.d0*d3h) &
            -8.d0*vtess(phi,flam,h+d3h)+6.d0*v -vtess(phi-d3phi,flam,h+2.d0*d3h) &
            +4.d0*vtess(phi-d3phi,flam,h+d3h)-3.d0*vtess(phi-d3phi,flam,h) &
            )/(2.d0*d3phi*d3phi*d3h)
! case 8
    elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggPhiPhiH=-(-vtess(phi+d3phi,flam,h-2.d0*d3h)+4.d0*vtess(phi+d3phi,flam,h-d3h) &
            -3.d0*vtess(phi+d3phi,flam,h)+2.d0*vtess(phi,flam,h-2.d0*d3h) &
            -8.d0*vtess(phi,flam,h-d3h)+6.d0*v -vtess(phi-d3phi,flam,h-2.d0*d3h) &
            +4.d0*vtess(phi-d3phi,flam,h-d3h)-3.d0*vtess(phi-d3phi,flam,h) &
            )/(2.d0*d3phi*d3phi*d3h)
! case 9
    else
        gggPhiPhiH=(vtess(phi+d3phi,flam,h+d3h)-vtess(phi+d3phi,flam,h-d3h) &
            -2.d0*vtess(phi,flam,h+d3h)+2.d0*vtess(phi,flam,h-d3h) &
            +vtess(phi-d3phi,flam,h+d3h)-vtess(phi-d3phi,flam,h-d3h) &
            )/(2.d0*d3phi*d3phi*d3h)
    endif
else
    gggPhiPhiH=(vtess(phi+d3phi,flam,h+d3h)-vtess(phi+d3phi,flam,h-d3h) &
            -2.d0*vtess(phi,flam,h+d3h)+2.d0*vtess(phi,flam,h-d3h) &
            +vtess(phi-d3phi,flam,h+d3h)-vtess(phi-d3phi,flam,h-d3h) &
            )/(2.d0*d3phi*d3phi*d3h)
endif
gggPhiPhiH=(gggPhiPhiH + gH + r*ggHH - 2.d0*r*ggPhiPhi)/(r*r)
!
!   Compute gggPhiLamH
!
! case 1
if(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
    gggPhiLamH=(-vtess(phi+2.d0*d3phi,flam+2.d0*d3lam,h+2.d0*d3h) &
        +4.d0*vtess(phi+2.d0*d3phi,flam+2.d0*d3lam,h+d3h)-3.d0*vtess(phi+2.d0*d3phi,flam+2.d0*d3lam,h) &
        +4.d0*vtess(phi+2.d0*d3phi,flam+d3lam,h+2.d0*d3h)-16.d0*vtess(phi+2.d0*d3phi,flam+d3lam,h+d3h) &
        +12.d0*vtess(phi+2.d0*d3phi,flam+d3lam,h)-3.d0*vtess(phi+2.d0*d3phi,flam,h+2.d0*d3h) &
        +12.d0*vtess(phi+2.d0*d3phi,flam,h+d3h)-9.d0*vtess(phi+2.d0*d3phi,flam,h) &
        +4.d0*vtess(phi+d3phi,flam+2.d0*d3lam,h+2.d0*d3h)-16.d0*vtess(phi+d3phi,flam+2.d0*d3lam,h+d3h) &
        +12.d0*vtess(phi+d3phi,flam+2.d0*d3lam,h)-16.d0*vtess(phi+d3phi,flam+d3lam,h+2.d0*d3h) &
        +64.d0*vtess(phi+d3phi,flam+d3lam,h+d3h)-48.d0*vtess(phi+d3phi,flam+d3lam,h) &
        +12.d0*vtess(phi+d3phi,flam,h+2.d0*d3h)-48.d0*vtess(phi+d3phi,flam,h+d3h) &
        +36.d0*vtess(phi+d3phi,flam,h)-3.d0*vtess(phi,flam+2.d0*d3lam,h+2.d0*d3h) &
        +12.d0*vtess(phi,flam+2.d0*d3lam,h+d3h)-9.d0*vtess(phi,flam+2.d0*d3lam,h) &
        +12.d0*vtess(phi,flam+d3lam,h+2.d0*d3h)-48.d0*vtess(phi,flam+d3lam,h+d3h) &
        +36.d0*vtess(phi,flam+d3lam,h)-9.d0*vtess(phi,flam,h+2.d0*d3h) &
        +36.d0*vtess(phi,flam,h+d3h)-27.d0*v &
        )/(8.d0*d3phi*d3lam*d3h)
! case 2
elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
    gggPhiLamH=-(-vtess(phi+2.d0*d3phi,flam+2.d0*d3lam,h-2.d0*d3h) &
        +4.d0*vtess(phi+2.d0*d3phi,flam+2.d0*d3lam,h-d3h)-3.d0*vtess(phi+2.d0*d3phi,flam+2.d0*d3lam,h) &
        +4.d0*vtess(phi+2.d0*d3phi,flam+d3lam,h-2.d0*d3h)-16.d0*vtess(phi+2.d0*d3phi,flam+d3lam,h-d3h) &
        +12.d0*vtess(phi+2.d0*d3phi,flam+d3lam,h)-3.d0*vtess(phi+2.d0*d3phi,flam,h-2.d0*d3h) &
        +12.d0*vtess(phi+2.d0*d3phi,flam,h-d3h)-9.d0*vtess(phi+2.d0*d3phi,flam,h) &
        +4.d0*vtess(phi+d3phi,flam+2.d0*d3lam,h-2.d0*d3h)-16.d0*vtess(phi+d3phi,flam+2.d0*d3lam,h-d3h) &
        +12.d0*vtess(phi+d3phi,flam+2.d0*d3lam,h)-16.d0*vtess(phi+d3phi,flam+d3lam,h-2.d0*d3h) &
        +64.d0*vtess(phi+d3phi,flam+d3lam,h-d3h)-48.d0*vtess(phi+d3phi,flam+d3lam,h) &
        +12.d0*vtess(phi+d3phi,flam,h-2.d0*d3h)-48.d0*vtess(phi+d3phi,flam,h-d3h) &
        +36.d0*vtess(phi+d3phi,flam,h)-3.d0*vtess(phi,flam+2.d0*d3lam,h-2.d0*d3h) &
        +12.d0*vtess(phi,flam+2.d0*d3lam,h-d3h)-9.d0*vtess(phi,flam+2.d0*d3lam,h) &
        +12.d0*vtess(phi,flam+d3lam,h-2.d0*d3h)-48.d0*vtess(phi,flam+d3lam,h-d3h) &
        +36.d0*vtess(phi,flam+d3lam,h)-9.d0*vtess(phi,flam,h-2.d0*d3h) &
        +36.d0*vtess(phi,flam,h-d3h)-27.d0*v &
        )/(8.d0*d3phi*d3lam*d3h)
! case 3
elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
    gggPhiLamH=-(-vtess(phi+2.d0*d3phi,flam-2.d0*d3lam,h+2.d0*d3h) &
        +4.d0*vtess(phi+2.d0*d3phi,flam-2.d0*d3lam,h+d3h)-3.d0*vtess(phi+2.d0*d3phi,flam-2.d0*d3lam,h) &
        +4.d0*vtess(phi+2.d0*d3phi,flam-d3lam,h+2.d0*d3h)-16.d0*vtess(phi+2.d0*d3phi,flam-d3lam,h+d3h) &
        +12.d0*vtess(phi+2.d0*d3phi,flam-d3lam,h)-3.d0*vtess(phi+2.d0*d3phi,flam,h+2.d0*d3h) &
        +12.d0*vtess(phi+2.d0*d3phi,flam,h+d3h)-9.d0*vtess(phi+2.d0*d3phi,flam,h) &
        +4.d0*vtess(phi+d3phi,flam-2.d0*d3lam,h+2.d0*d3h)-16.d0*vtess(phi+d3phi,flam-2.d0*d3lam,h+d3h) &
        +12.d0*vtess(phi+d3phi,flam-2.d0*d3lam,h)-16.d0*vtess(phi+d3phi,flam-d3lam,h+2.d0*d3h) &
        +64.d0*vtess(phi+d3phi,flam-d3lam,h+d3h)-48.d0*vtess(phi+d3phi,flam-d3lam,h) &
        +12.d0*vtess(phi+d3phi,flam,h+2.d0*d3h)-48.d0*vtess(phi+d3phi,flam,h+d3h) &
        +36.d0*vtess(phi+d3phi,flam,h)-3.d0*vtess(phi,flam-2.d0*d3lam,h+2.d0*d3h) &
        +12.d0*vtess(phi,flam-2.d0*d3lam,h+d3h)-9.d0*vtess(phi,flam-2.d0*d3lam,h) &
        +12.d0*vtess(phi,flam-d3lam,h+2.d0*d3h)-48.d0*vtess(phi,flam-d3lam,h+d3h) &
        +36.d0*vtess(phi,flam-d3lam,h)-9.d0*vtess(phi,flam,h+2.d0*d3h) &
        +36.d0*vtess(phi,flam,h+d3h)-27.d0*v &
        )/(8.d0*d3phi*d3lam*d3h)
! case 4
elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
    gggPhiLamH=(-vtess(phi+2.d0*d3phi,flam-2.d0*d3lam,h-2.d0*d3h) &
        +4.d0*vtess(phi+2.d0*d3phi,flam-2.d0*d3lam,h-d3h)-3.d0*vtess(phi+2.d0*d3phi,flam-2.d0*d3lam,h) &
        +4.d0*vtess(phi+2.d0*d3phi,flam-d3lam,h-2.d0*d3h)-16.d0*vtess(phi+2.d0*d3phi,flam-d3lam,h-d3h) &
        +12.d0*vtess(phi+2.d0*d3phi,flam-d3lam,h)-3.d0*vtess(phi+2.d0*d3phi,flam,h-2.d0*d3h) &
        +12.d0*vtess(phi+2.d0*d3phi,flam,h-d3h)-9.d0*vtess(phi+2.d0*d3phi,flam,h) &
        +4.d0*vtess(phi+d3phi,flam-2.d0*d3lam,h-2.d0*d3h)-16.d0*vtess(phi+d3phi,flam-2.d0*d3lam,h-d3h) &
        +12.d0*vtess(phi+d3phi,flam-2.d0*d3lam,h)-16.d0*vtess(phi+d3phi,flam-d3lam,h-2.d0*d3h) &
        +64.d0*vtess(phi+d3phi,flam-d3lam,h-d3h)-48.d0*vtess(phi+d3phi,flam-d3lam,h) &
        +12.d0*vtess(phi+d3phi,flam,h-2.d0*d3h)-48.d0*vtess(phi+d3phi,flam,h-d3h) &
        +36.d0*vtess(phi+d3phi,flam,h)-3.d0*vtess(phi,flam-2.d0*d3lam,h-2.d0*d3h) &
        +12.d0*vtess(phi,flam-2.d0*d3lam,h-d3h)-9.d0*vtess(phi,flam-2.d0*d3lam,h) &
        +12.d0*vtess(phi,flam-d3lam,h-2.d0*d3h)-48.d0*vtess(phi,flam-d3lam,h-d3h) &
        +36.d0*vtess(phi,flam-d3lam,h)-9.d0*vtess(phi,flam,h-2.d0*d3h) &
        +36.d0*vtess(phi,flam,h-d3h)-27.d0*v &
        )/(8.d0*d3phi*d3lam*d3h)
! case 5
elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
    gggPhiLamH=-(-vtess(phi-2.d0*d3phi,flam+2.d0*d3lam,h+2.d0*d3h) &
        +4.d0*vtess(phi-2.d0*d3phi,flam+2.d0*d3lam,h+d3h)-3.d0*vtess(phi-2.d0*d3phi,flam+2.d0*d3lam,h) &
        +4.d0*vtess(phi-2.d0*d3phi,flam+d3lam,h+2.d0*d3h)-16.d0*vtess(phi-2.d0*d3phi,flam+d3lam,h+d3h) &
        +12.d0*vtess(phi-2.d0*d3phi,flam+d3lam,h)-3.d0*vtess(phi-2.d0*d3phi,flam,h+2.d0*d3h) &
        +12.d0*vtess(phi-2.d0*d3phi,flam,h+d3h)-9.d0*vtess(phi-2.d0*d3phi,flam,h) &
        +4.d0*vtess(phi-d3phi,flam+2.d0*d3lam,h+2.d0*d3h)-16.d0*vtess(phi-d3phi,flam+2.d0*d3lam,h+d3h) &
        +12.d0*vtess(phi-d3phi,flam+2.d0*d3lam,h)-16.d0*vtess(phi-d3phi,flam+d3lam,h+2.d0*d3h) &
        +64.d0*vtess(phi-d3phi,flam+d3lam,h+d3h)-48.d0*vtess(phi-d3phi,flam+d3lam,h) &
        +12.d0*vtess(phi-d3phi,flam,h+2.d0*d3h)-48.d0*vtess(phi-d3phi,flam,h+d3h) &
        +36.d0*vtess(phi-d3phi,flam,h)-3.d0*vtess(phi,flam+2.d0*d3lam,h+2.d0*d3h) &
        +12.d0*vtess(phi,flam+2.d0*d3lam,h+d3h)-9.d0*vtess(phi,flam+2.d0*d3lam,h) &
        +12.d0*vtess(phi,flam+d3lam,h+2.d0*d3h)-48.d0*vtess(phi,flam+d3lam,h+d3h) &
        +36.d0*vtess(phi,flam+d3lam,h)-9.d0*vtess(phi,flam,h+2.d0*d3h) &
        +36.d0*vtess(phi,flam,h+d3h)-27.d0*v &
        )/(8.d0*d3phi*d3lam*d3h)
! case 6
elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
    gggPhiLamH=(-vtess(phi-2.d0*d3phi,flam+2.d0*d3lam,h-2.d0*d3h) &
        +4.d0*vtess(phi-2.d0*d3phi,flam+2.d0*d3lam,h-d3h)-3.d0*vtess(phi-2.d0*d3phi,flam+2.d0*d3lam,h) &
        +4.d0*vtess(phi-2.d0*d3phi,flam+d3lam,h-2.d0*d3h)-16.d0*vtess(phi-2.d0*d3phi,flam+d3lam,h-d3h) &
        +12.d0*vtess(phi-2.d0*d3phi,flam+d3lam,h)-3.d0*vtess(phi-2.d0*d3phi,flam,h-2.d0*d3h) &
        +12.d0*vtess(phi-2.d0*d3phi,flam,h-d3h)-9.d0*vtess(phi-2.d0*d3phi,flam,h) &
        +4.d0*vtess(phi-d3phi,flam+2.d0*d3lam,h-2.d0*d3h)-16.d0*vtess(phi-d3phi,flam+2.d0*d3lam,h-d3h) &
        +12.d0*vtess(phi-d3phi,flam+2.d0*d3lam,h)-16.d0*vtess(phi-d3phi,flam+d3lam,h-2.d0*d3h) &
        +64.d0*vtess(phi-d3phi,flam+d3lam,h-d3h)-48.d0*vtess(phi-d3phi,flam+d3lam,h) &
        +12.d0*vtess(phi-d3phi,flam,h-2.d0*d3h)-48.d0*vtess(phi-d3phi,flam,h-d3h) &
        +36.d0*vtess(phi-d3phi,flam,h)-3.d0*vtess(phi,flam+2.d0*d3lam,h-2.d0*d3h) &
        +12.d0*vtess(phi,flam+2.d0*d3lam,h-d3h)-9.d0*vtess(phi,flam+2.d0*d3lam,h) &
        +12.d0*vtess(phi,flam+d3lam,h-2.d0*d3h)-48.d0*vtess(phi,flam+d3lam,h-d3h) &
        +36.d0*vtess(phi,flam+d3lam,h)-9.d0*vtess(phi,flam,h-2.d0*d3h) &
        +36.d0*vtess(phi,flam,h-d3h)-27.d0*v &
        )/(8.d0*d3phi*d3lam*d3h)
! case 7
elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
    gggPhiLamH=(-vtess(phi-2.d0*d3phi,flam-2.d0*d3lam,h+2.d0*d3h) &
        +4.d0*vtess(phi-2.d0*d3phi,flam-2.d0*d3lam,h+d3h)-3.d0*vtess(phi-2.d0*d3phi,flam-2.d0*d3lam,h) &
        +4.d0*vtess(phi-2.d0*d3phi,flam-d3lam,h+2.d0*d3h)-16.d0*vtess(phi-2.d0*d3phi,flam-d3lam,h+d3h) &
        +12.d0*vtess(phi-2.d0*d3phi,flam-d3lam,h)-3.d0*vtess(phi-2.d0*d3phi,flam,h+2.d0*d3h) &
        +12.d0*vtess(phi-2.d0*d3phi,flam,h+d3h)-9.d0*vtess(phi-2.d0*d3phi,flam,h) &
        +4.d0*vtess(phi-d3phi,flam-2.d0*d3lam,h+2.d0*d3h)-16.d0*vtess(phi-d3phi,flam-2.d0*d3lam,h+d3h) &
        +12.d0*vtess(phi-d3phi,flam-2.d0*d3lam,h)-16.d0*vtess(phi-d3phi,flam-d3lam,h+2.d0*d3h) &
        +64.d0*vtess(phi-d3phi,flam-d3lam,h+d3h)-48.d0*vtess(phi-d3phi,flam-d3lam,h) &
        +12.d0*vtess(phi-d3phi,flam,h+2.d0*d3h)-48.d0*vtess(phi-d3phi,flam,h+d3h) &
        +36.d0*vtess(phi-d3phi,flam,h)-3.d0*vtess(phi,flam-2.d0*d3lam,h+2.d0*d3h) &
        +12.d0*vtess(phi,flam-2.d0*d3lam,h+d3h)-9.d0*vtess(phi,flam-2.d0*d3lam,h) &
        +12.d0*vtess(phi,flam-d3lam,h+2.d0*d3h)-48.d0*vtess(phi,flam-d3lam,h+d3h) &
        +36.d0*vtess(phi,flam-d3lam,h)-9.d0*vtess(phi,flam,h+2.d0*d3h) &
        +36.d0*vtess(phi,flam,h+d3h)-27.d0*v &
        )/(8.d0*d3phi*d3lam*d3h)
! case 8
elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
    gggPhiLamH=-(-vtess(phi-2.d0*d3phi,flam-2.d0*d3lam,h-2.d0*d3h) &
        +4.d0*vtess(phi-2.d0*d3phi,flam-2.d0*d3lam,h-d3h)-3.d0*vtess(phi-2.d0*d3phi,flam-2.d0*d3lam,h) &
        +4.d0*vtess(phi-2.d0*d3phi,flam-d3lam,h-2.d0*d3h)-16.d0*vtess(phi-2.d0*d3phi,flam-d3lam,h-d3h) &
        +12.d0*vtess(phi-2.d0*d3phi,flam-d3lam,h)-3.d0*vtess(phi-2.d0*d3phi,flam,h-2.d0*d3h) &
        +12.d0*vtess(phi-2.d0*d3phi,flam,h-d3h)-9.d0*vtess(phi-2.d0*d3phi,flam,h) &
        +4.d0*vtess(phi-d3phi,flam-2.d0*d3lam,h-2.d0*d3h)-16.d0*vtess(phi-d3phi,flam-2.d0*d3lam,h-d3h) &
        +12.d0*vtess(phi-d3phi,flam-2.d0*d3lam,h)-16.d0*vtess(phi-d3phi,flam-d3lam,h-2.d0*d3h) &
        +64.d0*vtess(phi-d3phi,flam-d3lam,h-d3h)-48.d0*vtess(phi-d3phi,flam-d3lam,h) &
        +12.d0*vtess(phi-d3phi,flam,h-2.d0*d3h)-48.d0*vtess(phi-d3phi,flam,h-d3h) &
        +36.d0*vtess(phi-d3phi,flam,h)-3.d0*vtess(phi,flam-2.d0*d3lam,h-2.d0*d3h) &
        +12.d0*vtess(phi,flam-2.d0*d3lam,h-d3h)-9.d0*vtess(phi,flam-2.d0*d3lam,h) &
        +12.d0*vtess(phi,flam-d3lam,h-2.d0*d3h)-48.d0*vtess(phi,flam-d3lam,h-d3h) &
        +36.d0*vtess(phi,flam-d3lam,h)-9.d0*vtess(phi,flam,h-2.d0*d3h) &
        +36.d0*vtess(phi,flam,h-d3h)-27.d0*v &
        )/(8.d0*d3phi*d3lam*d3h)
! case 9
elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
    gggPhiLamH=(vtess(phi+d3phi,flam+2.d0*d3lam,h+2.d0*d3h) &
        -4.d0*vtess(phi+d3phi,flam+2.d0*d3lam,h+d3h)+3.d0*vtess(phi+d3phi,flam+2.d0*d3lam,h) &
        -4.d0*vtess(phi+d3phi,flam+d3lam,h+2.d0*d3h)+16.d0*vtess(phi+d3phi,flam+d3lam,h+d3h) &
        -12.d0*vtess(phi+d3phi,flam+d3lam,h)+3.d0*vtess(phi+d3phi,flam,h+2.d0*d3h) &
        -12.d0*vtess(phi+d3phi,flam,h+d3h)+9.d0*vtess(phi+d3phi,flam,h) &
        -vtess(phi-d3phi,flam+2.d0*d3lam,h+2.d0*d3h) &
        +4.d0*vtess(phi-d3phi,flam+2.d0*d3lam,h+d3h)-3.d0*vtess(phi-d3phi,flam+2.d0*d3lam,h) &
        +4.d0*vtess(phi-d3phi,flam+d3lam,h+2.d0*d3h)-16.d0*vtess(phi-d3phi,flam+d3lam,h+d3h) &
        +12.d0*vtess(phi-d3phi,flam+d3lam,h)-3.d0*vtess(phi-d3phi,flam,h+2.d0*d3h) &
        +12.d0*vtess(phi-d3phi,flam,h+d3h)-9.d0*vtess(phi-d3phi,flam,h) &
        )/(8.d0*d3phi*d3lam*d3h)
! case 10
elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
    gggPhiLamH=-(vtess(phi+d3phi,flam+2.d0*d3lam,h-2.d0*d3h) &
        -4.d0*vtess(phi+d3phi,flam+2.d0*d3lam,h-d3h)+3.d0*vtess(phi+d3phi,flam+2.d0*d3lam,h) &
        -4.d0*vtess(phi+d3phi,flam+d3lam,h-2.d0*d3h)+16.d0*vtess(phi+d3phi,flam+d3lam,h-d3h) &
        -12.d0*vtess(phi+d3phi,flam+d3lam,h)+3.d0*vtess(phi+d3phi,flam,h-2.d0*d3h) &
        -12.d0*vtess(phi+d3phi,flam,h-d3h)+9.d0*vtess(phi+d3phi,flam,h) &
        -vtess(phi-d3phi,flam+2.d0*d3lam,h-2.d0*d3h) &
        +4.d0*vtess(phi-d3phi,flam+2.d0*d3lam,h-d3h)-3.d0*vtess(phi-d3phi,flam+2.d0*d3lam,h) &
        +4.d0*vtess(phi-d3phi,flam+d3lam,h-2.d0*d3h)-16.d0*vtess(phi-d3phi,flam+d3lam,h-d3h) &
        +12.d0*vtess(phi-d3phi,flam+d3lam,h)-3.d0*vtess(phi-d3phi,flam,h-2.d0*d3h) &
        +12.d0*vtess(phi-d3phi,flam,h-d3h)-9.d0*vtess(phi-d3phi,flam,h) &
        )/(8.d0*d3phi*d3lam*d3h)
! case 11
elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
    gggPhiLamH=-(vtess(phi+d3phi,flam-2.d0*d3lam,h+2.d0*d3h) &
        -4.d0*vtess(phi+d3phi,flam-2.d0*d3lam,h+d3h)+3.d0*vtess(phi+d3phi,flam-2.d0*d3lam,h) &
        -4.d0*vtess(phi+d3phi,flam-d3lam,h+2.d0*d3h)+16.d0*vtess(phi+d3phi,flam-d3lam,h+d3h) &
        -12.d0*vtess(phi+d3phi,flam-d3lam,h)+3.d0*vtess(phi+d3phi,flam,h+2.d0*d3h) &
        -12.d0*vtess(phi+d3phi,flam,h+d3h)+9.d0*vtess(phi+d3phi,flam,h) &
        -vtess(phi-d3phi,flam-2.d0*d3lam,h+2.d0*d3h) &
        +4.d0*vtess(phi-d3phi,flam-2.d0*d3lam,h+d3h)-3.d0*vtess(phi-d3phi,flam-2.d0*d3lam,h) &
        +4.d0*vtess(phi-d3phi,flam-d3lam,h+2.d0*d3h)-16.d0*vtess(phi-d3phi,flam-d3lam,h+d3h) &
        +12.d0*vtess(phi-d3phi,flam-d3lam,h)-3.d0*vtess(phi-d3phi,flam,h+2.d0*d3h) &
        +12.d0*vtess(phi-d3phi,flam,h+d3h)-9.d0*vtess(phi-d3phi,flam,h) &
        )/(8.d0*d3phi*d3lam*d3h)
! case 12
elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
    gggPhiLamH=(vtess(phi+d3phi,flam-2.d0*d3lam,h-2.d0*d3h) &
        -4.d0*vtess(phi+d3phi,flam-2.d0*d3lam,h-d3h)+3.d0*vtess(phi+d3phi,flam-2.d0*d3lam,h) &
        -4.d0*vtess(phi+d3phi,flam-d3lam,h-2.d0*d3h)+16.d0*vtess(phi+d3phi,flam-d3lam,h-d3h) &
        -12.d0*vtess(phi+d3phi,flam-d3lam,h)+3.d0*vtess(phi+d3phi,flam,h-2.d0*d3h) &
        -12.d0*vtess(phi+d3phi,flam,h-d3h)+9.d0*vtess(phi+d3phi,flam,h) &
        -vtess(phi-d3phi,flam-2.d0*d3lam,h-2.d0*d3h) &
        +4.d0*vtess(phi-d3phi,flam-2.d0*d3lam,h-d3h)-3.d0*vtess(phi-d3phi,flam-2.d0*d3lam,h) &
        +4.d0*vtess(phi-d3phi,flam-d3lam,h-2.d0*d3h)-16.d0*vtess(phi-d3phi,flam-d3lam,h-d3h) &
        +12.d0*vtess(phi-d3phi,flam-d3lam,h)-3.d0*vtess(phi-d3phi,flam,h-2.d0*d3h) &
        +12.d0*vtess(phi-d3phi,flam,h-d3h)-9.d0*vtess(phi-d3phi,flam,h) &
        )/(8.d0*d3phi*d3lam*d3h)
! case 13
elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
    gggPhiLamH=(vtess(phi+2.d0*d3phi,flam+d3lam,h+2.d0*d3h) &
        -4.d0*vtess(phi+2.d0*d3phi,flam+d3lam,h+d3h)+3.d0*vtess(phi+2.0d0*d3phi,flam+d3lam,h) &
        -4.d0*vtess(phi+d3phi,flam+d3lam,h+2.d0*d3h)+16.d0*vtess(phi+d3phi,flam+d3lam,h+d3h) &
        -12.d0*vtess(phi+d3phi,flam+d3lam,h)+3.d0*vtess(phi,flam+d3lam,h+2.d0*d3h) &
        -12.d0*vtess(phi,flam+d3lam,h+d3h)+9.d0*vtess(phi,flam+d3lam,h) &
        -vtess(phi+2.d0*d3phi,flam-d3lam,h+2.d0*d3h) &
        +4.d0*vtess(phi+2.d0*d3phi,flam-d3lam,h+d3h)-3.d0*vtess(phi+2.0d0*d3phi,flam-d3lam,h) &
        +4.d0*vtess(phi+d3phi,flam-d3lam,h+2.d0*d3h)-16.d0*vtess(phi+d3phi,flam-d3lam,h+d3h) &
        +12.d0*vtess(phi+d3phi,flam-d3lam,h)-3.d0*vtess(phi,flam-d3lam,h+2.d0*d3h) &
        +12.d0*vtess(phi,flam-d3lam,h+d3h)-9.d0*vtess(phi,flam-d3lam,h) &
        )/(8.d0*d3phi*d3lam*d3h)
! case 14
elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
    gggPhiLamH=-(vtess(phi+2.d0*d3phi,flam+d3lam,h-2.d0*d3h) &
        -4.d0*vtess(phi+2.d0*d3phi,flam+d3lam,h-d3h)+3.d0*vtess(phi+2.0d0*d3phi,flam+d3lam,h) &
        -4.d0*vtess(phi+d3phi,flam+d3lam,h-2.d0*d3h)+16.d0*vtess(phi+d3phi,flam+d3lam,h-d3h) &
        -12.d0*vtess(phi+d3phi,flam+d3lam,h)+3.d0*vtess(phi,flam+d3lam,h-2.d0*d3h) &
        -12.d0*vtess(phi,flam+d3lam,h-d3h)+9.d0*vtess(phi,flam+d3lam,h) &
        -vtess(phi+2.d0*d3phi,flam-d3lam,h-2.d0*d3h) &
        +4.d0*vtess(phi+2.d0*d3phi,flam-d3lam,h-d3h)-3.d0*vtess(phi+2.0d0*d3phi,flam-d3lam,h) &
        +4.d0*vtess(phi+d3phi,flam-d3lam,h-2.d0*d3h)-16.d0*vtess(phi+d3phi,flam-d3lam,h-d3h) &
        +12.d0*vtess(phi+d3phi,flam-d3lam,h)-3.d0*vtess(phi,flam-d3lam,h-2.d0*d3h) &
        +12.d0*vtess(phi,flam-d3lam,h-d3h)-9.d0*vtess(phi,flam-d3lam,h) &
        )/(8.d0*d3phi*d3lam*d3h)
! case 15
elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
    gggPhiLamH=-(vtess(phi-2.d0*d3phi,flam+d3lam,h+2.d0*d3h) &
        -4.d0*vtess(phi-2.d0*d3phi,flam+d3lam,h+d3h)+3.d0*vtess(phi-2.0d0*d3phi,flam+d3lam,h) &
        -4.d0*vtess(phi-d3phi,flam+d3lam,h+2.d0*d3h)+16.d0*vtess(phi-d3phi,flam+d3lam,h+d3h) &
        -12.d0*vtess(phi-d3phi,flam+d3lam,h)+3.d0*vtess(phi,flam+d3lam,h+2.d0*d3h) &
        -12.d0*vtess(phi,flam+d3lam,h+d3h)+9.d0*vtess(phi,flam+d3lam,h) &
        -vtess(phi-2.d0*d3phi,flam-d3lam,h+2.d0*d3h) &
        +4.d0*vtess(phi-2.d0*d3phi,flam-d3lam,h+d3h)-3.d0*vtess(phi-2.0d0*d3phi,flam-d3lam,h) &
        +4.d0*vtess(phi-d3phi,flam-d3lam,h+2.d0*d3h)-16.d0*vtess(phi-d3phi,flam-d3lam,h+d3h) &
        +12.d0*vtess(phi-d3phi,flam-d3lam,h)-3.d0*vtess(phi,flam-d3lam,h+2.d0*d3h) &
        +12.d0*vtess(phi,flam-d3lam,h+d3h)-9.d0*vtess(phi,flam-d3lam,h) &
        )/(8.d0*d3phi*d3lam*d3h)
! case 16
elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
    gggPhiLamH=(vtess(phi-2.d0*d3phi,flam+d3lam,h-2.d0*d3h) &
        -4.d0*vtess(phi-2.d0*d3phi,flam+d3lam,h-d3h)+3.d0*vtess(phi-2.0d0*d3phi,flam+d3lam,h) &
        -4.d0*vtess(phi-d3phi,flam+d3lam,h-2.d0*d3h)+16.d0*vtess(phi-d3phi,flam+d3lam,h-d3h) &
        -12.d0*vtess(phi-d3phi,flam+d3lam,h)+3.d0*vtess(phi,flam+d3lam,h-2.d0*d3h) &
        -12.d0*vtess(phi,flam+d3lam,h-d3h)+9.d0*vtess(phi,flam+d3lam,h) &
        -vtess(phi-2.d0*d3phi,flam-d3lam,h-2.d0*d3h) &
        +4.d0*vtess(phi-2.d0*d3phi,flam-d3lam,h-d3h)-3.d0*vtess(phi-2.0d0*d3phi,flam-d3lam,h) &
        +4.d0*vtess(phi-d3phi,flam-d3lam,h-2.d0*d3h)-16.d0*vtess(phi-d3phi,flam-d3lam,h-d3h) &
        +12.d0*vtess(phi-d3phi,flam-d3lam,h)-3.d0*vtess(phi,flam-d3lam,h-2.d0*d3h) &
        +12.d0*vtess(phi,flam-d3lam,h-d3h)-9.d0*vtess(phi,flam-d3lam,h) &
        )/(8.d0*d3phi*d3lam*d3h)
! case 17
elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
    gggPhiLamH=(vtess(phi+2.d0*d3phi,flam+2.d0*d3lam,h+d3h) &
        -4.d0*vtess(phi+2.d0*d3phi,flam+d3lam,h+d3h)+3.d0*vtess(phi+2.d0*d3phi,flam,h+d3h) &
        -4.d0*vtess(phi+d3phi,flam+2.d0*d3lam,h+d3h)+16.d0*vtess(phi+d3phi,flam+d3lam,h+d3h) &
        -12.d0*vtess(phi+d3phi,flam,h+d3h)+3.d0*vtess(phi,flam+2.d0*d3lam,h+d3h) &
        -12.d0*vtess(phi,flam+d3lam,h+d3h)+9.d0*vtess(phi,flam,h+d3h) &
        -vtess(phi+2.d0*d3phi,flam+2.d0*d3lam,h-d3h) &
        +4.d0*vtess(phi+2.d0*d3phi,flam+d3lam,h-d3h)-3.d0*vtess(phi+2.d0*d3phi,flam,h-d3h) &
        +4.d0*vtess(phi+d3phi,flam+2.d0*d3lam,h-d3h)-16.d0*vtess(phi+d3phi,flam+d3lam,h-d3h) &
        +12.d0*vtess(phi+d3phi,flam,h-d3h)-3.d0*vtess(phi,flam+2.d0*d3lam,h-d3h) &
        +12.d0*vtess(phi,flam+d3lam,h-d3h)-9.d0*vtess(phi,flam,h-d3h) &
        )/(8.d0*d3phi*d3lam*d3h)
! case 18
elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
    gggPhiLamH=-(vtess(phi+2.d0*d3phi,flam-2.d0*d3lam,h+d3h) &
        -4.d0*vtess(phi+2.d0*d3phi,flam-d3lam,h+d3h)+3.d0*vtess(phi+2.d0*d3phi,flam,h+d3h) &
        -4.d0*vtess(phi+d3phi,flam-2.d0*d3lam,h+d3h)+16.d0*vtess(phi+d3phi,flam-d3lam,h+d3h) &
        -12.d0*vtess(phi+d3phi,flam,h+d3h)+3.d0*vtess(phi,flam-2.d0*d3lam,h+d3h) &
        -12.d0*vtess(phi,flam-d3lam,h+d3h)+9.d0*vtess(phi,flam,h+d3h) &
        -vtess(phi+2.d0*d3phi,flam-2.d0*d3lam,h-d3h) &
        +4.d0*vtess(phi+2.d0*d3phi,flam-d3lam,h-d3h)-3.d0*vtess(phi+2.d0*d3phi,flam,h-d3h) &
        +4.d0*vtess(phi+d3phi,flam-2.d0*d3lam,h-d3h)-16.d0*vtess(phi+d3phi,flam-d3lam,h-d3h) &
        +12.d0*vtess(phi+d3phi,flam,h-d3h)-3.d0*vtess(phi,flam-2.d0*d3lam,h-d3h) &
        +12.d0*vtess(phi,flam-d3lam,h-d3h)-9.d0*vtess(phi,flam,h-d3h) &
        )/(8.d0*d3phi*d3lam*d3h)
! case 19
elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
    gggPhiLamH=-(vtess(phi-2.d0*d3phi,flam+2.d0*d3lam,h+d3h) &
        -4.d0*vtess(phi-2.d0*d3phi,flam+d3lam,h+d3h)+3.d0*vtess(phi-2.d0*d3phi,flam,h+d3h) &
        -4.d0*vtess(phi-d3phi,flam+2.d0*d3lam,h+d3h)+16.d0*vtess(phi-d3phi,flam+d3lam,h+d3h) &
        -12.d0*vtess(phi-d3phi,flam,h+d3h)+3.d0*vtess(phi,flam+2.d0*d3lam,h+d3h) &
        -12.d0*vtess(phi,flam+d3lam,h+d3h)+9.d0*vtess(phi,flam,h+d3h) &
        -vtess(phi-2.d0*d3phi,flam+2.d0*d3lam,h-d3h) &
        +4.d0*vtess(phi-2.d0*d3phi,flam+d3lam,h-d3h)-3.d0*vtess(phi-2.d0*d3phi,flam,h-d3h) &
        +4.d0*vtess(phi-d3phi,flam+2.d0*d3lam,h-d3h)-16.d0*vtess(phi-d3phi,flam+d3lam,h-d3h) &
        +12.d0*vtess(phi-d3phi,flam,h-d3h)-3.d0*vtess(phi,flam+2.d0*d3lam,h-d3h) &
        +12.d0*vtess(phi,flam+d3lam,h-d3h)-9.d0*vtess(phi,flam,h-d3h) &
        )/(8.d0*d3phi*d3lam*d3h)
! case 20
elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
    gggPhiLamH=(vtess(phi-2.d0*d3phi,flam-2.d0*d3lam,h+d3h) &
        -4.d0*vtess(phi-2.d0*d3phi,flam-d3lam,h+d3h)+3.d0*vtess(phi-2.d0*d3phi,flam,h+d3h) &
        -4.d0*vtess(phi-d3phi,flam-2.d0*d3lam,h+d3h)+16.d0*vtess(phi-d3phi,flam-d3lam,h+d3h) &
        -12.d0*vtess(phi-d3phi,flam,h+d3h)+3.d0*vtess(phi,flam-2.d0*d3lam,h+d3h) &
        -12.d0*vtess(phi,flam-d3lam,h+d3h)+9.d0*vtess(phi,flam,h+d3h) &
        -vtess(phi-2.d0*d3phi,flam-2.d0*d3lam,h-d3h) &
        +4.d0*vtess(phi-2.d0*d3phi,flam-d3lam,h-d3h)-3.d0*vtess(phi-2.d0*d3phi,flam,h-d3h) &
        +4.d0*vtess(phi-d3phi,flam-2.d0*d3lam,h-d3h)-16.d0*vtess(phi-d3phi,flam-d3lam,h-d3h) &
        +12.d0*vtess(phi-d3phi,flam,h-d3h)-3.d0*vtess(phi,flam-2.d0*d3lam,h-d3h) &
        +12.d0*vtess(phi,flam-d3lam,h-d3h)-9.d0*vtess(phi,flam,h-d3h) &
        )/(8.d0*d3phi*d3lam*d3h)
! case 21
elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
    gggPhiLamH=(-vtess(phi+d3phi,flam+d3lam,h+2.d0*d3h) &
        +4.d0*vtess(phi+d3phi,flam+d3lam,h+d3h)-3.d0*vtess(phi+d3phi,flam+d3lam,h) &
        +vtess(phi+d3phi,flam-d3lam,h+2.d0*d3h)-4.d0*vtess(phi+d3phi,flam-d3lam,h+d3h) &
        +3.d0*vtess(phi+d3phi,flam-d3lam,h)&
        +vtess(phi-d3phi,flam+d3lam,h+2.d0*d3h) &
        -4.d0*vtess(phi-d3phi,flam+d3lam,h+d3h)+3.d0*vtess(phi-d3phi,flam+d3lam,h) &
        -vtess(phi-d3phi,flam-d3lam,h+2.d0*d3h)+4.d0*vtess(phi-d3phi,flam-d3lam,h+d3h) &
        -3.d0*vtess(phi-d3phi,flam-d3lam,h) &
        )/(8.d0*d3phi*d3lam*d3h)
! case 22
elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
    gggPhiLamH=-(-vtess(phi+d3phi,flam+d3lam,h-2.d0*d3h) &
        +4.d0*vtess(phi+d3phi,flam+d3lam,h-d3h)-3.d0*vtess(phi+d3phi,flam+d3lam,h) &
        +vtess(phi+d3phi,flam-d3lam,h-2.d0*d3h)-4.d0*vtess(phi+d3phi,flam-d3lam,h-d3h) &
        +3.d0*vtess(phi+d3phi,flam-d3lam,h) &
        +vtess(phi-d3phi,flam+d3lam,h-2.d0*d3h) &
        -4.d0*vtess(phi-d3phi,flam+d3lam,h-d3h)+3.d0*vtess(phi-d3phi,flam+d3lam,h) &
        -vtess(phi-d3phi,flam-d3lam,h-2.d0*d3h)+4.d0*vtess(phi-d3phi,flam-d3lam,h-d3h) &
        -3.d0*vtess(phi-d3phi,flam-d3lam,h) &
        )/(8.d0*d3phi*d3lam*d3h)
! case 23
elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
    gggPhiLamH=(-vtess(phi+d3phi,flam+2.d0*d3lam,h+d3h) &
        +4.d0*vtess(phi+d3phi,flam+d3lam,h+d3h)-3.d0*vtess(phi+d3phi,flam,h+d3h) &
        +vtess(phi+d3phi,flam+2.d0*d3lam,h-d3h)-4.d0*vtess(phi+d3phi,flam+d3lam,h-d3h) &
        +3.d0*vtess(phi+d3phi,flam,h-d3h) &
        +vtess(phi-d3phi,flam+2.d0*d3lam,h+d3h) &
        -4.d0*vtess(phi-d3phi,flam+d3lam,h+d3h)+3.d0*vtess(phi-d3phi,flam,h+d3h) &
        -vtess(phi-d3phi,flam+2.d0*d3lam,h-d3h)+4.d0*vtess(phi-d3phi,flam+d3lam,h-d3h) &
        -3.d0*vtess(phi-d3phi,flam,h-d3h) &
        )/(8.d0*d3phi*d3lam*d3h)
! case 24
elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
    gggPhiLamH=-(-vtess(phi+d3phi,flam-2.d0*d3lam,h+d3h) &
        +4.d0*vtess(phi+d3phi,flam-d3lam,h+d3h)-3.d0*vtess(phi+d3phi,flam,h+d3h) &
        +vtess(phi+d3phi,flam-2.d0*d3lam,h-d3h)-4.d0*vtess(phi+d3phi,flam-d3lam,h-d3h) &
        +3.d0*vtess(phi+d3phi,flam,h-d3h) &
        +vtess(phi-d3phi,flam-2.d0*d3lam,h+d3h) &
        -4.d0*vtess(phi-d3phi,flam-d3lam,h+d3h)+3.d0*vtess(phi-d3phi,flam,h+d3h) &
        -vtess(phi-d3phi,flam-2.d0*d3lam,h-d3h)+4.d0*vtess(phi-d3phi,flam-d3lam,h-d3h) &
        -3.d0*vtess(phi-d3phi,flam,h-d3h) &
        )/(8.d0*d3phi*d3lam*d3h)
! case 25
elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
    gggPhiLamH=(-vtess(phi+2.d0*d3phi,flam+d3lam,h+d3h) &
        +4.d0*vtess(phi+d3phi,flam+d3lam,h+d3h)-3.d0*vtess(phi,flam+d3lam,h+d3h) &
        +vtess(phi+2.d0*d3phi,flam+d3lam,h-d3h)-4.d0*vtess(phi+d3phi,flam+d3lam,h-d3h) &
        +3.d0*vtess(phi,flam+d3lam,h-d3h) &
        +vtess(phi+2.d0*d3phi,flam-d3lam,h+d3h) &
        -4.d0*vtess(phi+d3phi,flam-d3lam,h+d3h)+3.d0*vtess(phi,flam-d3lam,h+d3h) &
        -vtess(phi+2.d0*d3phi,flam-d3lam,h-d3h)+4.d0*vtess(phi+d3phi,flam-d3lam,h-d3h) &
        -3.d0*vtess(phi,flam-d3lam,h-d3h) &
        )/(8.d0*d3phi*d3lam*d3h)
! case 26
elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
    gggPhiLamH=-(-vtess(phi-2.d0*d3phi,flam+d3lam,h+d3h) &
        +4.d0*vtess(phi-d3phi,flam+d3lam,h+d3h)-3.d0*vtess(phi,flam+d3lam,h+d3h) &
        +vtess(phi-2.d0*d3phi,flam+d3lam,h-d3h)-4.d0*vtess(phi-d3phi,flam+d3lam,h-d3h) &
        +3.d0*vtess(phi,flam+d3lam,h-d3h) &
        +vtess(phi-2.d0*d3phi,flam-d3lam,h+d3h) &
        -4.d0*vtess(phi-d3phi,flam-d3lam,h+d3h)+3.d0*vtess(phi,flam-d3lam,h+d3h) &
        -vtess(phi-2.d0*d3phi,flam-d3lam,h-d3h)+4.d0*vtess(phi-d3phi,flam-d3lam,h-d3h) &
        -3.d0*vtess(phi,flam-d3lam,h-d3h) &
        )/(8.d0*d3phi*d3lam*d3h)
! case 27
else
    gggPhiLamH=(vtess(phi+d3phi,flam+d3lam,h+d3h) &
        -vtess(phi+d3phi,flam+d3lam,h-d3h)-vtess(phi+d3phi,flam-d3lam,h+d3h) &
        +vtess(phi+d3phi,flam-d3lam,h-d3h)-vtess(phi-d3phi,flam+d3lam,h+d3h) &
        +vtess(phi-d3phi,flam+d3lam,h-d3h)+vtess(phi-d3phi,flam-d3lam,h+d3h) &
        -vtess(phi-d3phi,flam-d3lam,h-d3h))/(8.d0*d3phi*d3lam*d3h)
endif
gggPhiLamH=(gggPhiLamH/c-2.d0*r*ggPhiLam+t*(gLam+r*ggLamH))/(r*r)
!
!   Compute gggLamLamPhi
!
if(h.ge.HB.and.h.le.HT) then
! case 1
    if(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE))) then
        gggLamLamPhi=(vtess(phi+2.d0*d3phi,flam+3.d0*d3lam,h)-4.d0*vtess(phi+d3phi,flam+3.d0*d3lam,h) &
            +3.d0*vtess(phi,flam+3.d0*d3lam,h)-4.d0*vtess(phi+2.d0*d3phi,flam+2.d0*d3lam,h) &
            +16.d0*vtess(phi+d3phi,flam+2.d0*d3lam,h)-12.d0*vtess(phi,flam+2.d0*d3lam,h) &
            +5.d0*vtess(phi+2.d0*d3phi,flam+d3lam,h)-20.d0*vtess(phi+d3phi,flam+d3lam,h) &
            +15.d0*vtess(phi,flam+d3lam,h)-2.d0*vtess(phi+2.d0*d3phi,flam,h) &
            +8.d0*vtess(phi+d3phi,flam,h)-6.d0*v)/(2.d0*d3phi*d3lam*d3lam)
! case 2
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE))) then
        gggLamLamPhi=-(vtess(phi-2.d0*d3phi,flam+3.d0*d3lam,h)-4.d0*vtess(phi-d3phi,flam+3.d0*d3lam,h) &
            +3.d0*vtess(phi,flam+3.d0*d3lam,h)-4.d0*vtess(phi-2.d0*d3phi,flam+2.d0*d3lam,h) &
            +16.d0*vtess(phi-d3phi,flam+2.d0*d3lam,h)-12.d0*vtess(phi,flam+2.d0*d3lam,h) &
            +5.d0*vtess(phi-2.d0*d3phi,flam+d3lam,h)-20.d0*vtess(phi-d3phi,flam+d3lam,h) &
            +15.d0*vtess(phi,flam+d3lam,h)-2.d0*vtess(phi-2.d0*d3phi,flam,h) &
            +8.d0*vtess(phi-d3phi,flam,h)-6.d0*v)/(2.d0*d3phi*d3lam*d3lam)
! case 3
    elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE))) then
        gggLamLamPhi=(vtess(phi+2.d0*d3phi,flam-3.d0*d3lam,h)-4.d0*vtess(phi+d3phi,flam-3.d0*d3lam,h) &
            +3.d0*vtess(phi,flam-3.d0*d3lam,h)-4.d0*vtess(phi+2.d0*d3phi,flam-2.d0*d3lam,h) &
            +16.d0*vtess(phi+d3phi,flam-2.d0*d3lam,h)-12.d0*vtess(phi,flam-2.d0*d3lam,h) &
            +5.d0*vtess(phi+2.d0*d3phi,flam-d3lam,h)-20.d0*vtess(phi+d3phi,flam-d3lam,h) &
            +15.d0*vtess(phi,flam-d3lam,h)-2.d0*vtess(phi+2.d0*d3phi,flam,h) &
            +8.d0*vtess(phi+d3phi,flam,h)-6.d0*v)/(2.d0*d3phi*d3lam*d3lam)
! case 4
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE))) then
        gggLamLamPhi=-(vtess(phi-2.d0*d3phi,flam-3.d0*d3lam,h)-4.d0*vtess(phi-d3phi,flam-3.d0*d3lam,h) &
            +3.d0*vtess(phi,flam-3.d0*d3lam,h)-4.d0*vtess(phi-2.d0*d3phi,flam-2.d0*d3lam,h) &
            +16.d0*vtess(phi-d3phi,flam-2.d0*d3lam,h)-12.d0*vtess(phi,flam-2.d0*d3lam,h) &
            +5.d0*vtess(phi-2.d0*d3phi,flam-d3lam,h)-20.d0*vtess(phi-d3phi,flam-d3lam,h) &
            +15.d0*vtess(phi,flam-d3lam,h)-2.d0*vtess(phi-2.d0*d3phi,flam,h) &
            +8.d0*vtess(phi-d3phi,flam,h)-6.d0*v)/(2.d0*d3phi*d3lam*d3lam)
! case 5
    elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE))) then
        gggLamLamPhi=(-vtess(phi+d3phi,flam+3.d0*d3lam,h)+vtess(phi-d3phi,flam+3.d0*d3lam,h) &
            +4.d0*vtess(phi+d3phi,flam+2.d0*d3lam,h)-4.d0*vtess(phi-d3phi,flam+2.d0*d3lam,h) &
            -5.d0*vtess(phi+d3phi,flam+d3lam,h)+5.d0*vtess(phi-d3phi,flam+d3lam,h) &
            +2.d0*vtess(phi+d3phi,flam,h)-2.d0*vtess(phi-d3phi,flam,h))/(2.d0*d3phi*d3lam*d3lam)
! case 6
    elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE))) then
        gggLamLamPhi=(-vtess(phi+d3phi,flam-3.d0*d3lam,h)+vtess(phi-d3phi,flam-3.d0*d3lam,h) &
            +4.d0*vtess(phi+d3phi,flam-2.d0*d3lam,h)-4.d0*vtess(phi-d3phi,flam-2.d0*d3lam,h) &
            -5.d0*vtess(phi+d3phi,flam-d3lam,h)+5.d0*vtess(phi-d3phi,flam-d3lam,h) &
            +2.d0*vtess(phi+d3phi,flam,h)-2.d0*vtess(phi-d3phi,flam,h))/(2.d0*d3phi*d3lam*d3lam)
! case 7
    elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE))) then
        gggLamLamPhi=(-vtess(phi+2.d0*d3phi,flam+d3lam,h)+4.d0*vtess(phi+d3phi,flam+d3lam,h) &
            -3.d0*vtess(phi,flam+d3lam,h)+2.d0*vtess(phi+2.d0*d3phi,flam,h) &
            -8.d0*vtess(phi+d3phi,flam,h)+6.d0*v &
            -vtess(phi+2.d0*d3phi,flam-d3lam,h)+4.d0*vtess(phi+d3phi,flam-d3lam,h) &
            -3.d0*vtess(phi,flam-d3lam,h))/(2.d0*d3phi*d3lam*d3lam)
! case 8
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE))) then
        gggLamLamPhi=-(-vtess(phi-2.d0*d3phi,flam+d3lam,h)+4.d0*vtess(phi-d3phi,flam+d3lam,h) &
            -3.d0*vtess(phi,flam+d3lam,h)+2.d0*vtess(phi-2.d0*d3phi,flam,h) &
            -8.d0*vtess(phi-d3phi,flam,h)+6.d0*v &
            -vtess(phi-2.d0*d3phi,flam-d3lam,h)+4.d0*vtess(phi-d3phi,flam-d3lam,h) &
            -3.d0*vtess(phi,flam-d3lam,h))/(2.d0*d3phi*d3lam*d3lam)
! case 9
    else
        gggLamLamPhi=(vtess(phi+d3phi,flam+d3lam,h)-vtess(phi-d3phi,flam+d3lam,h) &
            -2.d0*vtess(phi+d3phi,flam,h)+2.d0*vtess(phi-d3phi,flam,h) &
            +vtess(phi+d3phi,flam-d3lam,h)-vtess(phi-d3phi,flam-d3lam,h))/(2.d0*d3phi*d3lam*d3lam)
    endif
else
    gggLamLamPhi=(vtess(phi+d3phi,flam+d3lam,h)-vtess(phi-d3phi,flam+d3lam,h) &
            -2.d0*vtess(phi+d3phi,flam,h)+2.d0*vtess(phi-d3phi,flam,h) &
            +vtess(phi+d3phi,flam-d3lam,h)-vtess(phi-d3phi,flam-d3lam,h))/(2.d0*d3phi*d3lam*d3lam)
endif
gggLamLamPhi=(gggLamLamPhi-r*gPhi)/(r*r*r*c*c) &
+ (gPhi+r*ggPhiH-t*(gH-2.d0*r*ggLamLam+r*ggPhiPhi)+2.d0*t*t*gPhi)/(r*r)
!
!   Compute gggLamLamLam
!
if((phi.ge.PhiS.and.phi.le.PhiN).and.(h.ge.HB.and.h.le.HT)) then
    if((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) then
        gggLamLamLam=(-3.d0*vtess(phi,flam+4.d0*d3lam,h)+14.d0*vtess(phi,flam+3.d0*d3lam,h) &
            -24.d0*vtess(phi,flam+2.d0*d3lam,h)+18.d0*vtess(phi,flam+d3lam,h) &
            -5.0d0*v)/(2.d0*d3lam*d3lam*d3lam)
    elseif((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) then
        gggLamLamLam=(3.d0*vtess(phi,flam-4.d0*d3lam,h)-14.d0*vtess(phi,flam-3.d0*d3lam,h) &
            +24.d0*vtess(phi,flam-2.d0*d3lam,h)-18.d0*vtess(phi,flam-d3lam,h) &
            +5.0d0*v)/(2.d0*d3lam*d3lam*d3lam)
    else
        gggLamLamLam=(vtess(phi,flam+2.d0*d3lam,h)-2.d0*vtess(phi,flam+d3lam,h) &
            +2.d0*vtess(phi,flam-d3lam,h)-vtess(phi,flam-2.d0*d3lam,h) &
            )/(2.d0*d3lam*d3lam*d3lam)
    endif
else
    gggLamLamLam=(vtess(phi,flam+2.d0*d3lam,h)-2.d0*vtess(phi,flam+d3lam,h) &
            +2.d0*vtess(phi,flam-d3lam,h)-vtess(phi,flam-2.d0*d3lam,h) &
            )/(2.d0*d3lam*d3lam*d3lam)
endif
gggLamLamLam=(r*c*gLam+gggLamLamLam)/(r*r*r*c*c*c)+3.d0*(ggLamH-t*ggPhiLam)/r
!
!   Compute gggLamLamH
!
if(phi.ge.PhiS.and.phi.le.PhiN) then
! case 1
    if(((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggLamLamH=(vtess(phi,flam+3.d0*d3lam,h+2.d0*d3h)-4.d0*vtess(phi,flam+3.d0*d3lam,h+d3h) &
            +3.d0*vtess(phi,flam+3.d0*d3lam,h)-4.d0*vtess(phi,flam+2.d0*d3lam,h+2.d0*d3h) &
            +16.d0*vtess(phi,flam+2.d0*d3lam,h+d3h)-12.d0*vtess(phi,flam+2.d0*d3lam,h) &
            +5.d0*vtess(phi,flam+d3lam,h+2.d0*d3h)-20.d0*vtess(phi,flam+d3lam,h+d3h) &
            +15.d0*vtess(phi,flam+d3lam,h)-2.d0*vtess(phi,flam,h+2.d0*d3h) &
            +8.d0*vtess(phi,flam,h+d3h)-6.d0*v &
            )/(2.d0*d3lam*d3lam*d3h)
! case 2
    elseif(((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggLamLamH=-(vtess(phi,flam+3.d0*d3lam,h-2.d0*d3h)-4.d0*vtess(phi,flam+3.d0*d3lam,h-d3h) &
            +3.d0*vtess(phi,flam+3.d0*d3lam,h)-4.d0*vtess(phi,flam+2.d0*d3lam,h-2.d0*d3h) &
            +16.d0*vtess(phi,flam+2.d0*d3lam,h-d3h)-12.d0*vtess(phi,flam+2.d0*d3lam,h) &
            +5.d0*vtess(phi,flam+d3lam,h-2.d0*d3h)-20.d0*vtess(phi,flam+d3lam,h-d3h) &
            +15.d0*vtess(phi,flam+d3lam,h)-2.d0*vtess(phi,flam,h-2.d0*d3h) &
            +8.d0*vtess(phi,flam,h-d3h)-6.d0*v &
            )/(2.d0*d3lam*d3lam*d3h)
! case 3
    elseif(((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggLamLamH=(vtess(phi,flam-3.d0*d3lam,h+2.d0*d3h)-4.d0*vtess(phi,flam-3.d0*d3lam,h+d3h) &
            +3.d0*vtess(phi,flam-3.d0*d3lam,h)-4.d0*vtess(phi,flam-2.d0*d3lam,h+2.d0*d3h) &
            +16.d0*vtess(phi,flam-2.d0*d3lam,h+d3h)-12.d0*vtess(phi,flam-2.d0*d3lam,h) &
            +5.d0*vtess(phi,flam-d3lam,h+2.d0*d3h)-20.d0*vtess(phi,flam-d3lam,h+d3h) &
            +15.d0*vtess(phi,flam-d3lam,h)-2.d0*vtess(phi,flam,h+2.d0*d3h) &
            +8.d0*vtess(phi,flam,h+d3h)-6.d0*v &
            )/(2.d0*d3lam*d3lam*d3h)
! case 4
    elseif(((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggLamLamH=-(vtess(phi,flam-3.d0*d3lam,h-2.d0*d3h)-4.d0*vtess(phi,flam-3.d0*d3lam,h-d3h) &
            +3.d0*vtess(phi,flam-3.d0*d3lam,h)-4.d0*vtess(phi,flam-2.d0*d3lam,h-2.d0*d3h) &
            +16.d0*vtess(phi,flam-2.d0*d3lam,h-d3h)-12.d0*vtess(phi,flam-2.d0*d3lam,h) &
            +5.d0*vtess(phi,flam-d3lam,h-2.d0*d3h)-20.d0*vtess(phi,flam-d3lam,h-d3h) &
            +15.d0*vtess(phi,flam-d3lam,h)-2.d0*vtess(phi,flam,h-2.d0*d3h) &
            +8.d0*vtess(phi,flam,h-d3h)-6.d0*v &
            )/(2.d0*d3lam*d3lam*d3h)
! case 5
    elseif(((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
        gggLamLamH=(-vtess(phi,flam+3.d0*d3lam,h+d3h)+vtess(phi,flam+3.d0*d3lam,h-d3h) &
            +4.d0*vtess(phi,flam+2.d0*d3lam,h+d3h)-4.d0*vtess(phi,flam+2.d0*d3lam,h-d3h) &
            -5.d0*vtess(phi,flam+d3lam,h+d3h)+5.d0*vtess(phi,flam+d3lam,h-d3h) &
            +2.d0*vtess(phi,flam,h+d3h)-2.d0*vtess(phi,flam,h-d3h) &
            )/(2.d0*d3lam*d3lam*d3h)
! case 6
    elseif(((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
        gggLamLamH=(-vtess(phi,flam-3.d0*d3lam,h+d3h)+vtess(phi,flam-3.d0*d3lam,h-d3h) &
            +4.d0*vtess(phi,flam-2.d0*d3lam,h+d3h)-4.d0*vtess(phi,flam-2.d0*d3lam,h-d3h) &
            -5.d0*vtess(phi,flam-d3lam,h+d3h)+5.d0*vtess(phi,flam-d3lam,h-d3h) &
            +2.d0*vtess(phi,flam,h+d3h)-2.d0*vtess(phi,flam,h-d3h) &
            )/(2.d0*d3lam*d3lam*d3h)
! case 7
    elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggLamLamH=(-vtess(phi,flam+d3lam,h+2.d0*d3h)+4.d0*vtess(phi,flam+d3lam,h+d3h) &
            -3.d0*vtess(phi,flam+d3lam,h)+2.d0*vtess(phi,flam,h+2.d0*d3h) &
            -8.d0*vtess(phi,flam,h+d3h)+6.d0*v -vtess(phi,flam-d3lam,h+2.d0*d3h) &
            +4.d0*vtess(phi,flam-d3lam,h+d3h)-3.d0*vtess(phi,flam-d3lam,h) &
            )/(2.d0*d3lam*d3lam*d3h)
! case 8
    elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggLamLamH=-(-vtess(phi,flam+d3lam,h-2.d0*d3h)+4.d0*vtess(phi,flam+d3lam,h-d3h) &
            -3.d0*vtess(phi,flam+d3lam,h)+2.d0*vtess(phi,flam,h-2.d0*d3h) &
            -8.d0*vtess(phi,flam,h-d3h)+6.d0*v -vtess(phi,flam-d3lam,h-2.d0*d3h) &
            +4.d0*vtess(phi,flam-d3lam,h-d3h)-3.d0*vtess(phi,flam-d3lam,h) &
            )/(2.d0*d3lam*d3lam*d3h)
! case 9
    else
        gggLamLamH=(vtess(phi,flam+d3lam,h+d3h)-vtess(phi,flam+d3lam,h-d3h) &
            -2.d0*vtess(phi,flam,h+d3h)+2.d0*vtess(phi,flam,h-d3h) &
            +vtess(phi,flam-d3lam,h+d3h)-vtess(phi,flam-d3lam,h-d3h) &
            )/(2.d0*d3lam*d3lam*d3h)
    endif
else
    gggLamLamH=(vtess(phi,flam+d3lam,h+d3h)-vtess(phi,flam+d3lam,h-d3h) &
            -2.d0*vtess(phi,flam,h+d3h)+2.d0*vtess(phi,flam,h-d3h) &
            +vtess(phi,flam-d3lam,h+d3h)-vtess(phi,flam-d3lam,h-d3h) &
            )/(2.d0*d3lam*d3lam*d3h)
endif
gggLamLamH=(gH+r*ggHH-2.d0*r*ggLamLam+gggLamLamH/(c*c)-t*(gPhi+r*ggPhiH))/(r*r)
!
!   Compute gggHHPhi
!
if(flam.ge.FlamW.and.flam.le.FlamE) then
! case 1
    if(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggHHPhi=(vtess(phi+2.d0*d3phi,flam,h+3.d0*d3h)-4.d0*vtess(phi+d3phi,flam,h+3.d0*d3h) &
            +3.d0*vtess(phi,flam,h+3.d0*d3h)-4.d0*vtess(phi+2.d0*d3phi,flam,h+2.d0*d3h) &
            +16.d0*vtess(phi+d3phi,flam,h+2.d0*d3h)-12.d0*vtess(phi,flam,h+2.d0*d3h) &
            +5.d0*vtess(phi+2.d0*d3phi,flam,h+d3h)-20.d0*vtess(phi+d3phi,flam,h+d3h) &
            +15.d0*vtess(phi,flam,h+d3h)-2.d0*vtess(phi+2.d0*d3phi,flam,h) &
            +8.d0*vtess(phi+d3phi,flam,h)-6.d0*v)/(2.d0*d3phi*d3h*d3h)
! case 2
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggHHPhi=-(vtess(phi-2.d0*d3phi,flam,h+3.d0*d3h)-4.d0*vtess(phi-d3phi,flam,h+3.d0*d3h) &
            +3.d0*vtess(phi,flam,h+3.d0*d3h)-4.d0*vtess(phi-2.d0*d3phi,flam,h+2.d0*d3h) &
            +16.d0*vtess(phi-d3phi,flam,h+2.d0*d3h)-12.d0*vtess(phi,flam,h+2.d0*d3h) &
            +5.d0*vtess(phi-2.d0*d3phi,flam,h+d3h)-20.d0*vtess(phi-d3phi,flam,h+d3h) &
            +15.d0*vtess(phi,flam,h+d3h)-2.d0*vtess(phi-2.d0*d3phi,flam,h) &
            +8.d0*vtess(phi-d3phi,flam,h)-6.d0*v)/(2.d0*d3phi*d3h*d3h)
! case 3
    elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggHHPhi=(vtess(phi+2.d0*d3phi,flam,h-3.d0*d3h)-4.d0*vtess(phi+d3phi,flam,h-3.d0*d3h) &
            +3.d0*vtess(phi,flam,h-3.d0*d3h)-4.d0*vtess(phi+2.d0*d3phi,flam,h-2.d0*d3h) &
            +16.d0*vtess(phi+d3phi,flam,h-2.d0*d3h)-12.d0*vtess(phi,flam,h-2.d0*d3h) &
            +5.d0*vtess(phi+2.d0*d3phi,flam,h-d3h)-20.d0*vtess(phi+d3phi,flam,h-d3h) &
            +15.d0*vtess(phi,flam,h-d3h)-2.d0*vtess(phi+2.d0*d3phi,flam,h) &
            +8.d0*vtess(phi+d3phi,flam,h)-6.d0*v)/(2.d0*d3phi*d3h*d3h)
! case 4
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggHHPhi=-(vtess(phi-2.d0*d3phi,flam,h-3.d0*d3h)-4.d0*vtess(phi-d3phi,flam,h-3.d0*d3h) &
            +3.d0*vtess(phi,flam,h-3.d0*d3h)-4.d0*vtess(phi-2.d0*d3phi,flam,h-2.d0*d3h) &
            +16.d0*vtess(phi-d3phi,flam,h-2.d0*d3h)-12.d0*vtess(phi,flam,h-2.d0*d3h) &
            +5.d0*vtess(phi-2.d0*d3phi,flam,h-d3h)-20.d0*vtess(phi-d3phi,flam,h-d3h) &
            +15.d0*vtess(phi,flam,h-d3h)-2.d0*vtess(phi-2.d0*d3phi,flam,h) &
            +8.d0*vtess(phi-d3phi,flam,h)-6.d0*v)/(2.d0*d3phi*d3h*d3h)
! case 5
    elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggHHPhi=(-vtess(phi+d3phi,flam,h+3.d0*d3h)+vtess(phi-d3phi,flam,h+3.d0*d3h) &
            +4.d0*vtess(phi+d3phi,flam,h+2.d0*d3h)-4.d0*vtess(phi-d3phi,flam,h+2.d0*d3h) &
            -5.d0*vtess(phi+d3phi,flam,h+d3h)+5.d0*vtess(phi-d3phi,flam,h+d3h) &
            +2.d0*vtess(phi+d3phi,flam,h)-2.d0*vtess(phi-d3phi,flam,h))/(2.d0*d3phi*d3h*d3h)
! case 6
    elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggHHPhi=(-vtess(phi+d3phi,flam,h-3.d0*d3h)+vtess(phi-d3phi,flam,h-3.d0*d3h) &
            +4.d0*vtess(phi+d3phi,flam,h-2.d0*d3h)-4.d0*vtess(phi-d3phi,flam,h-2.d0*d3h) &
            -5.d0*vtess(phi+d3phi,flam,h-d3h)+5.d0*vtess(phi-d3phi,flam,h-d3h) &
            +2.d0*vtess(phi+d3phi,flam,h)-2.d0*vtess(phi-d3phi,flam,h))/(2.d0*d3phi*d3h*d3h)
! case 7
    elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
        gggHHPhi=(-vtess(phi+2.d0*d3phi,flam,h+d3h)+4.d0*vtess(phi+d3phi,flam,h+d3h) &
            -3.d0*vtess(phi,flam,h+d3h)+2.d0*vtess(phi+2.d0*d3phi,flam,h) &
            -8.d0*vtess(phi+d3phi,flam,h)+6.d0*v &
            -vtess(phi+2.d0*d3phi,flam,h-d3h)+4.d0*vtess(phi+d3phi,flam,h-d3h) &
            -3.d0*vtess(phi,flam,h-d3h))/(2.d0*d3phi*d3h*d3h)
! case 8
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
        gggHHPhi=-(-vtess(phi-2.d0*d3phi,flam,h+d3h)+4.d0*vtess(phi-d3phi,flam,h+d3h) &
            -3.d0*vtess(phi,flam,h+d3h)+2.d0*vtess(phi-2.d0*d3phi,flam,h) &
            -8.d0*vtess(phi-d3phi,flam,h)+6.d0*v &
            -vtess(phi-2.d0*d3phi,flam,h-d3h)+4.d0*vtess(phi-d3phi,flam,h-d3h) &
            -3.d0*vtess(phi,flam,h-d3h))/(2.d0*d3phi*d3h*d3h)
! case 9
    else
        gggHHPhi=(vtess(phi+d3phi,flam,h+d3h)-vtess(phi-d3phi,flam,h+d3h) &
            -2.d0*vtess(phi+d3phi,flam,h)+2.d0*vtess(phi-d3phi,flam,h) &
            +vtess(phi+d3phi,flam,h-d3h)-vtess(phi-d3phi,flam,h-d3h))/(2.d0*d3phi*d3h*d3h)
    endif
else
    gggHHPhi=(vtess(phi+d3phi,flam,h+d3h)-vtess(phi-d3phi,flam,h+d3h) &
            -2.d0*vtess(phi+d3phi,flam,h)+2.d0*vtess(phi-d3phi,flam,h) &
            +vtess(phi+d3phi,flam,h-d3h)-vtess(phi-d3phi,flam,h-d3h))/(2.d0*d3phi*d3h*d3h)
endif
gggHHPhi=(gggHHPhi-2.d0*ggPhiH)/r
!
!   Compute gggHHLam
!
if(phi.ge.PhiS.and.phi.le.PhiN) then
! case 1
    if(((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
        .and.((h.gt.HB.and.h-d3h.lt.HB) &
        .or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggHHLam=(vtess(phi,flam+2.d0*d3lam,h+3.d0*d3h)-4.d0*vtess(phi,flam+d3lam,h+3.d0*d3h) &
            +3.d0*vtess(phi,flam,h+3.d0*d3h)-4.d0*vtess(phi,flam+2.d0*d3lam,h+2.d0*d3h) &
            +16.d0*vtess(phi,flam+d3lam,h+2.d0*d3h)-12.d0*vtess(phi,flam,h+2.d0*d3h) &
            +5.d0*vtess(phi,flam+2.d0*d3lam,h+d3h)-20.d0*vtess(phi,flam+d3lam,h+d3h) &
            +15.d0*vtess(phi,flam,h+d3h)-2.d0*vtess(phi,flam+2.d0*d3lam,h) &
            +8.d0*vtess(phi,flam+d3lam,h)-6.d0*v)/(2.d0*d3lam*d3h*d3h)
! case 2
    elseif(((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
        .and.((h.gt.HB.and.h-d3h.lt.HB) &
        .or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggHHLam=-(vtess(phi,flam-2.d0*d3lam,h+3.d0*d3h)-4.d0*vtess(phi,flam-d3lam,h+3.d0*d3h) &
            +3.d0*vtess(phi,flam,h+3.d0*d3h)-4.d0*vtess(phi,flam-2.d0*d3lam,h+2.d0*d3h) &
            +16.d0*vtess(phi,flam-d3lam,h+2.d0*d3h)-12.d0*vtess(phi,flam,h+2.d0*d3h) &
            +5.d0*vtess(phi,flam-2.d0*d3lam,h+d3h)-20.d0*vtess(phi,flam-d3lam,h+d3h) &
            +15.d0*vtess(phi,flam,h+d3h)-2.d0*vtess(phi,flam-2.d0*d3lam,h) &
            +8.d0*vtess(phi,flam-d3lam,h)-6.d0*v)/(2.d0*d3lam*d3h*d3h)
! case 3
    elseif(((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
        .and.((h.lt.HB.and.h+d3h.gt.HB) &
        .or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggHHLam=(vtess(phi,flam+2.d0*d3lam,h-3.d0*d3h)-4.d0*vtess(phi,flam+d3lam,h-3.d0*d3h) &
            +3.d0*vtess(phi,flam,h-3.d0*d3h)-4.d0*vtess(phi,flam+2.d0*d3lam,h-2.d0*d3h) &
            +16.d0*vtess(phi,flam+d3lam,h-2.d0*d3h)-12.d0*vtess(phi,flam,h-2.d0*d3h) &
            +5.d0*vtess(phi,flam+2.d0*d3lam,h-d3h)-20.d0*vtess(phi,flam+d3lam,h-d3h) &
            +15.d0*vtess(phi,flam,h-d3h)-2.d0*vtess(phi,flam+2.d0*d3lam,h) &
            +8.d0*vtess(phi,flam+d3lam,h)-6.d0*v)/(2.d0*d3lam*d3h*d3h)
! case 4
    elseif(((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
        .and.((h.lt.HB.and.h+d3h.gt.HB) &
        .or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggHHLam=-(vtess(phi,flam-2.d0*d3lam,h-3.d0*d3h)-4.d0*vtess(phi,flam-d3lam,h-3.d0*d3h) &
            +3.d0*vtess(phi,flam,h-3.d0*d3h)-4.d0*vtess(phi,flam-2.d0*d3lam,h-2.d0*d3h) &
            +16.d0*vtess(phi,flam-d3lam,h-2.d0*d3h)-12.d0*vtess(phi,flam,h-2.d0*d3h) &
            +5.d0*vtess(phi,flam-2.d0*d3lam,h-d3h)-20.d0*vtess(phi,flam-d3lam,h-d3h) &
            +15.d0*vtess(phi,flam,h-d3h)-2.d0*vtess(phi,flam-2.d0*d3lam,h) &
            +8.d0*vtess(phi,flam-d3lam,h)-6.d0*v)/(2.d0*d3lam*d3h*d3h)
! case 5
    elseif(((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE)) &
        .and.((h.gt.HB.and.h-d3h.lt.HB)&
        .or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggHHLam=(-vtess(phi,flam+d3lam,h+3.d0*d3h)+vtess(phi,flam-d3lam,h+3.d0*d3h) &
            +4.d0*vtess(phi,flam+d3lam,h+2.d0*d3h)-4.d0*vtess(phi,flam-d3lam,h+2.d0*d3h) &
            -5.d0*vtess(phi,flam+d3lam,h+d3h)+5.d0*vtess(phi,flam-d3lam,h+d3h) &
            +2.d0*vtess(phi,flam+d3lam,h)-2.d0*vtess(phi,flam-d3lam,h))/(2.d0*d3lam*d3h*d3h)
! case 6
    elseif(((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
        .and.((h.lt.HB.and.h+d3h.gt.HB)&
        .or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggHHLam=(-vtess(phi,flam+d3lam,h-3.d0*d3h)+vtess(phi,flam-d3lam,h-3.d0*d3h) &
            +4.d0*vtess(phi,flam+d3lam,h-2.d0*d3h)-4.d0*vtess(phi,flam-d3lam,h-2.d0*d3h) &
            -5.d0*vtess(phi,flam+d3lam,h-d3h)+5.d0*vtess(phi,flam-d3lam,h-d3h) &
            +2.d0*vtess(phi,flam+d3lam,h)-2.d0*vtess(phi,flam-d3lam,h))/(2.d0*d3lam*d3h*d3h)
! case 7
    elseif(((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
        .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
        gggHHLam=(-vtess(phi,flam+2.d0*d3lam,h+d3h)+4.d0*vtess(phi,flam+d3lam,h+d3h) &
            -3.d0*vtess(phi,flam,h+d3h)+2.d0*vtess(phi,flam+2.d0*d3lam,h) &
            -8.d0*vtess(phi,flam+d3lam,h)+6.d0*v &
            -vtess(phi,flam+2.d0*d3lam,h-d3h)+4.d0*vtess(phi,flam+d3lam,h-d3h) &
            -3.d0*vtess(phi,flam,h-d3h))/(2.d0*d3lam*d3h*d3h)
! case 8
    elseif(((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
        .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
        gggHHLam=-(-vtess(phi,flam-2.d0*d3lam,h+d3h)+4.d0*vtess(phi,flam-d3lam,h+d3h) &
            -3.d0*vtess(phi,flam,h+d3h)+2.d0*vtess(phi,flam-2.d0*d3lam,h) &
            -8.d0*vtess(phi,flam-d3lam,h)+6.d0*v &
            -vtess(phi,flam-2.d0*d3lam,h-d3h)+4.d0*vtess(phi,flam-d3lam,h-d3h) &
            -3.d0*vtess(phi,flam,h-d3h))/(2.d0*d3lam*d3h*d3h)
! case 9
    else
        gggHHLam=(vtess(phi,flam+d3lam,h+d3h)-vtess(phi,flam-d3lam,h+d3h) &
            -2.d0*vtess(phi,flam+d3lam,h)+2.d0*vtess(phi,flam-d3lam,h) &
            +vtess(phi,flam+d3lam,h-d3h)-vtess(phi,flam-d3lam,h-d3h))/(2.d0*d3lam*d3h*d3h)
    endif
else
    gggHHLam=(vtess(phi,flam+d3lam,h+d3h)-vtess(phi,flam-d3lam,h+d3h) &
            -2.d0*vtess(phi,flam+d3lam,h)+2.d0*vtess(phi,flam-d3lam,h) &
            +vtess(phi,flam+d3lam,h-d3h)-vtess(phi,flam-d3lam,h-d3h))/(2.d0*d3lam*d3h*d3h)
endif
gggHHLam=(gggHHLam/c-2.d0*ggLamH)/r
!
!
!   Compute gggHHH
!
if((phi.ge.PhiS.and.phi.le.PhiN).and.(flam.ge.FlamW.and.flam.le.FlamE)) then
    if((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT)) then
        gggHHH=(-3.d0*vtess(phi,flam,h+4.d0*d3h)+14.d0*vtess(phi,flam,h+3.d0*d3h) &
            -24.d0*vtess(phi,flam,h+2.d0*d3h)+18.d0*vtess(phi,flam,h+d3h) &
            -5.0d0*v)/(2.d0*d3h*d3h*d3h)
    elseif((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT)) then
        gggHHH=(3.d0*vtess(phi,flam,h-4.d0*d3h)-14.d0*vtess(phi,flam,h-3.d0*d3h) &
            +24.d0*vtess(phi,flam,h-2.d0*d3h)-18.d0*vtess(phi,flam,h-d3h) &
            +5.0d0*v)/(2.d0*d3h*d3h*d3h)
    else
        gggHHH=(vtess(phi,flam,h+2.d0*d3h)-2.d0*vtess(phi,flam,h+d3h) &
            +2.d0*vtess(phi,flam,h-d3h)-vtess(phi,flam,h-2.d0*d3h) &
            )/(2.d0*d3h*d3h*d3h)
    endif
else
    gggHHH=(vtess(phi,flam,h+2.d0*d3h)-2.d0*vtess(phi,flam,h+d3h) &
            +2.d0*vtess(phi,flam,h-d3h)-vtess(phi,flam,h-2.d0*d3h) &
            )/(2.d0*d3h*d3h*d3h)
endif
!
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dqde(f,a,b,delta0,s)
!
!   The subroutine to compute a given integral by the double precision
!    double exponential (DE) quadrature rule.
!
!    Reference: Fukushima (2017) Precise and fast computation of gravitational
!      field of general finite body and its application to gravitational study
!      of asteroid Eros. Astron J. in printing
!
!    Input function/variables:
!      f(t): a function to be integrated 
!      a: the lower end point of the integration interval
!      b: the upper end point of the integration interval
!      delta0: the relative error tolerance (typically 1d-15)
!    Output variable:
!      s: the definite integral, s=int_a^b f(t) dt
!        (Note) "errd" is the latest relative error estimate o "s"
!
real*8 f,a,b,delta0,s
integer MMAX
real*8 SAFETY,TWOPI
parameter (MMAX=1024,SAFETY=10.d0)
parameter (TWOPI=6.2831853071795865d0)
integer m
real*8 delta,hmax,deltax,factor,deltah,h0,eph,emh,deltat,apb,bma,sr,h 
real*8 sprev,srprev,t,ep,em,eppem,xw,xa,wg,fa,fb,fapfb,errt 
real*8 errh,errd
delta=max(1.d-16,min(0.01d0,delta0))
if(delta.gt.1.d-4) then
    hmax=2.75d0+(-log10(delta))*0.75d0
elseif(delta.gt.1.d-6) then
    hmax=3.75d0+(-log10(delta))*0.5d0
elseif(delta.gt.1.d-10) then
    hmax=5.25d0+(-log10(delta))*0.25d0
else
    hmax=6.5d0+(-log10(delta))*0.125d0
endif
deltax=SAFETY*delta
factor=1.d0-log(deltax)
deltah=sqrt(deltax)
h0=hmax/factor
eph=exp(h0)
emh=1.d0/eph
deltat=exp(-emh*factor)
apb=a+b
bma=b-a
sr=f(apb*0.5d0)*(bma*0.25d0)
s=sr*(2.d0*TWOPI)
err=abs(s)*deltat
h=2.d0*h0
m=1
1 continue
    sprev=s
    srprev=sr
    t=h*0.5d0
2 continue
        em=exp(t)
        ep=TWOPI*em
        em=TWOPI/em
3 continue
            eppem=ep+em
            xw=1.d0/(1.d0+exp(ep-em))
            xa=bma*xw
            wg=xa*(1.d0-xw)
            fa=f(a+xa)*wg
            fb=f(b-xa)*wg
            fapfb=fa+fb
            sr=sr+fapfb
            s=s+fapfb*eppem
            errt=(abs(fa)+abs(fb))*eppem
            if(m.eq.1) err=err+errt*deltat
            ep=ep*eph
            em=em*emh
            if(errt.gt.err.or.xw.gt.deltah) goto 3
        t=t+h
        if(t.lt.h0) goto 2
    if(m.eq.1) then
        errh=(err/deltat)*deltah*h0
        errd=1.d0+2.d0*errh
    else
        errd=h*(abs(s-2.d0*sprev)+4.d0*abs(sr-2.d0*srprev))
    endif
    h=h*0.5d0
    m=m*2
    if(errd.gt.errh.and.m.lt.MMAX) goto 1
!if(errd.gt.errh) write(*,*)"(dqde) Required too many grids: >", MMAX
s=s*h
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dqde1(f,a,b,delta0,s)
!
!    The first copy of "dqde" to be used in computing the internal
!    line integral called by "dqde"
!
real*8 f,a,b,delta0,s
integer MMAX
real*8 SAFETY,TWOPI
parameter (MMAX=1024,SAFETY=10.d0)
parameter (TWOPI=6.2831853071795865d0)
integer m
real*8 delta,hmax,deltax,factor,deltah,h0,eph,emh,deltat,apb,bma,sr,h 
real*8 sprev,srprev,t,ep,em,eppem,xw,xa,wg,fa,fb,fapfb,errt 
real*8 errh,errd
delta=max(1.d-16,min(0.01d0,delta0))
if(delta.gt.1.d-4) then
    hmax=2.75d0+(-log10(delta))*0.75d0
elseif(delta.gt.1.d-6) then
    hmax=3.75d0+(-log10(delta))*0.5d0
elseif(delta.gt.1.d-10) then
    hmax=5.25d0+(-log10(delta))*0.25d0
else
    hmax=6.5d0+(-log10(delta))*0.125d0
endif
deltax=SAFETY*delta
factor=1.d0-log(deltax)
deltah=sqrt(deltax)
h0=hmax/factor
eph=exp(h0)
emh=1.d0/eph
deltat=exp(-emh*factor)
apb=a+b
bma=b-a
sr=f(apb*0.5d0)*(bma*0.25d0)
s=sr*(2.d0*TWOPI)
err=abs(s)*deltat
h=2.d0*h0
m=1
1 continue
    sprev=s
    srprev=sr
    t=h*0.5d0
2 continue
        em=exp(t)
        ep=TWOPI*em
        em=TWOPI/em
3 continue
            eppem=ep+em
            xw=1.d0/(1.d0+exp(ep-em))
            xa=bma*xw
            wg=xa*(1.d0-xw)
            fa=f(a+xa)*wg
            fb=f(b-xa)*wg
            fapfb=fa+fb
            sr=sr+fapfb
            s=s+fapfb*eppem
            errt=(abs(fa)+abs(fb))*eppem
            if(m.eq.1) err=err+errt*deltat
            ep=ep*eph
            em=em*emh
            if(errt.gt.err.or.xw.gt.deltah) goto 3
        t=t+h
        if(t.lt.h0) goto 2
    if(m.eq.1) then
        errh=(err/deltat)*deltah*h0
        errd=1.d0+2.d0*errh
    else
        errd=h*(abs(s-2.d0*sprev)+4.d0*abs(sr-2.d0*srprev))
    endif
    h=h*0.5d0
    m=m*2
    if(errd.gt.errh.and.m.lt.MMAX) goto 1
!if(errd.gt.errh) write(*,*)"(dqde1) Required too many grids: >", MMAX
s=s*h
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8 function dlog1p(x0)
!
!    The double precision Fortran function to compute precisely log1p(x),
!    a special logarithm function defined as log1p(x) = ln(1+x),
!    by the (7,6) rational minimax approximation evaluated by Estrin's scheme
!
!    Reference: Fukushima (2017) Precise and fast computation of gravitational
!      field of general finite body and its application to gravitational study
!      of asteroid Eros. Astron J. in printing
!
real*8 x0,x1,x2,x4,a1,a2,a3,a4,a5,a6,a7,b1,b2,b3,b4,b5,b6
parameter (a1=0.99999999999999999405d0,a2=2.45235728562912886048d0)
parameter (a3=2.17053627298972253249d0,a4=0.83928994566440838378d0)
parameter (a5=0.13520496594993836479d0,a6=0.00682631751459270270d0)
parameter (a7=0.00002291289324181940d0)
parameter (b1=2.95235728562912599232d0,b2=3.31338158247117791600d0)
parameter (b3=1.76186164168333482938d0,b4=0.44976458082070468584d0)
parameter (b5=0.04896199808811261680d0,b6=0.00157389087429218809d0)
if(x0.lt.-0.5d0.or.x0.gt.1.d0) then
    dlog1p=log(1.d0+x0)
elseif(x0.lt.0.d0) then
    x1=-x0/(1.d0+x0)
    x2=x1*x1; x4=x2*x2
    dlog1p=-x1*(((a1+x1*a2)+x2*(a3+x1*a4))+x4*((a5+x1*a6)+x2*a7)) &
        /(((1.d0+x1*b1)+x2*(b2+x1*b3))+x4*((b4+x1*b5)+x2*b6))
else
    x1=x0
    x2=x1*x1; x4=x2*x2
    dlog1p=x1*(((a1+x1*a2)+x2*(a3+x1*a4))+x4*((a5+x1*a6)+x2*a7)) &
        /(((1.d0+x1*b1)+x2*(b2+x1*b3))+x4*((b4+x1*b5)+x2*b6))
endif
return; end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    Output of "xtessgc"
!      (Note) It took 11.750 s at a PC with an Intel Core i5-10400 running at 2.90 GHz
!      The ifort version is
!           Intel(R) Fortran Intel(R) 64 Compiler Classic for applications running on 
!            Intel(R) 64, Version 2021.4.0 Build 20210910_000000
!            Copyright (C) 1985-2021 Intel Corporation.  All rights reserved.
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        # delta=  1.0000000E-16
!     # R0,HT,HB=  6.3800000E+03  1.0000000E+01 -4.0000000E+01
!    # PhiN,PhiS=  2.8000000E+01  2.7000000E+01
! # LamdaE,LamdaW  8.7000000E+01  8.6000000E+01
!   # phi,lambda=  2.7500000E+01  8.6500000E+01
!             # H                        V
!               #                     gPhi                  gLambda                       gH
!               #              GammaPhiPhi           GammaPhiLambda                GammaPhiH
!               #        GammaLambdaLambda             GammaLambdaH                  GammaHH
!               #             gggPhiPhiPhi             gggPhiPhiLam               gggPhiPhiH
!               #               gggPhiLamH
!               #             gggLamLamPhi             gggLamLamLam               gggLamLamH
!               #                 gggHHPhi                 gggHHLam                   gggHHH
!  -9.5000000E+01    1.522385430859932E-04
!                   -6.913731817795242E-10    0.000000000000000E+00    1.564130608898577E-06
!                   -1.387341662716454E-08    0.000000000000000E+00   -1.613905285522335E-11
!                   -1.505917753693766E-08    0.000000000000000E+00    2.893259481144160E-08
!                    5.714368783126311E-13    0.000000000000000E+00   -3.280331145834972E-10
!                    0.000000000000000E+00
!                   -1.621242262241418E-13    0.000000000000000E+00   -3.787855739041423E-10
!                   -4.113184773586338E-13    0.000000000000000E+00    7.068188008468201E-10
!  -8.5000000E+01    1.694526871210555E-04
!                   -8.747066001775610E-10    0.000000000000000E+00    1.892198142232928E-06
!                   -1.760580950750948E-08    0.000000000000000E+00   -2.066069891693934E-11
!                   -1.943435108474880E-08    0.000000000000000E+00    3.704016020364891E-08
!                    8.350338006581876E-13    0.000000000000000E+00   -4.212082222628171E-10
!                    0.000000000000000E+00
!                   -3.427478758986157E-13    0.000000000000000E+00   -5.010627017808653E-10
!                   -4.893810382734468E-13    0.000000000000000E+00    9.222703533063007E-10
!  -7.5000000E+01    1.903905934236350E-04
!                   -1.106672579981757E-09    0.000000000000000E+00    2.312822513933103E-06
!                   -2.233996332630832E-08    0.000000000000000E+00   -2.578793931005678E-11
!                   -2.516840463620129E-08    0.000000000000000E+00    4.750836924785179E-08
!                    1.213314919246471E-12    0.000000000000000E+00   -5.266975081273169E-10
!                    0.000000000000000E+00
!                   -6.901140862846868E-13    0.000000000000000E+00   -6.492000960788878E-10
!                   -5.220997124390895E-13    0.000000000000000E+00    1.175897626205599E-09
!  -6.5000000E+01    2.161011857153642E-04
!                   -1.389893097215668E-09    0.000000000000000E+00    2.851074036420383E-06
!                   -2.812502336398109E-08    0.000000000000000E+00   -3.071419320558526E-11
!                   -3.244537518352263E-08    0.000000000000000E+00    6.057039583336358E-08
!                    1.737788591516231E-12    0.000000000000000E+00   -6.271702476584572E-10
!                    0.000000000000000E+00
!                   -1.304130314003106E-12    0.000000000000000E+00   -8.044453840992871E-10
!                   -4.366233028625275E-13    0.000000000000000E+00    1.431617406352432E-09
!  -5.5000000E+01    2.478880058505530E-04
!                   -1.715013458072008E-09    0.000000000000000E+00    3.531846611828726E-06
!                   -3.475500622664502E-08    0.000000000000000E+00   -3.383415931884386E-11
!                   -4.112692458673271E-08    0.000000000000000E+00    7.588192820201063E-08
!                    2.414714016685705E-12    0.000000000000000E+00   -6.886858755137172E-10
!                    0.000000000000000E+00
!                   -2.264756803821667E-12    0.000000000000000E+00   -9.195329056654826E-10
!                   -1.490810816320179E-13    0.000000000000000E+00    1.608217657409740E-09
!  -4.5000000E+01    2.872708844104648E-04
!                   -2.053259639166730E-09    0.000000000000000E+00    4.371789329686412E-06
!                   -4.161825122777698E-08    0.000000000000000E+00   -3.297657053417959E-11
!                   -5.045690322613585E-08    0.000000000000000E+00    9.207514074757049E-08
!                    3.182617982202230E-12    0.000000000000000E+00   -6.664205585018362E-10
!                    0.000000000000000E+00
!                   -3.535910435666732E-12    0.000000000000000E+00   -9.211681549242437E-10
!                    3.501670375520160E-13    0.000000000000000E+00    1.587587459139255E-09
!  -3.5000000E+01    3.319910086311147E-04
!                   -2.355567528758797E-09    0.000000000000000E+00    3.825615688316108E-06
!                   -4.793266758433300E-08    0.000000000000000E+00   -2.649377265963343E-11
!                   -5.920230924282540E-08    0.000000000000000E+00   -2.015876840700058E-07
!                    3.914133937066427E-12    0.000000000000000E+00   -5.759221069932876E-10
!                    0.000000000000000E+00
!                   -4.855887515236279E-12    0.000000000000000E+00   -7.963898741599008E-10
!                    9.441390727118867E-13    0.000000000000000E+00    1.372314254145308E-09
!  -2.5000000E+01    3.603747131612369E-04
!                   -2.564755213201881E-09    0.000000000000000E+00    1.869421744100308E-06
!                   -5.251553452556305E-08    0.000000000000000E+00   -1.457160314363856E-11
!                   -6.555716477277731E-08    0.000000000000000E+00   -1.906499679045349E-07
!                    4.434796757585294E-12    0.000000000000000E+00   -3.248997275396590E-10
!                    0.000000000000000E+00
!                   -5.839363027729824E-12    0.000000000000000E+00   -4.502407667635140E-10
!                    1.400413537246880E-12    0.000000000000000E+00    7.751405798112112E-10
!  -1.5000000E+01    3.696350462514114E-04
!                   -2.636222132600263E-09    0.000000000000000E+00   -1.064260798473822E-08
!                   -5.421371347135290E-08    0.000000000000000E+00    5.471010114351431E-13
!                   -6.787327496614387E-08    0.000000000000000E+00   -1.866356612175997E-07
!                    4.614951658308923E-12    0.000000000000000E+00   -9.558001225807770E-12
!                    0.000000000000000E+00
!                   -6.180386693519267E-12    0.000000000000000E+00   -4.978458100341342E-12
!                    1.563289723448127E-12    0.000000000000000E+00    1.453725341791131E-11
!  -5.0000000E+00    3.601669342163731E-04
!                   -2.554314977920888E-09    0.000000000000000E+00   -1.889139936024406E-06
!                   -5.272632315535408E-08    0.000000000000000E+00    1.551695685477210E-11
!                   -6.568213752784121E-08    0.000000000000000E+00   -1.903142199554053E-07
!                    4.407911359767593E-12    0.000000000000000E+00    2.999782786368029E-10
!                    0.000000000000000E+00
!                   -5.779662830340790E-12    0.000000000000000E+00    4.327759716468489E-10
!                    1.375699749661199E-12    0.000000000000000E+00   -7.327531982954847E-10
!   5.0000000E+00    3.316111352803727E-04
!                   -2.337440212306251E-09    0.000000000000000E+00   -3.839348173993401E-06
!                   -4.846251253622977E-08    0.000000000000000E+00    2.705975669803280E-11
!                   -5.959167367260831E-08    0.000000000000000E+00   -2.006684789352708E-07
!                    3.863525336150733E-12    0.000000000000000E+00    5.358970730199236E-10
!                    0.000000000000000E+00
!                   -4.762001588928537E-12    0.000000000000000E+00    7.595877798290157E-10
!                    8.981387149968035E-13    0.000000000000000E+00   -1.295484039516153E-09
!   1.5000000E+01    2.868101440205006E-04
!                   -2.031662182754117E-09    0.000000000000000E+00   -4.374218164463007E-06
!                   -4.215634275488205E-08    0.000000000000000E+00    3.311074245470626E-11
!                   -5.084513320252664E-08    0.000000000000000E+00    9.300147179106866E-08
!                    3.126922120187886E-12    0.000000000000000E+00    7.047555641314654E-10
!                    0.000000000000000E+00
!                   -3.437993268490088E-12    0.000000000000000E+00    9.583646775444268E-10
!                    3.123457973501095E-13    0.000000000000000E+00   -1.663117893407240E-09
!   2.5000000E+01    2.474382990881565E-04
!                   -1.693807030640827E-09    0.000000000000000E+00   -3.528159772641218E-06
!                   -3.499502159672540E-08    0.000000000000000E+00    3.364860021977298E-11
!                   -4.124208017518744E-08    0.000000000000000E+00    7.623710193663058E-08
!                    2.363492339594734E-12    0.000000000000000E+00    7.107103580901066E-10
!                    0.000000000000000E+00
!                   -2.191230522972409E-12    0.000000000000000E+00    9.381141168900972E-10
!                   -1.743259816939549E-13    0.000000000000000E+00   -1.648822777939899E-09
!   3.5000000E+01    2.157004300813908E-04
!                   -1.371418543151609E-09    0.000000000000000E+00   -2.845439875476928E-06
!                   -2.820499512163055E-08    0.000000000000000E+00    3.038032155719775E-11
!                   -3.243694585603805E-08    0.000000000000000E+00    6.064194241272785E-08
!                    1.699602397752741E-12    0.000000000000000E+00    6.380617164917233E-10
!                    0.000000000000000E+00
!                   -1.257776550844332E-12    0.000000000000000E+00    8.117647241191198E-10
!                   -4.441530324476631E-13    0.000000000000000E+00   -1.449822522358936E-09
!   4.5000000E+01    1.900473207915235E-04
!                   -1.091721162752677E-09    0.000000000000000E+00   -2.307151583407502E-06
!                   -2.234592666376857E-08    0.000000000000000E+00    2.543352974381340E-11
!                   -2.511765950393008E-08    0.000000000000000E+00    4.746358236367228E-08
!                    1.185399097008428E-12    0.000000000000000E+00    5.312915017796743E-10
!                    0.000000000000000E+00
!                   -6.660675120401457E-13    0.000000000000000E+00    6.510953942963361E-10
!                   -5.205807666824120E-13    0.000000000000000E+00   -1.182386766194161E-09
!   5.5000000E+01    1.691630798625040E-04
!                   -8.631050177663385E-10    0.000000000000000E+00   -1.887193033063033E-06
!                   -1.758322031949842E-08    0.000000000000000E+00    2.035117683351480E-11
!                   -1.937716736234116E-08    0.000000000000000E+00    3.696038632630143E-08
!                    8.153217323658647E-13    0.000000000000000E+00    4.227216114279704E-10
!                    0.000000000000000E+00
!                   -3.310206416580767E-13    0.000000000000000E+00    5.008213723779776E-10
!                   -4.845419942209920E-13    0.000000000000000E+00   -9.235433167211372E-10
!   6.5000000E+01    1.519948957924417E-04
!                   -6.825567987216375E-10    0.000000000000000E+00   -1.559945026069055E-06
!                   -1.384333972980420E-08    0.000000000000000E+00    1.589248141073156E-11
!                   -1.500835978344898E-08    0.000000000000000E+00    2.885169384761639E-08
!                    5.600892396096175E-13    0.000000000000000E+00    3.282168080415752E-10
!                    0.000000000000000E+00
!                   -1.548814465900948E-13    0.000000000000000E+00    3.779064407573439E-10
!                   -4.039985534480120E-13    0.000000000000000E+00   -7.061250968762420E-10
!   7.5000000E+01    1.377280082166444E-04
!                   -5.424089697527433E-10    0.000000000000000E+00   -1.303706863329740E-06
!                   -1.095668286747153E-08    0.000000000000000E+00    1.227954690632847E-11
!                   -1.172281568705910E-08    0.000000000000000E+00    2.267949564325045E-08
!                    3.869955618871021E-13    0.000000000000000E+00    2.521948138692616E-10
!                    0.000000000000000E+00
!                   -6.704393634927498E-14    0.000000000000000E+00    2.836831349150813E-10
!                   -3.192369915676482E-13    0.000000000000000E+00   -5.358791180978119E-10
!   8.5000000E+01    1.257414395776617E-04
!                   -4.342891730356939E-10    0.000000000000000E+00   -1.101402778814564E-06
!                   -8.741565399837410E-09    0.000000000000000E+00    9.468766808869749E-12
!                   -9.254182417354188E-09    0.000000000000000E+00    1.799575435499842E-08
!                    2.696192797525769E-13    0.000000000000000E+00    1.934773135408658E-10
!                    0.000000000000000E+00
!                   -2.523479310720294E-14    0.000000000000000E+00    2.135899251384668E-10
!                   -2.458282675013455E-13    0.000000000000000E+00   -4.070675856947556E-10
!   9.5000000E+01    1.155637184360952E-04
!                   -3.508161779288250E-10    0.000000000000000E+00   -9.400780632703064E-07
!                   -7.039685098356269E-09    0.000000000000000E+00    7.324942917395352E-12
!                   -7.389354992274144E-09    0.000000000000000E+00    1.442904760230579E-08
!                    1.905609856351106E-13    0.000000000000000E+00    1.489935749029994E-10
!                    0.000000000000000E+00
!                   -4.995125651117251E-15    0.000000000000000E+00    1.620256562135275E-10
!                   -1.862357340269648E-13    0.000000000000000E+00   -3.110190334571592E-10
!   1.0500000E+02    1.068357959386217E-04
!                   -2.860576345917093E-10    0.000000000000000E+00   -8.100601080971190E-07
!                   -5.725033012710587E-09    0.000000000000000E+00    5.702290748351356E-12
!                   -5.968285630016786E-09    0.000000000000000E+00    1.169332767098393E-08
!                    1.370779441080873E-13    0.000000000000000E+00    1.155263157076373E-10
!                    0.000000000000000E+00
!                    3.772347352689872E-15    0.000000000000000E+00    1.241184287375809E-10
!                   -1.404830112969229E-13    0.000000000000000E+00   -2.396458579951686E-10
!   1.1500000E+02    9.928233337092821E-05
!                   -2.354547369303172E-10    0.000000000000000E+00   -7.041574277909766E-07
!                   -4.701624419951508E-09    0.000000000000000E+00    4.475025105084251E-12
!                   -4.874108930517592E-09    0.000000000000000E+00    9.575736618672118E-09
!                    9.975250657444415E-14    0.000000000000000E+00    9.033981980784378E-11
!                    0.000000000000000E+00
!                    6.934436960512076E-15    0.000000000000000E+00    9.611068675786290E-11
!                   -1.066451874060181E-13    0.000000000000000E+00   -1.864497055911042E-10
!   1.2500000E+02    9.269027708860836E-05
!                   -1.955796981354374E-10    0.000000000000000E+00   -6.170099766597611E-07
!                   -3.897813433179675E-09    0.000000000000000E+00    3.542452393056092E-12
!                   -4.022351452862562E-09    0.000000000000000E+00    7.920162281693801E-09
!                    7.319034682818025E-14    0.000000000000000E+00    7.129839787092590E-11
!                    0.000000000000000E+00
!                    7.626737191403898E-15    0.000000000000000E+00    7.524675356449019E-11
!                   -8.140015604110879E-14    0.000000000000000E+00   -1.465461683360995E-10
!   1.3500000E+02    8.689312433613946E-05
!                   -1.638766556509691E-10    0.000000000000000E+00   -5.445976434900247E-07
!                   -3.260546696865931E-09    0.000000000000000E+00    2.829623365566593E-12
!                   -3.351998343528302E-09    0.000000000000000E+00    6.612545978082935E-09
!                    5.485303849002564E-14    0.000000000000000E+00    5.680329123124639E-11
!                    0.000000000000000E+00
!                    7.254173106525802E-15    0.000000000000000E+00    5.955367906126566E-11
!                   -6.239971207457662E-14    0.000000000000000E+00   -1.163569832114101E-10
!   1.4500000E+02    8.175942046629835E-05
!                   -1.384428188603235E-10    0.000000000000000E+00   -4.838801343670081E-07
!                   -2.750555989344125E-09    0.000000000000000E+00    2.280479153700626E-12
!                   -2.818773501162813E-09    0.000000000000000E+00    5.569327866738967E-09
!                    4.152360353159709E-14    0.000000000000000E+00    4.567707661133825E-11
!                    0.000000000000000E+00
!                    6.414746746602509E-15    0.000000000000000E+00    4.762614597356267E-11
!                   -4.828379148917547E-14    0.000000000000000E+00   -9.330274277201375E-11
!   1.5500000E+02    7.718433330381282E-05
!                   -1.178576950818463E-10    0.000000000000000E+00   -4.325366576174139E-07
!                   -2.338665048770633E-09    0.000000000000000E+00    1.854120904701192E-12
!                   -2.390292477094614E-09    0.000000000000000E+00    4.728959598012357E-09
!                    3.202573480499018E-14    0.000000000000000E+00    3.706119129960413E-11
!                    0.000000000000000E+00
!                    5.770928561895545E-15    0.000000000000000E+00    3.846441769207098E-11
!                   -3.751021779107088E-14    0.000000000000000E+00   -7.552620795589124E-11
!   1.6500000E+02    7.308344678779021E-05
!                   -1.010539761487662E-10    0.000000000000000E+00   -3.887784141257901E-07
!                   -2.003069470387724E-09    0.000000000000000E+00    1.520021361173847E-12
!                   -2.042663925029824E-09    0.000000000000000E+00    4.045729268493047E-09
!                    2.465985881253098E-14    0.000000000000000E+00    3.032928043604936E-11
!                    0.000000000000000E+00
!                    4.762886985917925E-15    0.000000000000000E+00    3.135493874552755E-11
!                   -2.957369910860912E-14    0.000000000000000E+00   -6.168573579240012E-11
!   1.7500000E+02    6.938815426986305E-05
!                   -8.722577072088447E-11    0.000000000000000E+00   -3.512132393063509E-07
!                   -1.727335015796498E-09    0.000000000000000E+00    1.256037829399525E-12
!                   -1.758078359898185E-09    0.000000000000000E+00    3.485412429609487E-09
!                    1.949472009928672E-14    0.000000000000000E+00    2.502298536962056E-11
!                    0.000000000000000E+00
!                    4.172803122110900E-15    0.000000000000000E+00    2.578314114053557E-11
!                   -2.339944012528545E-14    0.000000000000000E+00   -5.080663489028059E-11
!   1.8500000E+02    6.604220927845374E-05
!                   -7.575771705036650E-11    0.000000000000000E+00   -3.187473618801081E-07
!                   -1.498993012879632E-09    0.000000000000000E+00    1.045300997498653E-12
!                   -1.523135745644453E-09    0.000000000000000E+00    3.022127429493267E-09
!                    1.529965892808741E-14    0.000000000000000E+00    2.080359982958226E-11
!                    0.000000000000000E+00
!                    3.290517949740219E-15    0.000000000000000E+00    2.137432696302121E-11
!                   -1.889871154888758E-14    0.000000000000000E+00   -4.217803094748410E-11
!   1.9500000E+02    6.299911914505307E-05
!                   -6.617778348804906E-11    0.000000000000000E+00   -2.905135868753384E-07
!                   -1.308482993558236E-09    0.000000000000000E+00    8.764510796233346E-13
!                   -1.327639244998674E-09    0.000000000000000E+00    2.636123200298561E-09
!                    1.219580137919097E-14    0.000000000000000E+00    1.742040717045115E-11
!                    0.000000000000000E+00
!                    2.690777847345353E-15    0.000000000000000E+00    1.785402487941451E-11
!                   -1.519253544881670E-14    0.000000000000000E+00   -3.527463777999409E-11
!   2.0500000E+02    6.022015719417137E-05
!                   -5.812021753502276E-11    0.000000000000000E+00   -2.658183543094025E-07
!                   -1.148425904420230E-09    0.000000000000000E+00    7.395320754572204E-13
!                   -1.163773382969256E-09    0.000000000000000E+00    2.312197246119771E-09
!                    9.926170149513865E-15    0.000000000000000E+00    1.468656513873652E-11
!                    0.000000000000000E+00
!                    2.363681403612682E-15    0.000000000000000E+00    1.501979761797691E-11
!                   -1.238291510643198E-14    0.000000000000000E+00   -2.970582754594552E-11
!   2.1500000E+02    5.767283276742178E-05
!                   -5.129914737953770E-11    0.000000000000000E+00   -2.441023741763074E-07
!                   -1.013065267991374E-09    0.000000000000000E+00    6.282389313284449E-13
!                   -1.025472707678019E-09    0.000000000000000E+00    2.038537560927707E-09
!                    8.142391601546663E-15    0.000000000000000E+00    1.245998568221136E-11
!                    0.000000000000000E+00
!                    2.022184103719468E-15    0.000000000000000E+00    1.271865855574601E-11
!                   -1.009923020263910E-14    0.000000000000000E+00   -2.517918173188958E-11
!   2.2500000E+02    5.532970311619008E-05
!                   -4.549026728311057E-11    0.000000000000000E+00   -2.249110992018623E-07
!                   -8.978893823358672E-10    0.000000000000000E+00    5.364004826272384E-13
!                   -9.080023938892072E-10    0.000000000000000E+00    1.805888513698460E-09
!                    6.656018779915141E-15    0.000000000000000E+00    1.063442789393701E-11
!                    0.000000000000000E+00
!                    1.720823870361711E-15    0.000000000000000E+00    1.083715060460564E-11
!                   -8.218967870841289E-15    0.000000000000000E+00   -2.147140652451971E-11
!   2.3500000E+02    5.316744283951891E-05
!                   -4.051502936112462E-11    0.000000000000000E+00   -2.078723714455611E-07
!                   -7.993184764310249E-10    0.000000000000000E+00    4.608336911214137E-13
!                   -8.076268253254436E-10    0.000000000000000E+00    1.606946455241786E-09
!                    5.420416675961874E-15    0.000000000000000E+00    9.126811495364922E-12
!                    0.000000000000000E+00
!                    1.376064633648046E-15    0.000000000000000E+00    9.287249961671811E-12
!                   -6.820078301209406E-15    0.000000000000000E+00   -1.841401754957270E-11
!   2.4500000E+02    5.116610904202815E-05
!                   -3.623095682449497E-11    0.000000000000000E+00   -1.926793662886023E-07
!                   -7.145051646707031E-10    0.000000000000000E+00    3.978605955091317E-13
!                   -7.213773985284890E-10    0.000000000000000E+00    1.435881170388446E-09
!                    4.496913684464920E-15    0.000000000000000E+00    7.874216409261084E-12
!                    0.000000000000000E+00
!                    1.152310874454832E-15    0.000000000000000E+00    8.002139114981311E-12
!                   -5.792203299094049E-15    0.000000000000000E+00   -1.587618543086920E-11
!   2.5500000E+02    4.930855649651286E-05
!                   -3.252370700956010E-11    0.000000000000000E+00   -1.790774575837581E-07
!                   -6.411538484083760E-10    0.000000000000000E+00    3.450907210372193E-13
!                   -6.468772022143885E-10    0.000000000000000E+00    1.288030080178117E-09
!                    3.753021022247467E-15    0.000000000000000E+00    6.826899436994928E-12
!                    0.000000000000000E+00
!                    1.007065653483266E-15    0.000000000000000E+00    6.929780444914350E-12
!                   -4.820128543011578E-15    0.000000000000000E+00   -1.375660244063963E-11
!   2.6500000E+02    4.757996872078648E-05
!                   -2.930002442042743E-11    0.000000000000000E+00   -1.668540276287930E-07
!                   -5.774123234724967E-10    0.000000000000000E+00    3.007585040241532E-13
!                   -5.822099218211791E-10    0.000000000000000E+00    1.159621692240454E-09
!                    3.153909496531086E-15    0.000000000000000E+00    5.946281845758397E-12
!                    0.000000000000000E+00
!                    8.638063734903088E-16    0.000000000000000E+00    6.029646817567191E-12
!                   -4.113322675839351E-15    0.000000000000000E+00   -1.197619050584830E-11
!   2.7500000E+02    4.596747933454686E-05
!                   -2.648492557584876E-11    0.000000000000000E+00   -1.558305006175582E-07
!                   -5.217727389464347E-10    0.000000000000000E+00    2.632785455722527E-13
!                   -5.258166525598973E-10    0.000000000000000E+00    1.047589512757920E-09
!                    2.656595670811530E-15    0.000000000000000E+00    5.202038874125444E-12
!                    0.000000000000000E+00
!                    6.959028913124725E-16    0.000000000000000E+00    5.269956901972507E-12
!                   -3.468227289266320E-15    0.000000000000000E+00   -1.047209578802512E-11
!   2.8500000E+02    4.445986426369510E-05
!                   -2.401592365720716E-11    0.000000000000000E+00   -1.458560632822812E-07
!                   -4.729987203902264E-10    0.000000000000000E+00    2.313760957739202E-13
!                   -4.764253629564484E-10    0.000000000000000E+00    9.494266869605439E-10
!                    2.260748806421330E-15    0.000000000000000E+00    4.569558346837683E-12
!                    0.000000000000000E+00
!                    5.919060208512516E-16    0.000000000000000E+00    4.625385146717483E-12
!                   -2.983922776267628E-15    0.000000000000000E+00   -9.195145189322153E-12
!   2.9500000E+02    4.304728994459764E-05
!                   -2.184203249430225E-11    0.000000000000000E+00   -1.368026878325614E-07
!                   -4.300723360646307E-10    0.000000000000000E+00    2.040287185915461E-13
!                   -4.329928028963195E-10    0.000000000000000E+00    8.630637676541201E-10
!                    1.977890114731782E-15    0.000000000000000E+00    4.029689216062964E-12
!                    0.000000000000000E+00
!                    5.455264680405333E-16    0.000000000000000E+00    4.075724385315122E-12
!                   -2.518754009469029E-15    0.000000000000000E+00   -8.105410589257435E-12
!  The total time is   11.7500530000000      seconds.
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!