!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program xqtessgc
!
!   Test driver of "qtessgc", the quadruple precision Fortran subroutine to compute
!   the gravitational potential/acceleration/gradient tensor/curvatures of a tesseroid.
!
!    Reference: 
!           Fukushima T (2018) Accurate computation of gravitational field of a tesseroid.
!                Journal of Geodesy 92(12):1371–1386, doi:10.1007/s00190-018-1126-2
!           Deng XL (2023) Corrections to: “Accurate computation of gravitational field 
!                of a tesseroid” by Fukushima (2018) in J. Geod. 92(12):1371–1386, Journal of Geodesy,
!                doi: 10.1007/s00190-022-01673-2
!
implicit integer (i-n)
implicit real*16 (q)
implicit real*8 (a-h,o,p,r-z)
!
qPI=atan(1.q0)*4.q0
qraddeg=qPI/180.q0
qdelta=1.q-33
write (*,"(a15,1pe15.7)") "# qdelta=", qdelta
!
!   A sample tesseroid covering an eastern Himalaya:
!   (Note) Mt. Everest: h = 8.848 km, phi = 27 deg 59 min, lat = 86 deg 55 min
!
qR0=6380.q0
qHB=-40.q0
qHT=10.q0
qPhiS=27.q0*qraddeg
qPhiN=28.q0*qraddeg
qFlamW=86.q0*qraddeg
qFlamE=87.q0*qraddeg
!
!   Test evaluation points are located along a radial straight line passing through
!    the geometrical center of the tesseroid
!
qphi=27.5q0*qraddeg
qflam=86.5q0*qraddeg
!
write (*,"(a15,1p3e15.7)") "# R0,HT,HB=",qR0,qHT,qHB
write (*,"(a15,1p2e15.7)") "# PhiN,PhiS=",qPhiN/qraddeg,qPhiS/qraddeg
write (*,"(a15,1p2e15.7)") "# LamdaE,LamdaW=",qFlamE/qraddeg,qFlamW/qraddeg
write (*,"(a15,1p2e15.7)") "# phi,lambda=",qphi/qraddeg,qflam/qraddeg
write (*,"(a15,a45)") "# H","V"
write (*,"(a15,3a45)") "#","gPhi","gLambda","gH"
write (*,"(a15,3a45)") "#","GammaPhiPhi","GammaPhiLambda","GammaPhiH"
write (*,"(a15,3a45)") "#","GammaLambdaLambda","GammaLambdaH","GammaHH"
write (*,"(a15,3a45)") "#","gggPhiPhiPhi","gggPhiPhiLam","gggPhiPhiH"
write (*,"(a15,1a45)") "#","gggPhiLamH"
write (*,"(a15,3a45)") "#","gggLamLamPhi","gggLamLamLam","gggLamLamH"
write (*,"(a15,3a45)") "#","gggHHPhi","gggHHLam","gggHHH"
!
!    Compute the field at the reference sphere
!
qh=0.q0
!
call cpu_time(time_begin)
call qtess(qPhiN,qPhiS,qFlamE,qFlamW,qHT,qHB,qR0,qdelta, &
    qphi,qflam,qh,qv,qgPhi,qgLam,qgH, &
    qggPhiPhi,qggPhiLam,qggPhiH,qggLamLam,qggLamH,qggHH, &
    qgggPhiPhiPhi,qgggPhiPhiLam,qgggPhiPhiH,qgggPhiLamH,qgggLamLamPhi, &
    qgggLamLamLam,qgggLamLamH,qgggHHPhi,qgggHHLam,qgggHHH)
write (*,"(1pe15.7,1pe45.35)") qh,qv
write (*,"(15x,1p3e45.35)") qgPhi,qgLam,qgH
write (*,"(15x,1p3e45.35)") qggPhiPhi,qggPhiLam,qggPhiH
write (*,"(15x,1p3e45.35)") qggLamLam,qggLamH,qggHH
write (*,"(15x,1p3e45.35)") qgggPhiPhiPhi,qgggPhiPhiLam,qgggPhiPhiH
write (*,"(15x,1p1e45.35)") qgggPhiLamH
write (*,"(15x,1p3e45.35)") qgggLamLamPhi,qgggLamLamLam,qgggLamLamH
write (*,"(15x,1p3e45.35)") qgggHHPhi,qgggHHLam,qgggHHH
!
call cpu_time(time_end)
write(*,*) "The total time is",time_end-time_begin,"seconds."
stop
end program xqtessgc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine qtess(PhiN0,PhiS0,FlamE0,FlamW0,HT0,HB0,R00,delta0, &
    phi,flam,h, &
    v,gPhi,gLam,gH,ggPhiPhi,ggPhiLam,ggPhiH,ggLamLam,ggLamH,ggHH, &
    gggPhiPhiPhi,gggPhiPhiLam,gggPhiPhiH,gggPhiLamH,gggLamLamPhi, &
    gggLamLamLam,gggLamLamH,gggHHPhi,gggHHLam,gggHHH)
!
!   The quadruple precision Fortran subroutine to compute
!   the gravitational potential, V, the gravitational acceleration vector, g,
!   the gravity gradient tensor, Gamma, 
!   and the gravitational curvatures, GGG, of a tesseroid.
!
!   Input parameters:
!       PhiN0: Latitude of the north end point of the tesseroid
!       PhiS0: Latitude of the south end point of the tesseroid, PhiS0 < PhiN0
!       FlamE0: Longitude of the east end point of the tesseroid
!       FlamW0: Longitude of the west end point of the tesseroid, FlamW0 < FlamE0
!       HT0: Height of the top end point of the tesseroid
!       HB0: Height of the bottom end point of the tesseroid, HB0 < HT0
!       R00: Radius of the reference spherical surface from which H is counted
!       delta0: the relative error tolerance (typical value is 1q-24 or 1q-30)
!   Input variables:
!       phi: Latitude of the evaluation point
!       flam: Longitude of the evaluation point
!       h: Height from the reference sphere of the evaluation point
!   Output variables:
!       v: Normalized gravitational potential at the evaluation point
!       (Note) normalization constant: V0 = G rho R00^2
!       gPhi: Latitude component of normalized gravitational acceleration vector
!       gLam: Longitude component of normalized gravitational acceleration vector
!       gH: Height component of normalized gravitational acceleration vector
!       (Note) normalization constant: g0 = G rho R00^2
!       ggPhiPhi: Latitude-latitude component of normalized Gamma
!       ggPhiLam: Latitude-longitude component of normalized Gamma
!       ggPhiH: Latitude-height component of normalized Gamma
!       ggLamLam: Longitude-longitude component of normalized Gamma
!       ggLamH: Longitude-height component of normalized Gamma
!       ggHH: Height-height component of normalized Gamma
!       (Note) normalization constant: Gamma0 = G rho R00^2
! add the output variables of ggg
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
implicit real*16 (a-h,o-z)
!
external qvtess
!
real*16 eps
common /qepsV/eps
real*16 R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
common /qparamV/R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
real*16 delta,delta1,delta2,d1phi,d1lam,d1h,d2phi,d2lam,d2h,delta3,d3phi,d3lam,d3h
common /qdeltaV/delta,delta1,delta2,d1phi,d1lam,d1h,d2phi,d2lam,d2h,delta3,d3phi,d3lam,d3h
!
PhiN=PhiN0
PhiS=PhiS0
FlamE=FlamE0
FlamW=FlamW0
HT=HT0
HB=HB0
R0=R00
eps=delta0
delta=delta0
delta1=delta**(1.q0/3.q0)
delta2=sqrt(sqrt(delta))
delta3=delta**(1.q0/5.q0)
!
deltaH=HT-HB
deltaPhi=PhiN-PhiS
deltaFlam=FlamE-FlamW
beta=deltaH/R0
RC=R0+(HT+HB)*0.5q0
PhiC=(PhiN+PhiS)*0.5q0
FlamC=(FlamE+FlamW)*0.5q0
sinPhiC=sin(PhiC)
cosPhiC=cos(PhiC)
!
r=R0+h
sinphi=sin(phi)
cosphi=cos(phi)
cosdlam=cos(flam-FlamC)
scale0=0.5q0*sqrt(deltaH*deltaH+RC*RC*(deltaPhi*deltaPhi &
    +cosPhiC*cosPhiC*deltaFlam*deltaFlam))
scale1=sqrt(r*r+RC*RC-2.q0*RC*r*(sinphi*sinPhiC+cosphi*cosPhiC*cosdlam))
scale=max(scale0,scale1)
d1h=scale*delta1
d2h=scale*delta2
! add d1phi, d1lam
d1phi=d1h/r
d1lam=d1h/(r*cosphi)
! add d2phi, d2lam
d2phi=d2h/r
d2lam=d2h/(r*cosphi)
d3h=scale*delta3
d3phi=d3h/r
d3lam=d3h/(r*cosphi)
!
v=qvtess(phi,flam,h)
call qgtess(qvtess,phi,flam,h,v,gPhi,gLam,gH)
call qggtess(qvtess,phi,flam,h,v,gPhi,gLam,gH, &
    ggPhiPhi,ggPhiLam,ggPhiH,ggLamLam,ggLamH,ggHH)
call qgggtess(qvtess,phi,flam,h,v,gPhi,gLam,gH, &
    ggPhiPhi,ggPhiLam,ggPhiH,ggLamLam,ggLamH,ggHH, &
    gggPhiPhiPhi,gggPhiPhiLam,gggPhiPhiH,gggPhiLamH,gggLamLamPhi, &
    gggLamLamLam,gggLamLamH,gggHHPhi,gggHHLam,gggHHH)
!
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*16 function qvtess(phi,flam,h)
!
!   The quadruple precision Fortran function to compute
!   the gravitational potential by the conditional split quadrature method
!
implicit integer (i-n)
implicit real*16 (a-h,o-z)
!
external qstess
!
real*16 eps
common /qepsV/eps
real*16 R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
common /qparamV/R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
real*16 phi0,flam0,h0,xiN,xiS,etaE,etaW
common /qargV/phi0,flam0,h0,xiN,xiS,etaE,etaW
real*16 alpha,gamma,zetaT,zetaB,sigma,fmu
common /qparamV2/alpha,gamma,zetaT,zetaB,sigma,fmu
!
phi0=phi
flam0=flam
h0=h
!
xiN=PhiN-phi
xiS=PhiS-phi
etaE=FlamE-flam
etaW=FlamW-flam
!
alpha=1.q0+h/R0
gamma=cos(phi)
zetaB=(HB-h)/R0
zetaT=zetaB+beta
!
!   Split quadrature
!
if(xiS.lt.0.q0.and.xiN.gt.0.q0) then
    call qqde(qstess,xiS,0.q0,eps,v1)
    call qqde(qstess,0.q0,xiN,eps,v2)
    v=v1+v2
else
    call qqde(qstess,xiS,xiN,eps,v)
endif
qvtess=v
!
return; end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*16 function qstess(xi)
!
!   The quadruple precision Fortran function to compute
!   the latitude line integral requried in "vtess"
!
implicit integer (i-n)
implicit real*16 (a-h,o-z)
!
external qvkern
!
real*16 eps
common /qepsV/eps
real*16 phi0,flam0,h0,xiN,xiS,etaE,etaW
common /qargV/phi0,flam0,h0,xiN,xiS,etaE,etaW
real*16 alpha,gamma,zetaT,zetaB,sigma,fmu
common /qparamV2/alpha,gamma,zetaT,zetaB,sigma,fmu
!
sigma=cos(phi0+xi)
fmu=sin(xi*0.5q0)
!
!   Split quadrature
!
if(etaW.lt.0.q0.and.etaE.gt.0.q0) then
    call qqde1(qvkern,etaW,0.q0,eps,s1)
    call qqde1(qvkern,0.q0,etaE,eps,s2)
    s=s1+s2
else
    call qqde1(qvkern,etaW,etaE,eps,s)
endif
qstess=s
!write (*,"(a15,1p2e15.7)") "(qstess) xi,s=",xi,s
!
return; end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*16 function qvkern(eta)
!
!   The quadruple precision Fortran function to compute
!   the kernel function expressed in a cancellation-free form
!
implicit integer (i-n)
implicit real*16 (a-h,o-z)
!
real*16 R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
common /qparamV/R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
real*16 alpha,gamma,zetaT,zetaB,sigma,fmu
common /qparamV2/alpha,gamma,zetaT,zetaB,sigma,fmu
!
fnu=sin(eta*0.5q0)
B=2.q0*alpha*(fmu*fmu+gamma*sigma*fnu*fnu)
A=B*(2.q0*alpha-B)
C=(alpha-B)*(alpha-B)-A*0.5q0
D=0.5q0*(zetaT+4.q0*alpha-3.q0*B)
ST=sqrt(A+(B+zetaT)*(B+zetaT))
SB=sqrt(A+(B+zetaB)*(B+zetaB))
T=(zetaT+zetaB+2.q0*B)/(ST+SB)
if(zetaB+B.lt.0.q0) then
    fl=qlog1p((1.q0+T)*beta*(-zetaB-B+SB)/A)
else
    fl=qlog1p((1.q0+T)*beta/(zetaB+B+SB))
endif
qvkern=sigma*(C*fl+beta*(D*T+SB*0.5q0))
!write (*,"(a15,1p2e15.7)") "(qvkern) eta,vk=",eta,qvkern
!
return; end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine qgtess(qvtess,phi,flam,h,v,gPhi,gLam,gH)
!
!   Compute the gravitational acceleration vector by the conditional switch
!   of the central and single-sided second order finite difference formulas
!
implicit integer (i-n)
implicit real*16 (a-h,o-z)
real*16 R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
common /qparamV/R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
real*16 delta,delta1,delta2,d1phi,d1lam,d1h,d2phi,d2lam,d2h,delta3,d3phi,d3lam,d3h
common /qdeltaV/delta,delta1,delta2,d1phi,d1lam,d1h,d2phi,d2lam,d2h,delta3,d3phi,d3lam,d3h
r=R0+h
c=cos(phi)
if((flam.ge.FlamW.and.flam.le.FlamE).and.(h.ge.HB.and.h.le.HT)) then
    if((phi.gt.PhiS.and.phi-d1phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d1phi.lt.PhiN)) then
        gPhi=(-qvtess(phi+2.q0*d1phi,flam,h)+4.q0*qvtess(phi+d1phi,flam,h)-3.q0*v)/(2.q0*d1phi)
    elseif((phi.lt.PhiS.and.phi+d1phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d1phi.gt.PhiN)) then
        gPhi=(qvtess(phi-2.q0*d1phi,flam,h)-4.q0*qvtess(phi-d1phi,flam,h)+3.q0*v)/(2.q0*d1phi)
    else
        gPhi=(qvtess(phi+d1phi,flam,h)-qvtess(phi-d1phi,flam,h))/(2.q0*d1phi)
    endif
else
    gPhi=(qvtess(phi+d1phi,flam,h)-qvtess(phi-d1phi,flam,h))/(2.q0*d1phi)
endif
gPhi=gPhi/r
if((phi.ge.PhiS.and.phi.le.PhiN).and.(h.ge.HB.and.h.le.HT)) then
    if((flam.gt.FlamW.and.flam-d1lam.lt.FlamW).or.(flam.gt.FlamE.and.flam-d1lam.lt.FlamE)) then
        gLam=(-qvtess(phi,flam+2.q0*d1lam,h)+4.q0*qvtess(phi,flam+d1lam,h)-3.q0*v)/(2.q0*d1lam)
    elseif((flam.lt.FlamW.and.flam+d1lam.gt.FlamW).or.(flam.lt.FlamE.and.flam+d1lam.gt.FlamE)) then
        gLam=(qvtess(phi,flam-2.q0*d1lam,h)-4.q0*qvtess(phi,flam-d1lam,h)+3.q0*v)/(2.q0*d1lam)
    else
        gLam=(qvtess(phi,flam+d1lam,h)-qvtess(phi,flam-d1lam,h))/(2.q0*d1lam)
    endif
else
    gLam=(qvtess(phi,flam+d1lam,h)-qvtess(phi,flam-d1lam,h))/(2.q0*d1lam)
endif
gLam=gLam/(r*c)
if((phi.ge.PhiS.and.phi.le.PhiN).and.(flam.ge.FlamW.and.flam.le.FlamE)) then
    if((h.gt.HB.and.h-d1h.lt.HB).or.(h.gt.HT.and.h-d1h.lt.HT)) then
        gH=(-qvtess(phi,flam,h+2.q0*d1h)+4.q0*qvtess(phi,flam,h+d1h)-3.q0*v)/(2.q0*d1h)
    elseif((h.lt.HB.and.h+d1h.gt.HB).or.(h.lt.HT.and.h+d1h.gt.HT)) then
        gH=(qvtess(phi,flam,h-2.q0*d1h)-4.q0*qvtess(phi,flam,h-d1h)+3.q0*v)/(2.q0*d1h)
    else
        gH=(qvtess(phi,flam,h+d1h)-qvtess(phi,flam,h-d1h))/(2.q0*d1h)
    endif
else
    gH=(qvtess(phi,flam,h+d1h)-qvtess(phi,flam,h-d1h))/(2.q0*d1h)
endif
return; end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine qggtess(qvtess,phi,flam,h,v,gPhi,gLam,gH, &
    ggPhiPhi,ggPhiLam,ggPhiH,ggLamLam,ggLamH,ggHH)
!
!   Compute the gravity gradient tensor by the conditional switch
!   of the central and single-sided second order finite difference formulas
!
implicit integer (i-n)
implicit real*16 (a-h,o-z)
real*16 R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
common /qparamV/R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta
real*16 delta,delta1,delta2,d1phi,d1lam,d1h,d2phi,d2lam,d2h,delta3,d3phi,d3lam,d3h
common /qdeltaV/delta,delta1,delta2,d1phi,d1lam,d1h,d2phi,d2lam,d2h,delta3,d3phi,d3lam,d3h
r=R0+h
c=cos(phi)
t=tan(phi)
if((flam.ge.FlamW.and.flam.le.FlamE).and.(h.ge.HB.and.h.le.HT)) then
    if((phi.gt.PhiS.and.phi-d2phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d2phi.lt.PhiN)) then
        ggPhiPhi=(-qvtess(phi+3.q0*d2phi,flam,h)+4.q0*qvtess(phi+2.q0*d2phi,flam,h) &
            -5.q0*qvtess(phi+d2phi,flam,h)+2.q0*v)/(d2phi*d2phi)
    elseif((phi.lt.PhiS.and.phi+d2phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d2phi.gt.PhiN)) then
        ggPhiPhi=(-qvtess(phi-3.q0*d2phi,flam,h)+4.q0*qvtess(phi-2.q0*d2phi,flam,h) &
            -5.q0*qvtess(phi-d2phi,flam,h)+2.q0*v)/(d2phi*d2phi)
    else
        ggPhiPhi=(qvtess(phi+d2phi,flam,h)-2.q0*v+qvtess(phi-d2phi,flam,h))/(d2phi*d2phi)
    endif
else
    ggPhiPhi=(qvtess(phi+d2phi,flam,h)-2.q0*v+qvtess(phi-d2phi,flam,h))/(d2phi*d2phi)
endif
ggPhiPhi=ggPhiPhi/(r*r)+gH/r
if((phi.ge.PhiS.and.phi.le.PhiN).and.(h.ge.HB.and.h.le.HT)) then
    if((flam.gt.FlamW.and.flam-d2lam.lt.FlamW).or.(flam.gt.FlamE.and.flam-d2lam.lt.FlamE)) then
        ggLamLam=(-qvtess(phi,flam+3.q0*d2lam,h)+4.q0*qvtess(phi,flam+2.q0*d2lam,h) &
            -5.q0*qvtess(phi,flam+d2lam,h)+2.q0*v)/(d2lam*d2lam)
    elseif((flam.lt.FlamW.and.flam+d2lam.gt.FlamW).or.(flam.lt.FlamE.and.flam+d2lam.gt.FlamE)) then
        ggLamLam=(-qvtess(phi,flam-3.q0*d2lam,h)+4.q0*qvtess(phi,flam-2.q0*d2lam,h) &
            -5.q0*qvtess(phi,flam-d2lam,h)+2.q0*v)/(d2lam*d2lam)
    else
        ggLamLam=(qvtess(phi,flam+d2lam,h)-2.q0*v+qvtess(phi,flam-d2lam,h))/(d2lam*d2lam)
    endif
else
    ggLamLam=(qvtess(phi,flam+d2lam,h)-2.q0*v+qvtess(phi,flam-d2lam,h))/(d2lam*d2lam)
endif
ggLamLam=ggLamLam/(r*r*c*c)+(gH-gPhi*t)/r
if((phi.ge.PhiS.and.phi.le.PhiN).and.(flam.ge.FlamW.and.flam.le.FlamE)) then
    if((h.gt.HB.and.h-d2h.lt.HB).or.(h.gt.HT.and.h-d2h.lt.HT)) then
        ggHH=(-qvtess(phi,flam,h+3.q0*d2h)+4.q0*qvtess(phi,flam,h+2.q0*d2h) &
            -5.q0*qvtess(phi,flam,h+d2h)+2.q0*v)/(d2h*d2h)
    elseif((h.lt.HB.and.h+d2h.gt.HB).or.(h.lt.HT.and.h+d2h.gt.HT)) then
        ggHH=(-qvtess(phi,flam,h-3.q0*d2h)+4.q0*qvtess(phi,flam,h-2.q0*d2h) &
            -5.q0*qvtess(phi,flam,h-d2h)+2.q0*v)/(d2h*d2h)
    else
        ggHH=(qvtess(phi,flam,h+d2h)-2.q0*v+qvtess(phi,flam,h-d2h))/(d2h*d2h)
    endif
else
    ggHH=(qvtess(phi,flam,h+d2h)-2.q0*v+qvtess(phi,flam,h-d2h))/(d2h*d2h)
endif

if(h.ge.HB.and.h.le.HT) then
    if(((phi.gt.PhiS.and.phi-d2phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d2phi.lt.PhiN)).and. &
        ((flam.gt.FlamW.and.flam-d2lam.lt.FlamW).or.(flam.gt.FlamE.and.flam-d2lam.lt.FlamE))) then
        ggPhiLam=(qvtess(phi+2.q0*d2phi,flam+2.q0*d2lam,h) &
            -4.q0*qvtess(phi+2.q0*d2phi,flam+d2lam,h)+3.q0*qvtess(phi+2.q0*d2phi,flam,h) &
            -4.q0*qvtess(phi+d2phi,flam+2.q0*d2lam,h)+16.q0*qvtess(phi+d2phi,flam+d2lam,h) &
            -12.q0*qvtess(phi+d2phi,flam,h)+3.q0*qvtess(phi,flam+2.q0*d2lam,h) &
            -12.q0*qvtess(phi,flam+d2lam,h)+9.q0*v)/(4.q0*d2phi*d2lam)
    elseif(((phi.gt.PhiS.and.phi-d2phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d2phi.lt.PhiN)).and. &
        ((flam.lt.FlamW.and.flam+d2lam.gt.FlamW).or.(flam.lt.FlamE.and.flam+d2lam.gt.FlamE))) then
        ggPhiLam=-(qvtess(phi+2.q0*d2phi,flam-2.q0*d2lam,h) &
            -4.q0*qvtess(phi+2.q0*d2phi,flam-d2lam,h)+3.q0*qvtess(phi+2.q0*d2phi,flam,h) &
            -4.q0*qvtess(phi+d2phi,flam-2.q0*d2lam,h)+16.q0*qvtess(phi+d2phi,flam-d2lam,h) &
            -12.q0*qvtess(phi+d2phi,flam,h)+3.q0*qvtess(phi,flam-2.q0*d2lam,h) &
            -12.q0*qvtess(phi,flam-d2lam,h)+9.q0*v)/(4.q0*d2phi*d2lam)
    elseif(((phi.lt.PhiS.and.phi+d2phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d2phi.gt.PhiN)).and. &
        ((flam.gt.FlamW.and.flam-d2lam.lt.FlamW).or.(flam.gt.FlamE.and.flam-d2lam.lt.FlamE))) then
        ggPhiLam=-(qvtess(phi-2.q0*d2phi,flam+2.q0*d2lam,h) &
            -4.q0*qvtess(phi-2.q0*d2phi,flam+d2lam,h)+3.q0*qvtess(phi-2.q0*d2phi,flam,h) &
            -4.q0*qvtess(phi-d2phi,flam+2.q0*d2lam,h)+16.q0*qvtess(phi-d2phi,flam+d2lam,h) &
            -12.q0*qvtess(phi-d2phi,flam,h)+3.q0*qvtess(phi,flam+2.q0*d2lam,h) &
            -12.q0*qvtess(phi,flam+d2lam,h)+9.q0*v)/(4.q0*d2phi*d2lam)
    elseif(((phi.lt.PhiS.and.phi+d2phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d2phi.gt.PhiN)).and. &
        ((flam.lt.FlamW.and.flam+d2lam.gt.FlamW).or.(flam.lt.FlamE.and.flam+d2lam.gt.FlamE))) then
        ggPhiLam=(qvtess(phi-2.q0*d2phi,flam-2.q0*d2lam,h) &
            -4.q0*qvtess(phi-2.q0*d2phi,flam-d2lam,h)+3.q0*qvtess(phi-2.q0*d2phi,flam,h) &
            -4.q0*qvtess(phi-d2phi,flam-2.q0*d2lam,h)+16.q0*qvtess(phi-d2phi,flam-d2lam,h) &
            -12.q0*qvtess(phi-d2phi,flam,h)+3.q0*qvtess(phi,flam-2.q0*d2lam,h) &
            -12.q0*qvtess(phi,flam-d2lam,h)+9.q0*v)/(4.q0*d2phi*d2lam)
    elseif(((phi.gt.PhiS.and.phi-d2phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d2phi.lt.PhiN)).and. &
        ((flam+d2lam.lt.FlamW).or.(flam-d2lam.gt.FlamE).or. &
        (flam-d2lam.gt.FlamW.and.flam+d2lam.lt.FlamE))) then
        ggPhiLam=(-qvtess(phi+2.q0*d2phi,flam+d2lam,h)+qvtess(phi+2.q0*d2phi,flam-d2lam,h) &
            +4.q0*qvtess(phi+d2phi,flam+d2lam,h)-4.q0*qvtess(phi+d2phi,flam-d2lam,h) &
            -3.q0*qvtess(phi,flam+d2lam,h)+3.q0*qvtess(phi,flam-d2lam,h))/(4.q0*d2phi*d2lam)
    elseif(((phi.lt.PhiS.and.phi+d2phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d2phi.gt.PhiN)).and. &
        ((flam+d2lam.lt.FlamW).or.(flam-d2lam.gt.FlamE).or. &
        (flam-d2lam.gt.FlamW.and.flam+d2lam.lt.FlamE))) then
        ggPhiLam=(qvtess(phi-2.q0*d2phi,flam+d2lam,h)-qvtess(phi-2.q0*d2phi,flam-d2lam,h) &
            -4.q0*qvtess(phi-d2phi,flam+d2lam,h)+4.q0*qvtess(phi-d2phi,flam-d2lam,h) &
            +3.q0*qvtess(phi,flam+d2lam,h)-3.q0*qvtess(phi,flam-d2lam,h))/(4.q0*d2phi*d2lam)
    elseif(((flam.gt.FlamW.and.flam-d2lam.lt.FlamW).or. &
        (flam.gt.FlamE.and.flam-d2lam.lt.FlamE)).and. &
        ((phi+d2phi.lt.PhiS).or.(phi-d2phi.gt.PhiN).or.(phi-d2phi.gt.PhiS.and.phi+d2phi.lt.PhiN))) then
        ggPhiLam=(-qvtess(phi+d2phi,flam+2.q0*d2lam,h)+qvtess(phi-d2phi,flam+2.q0*d2lam,h) &
            +4.q0*qvtess(phi+d2phi,flam+d2lam,h)-4.q0*qvtess(phi-d2phi,flam+d2lam,h) &
            -3.q0*qvtess(phi+d2phi,flam,h)+3.q0*qvtess(phi-d2phi,flam,h))/(4.q0*d2phi*d2lam)
    elseif(((flam.lt.FlamW.and.flam+d2lam.gt.FlamW).or. &
        (flam.lt.FlamE.and.flam+d2lam.gt.FlamE)).and. &
        ((phi+d2phi.lt.PhiS).or.(phi-d2phi.gt.PhiN).or.(phi-d2phi.gt.PhiS.and.phi+d2phi.lt.PhiN))) then
        ggPhiLam=-(-qvtess(phi+d2phi,flam-2.q0*d2lam,h)+qvtess(phi-d2phi,flam-2.q0*d2lam,h) &
            +4.q0*qvtess(phi+d2phi,flam-d2lam,h)-4.q0*qvtess(phi-d2phi,flam-d2lam,h) &
            -3.q0*qvtess(phi+d2phi,flam,h)+3.q0*qvtess(phi-d2phi,flam,h))/(4.q0*d2phi*d2lam)
    else
        ggPhiLam=(qvtess(phi+d2phi,flam+d2lam,h)-qvtess(phi-d2phi,flam+d2lam,h) &
            -qvtess(phi+d2phi,flam-d2lam,h)+qvtess(phi-d2phi,flam-d2lam,h))/(4.q0*d2phi*d2lam)
    endif
else
    ggPhiLam=(qvtess(phi+d2phi,flam+d2lam,h)-qvtess(phi-d2phi,flam+d2lam,h) &
        -qvtess(phi+d2phi,flam-d2lam,h)+qvtess(phi-d2phi,flam-d2lam,h))/(4.q0*d2phi*d2lam)
endif
ggPhiLam=ggPhiLam/(r*r*c)+gLam*t/r

if(flam.ge.FlamW.and.flam.le.FlamE) then
    if(((phi.gt.PhiS.and.phi-d2phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d2phi.lt.PhiN)).and. &
        ((h.gt.HB.and.h-d2h.lt.HB).or.(h.gt.HT.and.h-d2h.lt.HT))) then
        ggPhiH=(qvtess(phi+2.q0*d2phi,flam,h+2.q0*d2h)-4.q0*qvtess(phi+2.q0*d2phi,flam,h+d2h) &
            +3.q0*qvtess(phi+2.q0*d2phi,flam,h)-4.q0*qvtess(phi+d2phi,flam,h+2.q0*d2h) &
            +16.q0*qvtess(phi+d2phi,flam,h+d2h)-12.q0*qvtess(phi+d2phi,flam,h) &
            +3.q0*qvtess(phi,flam,h+2.q0*d2h)-12.q0*qvtess(phi,flam,h+d2h)+9.q0*v)/(4.q0*d2phi*d2h)
    elseif(((phi.gt.PhiS.and.phi-d2phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d2phi.lt.PhiN)).and. &
        ((h.lt.HB.and.h+d2h.gt.HB).or.(h.lt.HT.and.h+d2h.gt.HT))) then
        ggPhiH=-(qvtess(phi+2.q0*d2phi,flam,h-2.q0*d2h)-4.q0*qvtess(phi+2.q0*d2phi,flam,h-d2h) &
            +3.q0*qvtess(phi+2.q0*d2phi,flam,h)-4.q0*qvtess(phi+d2phi,flam,h-2.q0*d2h) &
            +16.q0*qvtess(phi+d2phi,flam,h-d2h)-12.q0*qvtess(phi+d2phi,flam,h) &
            +3.q0*qvtess(phi,flam,h-2.q0*d2h)-12.q0*qvtess(phi,flam,h-d2h)+9.q0*v)/(4.q0*d2phi*d2h)
    elseif(((phi.lt.PhiS.and.phi+d2phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d2phi.gt.PhiN)).and. &
        ((h.gt.HB.and.h-d2h.lt.HB).or.(h.gt.HT.and.h-d2h.lt.HT))) then
        ggPhiH=-(qvtess(phi-2.q0*d2phi,flam,h+2.q0*d2h)-4.q0*qvtess(phi-2.q0*d2phi,flam,h+d2h) &
            +3.q0*qvtess(phi-2.q0*d2phi,flam,h)-4.q0*qvtess(phi-d2phi,flam,h+2.q0*d2h) &
            +16.q0*qvtess(phi-d2phi,flam,h+d2h)-12.q0*qvtess(phi-d2phi,flam,h) &
            +3.q0*qvtess(phi,flam,h+2.q0*d2h)-12.q0*qvtess(phi,flam,h+d2h)+9.q0*v)/(4.q0*d2phi*d2h)
    elseif(((phi.lt.PhiS.and.phi+d2phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d2phi.gt.PhiN)).and. &
        ((h.lt.HB.and.h+d2h.gt.HB).or.(h.lt.HT.and.h+d2h.gt.HT))) then
        ggPhiH=(qvtess(phi-2.q0*d2phi,flam,h-2.q0*d2h)-4.q0*qvtess(phi-2.q0*d2phi,flam,h-d2h) &
            +3.q0*qvtess(phi-2.q0*d2phi,flam,h)-4.q0*qvtess(phi-d2phi,flam,h-2.q0*d2h) &
            +16.q0*qvtess(phi-d2phi,flam,h-d2h)-12.q0*qvtess(phi-d2phi,flam,h) &
            +3.q0*qvtess(phi,flam,h-2.q0*d2h)-12.q0*qvtess(phi,flam,h-d2h)+9.q0*v)/(4.q0*d2phi*d2h)
    elseif(((phi.gt.PhiS.and.phi-d2phi.lt.PhiS).or.(phi.gt.PhiN.and.phi-d2phi.lt.PhiN)).and. &
        ((h+d2h.lt.HB).or.(h-d2h.gt.HT).or.(h-d2h.gt.HB.and.h+d2h.lt.HT))) then
        ggPhiH=(-qvtess(phi+2.q0*d2phi,flam,h+d2h)+qvtess(phi+2.q0*d2phi,flam,h-d2h) &
            +4.q0*qvtess(phi+d2phi,flam,h+d2h)-4.q0*qvtess(phi+d2phi,flam,h-d2h) &
            -3.q0*qvtess(phi,flam,h+d2h)+3.q0*qvtess(phi,flam,h-d2h))/(4.q0*d2phi*d2h)
    elseif(((phi.lt.PhiS.and.phi+d2phi.gt.PhiS).or.(phi.lt.PhiN.and.phi+d2phi.gt.PhiN)).and. &
        ((h+d2h.lt.HB).or.(h-d2h.gt.HT).or.(h-d2h.gt.HB.and.h+d2h.lt.HT))) then
        ggPhiH=-(-qvtess(phi-2.q0*d2phi,flam,h+d2h)+qvtess(phi-2.q0*d2phi,flam,h-d2h) &
            +4.q0*qvtess(phi-d2phi,flam,h+d2h)-4.q0*qvtess(phi-d2phi,flam,h-d2h) &
            -3.q0*qvtess(phi,flam,h+d2h)+3.q0*qvtess(phi,flam,h-d2h))/(4.q0*d2phi*d2h)
    elseif(((h.gt.HB.and.h-d2h.lt.HB).or.(h.gt.HT.and.h-d2h.lt.HT)).and.((phi+d2phi.lt.PhiS) &
        .or.(phi-d2phi.gt.PhiN).or.(phi-d2phi.gt.PhiS.and.phi+d2phi.lt.PhiN))) then
        ggPhiH=(-qvtess(phi+d2phi,flam,h+2.q0*d2h)+qvtess(phi-d2phi,flam,h+2.q0*d2h) &
            +4.q0*qvtess(phi+d2phi,flam,h+d2h)-4.q0*qvtess(phi-d2phi,flam,h+d2h) &
            -3.q0*qvtess(phi+d2phi,flam,h)+3.q0*qvtess(phi-d2phi,flam,h))/(4.q0*d2phi*d2h)
    elseif(((h.lt.HB.and.h+d2h.gt.HB).or.(h.lt.HT.and.h+d2h.gt.HT)).and.((phi+d2phi.lt.PhiS) &
        .or.(phi-d2phi.gt.PhiN).or.(phi-d2phi.gt.PhiS.and.phi+d2phi.lt.PhiN))) then
        ggPhiH=-(-qvtess(phi+d2phi,flam,h-2.q0*d2h)+qvtess(phi-d2phi,flam,h-2.q0*d2h) &
            +4.q0*qvtess(phi+d2phi,flam,h-d2h)-4.q0*qvtess(phi-d2phi,flam,h-d2h) &
            -3.q0*qvtess(phi+d2phi,flam,h)+3.q0*qvtess(phi-d2phi,flam,h))/(4.q0*d2phi*d2h)
    else
        ggPhiH=(qvtess(phi+d2phi,flam,h+d2h)-qvtess(phi-d2phi,flam,h+d2h) &
            -qvtess(phi+d2phi,flam,h-d2h)+qvtess(phi-d2phi,flam,h-d2h))/(4.q0*d2phi*d2h)
    endif
else
    ggPhiH=(qvtess(phi+d2phi,flam,h+d2h)-qvtess(phi-d2phi,flam,h+d2h) &
        -qvtess(phi+d2phi,flam,h-d2h)+qvtess(phi-d2phi,flam,h-d2h))/(4.q0*d2phi*d2h)
endif
ggPhiH=ggPhiH/r-gPhi/r

if(phi.ge.PhiS.and.phi.le.PhiN) then
    if(((flam.gt.FlamW.and.flam-d2flam.lt.FlamW).or.(flam.gt.FlamE.and.flam-d2flam.lt.FlamE)).and. &
        ((h.gt.HB.and.h-d2h.lt.HB).or.(h.gt.HT.and.h-d2h.lt.HT))) then
        ggLamH=(qvtess(phi,flam+2.q0*d2lam,h+2.q0*d2h)-4.q0*qvtess(phi,flam+2.q0*d2lam,h+d2h) &
            +3.q0*qvtess(phi,flam+2.q0*d2lam,h)-4.q0*qvtess(phi,flam+d2lam,h+2.q0*d2h) &
            +16.q0*qvtess(phi,flam+d2lam,h+d2h)-12.q0*qvtess(phi,flam+d2lam,h) &
            +3.q0*qvtess(phi,flam,h+2.q0*d2h)-12.q0*qvtess(phi,flam,h+d2h)+9.q0*v)/(4.q0*d2lam*d2h)
    elseif(((flam.gt.FlamW.and.flam-d2lam.lt.FlamW).or.(flam.gt.FlamE.and.flam-d2lam.lt.FlamE)).and. &
        ((h.lt.HB.and.h+d2h.gt.HB).or.(h.lt.HT.and.h+d2h.gt.HT))) then
        ggLamH=-(qvtess(phi,flam+2.q0*d2lam,h-2.q0*d2h)-4.q0*qvtess(phi,flam+2.q0*d2lam,h-d2h) &
            +3.q0*qvtess(phi,flam+2.q0*d2lam,h)-4.q0*qvtess(phi,flam+d2lam,h-2.q0*d2h) &
            +16.q0*qvtess(phi,flam+d2lam,h-d2h)-12.q0*qvtess(phi,flam+d2lam,h) &
            +3.q0*qvtess(phi,flam,h-2.q0*d2h)-12.q0*qvtess(phi,flam,h-d2h)+9.q0*v)/(4.q0*d2lam*d2h)
    elseif(((flam.lt.FlamW.and.flam+d2lam.gt.FlamW).or.(flam.lt.FlamE.and.flam+d2lam.gt.FlamE)).and. &
        ((h.gt.HB.and.h-d2h.lt.HB).or.(h.gt.HT.and.h-d2h.lt.HT))) then
        ggLamH=-(qvtess(phi,flam-2.q0*d2lam,h+2.q0*d2h)-4.q0*qvtess(phi,flam-2.q0*d2lam,h+d2h) &
            +3.q0*qvtess(phi,flam-2.q0*d2lam,h)-4.q0*qvtess(phi,flam-d2lam,h+2.q0*d2h) &
            +16.q0*qvtess(phi,flam-d2lam,h+d2h)-12.q0*qvtess(phi,flam-d2lam,h) &
            +3.q0*qvtess(phi,flam,h+2.q0*d2h)-12.q0*qvtess(phi,flam,h+d2h)+9.q0*v)/(4.q0*d2lam*d2h)
    elseif(((flam.lt.FlamW.and.flam+d2lam.gt.FlamW).or.(flam.lt.FlamE.and.flam+d2lam.gt.FlamE)).and. &
        ((h.lt.HB.and.h+d2h.gt.HB).or.(h.lt.HT.and.h+d2h.gt.HT))) then
        ggLamH=(qvtess(phi,flam-2.q0*d2lam,h-2.q0*d2h)-4.q0*qvtess(phi,flam-2.q0*d2lam,h-d2h) &
            +3.q0*qvtess(phi,flam-2.q0*d2lam,h)-4.q0*qvtess(phi,flam-d2lam,h-2.q0*d2h) &
            +16.q0*qvtess(phi,flam-d2lam,h-d2h)-12.q0*qvtess(phi,flam-d2lam,h) &
            +3.q0*qvtess(phi,flam,h-2.q0*d2h)-12.q0*qvtess(phi,flam,h-d2h)+9.q0*v)/(4.q0*d2lam*d2h)
    elseif(((flam.gt.FlamW.and.flam-d2lam.lt.FlamW).or.(flam.gt.FlamE.and.flam-d2lam.lt.FlamE)).and. &
        ((h+d2h.lt.HB).or.(h-d2h.gt.HT).or.(h-d2h.gt.HB.and.h+d2h.lt.HT))) then
        ggLamH=(-qvtess(phi,flam+2.q0*d2lam,h+d2h)+qvtess(phi,flam+2.q0*d2lam,h-d2h) &
            +4.q0*qvtess(phi,flam+d2lam,h+d2h)-4.q0*qvtess(phi,flam+d2lam,h-d2h) &
            -3.q0*qvtess(phi,flam,h+d2h)+3.q0*qvtess(phi,flam,h-d2h))/(4.q0*d2lam*d2h)
    elseif(((flam.lt.FlamW.and.flam+d2lam.gt.FlamW).or.(flam.lt.FlamE.and.flam+d2lam.gt.FlamE)).and. &
        ((h+d2h.lt.HB).or.(h-d2h.gt.HT).or.(h-d2h.gt.HB.and.h+d2h.lt.HT))) then
        ggLamH=-(-qvtess(phi,flam-2.q0*d2lam,h+d2h)+qvtess(phi,flam-2.q0*d2lam,h-d2h) &
            +4.q0*qvtess(phi,flam-d2lam,h+d2h)-4.q0*qvtess(phi,flam-d2lam,h-d2h) &
            -3.q0*qvtess(phi,flam,h+d2h)+3.q0*qvtess(phi,flam,h-d2h))/(4.q0*d2lam*d2h)
    elseif(((h.gt.HB.and.h-d2h.lt.HB).or.(h.gt.HT.and.h-d2h.lt.HT)).and.((flam+d2lam.lt.FlamW).or. &
        (flam-d2lam.gt.FlamE).or.(flam-d2lam.gt.FlamW.and.flam+d2lam.lt.FlamE))) then
        ggLamH=(-qvtess(phi,flam+d2lam,h+2.q0*d2h)+qvtess(phi,flam-d2lam,h+2.q0*d2h) &
            +4.q0*qvtess(phi,flam+d2lam,h+d2h)-4.q0*qvtess(phi,flam-d2lam,h+d2h) &
            -3.q0*qvtess(phi,flam+d2lam,h)+3.q0*qvtess(phi,flam-d2lam,h))/(4.q0*d2lam*d2h)
    elseif(((h.lt.HB.and.h+d2h.gt.HB).or.(h.lt.HT.and.h+d2h.gt.HT)).and.((flam+d2lam.lt.FlamW).or. &
        (flam-d2lam.gt.FlamE).or.(flam-d2lam.gt.FlamW.and.flam+d2lam.lt.FlamE))) then
        ggLamH=-(-qvtess(phi,flam+d2lam,h-2.q0*d2h)+qvtess(phi,flam-d2lam,h-2.q0*d2h) &
            +4.q0*qvtess(phi,flam+d2lam,h-d2h)-4.q0*qvtess(phi,flam-d2lam,h-d2h) &
            -3.q0*qvtess(phi,flam+d2lam,h)+3.q0*qvtess(phi,flam-d2lam,h))/(4.q0*d2lam*d2h)
    else
        ggLamH=(qvtess(phi,flam+d2lam,h+d2h)-qvtess(phi,flam-d2lam,h+d2h) &
            -qvtess(phi,flam+d2lam,h-d2h)+qvtess(phi,flam-d2lam,h-d2h))/(4.q0*d2lam*d2h)
    endif
else
    ggLamH=(qvtess(phi,flam+d2lam,h+d2h)-qvtess(phi,flam-d2lam,h+d2h) &
        -qvtess(phi,flam+d2lam,h-d2h)+qvtess(phi,flam-d2lam,h-d2h))/(4.q0*d2lam*d2h)
endif
ggLamH=ggLamH/(r*c)-gLam/r
return; end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine qgggtess(qvtess,phi,flam,h,v,gPhi,gLam,gH, &
    ggPhiPhi,ggPhiLam,ggPhiH,ggLamLam,ggLamH,ggHH,&
    gggPhiPhiPhi,gggPhiPhiLam,gggPhiPhiH,gggPhiLamH,gggLamLamPhi, &
    gggLamLamLam,gggLamLamH,gggHHPhi,gggHHLam,gggHHH)
!
!   Compute the gravitational curvatures by the conditional switch
!   of the third-order central and single-sided finite difference formulas
!
implicit integer (i-n)
implicit real*16 (a-h,o-z)
!
!    Common blocks for passing the parameters and indirect variables
!
real*16 R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta,r,s,c,t
common /qparamV/R0,HT,HB,PhiN,PhiS,FlamE,FlamW,beta,r,s,c,t
real*16 delta,delta1,delta2,d1phi,d1lam,d1h,d2phi,d2lam,d2h,delta3,d3phi,d3lam,d3h
common /qdeltaV/delta,delta1,delta2,d1phi,d1lam,d1h,d2phi,d2lam,d2h,delta3,d3phi,d3lam,d3h
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
        gggPhiPhiPhi=(-3.q0*qvtess(phi+4.q0*d3phi,flam,h)+14.q0*qvtess(phi+3.q0*d3phi,flam,h) &
            -24.q0*qvtess(phi+2.q0*d3phi,flam,h)+18.q0*qvtess(phi+d3phi,flam,h) &
            -5.q0*v)/(2.q0*d3phi*d3phi*d3phi)
    elseif((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) then
        gggPhiPhiPhi=(3.q0*qvtess(phi-4.q0*d3phi,flam,h)-14.q0*qvtess(phi-3.q0*d3phi,flam,h) &
            +24.q0*qvtess(phi-2.q0*d3phi,flam,h)-18.q0*qvtess(phi-d3phi,flam,h) &
            +5.q0*v)/(2.q0*d3phi*d3phi*d3phi)
    else
        gggPhiPhiPhi=(qvtess(phi+2.q0*d3phi,flam,h)-2.q0*qvtess(phi+d3phi,flam,h) &
            +2.q0*qvtess(phi-d3phi,flam,h)-qvtess(phi-2.q0*d3phi,flam,h) &
            )/(2.q0*d3phi*d3phi*d3phi)
    endif
else
    gggPhiPhiPhi=(qvtess(phi+2.q0*d3phi,flam,h)-2.q0*qvtess(phi+d3phi,flam,h) &
            +2.q0*qvtess(phi-d3phi,flam,h)-qvtess(phi-2.q0*d3phi,flam,h) &
            )/(2.q0*d3phi*d3phi*d3phi)
endif
gggPhiPhiPhi=(gggPhiPhiPhi+r*(gPhi + 3.q0*r*ggPhiH))/(r*r*r)
!
!   Compute gggPhiPhiLam
!
if(h.ge.HB.and.h.le.HT) then
! case 1
    if(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE))) then
        gggPhiPhiLam=(qvtess(phi+3.q0*d3phi,flam+2.q0*d3lam,h)-4.q0*qvtess(phi+3.q0*d3phi,flam+d3lam,h) &
            +3.q0*qvtess(phi+3.q0*d3phi,flam,h)-4.q0*qvtess(phi+2.q0*d3phi,flam+2.q0*d3lam,h) &
            +16.q0*qvtess(phi+2.q0*d3phi,flam+d3lam,h)-12.q0*qvtess(phi+2.q0*d3phi,flam,h) &
            +5.q0*qvtess(phi+d3phi,flam+2.q0*d3lam,h)-20.q0*qvtess(phi+d3phi,flam+d3lam,h) &
            +15.q0*qvtess(phi+d3phi,flam,h)-2.q0*qvtess(phi,flam+2.q0*d3lam,h) &
            +8.q0*qvtess(phi,flam+d3lam,h)-6.q0*v &
            )/(2.q0*d3phi*d3phi*d3lam)
! case 2
    elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE))) then
        gggPhiPhiLam=-(qvtess(phi+3.q0*d3phi,flam-2.q0*d3lam,h)-4.q0*qvtess(phi+3.q0*d3phi,flam-d3lam,h) &
            +3.q0*qvtess(phi+3.q0*d3phi,flam,h)-4.q0*qvtess(phi+2.q0*d3phi,flam-2.q0*d3lam,h) &
            +16.q0*qvtess(phi+2.q0*d3phi,flam-d3lam,h)-12.q0*qvtess(phi+2.q0*d3phi,flam,h) &
            +5.q0*qvtess(phi+d3phi,flam-2.q0*d3lam,h)-20.q0*qvtess(phi+d3phi,flam-d3lam,h) &
            +15.q0*qvtess(phi+d3phi,flam,h)-2.q0*qvtess(phi,flam-2.q0*d3lam,h) &
            +8.q0*qvtess(phi,flam-d3lam,h)-6.q0*v &
            )/(2.q0*d3phi*d3phi*d3lam)
! case 3
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE))) then
        gggPhiPhiLam=(qvtess(phi-3.q0*d3phi,flam+2.q0*d3lam,h)-4.q0*qvtess(phi-3.q0*d3phi,flam+d3lam,h) &
            +3.q0*qvtess(phi-3.q0*d3phi,flam,h)-4.q0*qvtess(phi-2.q0*d3phi,flam+2.q0*d3lam,h) &
            +16.q0*qvtess(phi-2.q0*d3phi,flam+d3lam,h)-12.q0*qvtess(phi-2.q0*d3phi,flam,h) &
            +5.q0*qvtess(phi-d3phi,flam+2.q0*d3lam,h)-20.q0*qvtess(phi-d3phi,flam+d3lam,h) &
            +15.q0*qvtess(phi-d3phi,flam,h)-2.q0*qvtess(phi,flam+2.q0*d3lam,h) &
            +8.q0*qvtess(phi,flam+d3lam,h)-6.q0*v &
            )/(2.q0*d3phi*d3phi*d3lam)
! case 4
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE))) then
        gggPhiPhiLam=-(qvtess(phi-3.q0*d3phi,flam-2.q0*d3lam,h)-4.q0*qvtess(phi-3.q0*d3phi,flam-d3lam,h) &
            +3.q0*qvtess(phi-3.q0*d3phi,flam,h)-4.q0*qvtess(phi-2.q0*d3phi,flam-2.q0*d3lam,h) &
            +16.q0*qvtess(phi-2.q0*d3phi,flam-d3lam,h)-12.q0*qvtess(phi-2.q0*d3phi,flam,h) &
            +5.q0*qvtess(phi-d3phi,flam-2.q0*d3lam,h)-20.q0*qvtess(phi-d3phi,flam-d3lam,h) &
            +15.q0*qvtess(phi-d3phi,flam,h)-2.q0*qvtess(phi,flam-2.q0*d3lam,h) &
            +8.q0*qvtess(phi,flam-d3lam,h)-6.q0*v &
            )/(2.q0*d3phi*d3phi*d3lam)
! case 5
    elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE))) then
        gggPhiPhiLam=(-qvtess(phi+3.q0*d3phi,flam+d3lam,h)+qvtess(phi+3.q0*d3phi,flam-d3lam,h) &
            +4.q0*qvtess(phi+2.q0*d3phi,flam+d3lam,h)-4.q0*qvtess(phi+2.q0*d3phi,flam-d3lam,h) &
            -5.q0*qvtess(phi+d3phi,flam+d3lam,h)+5.q0*qvtess(phi+d3phi,flam-d3lam,h) &
            +2.q0*qvtess(phi,flam+d3lam,h)-2.q0*qvtess(phi,flam-d3lam,h) &
            )/(2.q0*d3phi*d3phi*d3lam)
! case 6
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE))) then
        gggPhiPhiLam=(-qvtess(phi-3.q0*d3phi,flam+d3lam,h)+qvtess(phi-3.q0*d3phi,flam-d3lam,h) &
            +4.q0*qvtess(phi-2.q0*d3phi,flam+d3lam,h)-4.q0*qvtess(phi-2.q0*d3phi,flam-d3lam,h) &
            -5.q0*qvtess(phi-d3phi,flam+d3lam,h)+5.q0*qvtess(phi-d3phi,flam-d3lam,h) &
            +2.q0*qvtess(phi,flam+d3lam,h)-2.q0*qvtess(phi,flam-d3lam,h) &
            )/(2.q0*d3phi*d3phi*d3lam)
! case 7
    elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE))) then
        gggPhiPhiLam=(-qvtess(phi+d3phi,flam+2.q0*d3lam,h)+4.q0*qvtess(phi+d3phi,flam+d3lam,h) &
            -3.q0*qvtess(phi+d3phi,flam,h)+2.q0*qvtess(phi,flam+2.q0*d3lam,h) &
            -8.q0*qvtess(phi,flam+d3lam,h)+6.q0*v-qvtess(phi-d3phi,flam+2.q0*d3lam,h) &
            +4.q0*qvtess(phi-d3phi,flam+d3lam,h)-3.q0*qvtess(phi-d3phi,flam,h) &
            )/(2.q0*d3phi*d3phi*d3lam)
! case 8
    elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE))) then
        gggPhiPhiLam=-(-qvtess(phi+d3phi,flam-2.q0*d3lam,h)+4.q0*qvtess(phi+d3phi,flam-d3lam,h) &
            -3.q0*qvtess(phi+d3phi,flam,h)+2.q0*qvtess(phi,flam-2.q0*d3lam,h) &
            -8.q0*qvtess(phi,flam-d3lam,h)+6.q0*v-qvtess(phi-d3phi,flam-2.q0*d3lam,h) &
            +4.q0*qvtess(phi-d3phi,flam-d3lam,h)-3.q0*qvtess(phi-d3phi,flam,h) &
            )/(2.q0*d3phi*d3phi*d3lam)
! case 9
    else
        gggPhiPhiLam=(qvtess(phi+d3phi,flam+d3lam,h)-qvtess(phi+d3phi,flam-d3lam,h) &
            -2.q0*qvtess(phi,flam+d3lam,h)+2.q0*qvtess(phi,flam-d3lam,h) &
            +qvtess(phi-d3phi,flam+d3lam,h)-qvtess(phi-d3phi,flam-d3lam,h) &
            )/(2.q0*d3phi*d3phi*d3lam)
    endif
else
    gggPhiPhiLam=(qvtess(phi+d3phi,flam+d3lam,h)-qvtess(phi+d3phi,flam-d3lam,h) &
            -2.q0*qvtess(phi,flam+d3lam,h)+2.q0*qvtess(phi,flam-d3lam,h) &
            +qvtess(phi-d3phi,flam+d3lam,h)-qvtess(phi-d3phi,flam-d3lam,h) &
            )/(2.q0*d3phi*d3phi*d3lam)
endif
gggPhiPhiLam=gggPhiPhiLam/(r*r*r*c)+(gLam+r*ggLamH+2.q0*r*t*ggPhiLam)/(r*r)
!
!   Compute gggPhiPhiH
!
if(flam.ge.FlamW.and.flam.le.FlamE) then
! case 1
    if(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggPhiPhiH=(qvtess(phi+3.q0*d3phi,flam,h+2.q0*d3h)-4.q0*qvtess(phi+3.q0*d3phi,flam,h+d3h) &
            +3.q0*qvtess(phi+3.q0*d3phi,flam,h)-4.q0*qvtess(phi+2.q0*d3phi,flam,h+2.q0*d3h) &
            +16.q0*qvtess(phi+2.q0*d3phi,flam,h+d3h)-12.q0*qvtess(phi+2.q0*d3phi,flam,h) &
            +5.q0*qvtess(phi+d3phi,flam,h+2.q0*d3h)-20.q0*qvtess(phi+d3phi,flam,h+d3h) &
            +15.q0*qvtess(phi+d3phi,flam,h)-2.q0*qvtess(phi,flam,h+2.q0*d3h) &
            +8.q0*qvtess(phi,flam,h+d3h)-6.q0*v &
            )/(2.q0*d3phi*d3phi*d3h)
! case 2
    elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggPhiPhiH=-(qvtess(phi+3.q0*d3phi,flam,h-2.q0*d3h)-4.q0*qvtess(phi+3.q0*d3phi,flam,h-d3h) &
            +3.q0*qvtess(phi+3.q0*d3phi,flam,h)-4.q0*qvtess(phi+2.q0*d3phi,flam,h-2.q0*d3h) &
            +16.q0*qvtess(phi+2.q0*d3phi,flam,h-d3h)-12.q0*qvtess(phi+2.q0*d3phi,flam,h) &
            +5.q0*qvtess(phi+d3phi,flam,h-2.q0*d3h)-20.q0*qvtess(phi+d3phi,flam,h-d3h) &
            +15.q0*qvtess(phi+d3phi,flam,h)-2.q0*qvtess(phi,flam,h-2.q0*d3h) &
            +8.q0*qvtess(phi,flam,h-d3h)-6.q0*v &
            )/(2.q0*d3phi*d3phi*d3h)
! case 3
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggPhiPhiH=(qvtess(phi-3.q0*d3phi,flam,h+2.q0*d3h)-4.q0*qvtess(phi-3.q0*d3phi,flam,h+d3h) &
            +3.q0*qvtess(phi-3.q0*d3phi,flam,h)-4.q0*qvtess(phi-2.q0*d3phi,flam,h+2.q0*d3h) &
            +16.q0*qvtess(phi-2.q0*d3phi,flam,h+d3h)-12.q0*qvtess(phi-2.q0*d3phi,flam,h) &
            +5.q0*qvtess(phi-d3phi,flam,h+2.q0*d3h)-20.q0*qvtess(phi-d3phi,flam,h+d3h) &
            +15.q0*qvtess(phi-d3phi,flam,h)-2.q0*qvtess(phi,flam,h+2.q0*d3h) &
            +8.q0*qvtess(phi,flam,h+d3h)-6.q0*v &
            )/(2.q0*d3phi*d3phi*d3h)
! case 4
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggPhiPhiH=-(qvtess(phi-3.q0*d3phi,flam,h-2.q0*d3h)-4.q0*qvtess(phi-3.q0*d3phi,flam,h-d3h) &
            +3.q0*qvtess(phi-3.q0*d3phi,flam,h)-4.q0*qvtess(phi-2.q0*d3phi,flam,h-2.q0*d3h) &
            +16.q0*qvtess(phi-2.q0*d3phi,flam,h-d3h)-12.q0*qvtess(phi-2.q0*d3phi,flam,h) &
            +5.q0*qvtess(phi-d3phi,flam,h-2.q0*d3h)-20.q0*qvtess(phi-d3phi,flam,h-d3h) &
            +15.q0*qvtess(phi-d3phi,flam,h)-2.q0*qvtess(phi,flam,h-2.q0*d3h) &
            +8.q0*qvtess(phi,flam,h-d3h)-6.q0*v &
            )/(2.q0*d3phi*d3phi*d3h)
! case 5
    elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
        gggPhiPhiH=(-qvtess(phi+3.q0*d3phi,flam,h+d3h)+qvtess(phi+3.q0*d3phi,flam,h-d3h) &
            +4.q0*qvtess(phi+2.q0*d3phi,flam,h+d3h)-4.q0*qvtess(phi+2.q0*d3phi,flam,h-d3h) &
            -5.q0*qvtess(phi+d3phi,flam,h+d3h)+5.q0*qvtess(phi+d3phi,flam,h-d3h) &
            +2.q0*qvtess(phi,flam,h+d3h)-2.q0*qvtess(phi,flam,h-d3h) &
            )/(2.q0*d3phi*d3phi*d3h)
! case 6
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
        gggPhiPhiH=(-qvtess(phi-3.q0*d3phi,flam,h+d3h)+qvtess(phi-3.q0*d3phi,flam,h-d3h) &
            +4.q0*qvtess(phi-2.q0*d3phi,flam,h+d3h)-4.q0*qvtess(phi-2.q0*d3phi,flam,h-d3h) &
            -5.q0*qvtess(phi-d3phi,flam,h+d3h)+5.q0*qvtess(phi-d3phi,flam,h-d3h) &
            +2.q0*qvtess(phi,flam,h+d3h)-2.q0*qvtess(phi,flam,h-d3h) &
            )/(2.q0*d3phi*d3phi*d3h)
! case 7
    elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggPhiPhiH=(-qvtess(phi+d3phi,flam,h+2.q0*d3h)+4.q0*qvtess(phi+d3phi,flam,h+d3h) &
            -3.q0*qvtess(phi+d3phi,flam,h)+2.q0*qvtess(phi,flam,h+2.q0*d3h) &
            -8.q0*qvtess(phi,flam,h+d3h)+6.q0*v -qvtess(phi-d3phi,flam,h+2.q0*d3h) &
            +4.q0*qvtess(phi-d3phi,flam,h+d3h)-3.q0*qvtess(phi-d3phi,flam,h) &
            )/(2.q0*d3phi*d3phi*d3h)
! case 8
    elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggPhiPhiH=-(-qvtess(phi+d3phi,flam,h-2.q0*d3h)+4.q0*qvtess(phi+d3phi,flam,h-d3h) &
            -3.q0*qvtess(phi+d3phi,flam,h)+2.q0*qvtess(phi,flam,h-2.q0*d3h) &
            -8.q0*qvtess(phi,flam,h-d3h)+6.q0*v -qvtess(phi-d3phi,flam,h-2.q0*d3h) &
            +4.q0*qvtess(phi-d3phi,flam,h-d3h)-3.q0*qvtess(phi-d3phi,flam,h) &
            )/(2.q0*d3phi*d3phi*d3h)
! case 9
    else
        gggPhiPhiH=(qvtess(phi+d3phi,flam,h+d3h)-qvtess(phi+d3phi,flam,h-d3h) &
            -2.q0*qvtess(phi,flam,h+d3h)+2.q0*qvtess(phi,flam,h-d3h) &
            +qvtess(phi-d3phi,flam,h+d3h)-qvtess(phi-d3phi,flam,h-d3h) &
            )/(2.q0*d3phi*d3phi*d3h)
    endif
else
    gggPhiPhiH=(qvtess(phi+d3phi,flam,h+d3h)-qvtess(phi+d3phi,flam,h-d3h) &
            -2.q0*qvtess(phi,flam,h+d3h)+2.q0*qvtess(phi,flam,h-d3h) &
            +qvtess(phi-d3phi,flam,h+d3h)-qvtess(phi-d3phi,flam,h-d3h) &
            )/(2.q0*d3phi*d3phi*d3h)
endif
gggPhiPhiH=(gggPhiPhiH + gH + r*ggHH - 2.q0*r*ggPhiPhi)/(r*r)
!
!   Compute gggPhiLamH
!
! case 1
if(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
    gggPhiLamH=(-qvtess(phi+2.q0*d3phi,flam+2.q0*d3lam,h+2.q0*d3h) &
        +4.q0*qvtess(phi+2.q0*d3phi,flam+2.q0*d3lam,h+d3h)-3.q0*qvtess(phi+2.q0*d3phi,flam+2.q0*d3lam,h) &
        +4.q0*qvtess(phi+2.q0*d3phi,flam+d3lam,h+2.q0*d3h)-16.q0*qvtess(phi+2.q0*d3phi,flam+d3lam,h+d3h) &
        +12.q0*qvtess(phi+2.q0*d3phi,flam+d3lam,h)-3.q0*qvtess(phi+2.q0*d3phi,flam,h+2.q0*d3h) &
        +12.q0*qvtess(phi+2.q0*d3phi,flam,h+d3h)-9.q0*qvtess(phi+2.q0*d3phi,flam,h) &
        +4.q0*qvtess(phi+d3phi,flam+2.q0*d3lam,h+2.q0*d3h)-16.q0*qvtess(phi+d3phi,flam+2.q0*d3lam,h+d3h) &
        +12.q0*qvtess(phi+d3phi,flam+2.q0*d3lam,h)-16.q0*qvtess(phi+d3phi,flam+d3lam,h+2.q0*d3h) &
        +64.q0*qvtess(phi+d3phi,flam+d3lam,h+d3h)-48.q0*qvtess(phi+d3phi,flam+d3lam,h) &
        +12.q0*qvtess(phi+d3phi,flam,h+2.q0*d3h)-48.q0*qvtess(phi+d3phi,flam,h+d3h) &
        +36.q0*qvtess(phi+d3phi,flam,h)-3.q0*qvtess(phi,flam+2.q0*d3lam,h+2.q0*d3h) &
        +12.q0*qvtess(phi,flam+2.q0*d3lam,h+d3h)-9.q0*qvtess(phi,flam+2.q0*d3lam,h) &
        +12.q0*qvtess(phi,flam+d3lam,h+2.q0*d3h)-48.q0*qvtess(phi,flam+d3lam,h+d3h) &
        +36.q0*qvtess(phi,flam+d3lam,h)-9.q0*qvtess(phi,flam,h+2.q0*d3h) &
        +36.q0*qvtess(phi,flam,h+d3h)-27.q0*v &
        )/(8.q0*d3phi*d3lam*d3h)
! case 2
elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
    gggPhiLamH=-(-qvtess(phi+2.q0*d3phi,flam+2.q0*d3lam,h-2.q0*d3h) &
        +4.q0*qvtess(phi+2.q0*d3phi,flam+2.q0*d3lam,h-d3h)-3.q0*qvtess(phi+2.q0*d3phi,flam+2.q0*d3lam,h) &
        +4.q0*qvtess(phi+2.q0*d3phi,flam+d3lam,h-2.q0*d3h)-16.q0*qvtess(phi+2.q0*d3phi,flam+d3lam,h-d3h) &
        +12.q0*qvtess(phi+2.q0*d3phi,flam+d3lam,h)-3.q0*qvtess(phi+2.q0*d3phi,flam,h-2.q0*d3h) &
        +12.q0*qvtess(phi+2.q0*d3phi,flam,h-d3h)-9.q0*qvtess(phi+2.q0*d3phi,flam,h) &
        +4.q0*qvtess(phi+d3phi,flam+2.q0*d3lam,h-2.q0*d3h)-16.q0*qvtess(phi+d3phi,flam+2.q0*d3lam,h-d3h) &
        +12.q0*qvtess(phi+d3phi,flam+2.q0*d3lam,h)-16.q0*qvtess(phi+d3phi,flam+d3lam,h-2.q0*d3h) &
        +64.q0*qvtess(phi+d3phi,flam+d3lam,h-d3h)-48.q0*qvtess(phi+d3phi,flam+d3lam,h) &
        +12.q0*qvtess(phi+d3phi,flam,h-2.q0*d3h)-48.q0*qvtess(phi+d3phi,flam,h-d3h) &
        +36.q0*qvtess(phi+d3phi,flam,h)-3.q0*qvtess(phi,flam+2.q0*d3lam,h-2.q0*d3h) &
        +12.q0*qvtess(phi,flam+2.q0*d3lam,h-d3h)-9.q0*qvtess(phi,flam+2.q0*d3lam,h) &
        +12.q0*qvtess(phi,flam+d3lam,h-2.q0*d3h)-48.q0*qvtess(phi,flam+d3lam,h-d3h) &
        +36.q0*qvtess(phi,flam+d3lam,h)-9.q0*qvtess(phi,flam,h-2.q0*d3h) &
        +36.q0*qvtess(phi,flam,h-d3h)-27.q0*v &
        )/(8.q0*d3phi*d3lam*d3h)
! case 3
elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
    gggPhiLamH=-(-qvtess(phi+2.q0*d3phi,flam-2.q0*d3lam,h+2.q0*d3h) &
        +4.q0*qvtess(phi+2.q0*d3phi,flam-2.q0*d3lam,h+d3h)-3.q0*qvtess(phi+2.q0*d3phi,flam-2.q0*d3lam,h) &
        +4.q0*qvtess(phi+2.q0*d3phi,flam-d3lam,h+2.q0*d3h)-16.q0*qvtess(phi+2.q0*d3phi,flam-d3lam,h+d3h) &
        +12.q0*qvtess(phi+2.q0*d3phi,flam-d3lam,h)-3.q0*qvtess(phi+2.q0*d3phi,flam,h+2.q0*d3h) &
        +12.q0*qvtess(phi+2.q0*d3phi,flam,h+d3h)-9.q0*qvtess(phi+2.q0*d3phi,flam,h) &
        +4.q0*qvtess(phi+d3phi,flam-2.q0*d3lam,h+2.q0*d3h)-16.q0*qvtess(phi+d3phi,flam-2.q0*d3lam,h+d3h) &
        +12.q0*qvtess(phi+d3phi,flam-2.q0*d3lam,h)-16.q0*qvtess(phi+d3phi,flam-d3lam,h+2.q0*d3h) &
        +64.q0*qvtess(phi+d3phi,flam-d3lam,h+d3h)-48.q0*qvtess(phi+d3phi,flam-d3lam,h) &
        +12.q0*qvtess(phi+d3phi,flam,h+2.q0*d3h)-48.q0*qvtess(phi+d3phi,flam,h+d3h) &
        +36.q0*qvtess(phi+d3phi,flam,h)-3.q0*qvtess(phi,flam-2.q0*d3lam,h+2.q0*d3h) &
        +12.q0*qvtess(phi,flam-2.q0*d3lam,h+d3h)-9.q0*qvtess(phi,flam-2.q0*d3lam,h) &
        +12.q0*qvtess(phi,flam-d3lam,h+2.q0*d3h)-48.q0*qvtess(phi,flam-d3lam,h+d3h) &
        +36.q0*qvtess(phi,flam-d3lam,h)-9.q0*qvtess(phi,flam,h+2.q0*d3h) &
        +36.q0*qvtess(phi,flam,h+d3h)-27.q0*v &
        )/(8.q0*d3phi*d3lam*d3h)
! case 4
elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
    gggPhiLamH=(-qvtess(phi+2.q0*d3phi,flam-2.q0*d3lam,h-2.q0*d3h) &
        +4.q0*qvtess(phi+2.q0*d3phi,flam-2.q0*d3lam,h-d3h)-3.q0*qvtess(phi+2.q0*d3phi,flam-2.q0*d3lam,h) &
        +4.q0*qvtess(phi+2.q0*d3phi,flam-d3lam,h-2.q0*d3h)-16.q0*qvtess(phi+2.q0*d3phi,flam-d3lam,h-d3h) &
        +12.q0*qvtess(phi+2.q0*d3phi,flam-d3lam,h)-3.q0*qvtess(phi+2.q0*d3phi,flam,h-2.q0*d3h) &
        +12.q0*qvtess(phi+2.q0*d3phi,flam,h-d3h)-9.q0*qvtess(phi+2.q0*d3phi,flam,h) &
        +4.q0*qvtess(phi+d3phi,flam-2.q0*d3lam,h-2.q0*d3h)-16.q0*qvtess(phi+d3phi,flam-2.q0*d3lam,h-d3h) &
        +12.q0*qvtess(phi+d3phi,flam-2.q0*d3lam,h)-16.q0*qvtess(phi+d3phi,flam-d3lam,h-2.q0*d3h) &
        +64.q0*qvtess(phi+d3phi,flam-d3lam,h-d3h)-48.q0*qvtess(phi+d3phi,flam-d3lam,h) &
        +12.q0*qvtess(phi+d3phi,flam,h-2.q0*d3h)-48.q0*qvtess(phi+d3phi,flam,h-d3h) &
        +36.q0*qvtess(phi+d3phi,flam,h)-3.q0*qvtess(phi,flam-2.q0*d3lam,h-2.q0*d3h) &
        +12.q0*qvtess(phi,flam-2.q0*d3lam,h-d3h)-9.q0*qvtess(phi,flam-2.q0*d3lam,h) &
        +12.q0*qvtess(phi,flam-d3lam,h-2.q0*d3h)-48.q0*qvtess(phi,flam-d3lam,h-d3h) &
        +36.q0*qvtess(phi,flam-d3lam,h)-9.q0*qvtess(phi,flam,h-2.q0*d3h) &
        +36.q0*qvtess(phi,flam,h-d3h)-27.q0*v &
        )/(8.q0*d3phi*d3lam*d3h)
! case 5
elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
    gggPhiLamH=-(-qvtess(phi-2.q0*d3phi,flam+2.q0*d3lam,h+2.q0*d3h) &
        +4.q0*qvtess(phi-2.q0*d3phi,flam+2.q0*d3lam,h+d3h)-3.q0*qvtess(phi-2.q0*d3phi,flam+2.q0*d3lam,h) &
        +4.q0*qvtess(phi-2.q0*d3phi,flam+d3lam,h+2.q0*d3h)-16.q0*qvtess(phi-2.q0*d3phi,flam+d3lam,h+d3h) &
        +12.q0*qvtess(phi-2.q0*d3phi,flam+d3lam,h)-3.q0*qvtess(phi-2.q0*d3phi,flam,h+2.q0*d3h) &
        +12.q0*qvtess(phi-2.q0*d3phi,flam,h+d3h)-9.q0*qvtess(phi-2.q0*d3phi,flam,h) &
        +4.q0*qvtess(phi-d3phi,flam+2.q0*d3lam,h+2.q0*d3h)-16.q0*qvtess(phi-d3phi,flam+2.q0*d3lam,h+d3h) &
        +12.q0*qvtess(phi-d3phi,flam+2.q0*d3lam,h)-16.q0*qvtess(phi-d3phi,flam+d3lam,h+2.q0*d3h) &
        +64.q0*qvtess(phi-d3phi,flam+d3lam,h+d3h)-48.q0*qvtess(phi-d3phi,flam+d3lam,h) &
        +12.q0*qvtess(phi-d3phi,flam,h+2.q0*d3h)-48.q0*qvtess(phi-d3phi,flam,h+d3h) &
        +36.q0*qvtess(phi-d3phi,flam,h)-3.q0*qvtess(phi,flam+2.q0*d3lam,h+2.q0*d3h) &
        +12.q0*qvtess(phi,flam+2.q0*d3lam,h+d3h)-9.q0*qvtess(phi,flam+2.q0*d3lam,h) &
        +12.q0*qvtess(phi,flam+d3lam,h+2.q0*d3h)-48.q0*qvtess(phi,flam+d3lam,h+d3h) &
        +36.q0*qvtess(phi,flam+d3lam,h)-9.q0*qvtess(phi,flam,h+2.q0*d3h) &
        +36.q0*qvtess(phi,flam,h+d3h)-27.q0*v &
        )/(8.q0*d3phi*d3lam*d3h)
! case 6
elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
    gggPhiLamH=(-qvtess(phi-2.q0*d3phi,flam+2.q0*d3lam,h-2.q0*d3h) &
        +4.q0*qvtess(phi-2.q0*d3phi,flam+2.q0*d3lam,h-d3h)-3.q0*qvtess(phi-2.q0*d3phi,flam+2.q0*d3lam,h) &
        +4.q0*qvtess(phi-2.q0*d3phi,flam+d3lam,h-2.q0*d3h)-16.q0*qvtess(phi-2.q0*d3phi,flam+d3lam,h-d3h) &
        +12.q0*qvtess(phi-2.q0*d3phi,flam+d3lam,h)-3.q0*qvtess(phi-2.q0*d3phi,flam,h-2.q0*d3h) &
        +12.q0*qvtess(phi-2.q0*d3phi,flam,h-d3h)-9.q0*qvtess(phi-2.q0*d3phi,flam,h) &
        +4.q0*qvtess(phi-d3phi,flam+2.q0*d3lam,h-2.q0*d3h)-16.q0*qvtess(phi-d3phi,flam+2.q0*d3lam,h-d3h) &
        +12.q0*qvtess(phi-d3phi,flam+2.q0*d3lam,h)-16.q0*qvtess(phi-d3phi,flam+d3lam,h-2.q0*d3h) &
        +64.q0*qvtess(phi-d3phi,flam+d3lam,h-d3h)-48.q0*qvtess(phi-d3phi,flam+d3lam,h) &
        +12.q0*qvtess(phi-d3phi,flam,h-2.q0*d3h)-48.q0*qvtess(phi-d3phi,flam,h-d3h) &
        +36.q0*qvtess(phi-d3phi,flam,h)-3.q0*qvtess(phi,flam+2.q0*d3lam,h-2.q0*d3h) &
        +12.q0*qvtess(phi,flam+2.q0*d3lam,h-d3h)-9.q0*qvtess(phi,flam+2.q0*d3lam,h) &
        +12.q0*qvtess(phi,flam+d3lam,h-2.q0*d3h)-48.q0*qvtess(phi,flam+d3lam,h-d3h) &
        +36.q0*qvtess(phi,flam+d3lam,h)-9.q0*qvtess(phi,flam,h-2.q0*d3h) &
        +36.q0*qvtess(phi,flam,h-d3h)-27.q0*v &
        )/(8.q0*d3phi*d3lam*d3h)
! case 7
elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
    gggPhiLamH=(-qvtess(phi-2.q0*d3phi,flam-2.q0*d3lam,h+2.q0*d3h) &
        +4.q0*qvtess(phi-2.q0*d3phi,flam-2.q0*d3lam,h+d3h)-3.q0*qvtess(phi-2.q0*d3phi,flam-2.q0*d3lam,h) &
        +4.q0*qvtess(phi-2.q0*d3phi,flam-d3lam,h+2.q0*d3h)-16.q0*qvtess(phi-2.q0*d3phi,flam-d3lam,h+d3h) &
        +12.q0*qvtess(phi-2.q0*d3phi,flam-d3lam,h)-3.q0*qvtess(phi-2.q0*d3phi,flam,h+2.q0*d3h) &
        +12.q0*qvtess(phi-2.q0*d3phi,flam,h+d3h)-9.q0*qvtess(phi-2.q0*d3phi,flam,h) &
        +4.q0*qvtess(phi-d3phi,flam-2.q0*d3lam,h+2.q0*d3h)-16.q0*qvtess(phi-d3phi,flam-2.q0*d3lam,h+d3h) &
        +12.q0*qvtess(phi-d3phi,flam-2.q0*d3lam,h)-16.q0*qvtess(phi-d3phi,flam-d3lam,h+2.q0*d3h) &
        +64.q0*qvtess(phi-d3phi,flam-d3lam,h+d3h)-48.q0*qvtess(phi-d3phi,flam-d3lam,h) &
        +12.q0*qvtess(phi-d3phi,flam,h+2.q0*d3h)-48.q0*qvtess(phi-d3phi,flam,h+d3h) &
        +36.q0*qvtess(phi-d3phi,flam,h)-3.q0*qvtess(phi,flam-2.q0*d3lam,h+2.q0*d3h) &
        +12.q0*qvtess(phi,flam-2.q0*d3lam,h+d3h)-9.q0*qvtess(phi,flam-2.q0*d3lam,h) &
        +12.q0*qvtess(phi,flam-d3lam,h+2.q0*d3h)-48.q0*qvtess(phi,flam-d3lam,h+d3h) &
        +36.q0*qvtess(phi,flam-d3lam,h)-9.q0*qvtess(phi,flam,h+2.q0*d3h) &
        +36.q0*qvtess(phi,flam,h+d3h)-27.q0*v &
        )/(8.q0*d3phi*d3lam*d3h)
! case 8
elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
    gggPhiLamH=-(-qvtess(phi-2.q0*d3phi,flam-2.q0*d3lam,h-2.q0*d3h) &
        +4.q0*qvtess(phi-2.q0*d3phi,flam-2.q0*d3lam,h-d3h)-3.q0*qvtess(phi-2.q0*d3phi,flam-2.q0*d3lam,h) &
        +4.q0*qvtess(phi-2.q0*d3phi,flam-d3lam,h-2.q0*d3h)-16.q0*qvtess(phi-2.q0*d3phi,flam-d3lam,h-d3h) &
        +12.q0*qvtess(phi-2.q0*d3phi,flam-d3lam,h)-3.q0*qvtess(phi-2.q0*d3phi,flam,h-2.q0*d3h) &
        +12.q0*qvtess(phi-2.q0*d3phi,flam,h-d3h)-9.q0*qvtess(phi-2.q0*d3phi,flam,h) &
        +4.q0*qvtess(phi-d3phi,flam-2.q0*d3lam,h-2.q0*d3h)-16.q0*qvtess(phi-d3phi,flam-2.q0*d3lam,h-d3h) &
        +12.q0*qvtess(phi-d3phi,flam-2.q0*d3lam,h)-16.q0*qvtess(phi-d3phi,flam-d3lam,h-2.q0*d3h) &
        +64.q0*qvtess(phi-d3phi,flam-d3lam,h-d3h)-48.q0*qvtess(phi-d3phi,flam-d3lam,h) &
        +12.q0*qvtess(phi-d3phi,flam,h-2.q0*d3h)-48.q0*qvtess(phi-d3phi,flam,h-d3h) &
        +36.q0*qvtess(phi-d3phi,flam,h)-3.q0*qvtess(phi,flam-2.q0*d3lam,h-2.q0*d3h) &
        +12.q0*qvtess(phi,flam-2.q0*d3lam,h-d3h)-9.q0*qvtess(phi,flam-2.q0*d3lam,h) &
        +12.q0*qvtess(phi,flam-d3lam,h-2.q0*d3h)-48.q0*qvtess(phi,flam-d3lam,h-d3h) &
        +36.q0*qvtess(phi,flam-d3lam,h)-9.q0*qvtess(phi,flam,h-2.q0*d3h) &
        +36.q0*qvtess(phi,flam,h-d3h)-27.q0*v &
        )/(8.q0*d3phi*d3lam*d3h)
! case 9
elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
    gggPhiLamH=(qvtess(phi+d3phi,flam+2.q0*d3lam,h+2.q0*d3h) &
        -4.q0*qvtess(phi+d3phi,flam+2.q0*d3lam,h+d3h)+3.q0*qvtess(phi+d3phi,flam+2.q0*d3lam,h) &
        -4.q0*qvtess(phi+d3phi,flam+d3lam,h+2.q0*d3h)+16.q0*qvtess(phi+d3phi,flam+d3lam,h+d3h) &
        -12.q0*qvtess(phi+d3phi,flam+d3lam,h)+3.q0*qvtess(phi+d3phi,flam,h+2.q0*d3h) &
        -12.q0*qvtess(phi+d3phi,flam,h+d3h)+9.q0*qvtess(phi+d3phi,flam,h) &
        -qvtess(phi-d3phi,flam+2.q0*d3lam,h+2.q0*d3h) &
        +4.q0*qvtess(phi-d3phi,flam+2.q0*d3lam,h+d3h)-3.q0*qvtess(phi-d3phi,flam+2.q0*d3lam,h) &
        +4.q0*qvtess(phi-d3phi,flam+d3lam,h+2.q0*d3h)-16.q0*qvtess(phi-d3phi,flam+d3lam,h+d3h) &
        +12.q0*qvtess(phi-d3phi,flam+d3lam,h)-3.q0*qvtess(phi-d3phi,flam,h+2.q0*d3h) &
        +12.q0*qvtess(phi-d3phi,flam,h+d3h)-9.q0*qvtess(phi-d3phi,flam,h) &
        )/(8.q0*d3phi*d3lam*d3h)
! case 10
elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
    gggPhiLamH=-(qvtess(phi+d3phi,flam+2.q0*d3lam,h-2.q0*d3h) &
        -4.q0*qvtess(phi+d3phi,flam+2.q0*d3lam,h-d3h)+3.q0*qvtess(phi+d3phi,flam+2.q0*d3lam,h) &
        -4.q0*qvtess(phi+d3phi,flam+d3lam,h-2.q0*d3h)+16.q0*qvtess(phi+d3phi,flam+d3lam,h-d3h) &
        -12.q0*qvtess(phi+d3phi,flam+d3lam,h)+3.q0*qvtess(phi+d3phi,flam,h-2.q0*d3h) &
        -12.q0*qvtess(phi+d3phi,flam,h-d3h)+9.q0*qvtess(phi+d3phi,flam,h) &
        -qvtess(phi-d3phi,flam+2.q0*d3lam,h-2.q0*d3h) &
        +4.q0*qvtess(phi-d3phi,flam+2.q0*d3lam,h-d3h)-3.q0*qvtess(phi-d3phi,flam+2.q0*d3lam,h) &
        +4.q0*qvtess(phi-d3phi,flam+d3lam,h-2.q0*d3h)-16.q0*qvtess(phi-d3phi,flam+d3lam,h-d3h) &
        +12.q0*qvtess(phi-d3phi,flam+d3lam,h)-3.q0*qvtess(phi-d3phi,flam,h-2.q0*d3h) &
        +12.q0*qvtess(phi-d3phi,flam,h-d3h)-9.q0*qvtess(phi-d3phi,flam,h) &
        )/(8.q0*d3phi*d3lam*d3h)
! case 11
elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
    gggPhiLamH=-(qvtess(phi+d3phi,flam-2.q0*d3lam,h+2.q0*d3h) &
        -4.q0*qvtess(phi+d3phi,flam-2.q0*d3lam,h+d3h)+3.q0*qvtess(phi+d3phi,flam-2.q0*d3lam,h) &
        -4.q0*qvtess(phi+d3phi,flam-d3lam,h+2.q0*d3h)+16.q0*qvtess(phi+d3phi,flam-d3lam,h+d3h) &
        -12.q0*qvtess(phi+d3phi,flam-d3lam,h)+3.q0*qvtess(phi+d3phi,flam,h+2.q0*d3h) &
        -12.q0*qvtess(phi+d3phi,flam,h+d3h)+9.q0*qvtess(phi+d3phi,flam,h) &
        -qvtess(phi-d3phi,flam-2.q0*d3lam,h+2.q0*d3h) &
        +4.q0*qvtess(phi-d3phi,flam-2.q0*d3lam,h+d3h)-3.q0*qvtess(phi-d3phi,flam-2.q0*d3lam,h) &
        +4.q0*qvtess(phi-d3phi,flam-d3lam,h+2.q0*d3h)-16.q0*qvtess(phi-d3phi,flam-d3lam,h+d3h) &
        +12.q0*qvtess(phi-d3phi,flam-d3lam,h)-3.q0*qvtess(phi-d3phi,flam,h+2.q0*d3h) &
        +12.q0*qvtess(phi-d3phi,flam,h+d3h)-9.q0*qvtess(phi-d3phi,flam,h) &
        )/(8.q0*d3phi*d3lam*d3h)
! case 12
elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
    gggPhiLamH=(qvtess(phi+d3phi,flam-2.q0*d3lam,h-2.q0*d3h) &
        -4.q0*qvtess(phi+d3phi,flam-2.q0*d3lam,h-d3h)+3.q0*qvtess(phi+d3phi,flam-2.q0*d3lam,h) &
        -4.q0*qvtess(phi+d3phi,flam-d3lam,h-2.q0*d3h)+16.q0*qvtess(phi+d3phi,flam-d3lam,h-d3h) &
        -12.q0*qvtess(phi+d3phi,flam-d3lam,h)+3.q0*qvtess(phi+d3phi,flam,h-2.q0*d3h) &
        -12.q0*qvtess(phi+d3phi,flam,h-d3h)+9.q0*qvtess(phi+d3phi,flam,h) &
        -qvtess(phi-d3phi,flam-2.q0*d3lam,h-2.q0*d3h) &
        +4.q0*qvtess(phi-d3phi,flam-2.q0*d3lam,h-d3h)-3.q0*qvtess(phi-d3phi,flam-2.q0*d3lam,h) &
        +4.q0*qvtess(phi-d3phi,flam-d3lam,h-2.q0*d3h)-16.q0*qvtess(phi-d3phi,flam-d3lam,h-d3h) &
        +12.q0*qvtess(phi-d3phi,flam-d3lam,h)-3.q0*qvtess(phi-d3phi,flam,h-2.q0*d3h) &
        +12.q0*qvtess(phi-d3phi,flam,h-d3h)-9.q0*qvtess(phi-d3phi,flam,h) &
        )/(8.q0*d3phi*d3lam*d3h)
! case 13
elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
    gggPhiLamH=(qvtess(phi+2.q0*d3phi,flam+d3lam,h+2.q0*d3h) &
        -4.q0*qvtess(phi+2.q0*d3phi,flam+d3lam,h+d3h)+3.q0*qvtess(phi+2.q0*d3phi,flam+d3lam,h) &
        -4.q0*qvtess(phi+d3phi,flam+d3lam,h+2.q0*d3h)+16.q0*qvtess(phi+d3phi,flam+d3lam,h+d3h) &
        -12.q0*qvtess(phi+d3phi,flam+d3lam,h)+3.q0*qvtess(phi,flam+d3lam,h+2.q0*d3h) &
        -12.q0*qvtess(phi,flam+d3lam,h+d3h)+9.q0*qvtess(phi,flam+d3lam,h) &
        -qvtess(phi+2.q0*d3phi,flam-d3lam,h+2.q0*d3h) &
        +4.q0*qvtess(phi+2.q0*d3phi,flam-d3lam,h+d3h)-3.q0*qvtess(phi+2.q0*d3phi,flam-d3lam,h) &
        +4.q0*qvtess(phi+d3phi,flam-d3lam,h+2.q0*d3h)-16.q0*qvtess(phi+d3phi,flam-d3lam,h+d3h) &
        +12.q0*qvtess(phi+d3phi,flam-d3lam,h)-3.q0*qvtess(phi,flam-d3lam,h+2.q0*d3h) &
        +12.q0*qvtess(phi,flam-d3lam,h+d3h)-9.q0*qvtess(phi,flam-d3lam,h) &
        )/(8.q0*d3phi*d3lam*d3h)
! case 14
elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
    gggPhiLamH=-(qvtess(phi+2.q0*d3phi,flam+d3lam,h-2.q0*d3h) &
        -4.q0*qvtess(phi+2.q0*d3phi,flam+d3lam,h-d3h)+3.q0*qvtess(phi+2.q0*d3phi,flam+d3lam,h) &
        -4.q0*qvtess(phi+d3phi,flam+d3lam,h-2.q0*d3h)+16.q0*qvtess(phi+d3phi,flam+d3lam,h-d3h) &
        -12.q0*qvtess(phi+d3phi,flam+d3lam,h)+3.q0*qvtess(phi,flam+d3lam,h-2.q0*d3h) &
        -12.q0*qvtess(phi,flam+d3lam,h-d3h)+9.q0*qvtess(phi,flam+d3lam,h) &
        -qvtess(phi+2.q0*d3phi,flam-d3lam,h-2.q0*d3h) &
        +4.q0*qvtess(phi+2.q0*d3phi,flam-d3lam,h-d3h)-3.q0*qvtess(phi+2.q0*d3phi,flam-d3lam,h) &
        +4.q0*qvtess(phi+d3phi,flam-d3lam,h-2.q0*d3h)-16.q0*qvtess(phi+d3phi,flam-d3lam,h-d3h) &
        +12.q0*qvtess(phi+d3phi,flam-d3lam,h)-3.q0*qvtess(phi,flam-d3lam,h-2.q0*d3h) &
        +12.q0*qvtess(phi,flam-d3lam,h-d3h)-9.q0*qvtess(phi,flam-d3lam,h) &
        )/(8.q0*d3phi*d3lam*d3h)
! case 15
elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
    gggPhiLamH=-(qvtess(phi-2.q0*d3phi,flam+d3lam,h+2.q0*d3h) &
        -4.q0*qvtess(phi-2.q0*d3phi,flam+d3lam,h+d3h)+3.q0*qvtess(phi-2.q0*d3phi,flam+d3lam,h) &
        -4.q0*qvtess(phi-d3phi,flam+d3lam,h+2.q0*d3h)+16.q0*qvtess(phi-d3phi,flam+d3lam,h+d3h) &
        -12.q0*qvtess(phi-d3phi,flam+d3lam,h)+3.q0*qvtess(phi,flam+d3lam,h+2.q0*d3h) &
        -12.q0*qvtess(phi,flam+d3lam,h+d3h)+9.q0*qvtess(phi,flam+d3lam,h) &
        -qvtess(phi-2.q0*d3phi,flam-d3lam,h+2.q0*d3h) &
        +4.q0*qvtess(phi-2.q0*d3phi,flam-d3lam,h+d3h)-3.q0*qvtess(phi-2.q0*d3phi,flam-d3lam,h) &
        +4.q0*qvtess(phi-d3phi,flam-d3lam,h+2.q0*d3h)-16.q0*qvtess(phi-d3phi,flam-d3lam,h+d3h) &
        +12.q0*qvtess(phi-d3phi,flam-d3lam,h)-3.q0*qvtess(phi,flam-d3lam,h+2.q0*d3h) &
        +12.q0*qvtess(phi,flam-d3lam,h+d3h)-9.q0*qvtess(phi,flam-d3lam,h) &
        )/(8.q0*d3phi*d3lam*d3h)
! case 16
elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
    gggPhiLamH=(qvtess(phi-2.q0*d3phi,flam+d3lam,h-2.q0*d3h) &
        -4.q0*qvtess(phi-2.q0*d3phi,flam+d3lam,h-d3h)+3.q0*qvtess(phi-2.q0*d3phi,flam+d3lam,h) &
        -4.q0*qvtess(phi-d3phi,flam+d3lam,h-2.q0*d3h)+16.q0*qvtess(phi-d3phi,flam+d3lam,h-d3h) &
        -12.q0*qvtess(phi-d3phi,flam+d3lam,h)+3.q0*qvtess(phi,flam+d3lam,h-2.q0*d3h) &
        -12.q0*qvtess(phi,flam+d3lam,h-d3h)+9.q0*qvtess(phi,flam+d3lam,h) &
        -qvtess(phi-2.q0*d3phi,flam-d3lam,h-2.q0*d3h) &
        +4.q0*qvtess(phi-2.q0*d3phi,flam-d3lam,h-d3h)-3.q0*qvtess(phi-2.q0*d3phi,flam-d3lam,h) &
        +4.q0*qvtess(phi-d3phi,flam-d3lam,h-2.q0*d3h)-16.q0*qvtess(phi-d3phi,flam-d3lam,h-d3h) &
        +12.q0*qvtess(phi-d3phi,flam-d3lam,h)-3.q0*qvtess(phi,flam-d3lam,h-2.q0*d3h) &
        +12.q0*qvtess(phi,flam-d3lam,h-d3h)-9.q0*qvtess(phi,flam-d3lam,h) &
        )/(8.q0*d3phi*d3lam*d3h)
! case 17
elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
    gggPhiLamH=(qvtess(phi+2.q0*d3phi,flam+2.q0*d3lam,h+d3h) &
        -4.q0*qvtess(phi+2.q0*d3phi,flam+d3lam,h+d3h)+3.q0*qvtess(phi+2.q0*d3phi,flam,h+d3h) &
        -4.q0*qvtess(phi+d3phi,flam+2.q0*d3lam,h+d3h)+16.q0*qvtess(phi+d3phi,flam+d3lam,h+d3h) &
        -12.q0*qvtess(phi+d3phi,flam,h+d3h)+3.q0*qvtess(phi,flam+2.q0*d3lam,h+d3h) &
        -12.q0*qvtess(phi,flam+d3lam,h+d3h)+9.q0*qvtess(phi,flam,h+d3h) &
        -qvtess(phi+2.q0*d3phi,flam+2.q0*d3lam,h-d3h) &
        +4.q0*qvtess(phi+2.q0*d3phi,flam+d3lam,h-d3h)-3.q0*qvtess(phi+2.q0*d3phi,flam,h-d3h) &
        +4.q0*qvtess(phi+d3phi,flam+2.q0*d3lam,h-d3h)-16.q0*qvtess(phi+d3phi,flam+d3lam,h-d3h) &
        +12.q0*qvtess(phi+d3phi,flam,h-d3h)-3.q0*qvtess(phi,flam+2.q0*d3lam,h-d3h) &
        +12.q0*qvtess(phi,flam+d3lam,h-d3h)-9.q0*qvtess(phi,flam,h-d3h) &
        )/(8.q0*d3phi*d3lam*d3h)
! case 18
elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
    gggPhiLamH=-(qvtess(phi+2.q0*d3phi,flam-2.q0*d3lam,h+d3h) &
        -4.q0*qvtess(phi+2.q0*d3phi,flam-d3lam,h+d3h)+3.q0*qvtess(phi+2.q0*d3phi,flam,h+d3h) &
        -4.q0*qvtess(phi+d3phi,flam-2.q0*d3lam,h+d3h)+16.q0*qvtess(phi+d3phi,flam-d3lam,h+d3h) &
        -12.q0*qvtess(phi+d3phi,flam,h+d3h)+3.q0*qvtess(phi,flam-2.q0*d3lam,h+d3h) &
        -12.q0*qvtess(phi,flam-d3lam,h+d3h)+9.q0*qvtess(phi,flam,h+d3h) &
        -qvtess(phi+2.q0*d3phi,flam-2.q0*d3lam,h-d3h) &
        +4.q0*qvtess(phi+2.q0*d3phi,flam-d3lam,h-d3h)-3.q0*qvtess(phi+2.q0*d3phi,flam,h-d3h) &
        +4.q0*qvtess(phi+d3phi,flam-2.q0*d3lam,h-d3h)-16.q0*qvtess(phi+d3phi,flam-d3lam,h-d3h) &
        +12.q0*qvtess(phi+d3phi,flam,h-d3h)-3.q0*qvtess(phi,flam-2.q0*d3lam,h-d3h) &
        +12.q0*qvtess(phi,flam-d3lam,h-d3h)-9.q0*qvtess(phi,flam,h-d3h) &
        )/(8.q0*d3phi*d3lam*d3h)
! case 19
elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
    gggPhiLamH=-(qvtess(phi-2.q0*d3phi,flam+2.q0*d3lam,h+d3h) &
        -4.q0*qvtess(phi-2.q0*d3phi,flam+d3lam,h+d3h)+3.q0*qvtess(phi-2.q0*d3phi,flam,h+d3h) &
        -4.q0*qvtess(phi-d3phi,flam+2.q0*d3lam,h+d3h)+16.q0*qvtess(phi-d3phi,flam+d3lam,h+d3h) &
        -12.q0*qvtess(phi-d3phi,flam,h+d3h)+3.q0*qvtess(phi,flam+2.q0*d3lam,h+d3h) &
        -12.q0*qvtess(phi,flam+d3lam,h+d3h)+9.q0*qvtess(phi,flam,h+d3h) &
        -qvtess(phi-2.q0*d3phi,flam+2.q0*d3lam,h-d3h) &
        +4.q0*qvtess(phi-2.q0*d3phi,flam+d3lam,h-d3h)-3.q0*qvtess(phi-2.q0*d3phi,flam,h-d3h) &
        +4.q0*qvtess(phi-d3phi,flam+2.q0*d3lam,h-d3h)-16.q0*qvtess(phi-d3phi,flam+d3lam,h-d3h) &
        +12.q0*qvtess(phi-d3phi,flam,h-d3h)-3.q0*qvtess(phi,flam+2.q0*d3lam,h-d3h) &
        +12.q0*qvtess(phi,flam+d3lam,h-d3h)-9.q0*qvtess(phi,flam,h-d3h) &
        )/(8.q0*d3phi*d3lam*d3h)
! case 20
elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
    gggPhiLamH=(qvtess(phi-2.q0*d3phi,flam-2.q0*d3lam,h+d3h) &
        -4.q0*qvtess(phi-2.q0*d3phi,flam-d3lam,h+d3h)+3.q0*qvtess(phi-2.q0*d3phi,flam,h+d3h) &
        -4.q0*qvtess(phi-d3phi,flam-2.q0*d3lam,h+d3h)+16.q0*qvtess(phi-d3phi,flam-d3lam,h+d3h) &
        -12.q0*qvtess(phi-d3phi,flam,h+d3h)+3.q0*qvtess(phi,flam-2.q0*d3lam,h+d3h) &
        -12.q0*qvtess(phi,flam-d3lam,h+d3h)+9.q0*qvtess(phi,flam,h+d3h) &
        -qvtess(phi-2.q0*d3phi,flam-2.q0*d3lam,h-d3h) &
        +4.q0*qvtess(phi-2.q0*d3phi,flam-d3lam,h-d3h)-3.q0*qvtess(phi-2.q0*d3phi,flam,h-d3h) &
        +4.q0*qvtess(phi-d3phi,flam-2.q0*d3lam,h-d3h)-16.q0*qvtess(phi-d3phi,flam-d3lam,h-d3h) &
        +12.q0*qvtess(phi-d3phi,flam,h-d3h)-3.q0*qvtess(phi,flam-2.q0*d3lam,h-d3h) &
        +12.q0*qvtess(phi,flam-d3lam,h-d3h)-9.q0*qvtess(phi,flam,h-d3h) &
        )/(8.q0*d3phi*d3lam*d3h)
! case 21
elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
    gggPhiLamH=(-qvtess(phi+d3phi,flam+d3lam,h+2.q0*d3h) &
        +4.q0*qvtess(phi+d3phi,flam+d3lam,h+d3h)-3.q0*qvtess(phi+d3phi,flam+d3lam,h) &
        +qvtess(phi+d3phi,flam-d3lam,h+2.q0*d3h)-4.q0*qvtess(phi+d3phi,flam-d3lam,h+d3h) &
        +3.q0*qvtess(phi+d3phi,flam-d3lam,h)&
        +qvtess(phi-d3phi,flam+d3lam,h+2.q0*d3h) &
        -4.q0*qvtess(phi-d3phi,flam+d3lam,h+d3h)+3.q0*qvtess(phi-d3phi,flam+d3lam,h) &
        -qvtess(phi-d3phi,flam-d3lam,h+2.q0*d3h)+4.q0*qvtess(phi-d3phi,flam-d3lam,h+d3h) &
        -3.q0*qvtess(phi-d3phi,flam-d3lam,h) &
        )/(8.q0*d3phi*d3lam*d3h)
! case 22
elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
    gggPhiLamH=-(-qvtess(phi+d3phi,flam+d3lam,h-2.q0*d3h) &
        +4.q0*qvtess(phi+d3phi,flam+d3lam,h-d3h)-3.q0*qvtess(phi+d3phi,flam+d3lam,h) &
        +qvtess(phi+d3phi,flam-d3lam,h-2.q0*d3h)-4.q0*qvtess(phi+d3phi,flam-d3lam,h-d3h) &
        +3.q0*qvtess(phi+d3phi,flam-d3lam,h) &
        +qvtess(phi-d3phi,flam+d3lam,h-2.q0*d3h) &
        -4.q0*qvtess(phi-d3phi,flam+d3lam,h-d3h)+3.q0*qvtess(phi-d3phi,flam+d3lam,h) &
        -qvtess(phi-d3phi,flam-d3lam,h-2.q0*d3h)+4.q0*qvtess(phi-d3phi,flam-d3lam,h-d3h) &
        -3.q0*qvtess(phi-d3phi,flam-d3lam,h) &
        )/(8.q0*d3phi*d3lam*d3h)
! case 23
elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
    gggPhiLamH=(-qvtess(phi+d3phi,flam+2.q0*d3lam,h+d3h) &
        +4.q0*qvtess(phi+d3phi,flam+d3lam,h+d3h)-3.q0*qvtess(phi+d3phi,flam,h+d3h) &
        +qvtess(phi+d3phi,flam+2.q0*d3lam,h-d3h)-4.q0*qvtess(phi+d3phi,flam+d3lam,h-d3h) &
        +3.q0*qvtess(phi+d3phi,flam,h-d3h) &
        +qvtess(phi-d3phi,flam+2.q0*d3lam,h+d3h) &
        -4.q0*qvtess(phi-d3phi,flam+d3lam,h+d3h)+3.q0*qvtess(phi-d3phi,flam,h+d3h) &
        -qvtess(phi-d3phi,flam+2.q0*d3lam,h-d3h)+4.q0*qvtess(phi-d3phi,flam+d3lam,h-d3h) &
        -3.q0*qvtess(phi-d3phi,flam,h-d3h) &
        )/(8.q0*d3phi*d3lam*d3h)
! case 24
elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
    gggPhiLamH=-(-qvtess(phi+d3phi,flam-2.q0*d3lam,h+d3h) &
        +4.q0*qvtess(phi+d3phi,flam-d3lam,h+d3h)-3.q0*qvtess(phi+d3phi,flam,h+d3h) &
        +qvtess(phi+d3phi,flam-2.q0*d3lam,h-d3h)-4.q0*qvtess(phi+d3phi,flam-d3lam,h-d3h) &
        +3.q0*qvtess(phi+d3phi,flam,h-d3h) &
        +qvtess(phi-d3phi,flam-2.q0*d3lam,h+d3h) &
        -4.q0*qvtess(phi-d3phi,flam-d3lam,h+d3h)+3.q0*qvtess(phi-d3phi,flam,h+d3h) &
        -qvtess(phi-d3phi,flam-2.q0*d3lam,h-d3h)+4.q0*qvtess(phi-d3phi,flam-d3lam,h-d3h) &
        -3.q0*qvtess(phi-d3phi,flam,h-d3h) &
        )/(8.q0*d3phi*d3lam*d3h)
! case 25
elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
    gggPhiLamH=(-qvtess(phi+2.q0*d3phi,flam+d3lam,h+d3h) &
        +4.q0*qvtess(phi+d3phi,flam+d3lam,h+d3h)-3.q0*qvtess(phi,flam+d3lam,h+d3h) &
        +qvtess(phi+2.q0*d3phi,flam+d3lam,h-d3h)-4.q0*qvtess(phi+d3phi,flam+d3lam,h-d3h) &
        +3.q0*qvtess(phi,flam+d3lam,h-d3h) &
        +qvtess(phi+2.q0*d3phi,flam-d3lam,h+d3h) &
        -4.q0*qvtess(phi+d3phi,flam-d3lam,h+d3h)+3.q0*qvtess(phi,flam-d3lam,h+d3h) &
        -qvtess(phi+2.q0*d3phi,flam-d3lam,h-d3h)+4.q0*qvtess(phi+d3phi,flam-d3lam,h-d3h) &
        -3.q0*qvtess(phi,flam-d3lam,h-d3h) &
        )/(8.q0*d3phi*d3lam*d3h)
! case 26
elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
    gggPhiLamH=-(-qvtess(phi-2.q0*d3phi,flam+d3lam,h+d3h) &
        +4.q0*qvtess(phi-d3phi,flam+d3lam,h+d3h)-3.q0*qvtess(phi,flam+d3lam,h+d3h) &
        +qvtess(phi-2.q0*d3phi,flam+d3lam,h-d3h)-4.q0*qvtess(phi-d3phi,flam+d3lam,h-d3h) &
        +3.q0*qvtess(phi,flam+d3lam,h-d3h) &
        +qvtess(phi-2.q0*d3phi,flam-d3lam,h+d3h) &
        -4.q0*qvtess(phi-d3phi,flam-d3lam,h+d3h)+3.q0*qvtess(phi,flam-d3lam,h+d3h) &
        -qvtess(phi-2.q0*d3phi,flam-d3lam,h-d3h)+4.q0*qvtess(phi-d3phi,flam-d3lam,h-d3h) &
        -3.q0*qvtess(phi,flam-d3lam,h-d3h) &
        )/(8.q0*d3phi*d3lam*d3h)
! case 27
else
    gggPhiLamH=(qvtess(phi+d3phi,flam+d3lam,h+d3h) &
        -qvtess(phi+d3phi,flam+d3lam,h-d3h)-qvtess(phi+d3phi,flam-d3lam,h+d3h) &
        +qvtess(phi+d3phi,flam-d3lam,h-d3h)-qvtess(phi-d3phi,flam+d3lam,h+d3h) &
        +qvtess(phi-d3phi,flam+d3lam,h-d3h)+qvtess(phi-d3phi,flam-d3lam,h+d3h) &
        -qvtess(phi-d3phi,flam-d3lam,h-d3h))/(8.q0*d3phi*d3lam*d3h)
endif
gggPhiLamH=(gggPhiLamH/c-2.q0*r*ggPhiLam+t*(gLam+r*ggLamH))/(r*r)
!
!   Compute gggLamLamPhi
!
if(h.ge.HB.and.h.le.HT) then
! case 1
    if(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE))) then
        gggLamLamPhi=(qvtess(phi+2.q0*d3phi,flam+3.q0*d3lam,h)-4.q0*qvtess(phi+d3phi,flam+3.q0*d3lam,h) &
            +3.q0*qvtess(phi,flam+3.q0*d3lam,h)-4.q0*qvtess(phi+2.q0*d3phi,flam+2.q0*d3lam,h) &
            +16.q0*qvtess(phi+d3phi,flam+2.q0*d3lam,h)-12.q0*qvtess(phi,flam+2.q0*d3lam,h) &
            +5.q0*qvtess(phi+2.q0*d3phi,flam+d3lam,h)-20.q0*qvtess(phi+d3phi,flam+d3lam,h) &
            +15.q0*qvtess(phi,flam+d3lam,h)-2.q0*qvtess(phi+2.q0*d3phi,flam,h) &
            +8.q0*qvtess(phi+d3phi,flam,h)-6.q0*v)/(2.q0*d3phi*d3lam*d3lam)
! case 2
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE))) then
        gggLamLamPhi=-(qvtess(phi-2.q0*d3phi,flam+3.q0*d3lam,h)-4.q0*qvtess(phi-d3phi,flam+3.q0*d3lam,h) &
            +3.q0*qvtess(phi,flam+3.q0*d3lam,h)-4.q0*qvtess(phi-2.q0*d3phi,flam+2.q0*d3lam,h) &
            +16.q0*qvtess(phi-d3phi,flam+2.q0*d3lam,h)-12.q0*qvtess(phi,flam+2.q0*d3lam,h) &
            +5.q0*qvtess(phi-2.q0*d3phi,flam+d3lam,h)-20.q0*qvtess(phi-d3phi,flam+d3lam,h) &
            +15.q0*qvtess(phi,flam+d3lam,h)-2.q0*qvtess(phi-2.q0*d3phi,flam,h) &
            +8.q0*qvtess(phi-d3phi,flam,h)-6.q0*v)/(2.q0*d3phi*d3lam*d3lam)
! case 3
    elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE))) then
        gggLamLamPhi=(qvtess(phi+2.q0*d3phi,flam-3.q0*d3lam,h)-4.q0*qvtess(phi+d3phi,flam-3.q0*d3lam,h) &
            +3.q0*qvtess(phi,flam-3.q0*d3lam,h)-4.q0*qvtess(phi+2.q0*d3phi,flam-2.q0*d3lam,h) &
            +16.q0*qvtess(phi+d3phi,flam-2.q0*d3lam,h)-12.q0*qvtess(phi,flam-2.q0*d3lam,h) &
            +5.q0*qvtess(phi+2.q0*d3phi,flam-d3lam,h)-20.q0*qvtess(phi+d3phi,flam-d3lam,h) &
            +15.q0*qvtess(phi,flam-d3lam,h)-2.q0*qvtess(phi+2.q0*d3phi,flam,h) &
            +8.q0*qvtess(phi+d3phi,flam,h)-6.q0*v)/(2.q0*d3phi*d3lam*d3lam)
! case 4
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE))) then
        gggLamLamPhi=-(qvtess(phi-2.q0*d3phi,flam-3.q0*d3lam,h)-4.q0*qvtess(phi-d3phi,flam-3.q0*d3lam,h) &
            +3.q0*qvtess(phi,flam-3.q0*d3lam,h)-4.q0*qvtess(phi-2.q0*d3phi,flam-2.q0*d3lam,h) &
            +16.q0*qvtess(phi-d3phi,flam-2.q0*d3lam,h)-12.q0*qvtess(phi,flam-2.q0*d3lam,h) &
            +5.q0*qvtess(phi-2.q0*d3phi,flam-d3lam,h)-20.q0*qvtess(phi-d3phi,flam-d3lam,h) &
            +15.q0*qvtess(phi,flam-d3lam,h)-2.q0*qvtess(phi-2.q0*d3phi,flam,h) &
            +8.q0*qvtess(phi-d3phi,flam,h)-6.q0*v)/(2.q0*d3phi*d3lam*d3lam)
! case 5
    elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE))) then
        gggLamLamPhi=(-qvtess(phi+d3phi,flam+3.q0*d3lam,h)+qvtess(phi-d3phi,flam+3.q0*d3lam,h) &
            +4.q0*qvtess(phi+d3phi,flam+2.q0*d3lam,h)-4.q0*qvtess(phi-d3phi,flam+2.q0*d3lam,h) &
            -5.q0*qvtess(phi+d3phi,flam+d3lam,h)+5.q0*qvtess(phi-d3phi,flam+d3lam,h) &
            +2.q0*qvtess(phi+d3phi,flam,h)-2.q0*qvtess(phi-d3phi,flam,h))/(2.q0*d3phi*d3lam*d3lam)
! case 6
    elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE))) then
        gggLamLamPhi=(-qvtess(phi+d3phi,flam-3.q0*d3lam,h)+qvtess(phi-d3phi,flam-3.q0*d3lam,h) &
            +4.q0*qvtess(phi+d3phi,flam-2.q0*d3lam,h)-4.q0*qvtess(phi-d3phi,flam-2.q0*d3lam,h) &
            -5.q0*qvtess(phi+d3phi,flam-d3lam,h)+5.q0*qvtess(phi-d3phi,flam-d3lam,h) &
            +2.q0*qvtess(phi+d3phi,flam,h)-2.q0*qvtess(phi-d3phi,flam,h))/(2.q0*d3phi*d3lam*d3lam)
! case 7
    elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE))) then
        gggLamLamPhi=(-qvtess(phi+2.q0*d3phi,flam+d3lam,h)+4.q0*qvtess(phi+d3phi,flam+d3lam,h) &
            -3.q0*qvtess(phi,flam+d3lam,h)+2.q0*qvtess(phi+2.q0*d3phi,flam,h) &
            -8.q0*qvtess(phi+d3phi,flam,h)+6.q0*v &
            -qvtess(phi+2.q0*d3phi,flam-d3lam,h)+4.q0*qvtess(phi+d3phi,flam-d3lam,h) &
            -3.q0*qvtess(phi,flam-d3lam,h))/(2.q0*d3phi*d3lam*d3lam)
! case 8
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE))) then
        gggLamLamPhi=-(-qvtess(phi-2.q0*d3phi,flam+d3lam,h)+4.q0*qvtess(phi-d3phi,flam+d3lam,h) &
            -3.q0*qvtess(phi,flam+d3lam,h)+2.q0*qvtess(phi-2.q0*d3phi,flam,h) &
            -8.q0*qvtess(phi-d3phi,flam,h)+6.q0*v &
            -qvtess(phi-2.q0*d3phi,flam-d3lam,h)+4.q0*qvtess(phi-d3phi,flam-d3lam,h) &
            -3.q0*qvtess(phi,flam-d3lam,h))/(2.q0*d3phi*d3lam*d3lam)
! case 9
    else
        gggLamLamPhi=(qvtess(phi+d3phi,flam+d3lam,h)-qvtess(phi-d3phi,flam+d3lam,h) &
            -2.q0*qvtess(phi+d3phi,flam,h)+2.q0*qvtess(phi-d3phi,flam,h) &
            +qvtess(phi+d3phi,flam-d3lam,h)-qvtess(phi-d3phi,flam-d3lam,h))/(2.q0*d3phi*d3lam*d3lam)
    endif
else
    gggLamLamPhi=(qvtess(phi+d3phi,flam+d3lam,h)-qvtess(phi-d3phi,flam+d3lam,h) &
            -2.q0*qvtess(phi+d3phi,flam,h)+2.q0*qvtess(phi-d3phi,flam,h) &
            +qvtess(phi+d3phi,flam-d3lam,h)-qvtess(phi-d3phi,flam-d3lam,h))/(2.q0*d3phi*d3lam*d3lam)
endif
gggLamLamPhi=(gggLamLamPhi-r*gPhi)/(r*r*r*c*c) &
+ (gPhi+r*ggPhiH-t*(gH-2.q0*r*ggLamLam+r*ggPhiPhi)+2.q0*t*t*gPhi)/(r*r)
!
!   Compute gggLamLamLam
!
if((phi.ge.PhiS.and.phi.le.PhiN).and.(h.ge.HB.and.h.le.HT)) then
    if((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) then
        gggLamLamLam=(-3.q0*qvtess(phi,flam+4.q0*d3lam,h)+14.q0*qvtess(phi,flam+3.q0*d3lam,h) &
            -24.q0*qvtess(phi,flam+2.q0*d3lam,h)+18.q0*qvtess(phi,flam+d3lam,h) &
            -5.q0*v)/(2.q0*d3lam*d3lam*d3lam)
    elseif((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) then
        gggLamLamLam=(3.q0*qvtess(phi,flam-4.q0*d3lam,h)-14.q0*qvtess(phi,flam-3.q0*d3lam,h) &
            +24.q0*qvtess(phi,flam-2.q0*d3lam,h)-18.q0*qvtess(phi,flam-d3lam,h) &
            +5.q0*v)/(2.q0*d3lam*d3lam*d3lam)
    else
        gggLamLamLam=(qvtess(phi,flam+2.q0*d3lam,h)-2.q0*qvtess(phi,flam+d3lam,h) &
            +2.q0*qvtess(phi,flam-d3lam,h)-qvtess(phi,flam-2.q0*d3lam,h) &
            )/(2.q0*d3lam*d3lam*d3lam)
    endif
else
    gggLamLamLam=(qvtess(phi,flam+2.q0*d3lam,h)-2.q0*qvtess(phi,flam+d3lam,h) &
            +2.q0*qvtess(phi,flam-d3lam,h)-qvtess(phi,flam-2.q0*d3lam,h) &
            )/(2.q0*d3lam*d3lam*d3lam)
endif
gggLamLamLam=(r*c*gLam+gggLamLamLam)/(r*r*r*c*c*c)+3.q0*(ggLamH-t*ggPhiLam)/r
!
!   Compute gggLamLamH
!
if(phi.ge.PhiS.and.phi.le.PhiN) then
! case 1
    if(((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggLamLamH=(qvtess(phi,flam+3.q0*d3lam,h+2.q0*d3h)-4.q0*qvtess(phi,flam+3.q0*d3lam,h+d3h) &
            +3.q0*qvtess(phi,flam+3.q0*d3lam,h)-4.q0*qvtess(phi,flam+2.q0*d3lam,h+2.q0*d3h) &
            +16.q0*qvtess(phi,flam+2.q0*d3lam,h+d3h)-12.q0*qvtess(phi,flam+2.q0*d3lam,h) &
            +5.q0*qvtess(phi,flam+d3lam,h+2.q0*d3h)-20.q0*qvtess(phi,flam+d3lam,h+d3h) &
            +15.q0*qvtess(phi,flam+d3lam,h)-2.q0*qvtess(phi,flam,h+2.q0*d3h) &
            +8.q0*qvtess(phi,flam,h+d3h)-6.q0*v &
            )/(2.q0*d3lam*d3lam*d3h)
! case 2
    elseif(((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggLamLamH=-(qvtess(phi,flam+3.q0*d3lam,h-2.q0*d3h)-4.q0*qvtess(phi,flam+3.q0*d3lam,h-d3h) &
            +3.q0*qvtess(phi,flam+3.q0*d3lam,h)-4.q0*qvtess(phi,flam+2.q0*d3lam,h-2.q0*d3h) &
            +16.q0*qvtess(phi,flam+2.q0*d3lam,h-d3h)-12.q0*qvtess(phi,flam+2.q0*d3lam,h) &
            +5.q0*qvtess(phi,flam+d3lam,h-2.q0*d3h)-20.q0*qvtess(phi,flam+d3lam,h-d3h) &
            +15.q0*qvtess(phi,flam+d3lam,h)-2.q0*qvtess(phi,flam,h-2.q0*d3h) &
            +8.q0*qvtess(phi,flam,h-d3h)-6.q0*v &
            )/(2.q0*d3lam*d3lam*d3h)
! case 3
    elseif(((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggLamLamH=(qvtess(phi,flam-3.q0*d3lam,h+2.q0*d3h)-4.q0*qvtess(phi,flam-3.q0*d3lam,h+d3h) &
            +3.q0*qvtess(phi,flam-3.q0*d3lam,h)-4.q0*qvtess(phi,flam-2.q0*d3lam,h+2.q0*d3h) &
            +16.q0*qvtess(phi,flam-2.q0*d3lam,h+d3h)-12.q0*qvtess(phi,flam-2.q0*d3lam,h) &
            +5.q0*qvtess(phi,flam-d3lam,h+2.q0*d3h)-20.q0*qvtess(phi,flam-d3lam,h+d3h) &
            +15.q0*qvtess(phi,flam-d3lam,h)-2.q0*qvtess(phi,flam,h+2.q0*d3h) &
            +8.q0*qvtess(phi,flam,h+d3h)-6.q0*v &
            )/(2.q0*d3lam*d3lam*d3h)
! case 4
    elseif(((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggLamLamH=-(qvtess(phi,flam-3.q0*d3lam,h-2.q0*d3h)-4.q0*qvtess(phi,flam-3.q0*d3lam,h-d3h) &
            +3.q0*qvtess(phi,flam-3.q0*d3lam,h)-4.q0*qvtess(phi,flam-2.q0*d3lam,h-2.q0*d3h) &
            +16.q0*qvtess(phi,flam-2.q0*d3lam,h-d3h)-12.q0*qvtess(phi,flam-2.q0*d3lam,h) &
            +5.q0*qvtess(phi,flam-d3lam,h-2.q0*d3h)-20.q0*qvtess(phi,flam-d3lam,h-d3h) &
            +15.q0*qvtess(phi,flam-d3lam,h)-2.q0*qvtess(phi,flam,h-2.q0*d3h) &
            +8.q0*qvtess(phi,flam,h-d3h)-6.q0*v &
            )/(2.q0*d3lam*d3lam*d3h)
! case 5
    elseif(((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
        gggLamLamH=(-qvtess(phi,flam+3.q0*d3lam,h+d3h)+qvtess(phi,flam+3.q0*d3lam,h-d3h) &
            +4.q0*qvtess(phi,flam+2.q0*d3lam,h+d3h)-4.q0*qvtess(phi,flam+2.q0*d3lam,h-d3h) &
            -5.q0*qvtess(phi,flam+d3lam,h+d3h)+5.q0*qvtess(phi,flam+d3lam,h-d3h) &
            +2.q0*qvtess(phi,flam,h+d3h)-2.q0*qvtess(phi,flam,h-d3h) &
            )/(2.q0*d3lam*d3lam*d3h)
! case 6
    elseif(((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
        gggLamLamH=(-qvtess(phi,flam-3.q0*d3lam,h+d3h)+qvtess(phi,flam-3.q0*d3lam,h-d3h) &
            +4.q0*qvtess(phi,flam-2.q0*d3lam,h+d3h)-4.q0*qvtess(phi,flam-2.q0*d3lam,h-d3h) &
            -5.q0*qvtess(phi,flam-d3lam,h+d3h)+5.q0*qvtess(phi,flam-d3lam,h-d3h) &
            +2.q0*qvtess(phi,flam,h+d3h)-2.q0*qvtess(phi,flam,h-d3h) &
            )/(2.q0*d3lam*d3lam*d3h)
! case 7
    elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggLamLamH=(-qvtess(phi,flam+d3lam,h+2.q0*d3h)+4.q0*qvtess(phi,flam+d3lam,h+d3h) &
            -3.q0*qvtess(phi,flam+d3lam,h)+2.q0*qvtess(phi,flam,h+2.q0*d3h) &
            -8.q0*qvtess(phi,flam,h+d3h)+6.q0*v -qvtess(phi,flam-d3lam,h+2.q0*d3h) &
            +4.q0*qvtess(phi,flam-d3lam,h+d3h)-3.q0*qvtess(phi,flam-d3lam,h) &
            )/(2.q0*d3lam*d3lam*d3h)
! case 8
    elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggLamLamH=-(-qvtess(phi,flam+d3lam,h-2.q0*d3h)+4.q0*qvtess(phi,flam+d3lam,h-d3h) &
            -3.q0*qvtess(phi,flam+d3lam,h)+2.q0*qvtess(phi,flam,h-2.q0*d3h) &
            -8.q0*qvtess(phi,flam,h-d3h)+6.q0*v -qvtess(phi,flam-d3lam,h-2.q0*d3h) &
            +4.q0*qvtess(phi,flam-d3lam,h-d3h)-3.q0*qvtess(phi,flam-d3lam,h) &
            )/(2.q0*d3lam*d3lam*d3h)
! case 9
    else
        gggLamLamH=(qvtess(phi,flam+d3lam,h+d3h)-qvtess(phi,flam+d3lam,h-d3h) &
            -2.q0*qvtess(phi,flam,h+d3h)+2.q0*qvtess(phi,flam,h-d3h) &
            +qvtess(phi,flam-d3lam,h+d3h)-qvtess(phi,flam-d3lam,h-d3h) &
            )/(2.q0*d3lam*d3lam*d3h)
    endif
else
    gggLamLamH=(qvtess(phi,flam+d3lam,h+d3h)-qvtess(phi,flam+d3lam,h-d3h) &
            -2.q0*qvtess(phi,flam,h+d3h)+2.q0*qvtess(phi,flam,h-d3h) &
            +qvtess(phi,flam-d3lam,h+d3h)-qvtess(phi,flam-d3lam,h-d3h) &
            )/(2.q0*d3lam*d3lam*d3h)
endif
gggLamLamH=(gH+r*ggHH-2.q0*r*ggLamLam+gggLamLamH/(c*c)-t*(gPhi+r*ggPhiH))/(r*r)
!
!   Compute gggHHPhi
!
if(flam.ge.FlamW.and.flam.le.FlamE) then
! case 1
    if(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggHHPhi=(qvtess(phi+2.q0*d3phi,flam,h+3.q0*d3h)-4.q0*qvtess(phi+d3phi,flam,h+3.q0*d3h) &
            +3.q0*qvtess(phi,flam,h+3.q0*d3h)-4.q0*qvtess(phi+2.q0*d3phi,flam,h+2.q0*d3h) &
            +16.q0*qvtess(phi+d3phi,flam,h+2.q0*d3h)-12.q0*qvtess(phi,flam,h+2.q0*d3h) &
            +5.q0*qvtess(phi+2.q0*d3phi,flam,h+d3h)-20.q0*qvtess(phi+d3phi,flam,h+d3h) &
            +15.q0*qvtess(phi,flam,h+d3h)-2.q0*qvtess(phi+2.q0*d3phi,flam,h) &
            +8.q0*qvtess(phi+d3phi,flam,h)-6.q0*v)/(2.q0*d3phi*d3h*d3h)
! case 2
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggHHPhi=-(qvtess(phi-2.q0*d3phi,flam,h+3.q0*d3h)-4.q0*qvtess(phi-d3phi,flam,h+3.q0*d3h) &
            +3.q0*qvtess(phi,flam,h+3.q0*d3h)-4.q0*qvtess(phi-2.q0*d3phi,flam,h+2.q0*d3h) &
            +16.q0*qvtess(phi-d3phi,flam,h+2.q0*d3h)-12.q0*qvtess(phi,flam,h+2.q0*d3h) &
            +5.q0*qvtess(phi-2.q0*d3phi,flam,h+d3h)-20.q0*qvtess(phi-d3phi,flam,h+d3h) &
            +15.q0*qvtess(phi,flam,h+d3h)-2.q0*qvtess(phi-2.q0*d3phi,flam,h) &
            +8.q0*qvtess(phi-d3phi,flam,h)-6.q0*v)/(2.q0*d3phi*d3h*d3h)
! case 3
    elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggHHPhi=(qvtess(phi+2.q0*d3phi,flam,h-3.q0*d3h)-4.q0*qvtess(phi+d3phi,flam,h-3.q0*d3h) &
            +3.q0*qvtess(phi,flam,h-3.q0*d3h)-4.q0*qvtess(phi+2.q0*d3phi,flam,h-2.q0*d3h) &
            +16.q0*qvtess(phi+d3phi,flam,h-2.q0*d3h)-12.q0*qvtess(phi,flam,h-2.q0*d3h) &
            +5.q0*qvtess(phi+2.q0*d3phi,flam,h-d3h)-20.q0*qvtess(phi+d3phi,flam,h-d3h) &
            +15.q0*qvtess(phi,flam,h-d3h)-2.q0*qvtess(phi+2.q0*d3phi,flam,h) &
            +8.q0*qvtess(phi+d3phi,flam,h)-6.q0*v)/(2.q0*d3phi*d3h*d3h)
! case 4
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggHHPhi=-(qvtess(phi-2.q0*d3phi,flam,h-3.q0*d3h)-4.q0*qvtess(phi-d3phi,flam,h-3.q0*d3h) &
            +3.q0*qvtess(phi,flam,h-3.q0*d3h)-4.q0*qvtess(phi-2.q0*d3phi,flam,h-2.q0*d3h) &
            +16.q0*qvtess(phi-d3phi,flam,h-2.q0*d3h)-12.q0*qvtess(phi,flam,h-2.q0*d3h) &
            +5.q0*qvtess(phi-2.q0*d3phi,flam,h-d3h)-20.q0*qvtess(phi-d3phi,flam,h-d3h) &
            +15.q0*qvtess(phi,flam,h-d3h)-2.q0*qvtess(phi-2.q0*d3phi,flam,h) &
            +8.q0*qvtess(phi-d3phi,flam,h)-6.q0*v)/(2.q0*d3phi*d3h*d3h)
! case 5
    elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggHHPhi=(-qvtess(phi+d3phi,flam,h+3.q0*d3h)+qvtess(phi-d3phi,flam,h+3.q0*d3h) &
            +4.q0*qvtess(phi+d3phi,flam,h+2.q0*d3h)-4.q0*qvtess(phi-d3phi,flam,h+2.q0*d3h) &
            -5.q0*qvtess(phi+d3phi,flam,h+d3h)+5.q0*qvtess(phi-d3phi,flam,h+d3h) &
            +2.q0*qvtess(phi+d3phi,flam,h)-2.q0*qvtess(phi-d3phi,flam,h))/(2.q0*d3phi*d3h*d3h)
! case 6
    elseif(((phi+d3phi.lt.PhiS).or.(phi-d3phi.gt.PhiN) &
        .or.(phi-d3phi.gt.PhiS.and.phi+d3phi.lt.PhiN)) &
    .and.((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggHHPhi=(-qvtess(phi+d3phi,flam,h-3.q0*d3h)+qvtess(phi-d3phi,flam,h-3.q0*d3h) &
            +4.q0*qvtess(phi+d3phi,flam,h-2.q0*d3h)-4.q0*qvtess(phi-d3phi,flam,h-2.q0*d3h) &
            -5.q0*qvtess(phi+d3phi,flam,h-d3h)+5.q0*qvtess(phi-d3phi,flam,h-d3h) &
            +2.q0*qvtess(phi+d3phi,flam,h)-2.q0*qvtess(phi-d3phi,flam,h))/(2.q0*d3phi*d3h*d3h)
! case 7
    elseif(((phi.gt.PhiS.and.phi-d3phi.lt.PhiS) &
        .or.(phi.gt.PhiN.and.phi-d3phi.lt.PhiN)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
        gggHHPhi=(-qvtess(phi+2.q0*d3phi,flam,h+d3h)+4.q0*qvtess(phi+d3phi,flam,h+d3h) &
            -3.q0*qvtess(phi,flam,h+d3h)+2.q0*qvtess(phi+2.q0*d3phi,flam,h) &
            -8.q0*qvtess(phi+d3phi,flam,h)+6.q0*v &
            -qvtess(phi+2.q0*d3phi,flam,h-d3h)+4.q0*qvtess(phi+d3phi,flam,h-d3h) &
            -3.q0*qvtess(phi,flam,h-d3h))/(2.q0*d3phi*d3h*d3h)
! case 8
    elseif(((phi.lt.PhiS.and.phi+d3phi.gt.PhiS) &
        .or.(phi.lt.PhiN.and.phi+d3phi.gt.PhiN)) &
    .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
        gggHHPhi=-(-qvtess(phi-2.q0*d3phi,flam,h+d3h)+4.q0*qvtess(phi-d3phi,flam,h+d3h) &
            -3.q0*qvtess(phi,flam,h+d3h)+2.q0*qvtess(phi-2.q0*d3phi,flam,h) &
            -8.q0*qvtess(phi-d3phi,flam,h)+6.q0*v &
            -qvtess(phi-2.q0*d3phi,flam,h-d3h)+4.q0*qvtess(phi-d3phi,flam,h-d3h) &
            -3.q0*qvtess(phi,flam,h-d3h))/(2.q0*d3phi*d3h*d3h)
! case 9
    else
        gggHHPhi=(qvtess(phi+d3phi,flam,h+d3h)-qvtess(phi-d3phi,flam,h+d3h) &
            -2.q0*qvtess(phi+d3phi,flam,h)+2.q0*qvtess(phi-d3phi,flam,h) &
            +qvtess(phi+d3phi,flam,h-d3h)-qvtess(phi-d3phi,flam,h-d3h))/(2.q0*d3phi*d3h*d3h)
    endif
else
    gggHHPhi=(qvtess(phi+d3phi,flam,h+d3h)-qvtess(phi-d3phi,flam,h+d3h) &
            -2.q0*qvtess(phi+d3phi,flam,h)+2.q0*qvtess(phi-d3phi,flam,h) &
            +qvtess(phi+d3phi,flam,h-d3h)-qvtess(phi-d3phi,flam,h-d3h))/(2.q0*d3phi*d3h*d3h)
endif
gggHHPhi=(gggHHPhi-2.q0*ggPhiH)/r
!
!   Compute gggHHLam
!
if(phi.ge.PhiS.and.phi.le.PhiN) then
! case 1
    if(((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
        .and.((h.gt.HB.and.h-d3h.lt.HB) &
        .or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggHHLam=(qvtess(phi,flam+2.q0*d3lam,h+3.q0*d3h)-4.q0*qvtess(phi,flam+d3lam,h+3.q0*d3h) &
            +3.q0*qvtess(phi,flam,h+3.q0*d3h)-4.q0*qvtess(phi,flam+2.q0*d3lam,h+2.q0*d3h) &
            +16.q0*qvtess(phi,flam+d3lam,h+2.q0*d3h)-12.q0*qvtess(phi,flam,h+2.q0*d3h) &
            +5.q0*qvtess(phi,flam+2.q0*d3lam,h+d3h)-20.q0*qvtess(phi,flam+d3lam,h+d3h) &
            +15.q0*qvtess(phi,flam,h+d3h)-2.q0*qvtess(phi,flam+2.q0*d3lam,h) &
            +8.q0*qvtess(phi,flam+d3lam,h)-6.q0*v)/(2.q0*d3lam*d3h*d3h)
! case 2
    elseif(((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
        .and.((h.gt.HB.and.h-d3h.lt.HB) &
        .or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggHHLam=-(qvtess(phi,flam-2.q0*d3lam,h+3.q0*d3h)-4.q0*qvtess(phi,flam-d3lam,h+3.q0*d3h) &
            +3.q0*qvtess(phi,flam,h+3.q0*d3h)-4.q0*qvtess(phi,flam-2.q0*d3lam,h+2.q0*d3h) &
            +16.q0*qvtess(phi,flam-d3lam,h+2.q0*d3h)-12.q0*qvtess(phi,flam,h+2.q0*d3h) &
            +5.q0*qvtess(phi,flam-2.q0*d3lam,h+d3h)-20.q0*qvtess(phi,flam-d3lam,h+d3h) &
            +15.q0*qvtess(phi,flam,h+d3h)-2.q0*qvtess(phi,flam-2.q0*d3lam,h) &
            +8.q0*qvtess(phi,flam-d3lam,h)-6.q0*v)/(2.q0*d3lam*d3h*d3h)
! case 3
    elseif(((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
        .and.((h.lt.HB.and.h+d3h.gt.HB) &
        .or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggHHLam=(qvtess(phi,flam+2.q0*d3lam,h-3.q0*d3h)-4.q0*qvtess(phi,flam+d3lam,h-3.q0*d3h) &
            +3.q0*qvtess(phi,flam,h-3.q0*d3h)-4.q0*qvtess(phi,flam+2.q0*d3lam,h-2.q0*d3h) &
            +16.q0*qvtess(phi,flam+d3lam,h-2.q0*d3h)-12.q0*qvtess(phi,flam,h-2.q0*d3h) &
            +5.q0*qvtess(phi,flam+2.q0*d3lam,h-d3h)-20.q0*qvtess(phi,flam+d3lam,h-d3h) &
            +15.q0*qvtess(phi,flam,h-d3h)-2.q0*qvtess(phi,flam+2.q0*d3lam,h) &
            +8.q0*qvtess(phi,flam+d3lam,h)-6.q0*v)/(2.q0*d3lam*d3h*d3h)
! case 4
    elseif(((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
        .and.((h.lt.HB.and.h+d3h.gt.HB) &
        .or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggHHLam=-(qvtess(phi,flam-2.q0*d3lam,h-3.q0*d3h)-4.q0*qvtess(phi,flam-d3lam,h-3.q0*d3h) &
            +3.q0*qvtess(phi,flam,h-3.q0*d3h)-4.q0*qvtess(phi,flam-2.q0*d3lam,h-2.q0*d3h) &
            +16.q0*qvtess(phi,flam-d3lam,h-2.q0*d3h)-12.q0*qvtess(phi,flam,h-2.q0*d3h) &
            +5.q0*qvtess(phi,flam-2.q0*d3lam,h-d3h)-20.q0*qvtess(phi,flam-d3lam,h-d3h) &
            +15.q0*qvtess(phi,flam,h-d3h)-2.q0*qvtess(phi,flam-2.q0*d3lam,h) &
            +8.q0*qvtess(phi,flam-d3lam,h)-6.q0*v)/(2.q0*d3lam*d3h*d3h)
! case 5
    elseif(((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE)) &
        .and.((h.gt.HB.and.h-d3h.lt.HB)&
        .or.(h.gt.HT.and.h-d3h.lt.HT))) then
        gggHHLam=(-qvtess(phi,flam+d3lam,h+3.q0*d3h)+qvtess(phi,flam-d3lam,h+3.q0*d3h) &
            +4.q0*qvtess(phi,flam+d3lam,h+2.q0*d3h)-4.q0*qvtess(phi,flam-d3lam,h+2.q0*d3h) &
            -5.q0*qvtess(phi,flam+d3lam,h+d3h)+5.q0*qvtess(phi,flam-d3lam,h+d3h) &
            +2.q0*qvtess(phi,flam+d3lam,h)-2.q0*qvtess(phi,flam-d3lam,h))/(2.q0*d3lam*d3h*d3h)
! case 6
    elseif(((flam+d3lam.lt.FlamW).or.(flam-d3lam.gt.FlamE) &
        .or.(flam-d3lam.gt.FlamW.and.flam+d3lam.lt.FlamE)) &
        .and.((h.lt.HB.and.h+d3h.gt.HB)&
        .or.(h.lt.HT.and.h+d3h.gt.HT))) then
        gggHHLam=(-qvtess(phi,flam+d3lam,h-3.q0*d3h)+qvtess(phi,flam-d3lam,h-3.q0*d3h) &
            +4.q0*qvtess(phi,flam+d3lam,h-2.q0*d3h)-4.q0*qvtess(phi,flam-d3lam,h-2.q0*d3h) &
            -5.q0*qvtess(phi,flam+d3lam,h-d3h)+5.q0*qvtess(phi,flam-d3lam,h-d3h) &
            +2.q0*qvtess(phi,flam+d3lam,h)-2.q0*qvtess(phi,flam-d3lam,h))/(2.q0*d3lam*d3h*d3h)
! case 7
    elseif(((flam.gt.FlamW.and.flam-d3lam.lt.FlamW) &
        .or.(flam.gt.FlamE.and.flam-d3lam.lt.FlamE)) &
        .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
        gggHHLam=(-qvtess(phi,flam+2.q0*d3lam,h+d3h)+4.q0*qvtess(phi,flam+d3lam,h+d3h) &
            -3.q0*qvtess(phi,flam,h+d3h)+2.q0*qvtess(phi,flam+2.q0*d3lam,h) &
            -8.q0*qvtess(phi,flam+d3lam,h)+6.q0*v &
            -qvtess(phi,flam+2.q0*d3lam,h-d3h)+4.q0*qvtess(phi,flam+d3lam,h-d3h) &
            -3.q0*qvtess(phi,flam,h-d3h))/(2.q0*d3lam*d3h*d3h)
! case 8
    elseif(((flam.lt.FlamW.and.flam+d3lam.gt.FlamW) &
        .or.(flam.lt.FlamE.and.flam+d3lam.gt.FlamE)) &
        .and.((h+d3h.lt.HB).or.(h-d3h.gt.HT) &
        .or.(h-d3h.gt.HB.and.h+d3h.lt.HT))) then
        gggHHLam=-(-qvtess(phi,flam-2.q0*d3lam,h+d3h)+4.q0*qvtess(phi,flam-d3lam,h+d3h) &
            -3.q0*qvtess(phi,flam,h+d3h)+2.q0*qvtess(phi,flam-2.q0*d3lam,h) &
            -8.q0*qvtess(phi,flam-d3lam,h)+6.q0*v &
            -qvtess(phi,flam-2.q0*d3lam,h-d3h)+4.q0*qvtess(phi,flam-d3lam,h-d3h) &
            -3.q0*qvtess(phi,flam,h-d3h))/(2.q0*d3lam*d3h*d3h)
! case 9
    else
        gggHHLam=(qvtess(phi,flam+d3lam,h+d3h)-qvtess(phi,flam-d3lam,h+d3h) &
            -2.q0*qvtess(phi,flam+d3lam,h)+2.q0*qvtess(phi,flam-d3lam,h) &
            +qvtess(phi,flam+d3lam,h-d3h)-qvtess(phi,flam-d3lam,h-d3h))/(2.q0*d3lam*d3h*d3h)
    endif
else
    gggHHLam=(qvtess(phi,flam+d3lam,h+d3h)-qvtess(phi,flam-d3lam,h+d3h) &
            -2.q0*qvtess(phi,flam+d3lam,h)+2.q0*qvtess(phi,flam-d3lam,h) &
            +qvtess(phi,flam+d3lam,h-d3h)-qvtess(phi,flam-d3lam,h-d3h))/(2.q0*d3lam*d3h*d3h)
endif
gggHHLam=(gggHHLam/c-2.q0*ggLamH)/r
!
!
!   Compute gggHHH
!
if((phi.ge.PhiS.and.phi.le.PhiN).and.(flam.ge.FlamW.and.flam.le.FlamE)) then
    if((h.gt.HB.and.h-d3h.lt.HB).or.(h.gt.HT.and.h-d3h.lt.HT)) then
        gggHHH=(-3.q0*qvtess(phi,flam,h+4.q0*d3h)+14.q0*qvtess(phi,flam,h+3.q0*d3h) &
            -24.q0*qvtess(phi,flam,h+2.q0*d3h)+18.q0*qvtess(phi,flam,h+d3h) &
            -5.q0*v)/(2.q0*d3h*d3h*d3h)
    elseif((h.lt.HB.and.h+d3h.gt.HB).or.(h.lt.HT.and.h+d3h.gt.HT)) then
        gggHHH=(3.q0*qvtess(phi,flam,h-4.q0*d3h)-14.q0*qvtess(phi,flam,h-3.q0*d3h) &
            +24.q0*qvtess(phi,flam,h-2.q0*d3h)-18.q0*qvtess(phi,flam,h-d3h) &
            +5.q0*v)/(2.q0*d3h*d3h*d3h)
    else
        gggHHH=(qvtess(phi,flam,h+2.q0*d3h)-2.q0*qvtess(phi,flam,h+d3h) &
            +2.q0*qvtess(phi,flam,h-d3h)-qvtess(phi,flam,h-2.q0*d3h) &
            )/(2.q0*d3h*d3h*d3h)
    endif
else
    gggHHH=(qvtess(phi,flam,h+2.q0*d3h)-2.q0*qvtess(phi,flam,h+d3h) &
            +2.q0*qvtess(phi,flam,h-d3h)-qvtess(phi,flam,h-2.q0*d3h) &
            )/(2.q0*d3h*d3h*d3h)
endif
!
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine qqde(f,a,b,delta0,s)
!
!   The subroutine to compute a given integral by the quadruple precision
!    double exponential (DE) quadrature rule.
!
!   Reference: "Accurate computation of gravitational field of a tesseroid"
!   authored by T. Fukushima and submitted to J. Geodesy (2017)
!
!    Input function/variables:
!         f(t): a function to be integrated 
!         a: the lower end point of the integration interval
!         b: the upper end point of the integration interval
!         delta0: the relative error tolerance (typically 1d-15)
!    Output variable:
!         s: the definite integral, s=int_a^b f(t) dt
!         (Note) "errd" is the latest relative error estimate o "s"
!
real*16 f,a,b,delta0,s
integer MMAX
real*16 SAFETY,TWOPI
parameter (MMAX=4092,SAFETY=10.q0)
parameter (TWOPI=6.283185307179586476925286766559q0)
integer m
real*16 delta,hmax,deltax,factor,deltah,h0,eph,emh,deltat,apb,bma,sr,h 
real*16 sprev,srprev,t,ep,em,eppem,xw,xa,wg,fa,fb,fapfb,errt 
real*16 errh,errd
delta=max(1.q-33,min(0.01q0,delta0))
if(delta.gt.1.q-4) then
    hmax=2.75q0+(-log10(delta))*0.75q0
elseif(delta.gt.1.q-6) then
    hmax=3.75q0+(-log10(delta))*0.5q0
elseif(delta.gt.1.q-10) then
    hmax=5.25q0+(-log10(delta))*0.25q0
else
    hmax=6.5q0+(-log10(delta))*0.125q0
endif
!hmax=sqrt(-2.q0*log10(delta))+2.5q0
!write(*,*) "(qqde) hmax=",hmax
deltax=SAFETY*delta
factor=1.q0-log(deltax)
deltah=sqrt(deltax)
h0=hmax/factor
eph=exp(h0)
emh=1.q0/eph
deltat=exp(-emh*factor)
apb=a+b
bma=b-a
sr=f(apb*0.5q0)*(bma*0.25q0)
s=sr*(2.q0*TWOPI)
err=abs(s)*deltat
h=2.q0*h0
m=1
1 continue
    sprev=s
    srprev=sr
    t=h*0.5q0
2 continue
        em=exp(t)
        ep=TWOPI*em
        em=TWOPI/em
3 continue
            eppem=ep+em
            xw=1.q0/(1.q0+exp(ep-em))
            xa=bma*xw
            wg=xa*(1.q0-xw)
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
        errd=1.q0+2.q0*errh
    else
        errd=h*(abs(s-2.q0*sprev)+4.q0*abs(sr-2.q0*srprev))
    endif
    h=h*0.5q0
    m=m*2
    if(errd.gt.errh.and.m.lt.MMAX) goto 1
!if(errd.gt.errh) write(*,*)"(qqde) Required too many grids: >", MMAX
s=s*h
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine qqde1(f,a,b,delta0,s)
!
!   The first copy of "qqde" to be used in computing the internal
!    line integral called by "qqde"
!
real*16 f,a,b,delta0,s
integer MMAX
real*16 SAFETY,TWOPI
parameter (MMAX=4092,SAFETY=10.q0)
parameter (TWOPI=6.283185307179586476925286766559q0)
integer m
real*16 delta,hmax,deltax,factor,deltah,h0,eph,emh,deltat,apb,bma,sr,h 
real*16 sprev,srprev,t,ep,em,eppem,xw,xa,wg,fa,fb,fapfb,errt 
real*16 errh,errd
delta=max(1.q-33,min(0.01q0,delta0))
if(delta.gt.1.q-4) then
    hmax=2.75q0+(-log10(delta))*0.75q0
elseif(delta.gt.1.q-6) then
    hmax=3.75q0+(-log10(delta))*0.5q0
elseif(delta.gt.1.q-10) then
    hmax=5.25q0+(-log10(delta))*0.25q0
else
    hmax=6.5q0+(-log10(delta))*0.125q0
endif
!hmax=sqrt(-2.q0*log10(delta))+2.5q0
!write(*,*) "(qqde1) hmax=",hmax
deltax=SAFETY*delta
factor=1.q0-log(deltax)
deltah=sqrt(deltax)
h0=hmax/factor
eph=exp(h0)
emh=1.q0/eph
deltat=exp(-emh*factor)
apb=a+b
bma=b-a
sr=f(apb*0.5q0)*(bma*0.25q0)
s=sr*(2.q0*TWOPI)
err=abs(s)*deltat
h=2.q0*h0
m=1
1 continue
    sprev=s
    srprev=sr
    t=h*0.5q0
2 continue
        em=exp(t)
        ep=TWOPI*em
        em=TWOPI/em
3 continue
            eppem=ep+em
            xw=1.q0/(1.q0+exp(ep-em))
            xa=bma*xw
            wg=xa*(1.q0-xw)
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
        errd=1.q0+2.q0*errh
    else
        errd=h*(abs(s-2.q0*sprev)+4.q0*abs(sr-2.q0*srprev))
    endif
    h=h*0.5q0
    m=m*2
    if(errd.gt.errh.and.m.lt.MMAX) goto 1
!if(errd.gt.errh) write(*,*)"(qqde1) Required too many grids: >", MMAX
s=s*h
return;end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*16 function qlog1p(x0)
!
!   The quadruple precision Fortran function to compute precisely log1p(x),
!    a special logarithm function defined as log1p(x) = ln(1+x),
!    by the (13,13) rational minimax approximation evaluated by Estrin's scheme
!
!   Reference: "Accurate computation of gravitational field of a tesseroid"
!   authored by T. Fukushima and submitted to J. Geodesy (2017)
!
real*16 x0,x1,x2,x4,x8
real*16 a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13
real*16 b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13
parameter (a1=0.999999999999999999999999999999999989679405409q0)
parameter (a2=5.43952695470870990343440472736508508440853974q0)
parameter (a3=12.8687954156331531486576454765577285904854589q0)
parameter (a4=17.3813399893763161892491716455513072747241961q0)
parameter (a5=14.7937221247975315812946713068099324296210500q0)
parameter (a6=8.26554125834405953536667781138761190677622028q0)
parameter (a7=3.06422110888383600247671821221533143683683084q0)
parameter (a8=0.745288118102087966956230766831419822309727243q0)
parameter (a9=0.115047592893588673011032562622446366266055797q0)
parameter (a10=0.0105968202344748100801276898051500664263281933q0)
parameter (a11=0.000522466799810767549674343323841954134608619668q0)
parameter (a12=0.0000112217383140869479740242371068976645949743233q0)
parameter (a13=6.39672301036246645989908315040681256718633817q-8)
parameter (b1=5.93952695470870990343440472736506510635091212q0)
parameter (b2=15.5052255596541747670415145069133006422898887q0)
parameter (b3=23.4041104509671669382917939890824601741965693q0)
parameter (b4=22.6122505690735676039519980345629102098426624q0)
parameter (b5=14.6253640581969227356522409769061106689586552q0)
parameter (b6=6.43653279869637374357975617245086508089670491q0)
parameter (b7=1.92137412606261539134961029284884346128770878q0)
parameter (b8=0.381097234341555148065443000297505125713856247q0)
parameter (b9=0.0482176170045436717981339396150227533590537972q0)
parameter (b10=0.00362957310087890393686248940312126067485435786q0)
parameter (b11=0.000144194286147019122572434476322438425343273856q0)
parameter (b12=2.40586135243166027564864363686803405251652763q-6)
parameter (b13=9.53384644636009246025258296194840899175886989q-9)
if(x0.lt.-0.5q0.or.x0.gt.1.q0) then
    qlog1p=log(1.q0+x0)
elseif(x0.lt.0.q0) then
    x1=-x0/(1.q0+x0)
    x2=x1*x1; x4=x2*x2; x8=x4*x4
    qlog1p=-x1*((((a1+x1*a2)+x2*(a3+x1*a4)) &
        +x4*((a5+x1*a6)+x2*(a7+x1*a8))) &
        +x8*(((a9+x1*a10)+x2*(a11+x1*a12))+x4*a13)) &
        /((((1.q0+x1*b1)+x2*(b2+x1*b3)) &
        +x4*((b4+x1*b5)+x2*(b6+x1*b7))) &
        +x8*(((b8+x1*b9)+x2*(b10+x1*b11))+x4*(b12+x1*b13)))
else
    x1=x0
    x2=x1*x1; x4=x2*x2; x8=x4*x4
    qlog1p=x1*((((a1+x1*a2)+x2*(a3+x1*a4)) &
        +x4*((a5+x1*a6)+x2*(a7+x1*a8))) &
        +x8*(((a9+x1*a10)+x2*(a11+x1*a12))+x4*a13)) &
        /((((1.q0+x1*b1)+x2*(b2+x1*b3)) &
        +x4*((b4+x1*b5)+x2*(b6+x1*b7))) &
        +x8*(((b8+x1*b9)+x2*(b10+x1*b11))+x4*(b12+x1*b13)))
endif
return; end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    output of "xqtessgc"
!      (Note) It took 167.150 s at a PC with an Intel Core i5-10400 running at 2.90 GHz
!      The ifort version is
!           Intel(R) Fortran Intel(R) 64 Compiler Classic for applications running on 
!            Intel(R) 64, Version 2021.4.0 Build 20210910_000000
!            Copyright (C) 1985-2021 Intel Corporation.  All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       # qdelta=  1.0000000E-33
!     # R0,HT,HB=  6.3800000E+03  1.0000000E+01 -4.0000000E+01
!    # PhiN,PhiS=  2.8000000E+01  2.7000000E+01
! # LamdaE,LamdaW  8.7000000E+01  8.6000000E+01
!   # phi,lambda=  2.7500000E+01  8.6500000E+01
!             # H                                            V
!               #                                         gPhi                                      gLambda                                           gH
!               #                                  GammaPhiPhi                               GammaPhiLambda                                    GammaPhiH
!               #                            GammaLambdaLambda                                 GammaLambdaH                                      GammaHH
!               #                                 gggPhiPhiPhi                                 gggPhiPhiLam                                   gggPhiPhiH
!               #                                   gggPhiLamH
!               #                                 gggLamLamPhi                                 gggLamLamLam                                   gggLamLamH
!               #                                     gggHHPhi                                     gggHHLam                                       gggHHH
!   0.0000000E+00    3.48325322193870636629659255498939089E-04
!                   -2.46036392570567795088050350602583788E-09    0.00000000000000000000000000000000000E+00   -2.85123442878667527113560010904585051E-06
!                   -5.08904505133959948029056972137904267E-08    0.00000000000000000000000000000000000E+00    2.18878445224545712153812579729099966E-11
!                   -6.30472396588766653064121274000352571E-08    0.00000000000000000000000000000000000E+00   -1.94784969155936899128098553983335210E-07
!                    4.16924272027821873588140368362987521E-12    0.00000000000000000000000000000000000E+00    4.30622414210216864580768845473234987E-10
!                    0.00000000000000000000000000000000000E+00
!                   -5.33271611569043705384911273479679979E-12    0.00000000000000000000000000000000000E+00    6.15508712693833669794812126628998345E-10
!                    1.16347339531572460448684289969744009E-12    0.00000000000000000000000000000000000E+00   -1.04613112690392038979335927185562668E-09
!  The total time is   167.149618000000      seconds.
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!