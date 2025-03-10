clear all
clc
U0=0.86;
D=0.006;
rho0=997.2;
alpha=0.0806;
rhoA=997.2;
phi0=0;
zsea=0.45;
cp0=3840;
ws00=0.0266;
n=59982;
g=9.81;
m=zeros(1,n);
deltam=zeros(1,n);
U=zeros(1,n);
h=zeros(1,n);
b=zeros(1,n);
x=zeros(1,n);
deltaS=zeros(1,n);
Ealpha=zeros(1,n);
Q=zeros(1,n);
Pl=zeros(1,n);
P=zeros(1,n);
S=zeros(1,n);
Ss=zeros(1,n);
xj=zeros(1,n);
St=zeros(1,n);
rho=zeros(1,n);
u=zeros(1,n);
w=zeros(1,n);
phi=zeros(1,n);
z=zeros(1,n);
Pll=zeros(1,n);
phic=zeros(1,n);
DELTAS=zeros(1,n);
time=zeros(1,n);
bz1= zeros(1,n);
bz2=zeros(1,n);
x3=zeros(1,n/10);
bz3=zeros(1,n/10);
Q0=pi*D^2/4*U0
P0=Q0*cp0;
b0=0.5*D;
h0=0.5*D;
deltat=0.1*h0/U0
x0=0;
z0=0;
u0=U0*cos(phi0);
w0=U0*sin(phi0);
gj=((rhoA-rho0)/rho0)*9.81;
Ealpha0=alpha*rhoA*2*pi*b0*h0*U0;
m0= rho0*pi*b0^2*h0;
deltam0=Ealpha0*deltat;
m1=m0+deltam0;
rho1=(m0*rho0+deltam0*rhoA)/m1;
u1=m0*u0/m1;
w1=m0*w0/m1+((rhoA-rho1)/rho1)*9.81*deltat; 
U1=(u1^2+w1^2)^0.5;
h1=U1/U0*h0;
b1=(m1/(rho1*pi*h1))^0.5;
phi1=atan(w1/u1);
Q1=pi*b1^2/4*U1;
x1=x0+u0*deltat;
z1=z0+w0*deltat;
deltaS0=U0*deltat;
deltaS1=U1*deltat;
Ealpha(1)=Ealpha0;
m(1)=m0+Ealpha0*deltat;
rho(1)=rho1;
u(1)=u1;
w(1)=w1;
h(1)=h1;
x(1)=x1;
z(1)=z1;
b(1)=b1;
U(1)=U1;
phi(1)=phi1;
deltaS(1)= deltaS1;
deltam(1)=deltam0;
P(1)=P0;
Pl(1)=log(P0);
S(1)=0;
Ss(1)=0;
time(1)=deltat;
DELTAS(1)=deltaS1;
rb(1)=0;
wp(1)=w1;
 bz1(1)=0.15+b1;
 bz2(1)=0.15-b1;
Q=U0*pi*D^2/4;
H0=Q*u0;
M0=Q*U0;
B0=Q*gj;
LM=M0^0.75/B0^0.5;
Lj=(2*alpha)^(-0.5)*(pi)^(-0.25)*(B0)^(-0.5)*H0^0.75;
Qj=(pi)^(0.25)*(2*alpha)^0.5*(B0)^(-0.5)*H0^(5/4);
Vj=(2*alpha)^(0.5)*(pi)^(-0.25)*(B0)^(0.5)*H0^(-0.25);
St(1)=0;
bz3(1)=bz2(1);
 x3(1)=x1;
k=1;
while x(k)<0.8
m(k+1)=m(k)+Ealpha(k)*deltat;
deltam(k)=Ealpha(k)*deltat;
rho(k+1)=(m(k)*rho(k)+deltam(k)*rhoA)/m(k+1);
u(k+1)=m(k)*u(k)/m(k+1);
w(k+1)=m(k)*w(k)/m(k+1)+((rhoA-rho(k+1))/rho(k+1))*9.81*deltat;
U(k+1)=(u(k+1)^2)^0.5;
h(k+1)=(U(k+1)/U(k))*h(k);
b(k+1)=(m(k+1)/(rho(k+1)*pi*h(k+1)))^0.5;
phi(k+1)=atan(w(k+1)/u(k+1));
x(k+1)=x(k)+u(k)*deltat;
xj(k+1)=x(k+1)/Lj;
z(k+1)=z(k)+w(k)*deltat;
deltaS(k+1)=U(k+1)*deltat;
bz1(k+1)=0.15+b(k+1);
bz2(k+1)=0.15-b(k+1);
if bz2(k+1)>0
    bz3(k+1)=bz2(k+1);
    x3(k+1)=x(k+1);
end
Q(k)=U(k)*pi*b(k)^2;
 Ealpha(k+1)=alpha*rhoA*2*pi*b(k+1)*h(k+1)*U(k+1);

if ws00>alpha*U(k)
    Pl(k+1)=Pl(k)-2*(2*b(k))/Q(k)*(ws00-alpha*U(k))*((1-alpha*U(k)/ws00))*deltaS(k);
  
   bu=b(k);
end
P(k+1)=exp(Pl(k+1));
if P(k)>P(k+1)
Ss(k+1)=P(k)-P(k+1);
S(k+1)=Ss(k+1)/(deltaS(k+1)*cos(phi(k+1)));
end

if ws00*cos(phi(k))<alpha*U(k)
    Pl(k+1)=Pl(k);
    P(k+1)=P(k);  
    Ss(k+1)=0;
    S(k+1)=0;
end
St(k+1)=St(k)+Ss(k);
DELTAS(k+1)=DELTAS(k)+deltaS(k);
time(k+1)=time(k)+deltat;
if abs(alpha*cos(phi(k))*U(k)-ws00)<0.00002
dist=deltaS(k)
entvel=alpha*U(k)
slice=k
axvelo=U(k)
distD=deltaS(k)/D
axveloperc=U(k)/U0*100
angle=phi(k)*180/pi
end
k=k+1;
end

Smax=max(S)
xs=x(k)
zs=z(k)
Us=U(k)
phis=phi(k)
Rs=b(k) 
Qs=Us*pi*Rs^2
T=deltat*k
Ps=P(k)
Cs=Ps/Qs
stt=St(k);
ST=St(k)/P0*100
steps=k
Distance=DELTAS(k)
xe=[0.04 0.05 0.06 0.09 0.115 0.145 0.18 0.19 0.21 0.24 0.26 0.285 0.315 0.345 0.36 0.38 0.415 0.44 0.46 0.48 0.53 0.585 0.64 0.69 0.74 0.785];
Se=[0 0.03 0.175 0.33 0.415 0.415 0.375 0.34 0.305 0.24 0.215 0.2 0.14 0.115 0.09 0.075 0.05 0.045 0.035 0.02 0.015 0.005 0 0 0 0];

figure(1)
plot(x,S)
hold on
plot(xe,Se,'m*')
hold off
xlabel('{\it x} (m)')
ylabel('Dépôt (gr/m.s)')
title ('Cas 10')
legend('Modèle proposé',' Données expérimentales (cas 10) [Lee 2010]',1)
axis([0 0.8 0 0.8])
 
figure(2)
plot(x,bz1)
axis([-0.05 0.9 -0.05 0.6])
xlabel('{\it x} (m)')
ylabel('{\it z} (m)')
hold on
line([0;0.8],[0.15;0.15],'linestyle','-.')
line([-0.05;0.9],[0.45;0.45],'linestyle',':')
line([-0.05;0.9],[0;0],'linestyle','--','color','k')
title ('Cas 10')
plot(x3,bz3)
legend('Limites du jet','Axe du jet','Surface de l'' eau', 'Fond du bassin',1)
 hold off


figure(3)
plot(DELTAS,P/P0*100)
hold on
plot(DELTAS,St/P0*100,'--')
legend('Flux solide dans le jet {\it P/P_0 } (%)','Déposition cumulée (%)')
title('Cas 10')
hold off
xlabel('Axe du jet (m)')
ylabel('Flux solide normalisé (%) ')
axis([-0.05 0.9 -5 105])


figure(4)
plot(DELTAS,U)
xlabel('Axe du jet (m)')
ylabel('Vitesse du jet (m/s)')
legend('Vitesse axiale horizontale du jet ({\itV})' )
title('Cas 10')
axis([-0.05 0.9 0 1])


