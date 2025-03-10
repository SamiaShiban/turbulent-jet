clear all
clc
U0=0.194;
D=0.006;
rho0=998.3;
alpha=0.0806;
rhoA=998.3;
phi0=0;
zsea=0.192;
cp0=10000;
ws00=0.0042;
n=47511;
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
x3=zeros(1,n/6);
bz3=zeros(1,n/6);
Q0=pi*D^2/4*U0
P0=3.1569*10^(-6);
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
 bz1(1)=0.05+b1;
 bz2(1)=0.05-b1;
 bz3(1)=bz2(1);
 x3(1)=x1;
Q=U0*pi*D^2/4;
H0=Q*u0;
M0=Q*U0;
B0=Q*gj;
LM=M0^0.75/B0^0.5;
Lj=(2*alpha)^(-0.5)*(pi)^(-0.25)*(B0)^(-0.5)*H0^0.75;
Qj=(pi)^(0.25)*(2*alpha)^0.5*(B0)^(-0.5)*H0^(5/4);
Vj=(2*alpha)^(0.5)*(pi)^(-0.25)*(B0)^(0.5)*H0^(-0.25);
St(1)=0;

k=1;
while x(k)<0.71
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
bz1(k+1)=0.05+b(k+1);
bz2(k+1)=0.05-b(k+1);
if bz2(k+1)>0
    bz3(k+1)=bz2(k+1);
    x3(k+1)=x(k+1);
end
Q(k)=U(k)*pi*b(k)^2;
 Ealpha(k+1)=alpha*rhoA*2*pi*b(k+1)*h(k+1)*U(k+1);

if ws00>alpha*U(k)
    Pl(k+1)=Pl(k)-2*(2*b(k))/Q(k)*(ws00-alpha*U(k))*((1-alpha*U(k)/ws00))*deltaS(k);
   ru=4*(b(k)-D)/Q(k)*(ws00*cos(phi(k))-alpha*U(k))*((1-alpha*U(k)/(ws00*cos(phi(k)))))*deltaS(k);
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
if abs(alpha*cos(phi(k))*U(k)-ws00)<0.000002
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
xe=[0.008 0.024 0.04 0.056 0.072 0.088 0.104 0.120 0.136 0.152 0.168 0.184 0.200 0.216 0.232 0.248 0.264 0.280 0.296 0.312 0.328 0.344 0.360 0.376 0.392 0.408 0.424 0.44 0.456 0.472 0.488 0.504 0.520 0.536 0.552 0.568 0.584 0.600 0.616 0.632 0.648 0.664 0.68];
Se=[0 0 0 0 4 9 27 23 46 77 89 103 134 170 150 192 142 138 123 140 162 174 145 128 61 63 80 69 68 56 62 41 47 40 47 23 22 33 20 11 9 10 5];

figure(1)
 plot(x,S*100000)
 hold on
 plot(xe,Se*5.7332*10^(-6)/(0.016*53.42)*1000,'m*')
 hold off
xlabel('{\it x} (m)')
ylabel('Dépôt (gr/m.s)')
title ('Cas 9')
legend('Modèle proposé',' Données expérimentales (cas 9) [Bleninger et al. 2002]',1)

figure(2)
plot(x,bz1)
axis([-0.05 0.8 -0.05 0.25])
xlabel('{\it x} (m)')
ylabel('{\it z} (m)')
hold on
line([0;0.72],[0.05;0.05],'linestyle','-.')
line([-0.05;0.8],[0.197;0.197],'linestyle',':')
line([-0.05;0.8],[0;0],'linestyle','--','color','k')
title ('Cas 9')
plot(x3,bz3)
legend('Limites du jet','Axe du jet','Surface de l'' eau', 'Fond du bassin',1)

 hold off

figure(3)
plot(DELTAS,P/P0*100)
hold on
plot(DELTAS,St/P0*100,'--')
legend('Flux solide dans le jet {\it P/P_0 } (%)','Déposition cumulée (%)')
title('Cas 9')
hold off
xlabel('Axe du jet (m)')
ylabel('Flux solide normalisé (%) ')
axis([-0.05 0.8 -5 105])


figure(4)
plot(DELTAS,U)
xlabel('Axe du jet (m)')
ylabel('Vitesse du jet (m/s)')
legend('Vitesse axiale horizontale du jet ({\itV})' )
axis([-0.05 0.8 0 0.2])
title ('Cas 9')


