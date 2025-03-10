clear all
clc
U0=1.04;
D=0.0115;
rho0=1000;
rhoA=1020;
phi0=0;
zsea=0.747;
ws00=0.0299;
n=24906;
alpha=0.0806;
cp0=10000;
con0=100;
g=9.81;
m=zeros(1,n);
deltam=zeros(1,n);
U=zeros(1,n);
h=zeros(1,n);
b=zeros(1,n);
x=zeros(1,n);
AX=zeros(1,n);
H1=zeros(1,n);
H2=zeros(1,n);
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
DELTAS=zeros(1,n);
nn=zeros(1,n);
time=zeros(1,n);
bx1=zeros(1,n);
bx2=zeros(1,n);
bz1=zeros(1,n);
bz2=zeros(1,n);
r=zeros(1,n);
P0=1;
b0=0.5*D;
h0=0.5*D;
deltat=0.1*h0/U0
x0=0;
z0=0;
u0=U0*cos(phi0);
w0=U0*sin(phi0);
gj=((rhoA-rho0)/rhoA)*9.81;
fr=U0/(gj*D)^.5
Re=U0*D*1000000
Ealpha0=alpha*rhoA*2*pi*b0*h0*U0;
m0= rho0*pi*b0^2*h0;
deltam0=Ealpha0*deltat;
m1=m0+deltam0;
rho1=(m0*rho0+deltam0*rhoA)/m1;
u1=m0*u0/m1;
w1=m0*w0/m1+((rhoA-rho1)/rhoA)*9.81*deltat; 
U1=(u1^2+w1^2)^0.5;
h1=U1/U0*h0;
b1=(m1/(rho1*pi*h1))^0.5;
phi1=atan(w1/u1);
x1=x0+u0*deltat;
z1=z0+w0*deltat;
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
alphav(1)=alpha;
nn(1)=U1-U0;
bx1(1)=x1+b1*sin(phi1);
bz1(1)=z1-b1*cos(phi1);
bx2(1)=x1-b1*sin(phi1);
bz2(1)=z1+b1*cos(phi1);
Q(1)=U0*pi*D^2/4;
H0=Q(1)*u0;
M0=Q(1)*U0;
B0=Q(1)*gj;
LM=M0^0.75/B0^0.5
Lj=(2*alpha)^(-0.5)*(pi)^(-0.25)*(B0)^(-0.5)*H0^0.75;
St(1)=0;
r(1)=1-alpha*U1*cos(phi1)/ws00;
k=1;
while bz2(k)<zsea
m(k+1)=m(k)+Ealpha(k)*deltat;
deltam(k)=Ealpha(k)*deltat;
rho(k+1)=(m(k)*rho(k)+deltam(k)*rhoA)/m(k+1);
u(k+1)=m(k)*u(k)/m(k+1);
w(k+1)=m(k)*w(k)/m(k+1)+((rhoA-rho(k+1))/rhoA)*9.81*deltat;
U(k+1)=(u(k+1)^2+w(k+1)^2)^0.5;
h(k+1)=(U(k+1)/U(k))*h(k);
b(k+1)=(m(k+1)/(rho(k+1)*pi*h(k+1)))^0.5;
phi(k+1)=atan(w(k+1)/u(k+1));
x(k+1)=x(k)+u(k)*deltat;
xj(k+1)=x(k+1)/Lj;
z(k+1)=z(k)+w(k)*deltat;
bx1(k+1)=x(k+1)+b(k+1)*sin(phi(k+1));
bz1(k+1)=z(k+1)-b(k+1)*cos(phi(k+1));
bx2(k+1)=x(k+1)-b(k+1)*sin(phi(k+1));
bz2((k+1))=z(k+1)+b(k+1)*cos(phi(k+1));
r(k+1)=1-alpha*U(k+1)*cos(phi(k+1))/ws00;
deltaS(k+1)=U(k+1)*deltat;
Q(k+1)=U(k+1)*pi*b(k+1)^2;
Ealpha(k+1)=alpha*rhoA*2*pi*b(k+1)*h(k+1)*U(k+1);
if ws00>alpha*cos(phi(k))*U(k);
     Pl(k+1)=Pl(k)-1.5*(ws00-alpha*U(k)*cos(phi(k)))*2*b(k)/Q(k)*(1-alpha*U(k)*cos(phi(k))/ws00)*deltaS(k);
    ru=(ws00-alpha*U(k)*cos(phi(k)))*(2*b(k)-D)/Q(k)*(1-alpha*U(k)*cos(phi(k))/ws00)*deltaS(k);
    r(k+1)=1-alpha*U(k+1)*cos(phi(k+1))/ws00;
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
    r(k+1)=0;
end
St(k+1)=St(k)+Ss(k);
DELTAS(k+1)=DELTAS(k)+deltaS(k);
time(k+1)=time(k)+deltat;
if abs(alpha*cos(phi(k))*U(k)-ws00)<0.00001
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
Rs=b(k) 
Qs=Us*pi*Rs^2
T=deltat*k/60
Ps=P(k)
Cs=Ps/Qs
stt=St(k);
ST=St(k)/P0*100
steps=k
Distance=DELTAS(k)
xe=[0.015 0.075 0.125 0.175 0.225 0.275 0.325 0.375 0.4 0.425 0.525 0.550 0.625 0.7 0.725 0.8 0.81 0.875 0.92 0.98 1.03 1.1 1.14 1.18];
Se=[0 0.25 0.7 0.9 1.2 1.4 1.6 1.8 1.9 1.8 1.5 1.4 1.1 1 1 0.75 0.55 0.45 0.4 0.25 0.15 0.05 0 0];

figure(1)
plot(x,S/stt,'LineWidth',1.01)
hold on
plot(xe,Se,'m*')
hold off
xlabel('{\it x} (m)')
ylabel('PDF f_L(x) (m^{-1})')
axis([-0.05 1.2 0 2.5])
title ('Cas 1')
legend('Modèle proposé',' Données expérimentales (cas 1) [Cuthbertson et Davies 2008]',1)

figure(2)
plot(x,z,'-.','LineWidth',1.01)
axis([-0.05 1.2 -0.1 1])
xlabel('{\it x} (m)')
ylabel('{\it z} (m)')
hold on
plot(bx1,bz1,'LineWidth',1.01)
line([-0.05;1.5],[zsea;zsea],'LineStyle',':','LineWidth',1.01)
plot(bx2,bz2,'LineWidth',1.01)
legend('Axe du jet','Limites du jet','Surface de l'' eau',1)
title ('Cas 1')
hold off


figure(3)
plot(DELTAS,U,'LineWidth',1.01)
xlabel('Axe du jet (m)')
axis([-0.05 1.5 -0.05 1.2])
hold on
plot(DELTAS,u,'--','LineWidth',1.01)
plot(DELTAS,w,':','LineWidth',1.01)
ylabel('Vitesse du jet (m/s)')
legend('Vitesse axiale du jet ({\it V}  )',' Vitesse horizontale du jet ({\it u} )', ' Vitesse verticale du jet ({\it w} )',1)
title('Cas 1')
hold off

figure(4)
plot(DELTAS,P/P0*100,'LineWidth',1.01)
hold on
plot(DELTAS,St/P0*100,'--','LineWidth',1.01)
legend('Flux solide dans le jet {\it P/P_0 } (%)','Déposition cumulée (%)')
title('Cas 1')
hold off
xlabel('Axe du jet (m)')
ylabel('Flux solide normalisé (%)')
axis([-0.05 1.5 -5 105])

figure(5)
plot(x,S/stt,'LineWidth',1.01)
hold on
plot(x,Q*100,'r')
plot(x,U,'g')
plot(x,b*10,'y')
hold off
legend('depo','Q','U','b')
