
clear all
clc
U0=1.29;
D=0.0115;
rho0=1000;
alpha=0.0806;
rhoA=1040;
phi0=0;
zsea=0.749;
cp0=10000;
ws00=0.0285;
n=21911;
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
phic=zeros(1,n);
DELTAS=zeros(1,n);
nn=zeros(1,n);
time=zeros(1,n);
bx1=zeros(1,n);
bx2=zeros(1,n);
bz1=zeros(1,n);
bz2=zeros(1,n);
Q0=pi*D^2/4*U0
P0=1;
b0=0.5*D;
h0=0.5*D;
deltat=0.1*h0/U0
x0=0;
z0=0;
u0=U0*cos(phi0);
w0=U0*sin(phi0);
gj=((rhoA-rho0)/rhoA)*9.81;
Ealpha0=alpha*rhoA*2*pi*b0*h0*U0;
m0= rho0*pi*b0^2*h0;
deltam0=Ealpha0*deltat;
fr=U0/(gj*D)^.5
Re=U0*D*1000000
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
Q(1)=Q1;
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
Q=U0*pi*D^2/4;
H0=Q*u0;
M0=Q*U0;
B0=Q*gj;
LM=M0^0.75/B0^0.5
Lj=(2*alpha)^(-0.5)*(pi)^(-0.25)*(B0)^(-0.5)*H0^0.75;
Qj=(pi)^(0.25)*(2*alpha)^0.5*(B0)^(-0.5)*H0^(5/4);
Vj=(2*alpha)^(0.5)*(pi)^(-0.25)*(B0)^(0.5)*H0^(-0.25);
St(1)=0;

k=1;
while bz2(k)<zsea
m(k+1)=m(k)+Ealpha(k)*deltat;
deltam(k)=Ealpha(k)*deltat;
rho(k+1)=(m(k)*rho(k)+deltam(k)*rhoA)/m(k+1);
u(k+1)=m(k)*u(k)/m(k+1);
w(k+1)=m(k)*w(k)/m(k+1)+((rhoA-rho(k+1))/rho(k+1))*9.81*deltat;
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
 deltaS(k+1)=U(k+1)*deltat;

Q(k)=U(k)*pi*b(k)^2;
Ealpha(k+1)=alpha*rhoA*2*pi*b(k+1)*h(k+1)*U(k+1);

if ws00>alpha*cos(phi(k))*U(k);
    Pl(k+1)=Pl(k)-(ws00-alpha*U(k)*cos(phi(k)))*(2*b(k))/Q(k)*(1-alpha*U(k)*cos(phi(k))/ws00)*deltaS(k);
    ru=(ws00-alpha*U(k)*cos(phi(k)))*(2*b(k))/Q(k)*(1-alpha*U(k)*cos(phi(k))/ws00)*deltaS(k);
    bu=b(k);
end
%if ws00*cos(phi(k))>alpha*U(k)*cos(phi(k))&&DELTAS(k)>LM;
%Pl(k+1)=Pl(k)-(ws00-alpha*U(k)*cos(phi(k)))*2*(b(k)-bu)/Q(k)*(1-alpha*U(k)*cos(phi(k))/ws00)*deltaS(k)-ru;
%end
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
Rs=b(k) 
Qs=Us*pi*Rs^2
T=deltat*k/60
Ps=P(k)
Cs=Ps/Qs
stt=St(k);
phis=phi(k);
ST=St(k)/P0*100
steps=k
Distance=DELTAS(k)
xe=[0.025 0.08 0.14 0.19 0.23 0.3 0.34 0.38 0.42 0.48 0.52 0.57 0.7 0.785 0.83 0.895 0.915 0.99 1.020 1.080 1.120 1.180];
Se=[0 0.17 0.35 0.6 0.95 1.08 1.2 1.46 1.6 1.55 1.45 1.45 1.5 1.34 1.1 0.67 0.34 0.08 0 0 0 0];

figure(1)
plot(x,S/stt,'LineWidth',1.01)
hold on
plot(xe,Se,'m*')
hold off
xlabel('{\it x} (m)')
ylabel('PDF f_L(x) (m^{-1})')
axis([-0.05 1.2 0 2.5])
title ('Cas 2')
legend('Modèle proposé',' Données expérimentales (cas 2) [Cuthbertson et Davies 2008]',1)

figure(2)
plot(x,z,'-.','LineWidth',1.01)
axis([-0.05 1.2 -0.1 1])
xlabel('{\it x} (m)')
ylabel('{\it z} (m)')
hold on
title ('Cas 2')
plot(bx1,bz1,'LineWidth',1.01)
line([-0.05;1.5],[zsea;zsea],'LineStyle',':','LineWidth',1.01)
plot(bx2,bz2,'LineWidth',1.01)
legend('Axe du jet','Limites du jet','Surface de l'' eau',1)
hold off

figure(3)
plot(DELTAS,U,'LineWidth',1.01)
xlabel('Axe du jet (m)')
ylabel('Vitesse du jet (m/s)')
axis([-0.05 1.5 -0.05 1.4])
title ('Cas 2')
hold on
plot(DELTAS,u,'--','LineWidth',1.01)
plot(DELTAS,w,':','LineWidth',1.01)
legend('Vitesse axiale du jet ({\it V}  )',' Vitesse horizontale du jet ({\it u} )', ' Vitesse verticale du jet ({\it w} )',1)
hold off

figure(4)
plot(DELTAS,P/P0*100,'LineWidth',1.01)
hold on
plot(DELTAS,St/P0*100,'--','LineWidth',1.01)
legend('Flux solide dans le jet {\it P/P_0 } (%)','Déposition cumulée (%)')
hold off
title ('Cas 2')
xlabel('Axe du jet (m)')
ylabel('Flux solide normalisé (%) ')
axis([-0.05 1.5 -5 105])






