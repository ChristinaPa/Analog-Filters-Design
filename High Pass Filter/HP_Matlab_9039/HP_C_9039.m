%% Pardali Christina AEM 9039/September 2020 
%% Highpass Chebyshev Filter 

clear;
clc;

%% AEM 

a1 = 9;
a2 = 0;
a3 = 3;
a4 = 9;

%% Statheres filtrou me basi thn ekfwnisi

m=0;
fp=(3+m)*1000;
fs=(fp/1.8);

amin=25+(a3*(4/9));
amax=0.5+(a4*(0.25/9));

wp=2*pi*fp;
ws=2*pi*fs;

%% Transformation of specifications

Wp=1;
Ws=wp/ws;

%% Sintelestes HP filtrou kai taksi filtrou

e = sqrt(10^(amax/10)-1);
n_decimal = (acosh(((10^(amin/10)-1)/(10^(amax/10)-1))^(1/2)))/acosh(Ws);

n= ceil(n_decimal);
a = (asinh(1/e))/n;

%% Sixnotita 3dB

Whp = cosh((acosh((10^(amax/10)-1)^(-1/2)))/n);

%% Butterworth angles

ps_k1 = 22.5;
ps_k2 = - 22.5;
ps_k3 = 67.5;
ps_k4 = -67.5;

%% Poloi protipis sinartisis

p_1 = -sinh(a)*cosd(ps_k1) + j*cosh(a)*sind(ps_k1);
p_2 = -sinh(a)*cosd(ps_k2) + j*cosh(a)*sind(ps_k2);
p_3 = -sinh(a)*cosd(ps_k3) + j*cosh(a)*sind(ps_k3);
p_4 = -sinh(a)*cosd(ps_k4) + j*cosh(a)*sind(ps_k4);

%% w0, Q polwn chebyshev

w0_12 = sqrt((real(p_1))^2+(imag(p_1))^2);
Q_12 = w0_12 /(2 * abs(real(p_1)));

w0_34 = sqrt((real(p_3))^2+(imag(p_3))^2);
Q_34 = w0_34/(2*abs(real(p_3)));

%% Antistrofi polwn

w_hp = wp/Whp;

w12 = wp/w0_12;
w34 = wp/w0_34;

%% Piknotis

C=10^(-8);

%/-------------------Prwti monada-----------------------/

%% Ilopoiisi monadas 1
R11=1;
R12=R11;
C11=1;
C12=C11;
r11=1;
r12=2-(1/Q_12);
k1=3-(1/Q_12);

%% Klimakopoisi monadas 1 

kf1=w12;
C11new=C;
C12new=C;
km1=C11/(C11new*kf1);
r11new=r11*km1;
r12new=r12*km1;
R11new=R11*km1;
R12new=R12*km1;

%/-----------------Deuteri monada---------------------/

%% Ilopoiisi monadas 2 

R21=1;
R22=R21;
C21=1;
C22=C21;
r21=1;
r22=2-(1/Q_34);
k2=3-(1/Q_34);

%% Klimakopoisi monadas 2

kf2=w34;
C21new=C;
C22new=C;
km2=C21/(C21new*kf2);
r21new=r21*km2;
r22new=r22*km2;
R21new=R21*km2;
R22new=R22*km2;

%% Rithmisi kerdous

ktotal=k1*k2;
a_kerdos=(10^(0.5))/ktotal;

%% Ypologismos epimerous sunartisewn metaforas 

T1=tf([k1 0 0 ], [1 w12/Q_12 w12^2]);
T2=tf([k2 0 0 ], [1 w34/Q_34 w34^2]);
T_total=a_kerdos*T1*T2;
Invsys=inv(T_total);

%% Diaforetikes T gia kerdos 0dB

T1_1=tf([1 0 0 ], [1 w12/Q_12 w12^2]);
T2_2=tf([1 0 0 ], [1 w34/Q_34 w34^2]);

Invsys2=inv(T1_1*T2_2);

%% Diagrammata sinartisewn metaforas

plot_transfer_function( T1, [fp fs] );
plot_transfer_function( T2, [fp fs] );
plot_transfer_function(T_total, [fp fs]);

plot_transfer_function(Invsys,[fp fs]);

plot_transfer_function(Invsys2,[fp fs]);

%% Sxediasmos diagrammatos Bode

ltiview({'bodemag'}, T1, T2 , T_total);

%% Sixnotites simatos eisodou

f11=(0.4*ws)/(2*pi);
f12=(0.9*ws)/(2*pi);
f13=(1.4*wp)/(2*pi);
f14=(2.4*wp)/(2*pi);
f15=(4.5*wp)/(2*pi);

%% Fourier analysis

T =0.002 ;
f_s=200*(10^3);
dt = 1/f_s;
t = 0:dt:(T);

%% Dimiourgia kai emfanisi simatos eisodou

u1=cos(2*pi*f11*t)+0.5*cos(2*pi*f12*t)+cos(2*pi*f13*t)+0.7*cos(2*pi*f14*t)+0.5*cos(2*pi*f15*t);
figure
plot(u1)

%% Dimiourgia fasmatwn simatwn eisodou kai eksodou

N=T/dt;
figure
lsim(T_total,u1,t)
xt=lsim(T_total,u1,t);
figure
plot(t,xt)


n1=2^nextpow2(N);
xfourier= fft(xt,n1);
p2=abs(xfourier/n1);
p1=p2(1:n1/2+1);
p1(2:end-1)=2*p1(2:end-1);
f=200*10^3*(0:(n1/2))/n1;
figure
plot(f,p1)

nfft=n1;
y=fft(u1,nfft);
y = y(1:nfft/2); 
y_mag=abs(y);
f = (0:nfft/2-1)*f_s/nfft; 
figure
plot(f,y_mag)








