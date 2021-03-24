%% Pardali Christina AEM 9039/September 2020 
%% Lowpass Butterworth Filter  

clear;
clc;

%% AEM 

a1 = 9 ;
a2 = 0;
a3 = 3;
a4 = 9;

%% Statheres filtrou me basi thn ekfwnisi

m = 2; 

fp = 0.6 *(3+m) * 1000;
fs = 2 * fp;
wp = 2 * pi * fp;
ws = 2 * pi * fs;

amin = 17.5 + (max(1, a4)-5)*0.5;
amax = 0.6 +(max(1, a3)-3)/10 ;

%% Taksi filtrou 
n = (log((10^(amin/10)-1)/((10^(amax/10)-1))))/(2*log(ws/wp));
n_round = ceil (n);

%% Ypologismos sixnotitas wo

wo = wp /(((10^(amax/10)-1))^(1/(2*n_round)));

%% Gia n=5 pou prokyptei, ta Q me basi tis shmeiwseis einai 

Q1 = 0.5;
Q2 = 0.62;
Q3 = 1.62;

%/-------------------Prwti monada-----------------------/

%% 1 monada-Prwti taksi-katwdiabato filtro

R11 = 1;
C11 = 1;

%% Klimakopoiisi  prwtis monadas

kf = wo;
Cn = 10^(-8);
km1 = C11/(kf*Cn);
C11_new = Cn;
R11_new = km1*R11;


%/-------------------Deuteri monada-----------------------/

%% 2 monada stratigiki 2 sallen key 

R21=1;
R22=1;
k2=1;
C21=2*Q2;
C22=1/(2*Q2);

%% klimakopoiisi deuteris monadas

kf2 = wo;
C21_NEW = Cn;
km2 = (C21/(C21_NEW*kf2));
C22_NEW = (1/(kf*km2))*C22;
R21_NEW = km2*R21;
R22_NEW = km2*R22;


%/-------------------Triti monada-----------------------/

%% 3 monada stratigiki 2 sallen key 

R31=1;
R32=1;
k3=1;
C31=2*Q3;
C32=1/(2*Q3);

%% klimakopoiisi tritis monadas

kf3 = wo;
C31_NEW = Cn;
km3 = (C31/(C31_NEW*kf3));
C32_NEW = (1/(kf*km3))*C32;
R31_NEW = km3*R31;
R32_NEW = km3*R32;


%% Ypologismos epimerous sunartisewn metaforas
T1 = tf(wo, [1, wo]);
T2 = tf([k2*wo^2], [1, wo/Q2, wo^2]);
T3 = tf([k3*wo^2], [1, wo/Q3, wo^2]);

Tlp = T1 * T2 *T3;

%% Diagrammata sinartisewn metaforas
plot_transfer_function(T1, [fp fs]);
plot_transfer_function(T2, [fp fs]);
plot_transfer_function(T3, [fp fs]);
plot_transfer_function(Tlp, [fp fs]); 

InvSys_new = inv(Tlp)
plot_transfer_function(InvSys_new, [fp fs]);

%% Sxediasmos digrammatos Bode
ltiview({'bodemag'}, T1, T2, T3, Tlp);

%% Fourier analysis

f_s = 200*10^3;
T = 0.002     ;
dt = 1/f_s;
t = 0:dt:(T);
x = sawtooth(2*pi*2000*t,0.5);

%% Dimiourgia kai emfanisi trigonikou simatos
figure 
plot(t,x)
title('triangle wave') 
xlabel('t (sec)')
ylabel('Amplitude') 

%% Dimiourgia fasmatwn simatwn eisodou kai eksodou
N = T/dt;
xt = lsim(Tlp,x,t);
figure
plot(t,xt)
n = 2^nextpow2(N);
xfourier = fft(xt,n);
p2 = abs(xfourier/n);
p1 = p2(1:n/2+1);
p1(2:end-1) = 2*p1(2:end-1);
f = 200*10^3*(0:(n/2))/n;
figure
plot(f,p1)
nfft = n;
y = fft(x,nfft);
figure
plot(xt)
y = y(1:nfft/2); 
y_mag = abs(y);
f = (0:nfft/2-1)*f_s/nfft; 
figure
plot(f,y_mag)


