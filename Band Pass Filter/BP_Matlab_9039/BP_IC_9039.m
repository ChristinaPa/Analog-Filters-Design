%% Pardali Christina AEM 9039/September 2020 
%% Bandpass Inverse Chebyshev Filter 

clear;
clc;

%% AEM 

a1 = 9;
a2 = 0;
a3 = 3;
a4 = 9;

%% Prodiagrafes filtrou me basi thn ekfwnisi

f0 = 1000;
f1 = 650 + (25*a4);
f2 = (f0^2)/f1;
D  = 2.1*((f0^2)-(f1^2))/f1;
f3 = (-D + sqrt(D^2 + 4*f0^2))/2;
f4 = f0^2/f3;

amin = 35 - a3;
amax = 0.4 + (a4/36);

%% Gwniakes sixnotites

w0 = 2*pi*f0;
w1 = 2*pi*f1;
w2 = 2*pi*f2;
w3 = 2*pi*f3;
w4 = 2*pi*f4;

Wp=1;
Ws = (w4 - w3)/(w2 - w1);

BW = 2*pi*(f2 - f1);

qc = w0/BW;

%% Taksi filtrou 

n_decimal= acosh(sqrt((10^(amin/10)-1)/(10^(amax/10)-1)))/acosh(Ws);
n = ceil(n_decimal);

e = 1/sqrt(10^(amin/10)-1);
a = (1/n)*( asinh(1/e));

%% Sixnotita imiseias isxios

whp = 1/cosh(acosh(1/e)/n);

%% Butterworth angles

y1 = 22.5;
y2 = -22.5;
y3 = 67.5;
y4 = -67.5;

p1 = -sinh(a) * cosd(y1)+ (i)*cosh(a) * sind(y1);
p2 = -sinh(a) * cosd(y2) + (i)*cosh(a) * sind(y2);
p3 = -sinh(a) * cosd(y3) + (cosh(a) * sind (y3))*i;
p4 = -sinh(a) * cosd(y4) + (cosh(a) * sind (y4))*i;

W_12 = sqrt((real(p1))^2+(imag(p1))^2);
W_34 = sqrt((real(p3))^2+(imag(p3))^2);

Q12 = W_12/(2*abs(real(p1)));
Q34 = W_34/(2*abs(real(p3)));

%% antistrofi polon

InvW1_first= 1 / W_12;
InvW2_first = 1 / W_34;

%% klimakopoihsh sixnotitas

InvW1= InvW1_first * Ws ;
InvW2= InvW2_first * Ws;

S12 = - InvW1 / ( 2 * Q12 );
S34=  - InvW2 / ( 2 * Q34);
W12 = sqrt( InvW1^2 - S12^2 );
W34 = sqrt( InvW2^2 - S34^2 );

%% midenika gia k=1 kai k=3

Z1 = sec( pi /(2*n));
Z2 = sec( 3*pi /(2*n));

%% klimakopoihsh mhdenikwn

Z1new = Z1 * Ws;
Z2new = Z2 * Ws;

%% metasximatismos 1ou migadikou polou 

C1= S12^2 + W12^2;
D1= (-2*S12) / qc;
E1= 4 + (C1 / ( qc^2));
G1= sqrt (E1^2 - 4*(D1^2));
Q1_2= (1/D1)* (sqrt( 1/2 *( E1+G1)));
k1= (-S12 * Q1_2) / qc;

W1= k1 + sqrt((k1^2) -1);
w01 = (1/W1)*w0;
w02 = W1*w0;

%% metasximatismos 2ou migadikou polou 

C2= S34^2 + W34^2;
D2= (-2* S34)/ qc;
E2= 4 + (C2 / qc^2);
G2= sqrt (E2^2 - 4* (D2^2));
Q3_4= (1/D2) * sqrt((1/2) * ( E2+ G2));
k2= ((-S34)*Q3_4 ) / qc;
W2= k2 + sqrt( (k2^2) -1);
w03 = (1/W2)*w0;
w04 = W2*w0;

%% metasximatismos 1ou fantastikou midenikou

Kzero1 = 2 + ((Z1new^2) / (qc^2));
x1 = ( Kzero1 + sqrt( Kzero1^2 - 4 ) ) / 2;
wz1 = w0 * (sqrt(x1));
wz2 = w0 / (sqrt(x1));

%% metasximatismos 2ou fantastikou midenikou 

Kzero2 = 2 + ((Z2new^2) / (qc^2));
x2 = ( Kzero2+ sqrt( Kzero2^2 - 4 ) ) / 2;
wz3 = w0*(sqrt(x2));
wz4 = w0/(sqrt(x2));

%% Piknwtis
C1new = 10^(-7);

%/-------------------Prwti monada-LPN-----------------------/

%% Ilopoiisi monadas 1 - LPN

wzo1 = wz1 / w01 ;

R11 = 1;
R12 = 4*(Q1_2^2);
R13 = wzo1^2 / ( 2 * (Q1_2^2));
R14 = 1;
R15 = (4 * (Q1_2^2)) / ( wzo1^2 -1);

C1_lpn = 1 / (2 * Q1_2);

k1_lpn = 1 / ( 1 +(wzo1^2 / (2 * (Q1_2^2))));

%% klimakopoihsh 1is monadas

kf_1 = w01;

km1 = C1_lpn / (kf_1 * C1new);

R11new = R11 * km1;
R12new = R12 * km1;
R13new = R13 * km1;
R14new = R14 * km1;
R15new = R15 * km1;

%/-------------------Deuteri monada-HPN-----------------------/

%% Ilopoiisi monadas 2 - HPN

wzo2 = wz2 / w02;

R21 = 1 ;
R23 = 1;

k21 = ( 1 / (wzo2))^2  - 1;
k22 = ((2 + k21)^2) * Q1_2^2 / ((((2+k21)^2) * Q1_2^2)+1) ;

R22 = ((2 + k21) ^2) * (Q1_2^2);
R24 = (2 + k21) * (Q1_2 ^2);

C22 = 1/ ( (2+k21) * Q1_2 );
C21 = k21 * C22;

k2_hpn = k22 * ((1/ wzo2)^2);

%% klimakopoihsh monadas 2

kf_2 = w02;
km2 = C22/ ( kf_2 * C1new);

R21new = R21 * km2;
R22new = R22 * km2;
R23new = R23 * km2;
R24new = R24 * km2;

C21new = C21/ (km2* kf_2);
C22new = C22/ (km2* kf_2);

%/-------------------Triti monada-LPN-----------------------/

%% Ilopoiisi monadas 3 - LPN

wzo3 = wz3 / w03 ;

R31 = 1;
R32 = 4*(Q3_4^2);
R33 = wzo3^2 / ( 2 * (Q3_4^2));
R34 = 1;
R35 = 4 * (Q3_4^2) / ( wzo3^2 -1);

C3_lpn = 1 / (2 * Q3_4);

k3_lpn = 1 / ( 1 + (wzo3^2 / (2  * (Q3_4^2)) ));

%% klimakopoihsh monadas 3

kf_3 = w03;

km3 = C3_lpn / (kf_3 * C1new);

R31new = R31 * km3;
R32new = R32 * km3;
R33new = R33 * km3;
R34new = R34 * km3;
R35new = R35 * km3;

C1new = 10^(-7);

%/-------------------Tetarti monada-HPN-----------------------/

%% Ilopoiisi monadas 4 - HPN

wzo4= wz4/ w04;

R41=1 ;
R43=1;

k41 = (1/ wzo4)^2- 1;
k42 = (((2 + k41)^2) * (Q3_4^2)) / (((2+k41)^2)* (Q3_4^2) +1) ;

R42= ((2 + k41) ^2) * (Q3_4^2);
R44= (2 + k41) * (Q3_4^2);

C42 = 1 / ( (2+k41) * (Q3_4));
C41 = k41 * C42;

k4_hpn = k42 * (1 / wzo4^2);

%% klimakopoihsh monadas 4 

kf4 = w04;
km4 = C42 / ( kf4 * C1new);

R41new = R41 * km4;
R42new = R42 * km4;
R43new = R43 * km4;
R44new = R44 * km4;

C41new = C41 / (km4 * kf4);
C42new = C42 / (km4 * kf4);

%% synartiseis metaforas

T1 = tf( [k1_lpn 0 ( k1_lpn * wz1^2 ) ], [ 1 ( w01 / Q1_2 ) w01^2 ] );
T2 = tf( [k2_hpn 0 ( k2_hpn * wz2^2 ) ], [ 1 ( w02 / Q1_2 ) w02^2 ] );
T3 = tf( [k3_lpn 0 ( k3_lpn * wz3^2 ) ], [ 1 ( w03 / Q3_4 ) w03^2 ] );
T4 = tf( [k4_hpn 0 ( k4_hpn * wz4^2 ) ], [ 1 ( w04 / Q3_4 ) w04^2 ] );

%% Sinartiseis metaforas gia kerdos 0dB

T1_1 = tf( [1 0 ( 1 * wz1^2 ) ], [ 1 ( w01 / Q1_2 ) w01^2 ] );
T2_2 = tf( [1 0 ( 1 * wz2^2 ) ], [ 1 ( w02 / Q1_2 ) w02^2 ] );
T3_3 = tf( [1 0 ( 1 * wz3^2 ) ], [ 1 ( w03 / Q3_4 ) w03^2 ] );
T4_4 = tf( [1 0 ( 1 * wz4^2 ) ], [ 1 ( w04 / Q3_4 ) w04^2 ] );

K_total = k1_lpn * k2_hpn * k3_lpn * k4_hpn;

T_total_before_gain = T1*T2*T3*T4;

gain = abs(evalfr(T_total_before_gain, w0 * 1i));

a_kerdos= (10^(0.5))/gain;

T_total = a_kerdos * T_total_before_gain;

%% Dimiourgia diagrammatwn

plot_transfer_function(T1, [f1 f2 f3 f4])

plot_transfer_function(T2, [f1 f2 f3 f4])

plot_transfer_function(T3, [f1 f2 f3 f4])

plot_transfer_function(T4, [f1 f2 f3 f4])

ltiview({'bodemag'}, T1, T2, T3, T4, T_total)

plot_transfer_function(T_total, [ f0 f1 f2 f3 f4])

InvSys_new = inv (T_total);

InvSys_new_kerdos0 = inv (T1_1*T2_2*T3_3*T4_4);

plot_transfer_function(InvSys_new, [ f0 f1 f2 f3 f4]);

plot_transfer_function(InvSys_new_kerdos0, [ f0 f1 f2 f3 f4]);

%% Fourier Analysis

Time = (1/100);
Fs = 1000000;
dt = 1/Fs;
t = 0:dt:Time-dt;

fin_1=(w0-((w0-w1)/2)) /(2*pi);
fin_2 = (w0+((w0+w1)/3)) /(2*pi);
fin_3 = (0.4*w3)/(2*pi);
fin_4 = (2.5*w4)/(2*pi);
fin_5 = (3*w4)/(2*pi);

%% input signal
x = cos(2*pi*fin_1*t)+0.8*cos(2*pi*fin_2*t)+0.8*cos(2*pi*fin_3*t)+0.6*cos(2*pi*fin_4*t)+0.5*cos(2*pi*fin_5*t);

figure;
plot(t,x);
title('Input signal');
y=lsim(T_total,x,t);
figure;
plot(t,y);
title('Output signal');
figure;
plot(t,x);
hold on;
plot(t,y);
hold off;
title('Input and Output signals');
Xf=fft(x);
L=length(x);
P2 = abs(Xf/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure;
plot(f,P1);
axis([0.01 20000 0 inf]);
title('Single-Sided Amplitude Spectrum of Input signal');
Yf=fft(y);
L=length(y);
P2 = abs(Yf/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure;
plot(f,P1);
axis([0.01 20000 0 inf]);
title('Single-Sided Amplitude Spectrum of Output signal');




