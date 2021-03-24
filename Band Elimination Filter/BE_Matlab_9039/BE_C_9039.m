%% Pardali Christina AEM 9039/September 2020 
%% Band Elimination Chebyshev Filter 

clear;
clc;
 
%% AEM

a1=9;
a2=0;
a3=3;
a4=9;

%% Prodiagrafes filtrou me basi thn ekfwnisi

f0=1.05*1000;
f1=725+(25*a4);
f2=(f0^2)/(f1);
D=(1/(1.8))*(((f0^2)-(f1^2))/f1);
f3=(-D+sqrt((D^2)+4*(f0^2)))/2;
f4=(f0^2)/f3;

amin=25+((a3-5)/10);
amax=0.5+(a4/10);


%% Gwniakes sixnotites

w0=2*pi*f0;
w1=2*pi*f1;
w2=2*pi*f2;
w3=2*pi*f3;
w4=2*pi*f4;

bw=w2-w1;
qc= w0/bw;

Wp=1;
Ws=(w2-w1)/(w4-w3);

%% Taksi fliltrou 

n_decimal=acosh(sqrt((10^(amin/10)-1)/((10^(amax/10)-1))))/acosh(Ws);

n = ceil(n_decimal);

e = sqrt(10^(amax/10)-1);
a = (1/n)*( asinh(1/e));

%% Sixnotita imiseias isxios

whp = 1/(cosh(acosh(1/e))/n);

%% Butterworth angles

y1 = 22.5;
y2 = -22.5;
y3 = 67.5;
y4 = -67.5;

%% Poloi Chebyshev

p1 = -sinh(a) * cosd(y1)+ (1i)*cosh(a) * sind(y1);
p2 = -sinh(a) * cosd(y2) + (1i)*cosh(a) * sind(y2);
p3 = -sinh(a) * cosd(y3) + (cosh(a) * sind (y3))*(1i);
p4 = -sinh(a) * cosd(y4) + (cosh(a) * sind (y4))*(1i);


W_12=sqrt((real(p1))^2+(imag(p1))^2);
W_34=sqrt((real(p3))^2+(imag(p3))^2);


Q12=W_12/(2*abs(real(p1)));
Q34=W_34/(2*abs(real(p3)));


%% Antistrofi polwn 

InvW1_first= 1 / W_12;
InvW2_first = 1 / W_34;

%% Ypologismos kainouriwn gwniwn

y1_new= acosd(1/(2*Q12));
y2_new= - acosd(1/(2*Q12));
y3_new= acosd(1/(2*Q34));
y4_new= - acosd(1/(2*Q34));


%%  Kainourioi poloi 

p1_new = InvW1_first*(-cosd(y1_new)+(1i)* sind(y1_new));
p2_new = InvW1_first*(-cosd(y1_new)-(1i) * sind(y1_new));
p3_new = InvW2_first*(-cosd(y3_new)+ sind(y3_new)*(1i));
p4_new = InvW2_first*(-cosd(y3_new)- sind(y3_new)*(1i));


W_12_new=sqrt((real(p1_new))^2+(imag(p1_new))^2);
W_34_new=sqrt((real(p3_new))^2+(imag(p3_new))^2);


Q_12_new=W_12_new/(2*abs(real(p1_new)));
Q_34_new=W_34_new/(2*abs(real(p3_new)));


S12 = abs(real(p1_new));
S34 = abs(real(p3_new));


%% Piknwtis

Cnew = 10^(-8);

%% Metasximatismos 1ou migadikou polou 

C1= S12^2 + ((imag(p1_new))^2);
D1= (2*S12) / qc;
E1= 4 + (C1 / ( qc^2));
G1= sqrt (E1^2 - 4*(D1^2));
Q1_2= (1/D1)* (sqrt( 1/2 *( E1+G1)));
k1= (S12 * Q1_2) / qc;

W1= k1 + sqrt((k1^2) -1);

w01 = (1/W1)*w0;
w02 = W1*w0;

%% metasximatismos 2ou migadikou polou 

C2= S34^2 + ((imag(p3_new))^2);
D2= (2* S34)/ qc;
E2= 4 + (C2 / qc^2);
G2= sqrt (E2^2 - 4* (D2^2));
Q3_4= (1/D2) * sqrt((1/2) * ( E2+ G2));
k2= ((S34)*Q3_4 ) / qc;
W2= k2 + sqrt( (k2^2) -1);

w03 = (1/W2)*w0;
w04 = W2*w0;


%/-------------------Prwti monada-LPN-----------------------/

%% Ilopoiisi monadas 1 - LPN

wzo1 = w0 / w01 ;

R11 = 1;
R12 = 4*(Q1_2^2);
R13 = wzo1^2 / ( 2 * (Q1_2^2));
R14 = 1;
R15 = (4 * (Q1_2^2)) / ( wzo1^2 -1);

C1_lpn = 1 / (2 * Q1_2);


k1_lpn = (1/(R13+1));

H1= k1_lpn * (wzo1^2);

%% klimakopoihsh 1is monadas

kf_1 = w01;

km1 = C1_lpn / (kf_1 * Cnew);

R11new = R11 * km1;
R12new = R12 * km1;
R13new = R13 * km1;
R14new = R14 * km1;
R15new = R15 * km1;

%/-------------------Deuteri monada-HPN-----------------------/

%% Ilopoiisi monadas 2 - HPN

wzo2 = w0 / w02;

R21 = 1 ;
R23 = 1;

k21 = ( 1 / (wzo2))^2  - 1;

k22 = ((2 + k21)) * Q1_2^2 / ((((2+k21)) * Q1_2^2)+1) ;

R22 = ((2 + k21) ^2) * (Q1_2^2);
R24 = (2 + k21) * (Q1_2 ^2);

C22 = 1/ ( (2+k21) * Q1_2 );
C21 = k21 * C22;

k2_hpn = k22 * ((1/ wzo2)^2);

H2=k2_hpn*(wzo2^2);

%% klimakopoihsh monadas 2

kf_2 = w02;
km2 = C22/ ( kf_2 * Cnew);

R21new = R21 * km2;
R22new = R22 * km2;
R23new = R23 * km2;
R24new = R24 * km2;

C21new = C21/ (km2* kf_2);
C22new = C22/ (km2* kf_2);

%/-------------------Triti monada-LPN-----------------------/

%% Ilopoiisi monadas 3 - LPN

wzo3 = w0 / w03 ;

R31 = 1;
R32 = 4*(Q3_4^2);
R33 = wzo3^2 / ( 2 * (Q3_4^2));
R34 = 1;
R35 = 4 * (Q3_4^2) / ( wzo3^2 -1);

C3_lpn = 1 / (2 * Q3_4);

k3_lpn = (1/(R33+1));
 
H3=k3_lpn*(wzo3^2);

%% klimakopoihsh monadas 3

kf_3 = w03;

km3 = C3_lpn / (kf_3 * Cnew);

R31new = R31 * km3;
R32new = R32 * km3;
R33new = R33 * km3;
R34new = R34 * km3;
R35new = R35 * km3;


%/-------------------Tetarti monada-HPN-----------------------/

%% Ilopoiisi monadas 4 - HPN

wzo4= w0/ w04;

R41=1 ;
R43=1;

k41 = (1/ wzo4)^2- 1;
k42 = (((2 + k41)) * (Q3_4^2)) / (((2+k41))* (Q3_4^2) +1) ;

R42= ((2 + k41) ^2) * (Q3_4^2);
R44= (2 + k41) * (Q3_4^2);

C42 = 1 / ( (2+k41) * (Q3_4));
C41 = k41 * C42;

k4_hpn = k42 * (1 / wzo4^2);

H4=k4_hpn*(wzo4^2);

%% klimakopoihsh monadas 4 

kf4 = w04;
km4 = C42 / ( kf4 * Cnew);

R41new = R41 * km4;
R42new = R42 * km4;
R43new = R43 * km4;
R44new = R44 * km4;

C41new = C41 / (km4 * kf4);
C42new = C42 / (km4 * kf4);
 
%% synartiseis metaforas
T1 = tf( [k1_lpn 0 (k1_lpn * w0^2 ) ], [ 1 ( w01 / Q1_2 ) w01^2 ] );
T2 = tf( [k2_hpn 0 (k2_hpn * w0^2 ) ], [ 1 ( w02 / Q1_2 ) w02^2 ] );
T3 = tf( [k3_lpn 0 (k3_lpn * w0^2 ) ], [ 1 ( w03 / Q3_4 ) w03^2 ] );
T4 = tf( [k4_hpn 0 (k4_hpn * w0^2 ) ], [ 1 ( w04 / Q3_4 ) w04^2 ] );

K_total = k1_lpn * k2_hpn * k3_lpn * k4_hpn;

H_total=H1*H2*H3*H4;

T_total_before_gain = T1*T2*T3*T4;

a_kerdos= 1/(H_total);

T_total = a_kerdos * T_total_before_gain;

%% Dimiourgia diagrammatwn

plot_transfer_function(T1, [f1 f2 f3 f4])

plot_transfer_function(T2, [f1 f2 f3 f4])

plot_transfer_function(T3, [f1 f2 f3 f4])

plot_transfer_function(T4, [f1 f2 f3 f4])

ltiview({'bodemag'}, T1, T2, T3, T4, T_total)

plot_transfer_function(T_total, [f1 f2 f3 f4])

InvSys_new = inv (T_total);

plot_transfer_function(InvSys_new, [10 f1 f2 f3 f4])

%% Fourier analysis

f11=(w0-((w0-w3)/2)) /(2*pi);
f12 = (w0+((w0+w3)/3)) /(2*pi);
f13 = (0.4*w1)/(2*pi);
f14 = (2.5*w2)/(2*pi);
f15 = (3*w2)/(2*pi);

T =0.002 ;
%% Sixnotita deigmatolipsias

f_s=200*(10^3);
dt = 1/f_s;
t = 0:dt:(T);

%% Dimiourgia kai emfanisi simatos eisodou

u1= 0.5*cos((w0-((w0-w3)/2))*t)+0.8*cos((w0+((w0+w3)/3))*t)+0.8*cos(0.4*w1*t)+0.6*cos(2.5*w2*t)+1.2*cos(3*w2*t);
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
p2_f=abs(xfourier/n1);
p1_f=p2_f(1:n1/2+1);
p1_f(2:end-1)=2*p1_f(2:end-1);
f=200*10^3*(0:(n1/2))/n1;
figure
plot(f,p1_f)

nfft=n1;
y=fft(u1,nfft);
y = y(1:nfft/2); 
y_mag=abs(y);
f = (0:nfft/2-1)*f_s/nfft; 
figure
plot(f,y_mag) 

