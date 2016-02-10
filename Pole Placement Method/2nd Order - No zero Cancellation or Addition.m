% 2nd order Pole Placement method example
% no Zeros added \ Cancelled
clc;

format long;

s = tf('s');
Ts = 0.01 % given
z = tf('z',Ts);

% Required Transfer Function in Continuous Time
F = 10/(s*(s-10));

% Gain Correction due to ADC \ DAC Block
GainCorrection = 0.5;
% Transfer function in Discrete Time
G = c2d(F,Ts);
G = G*GainCorrection

% Extract the Numerator and Denominator from The Discrete time Transfer Function
[num,den] = tfdata(G,'v');
B = tf(num,1,Ts); % Numerator - 1st Order
A = tf(den,1,Ts); % Denominator - 2nd Order

% Required Denominator By Specifications 
Am = z^3-1.89*z^2+1.135*z-0.2147; %given

% From Degree Analysis of Am = A*R+B*S
% R = z + r0 ; 1st degree monic
% S = s1*z+s0 ; 1st degree not monic

% Solving System Of Equations from : Am = A*R+B*S to get coefficients of R & S
AA = [ 1 num(2) 0 den(2) ; den(2) num(3) num(2) den(3) ; den(3) 0 num(3) 0 ; 0 0 0 1 ]
v = [ -1.89 1.135 -0.2147 1 ]
th = inv(AA)*v';

r0 = th(1);
s1 = th(2);
s0 = th(3);

R = z + r0
S = s1*z+s0

% Resultant Closed loop Transfer Function
W = S*B/(A*R+S*B);

Gain = dcgain(W);
% Apply Gain Correction to meet specifications of Steady State Gain = 1
RequiredGain = 1; %given
alpha = RequiredGain/Gain;
W = alpha*W;

minreal(zpk(W))

step(W);