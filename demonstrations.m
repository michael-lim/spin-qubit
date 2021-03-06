%%
%Demonstrations of quBit class
clear
qb = quBit([1;0]) %Generates a quBit with psi=[1;0]
qb.rho %Displays the density matrix
isPure(qb) %Shows whether the quBit is pure

H = [0,-1i;1i,0]; %sy - this will rotate psi with respect to y-axis
t = linspace(0,pi,100);
evolve(qb,H,t) %This will evolve quBit for pi duration
%plotev(qb,H,t) %This shows an animated plot of quBit evolution

%%
%Demonstration of repeated experiments with normally distributed rotation rates
sx=[0,1;1,0];
sy=[0,-1i;1i,0];
sz=[1,0;0,-1];
mu = 1;
sigma = .03;

N = 50;

rot_rate = mu + sigma*randn(1,N);

t = linspace(0,10*pi,100);
dt = diff(t);

Z = zeros(length(t),length(rot_rate));
for i = 1:length(rot_rate)
   qb = quBit([1;0]);
   for j = 1:length(t)-1
       qb.evolve([0,1;1,0]*rot_rate(i),dt(j));
       tmpPsi = qb.psi;
       Z(j,i) = tmpPsi'*sz*tmpPsi;
   end
end

figure(1); clf;
plot(t,mean(Z,2))

%%
%Demonstrations of N experimental runs
clear;clc;
qb = quBit([1;0]) %Generates a quBit with psi=[1;0]

sx=[0,1;1,0];
sy=[0,-1i;1i,0];
sz=[1,0;0,-1];

H0 = sx; %this will rotate psi 
H1 = sx; %there is a random noise
t = linspace(0,10*pi,100);
steps = length(t)-1;
m = 1;
sg = .3;

N = 50;

H=qb.hamrunnoise(H0,H1,t,N,m,sg);

qb.nevolve(H,t,N,steps);
%%
X = measureSiN(qb,'x',N);
Y = measureSiN(qb,'y',N);
Z = measureSiN(qb,'z',N);
figure(1); clf;
hold on;
plot(t,mean(X,2),'r')
plot(t,mean(Y,2),'b')
plot(t,mean(Z,2),'k')
hold off;

%%
%expsine = @(b,x) (b(1).*exp(-b(2).*x).*sin(b(3).*x+b(4)))';
beta = qb.fitesin(t,mean(sn,2),[1,0.01,pi/2,0]);

figure(2); clf;
hold on;
plot(t,mean(sn,2), 'red');
expsine = @(b,x) (b(1).*exp(-b(2).*x).*sin(b(3).*x+b(4)))';
plot(t,expsine(beta,t), 'blue');
hold off;
