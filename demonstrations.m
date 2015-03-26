%%
%Demonstrations of qBit class
qb = qBit([1;0]) %Generates a qBit with psi=[1;0]
qb.rho %Displays the density matrix
isPure(qb) %Shows whether the qBit is pure

H = [0,-1i;1i,0]; %sy - this will rotate psi with respect to y-axis
t = linspace(0,pi,100);
evolve(qb,H,t) %This will evolve qBit for pi duration
plotev(qb,H,t) %This shows an animated plot of qBit evolution

%%
%Demonstration of repeated experiments with normally distributed rotation rates
sx=[0,1;1,0];
sy=[0,-1i;1i,0];
sz=[1,0;0,-1];
mu = 1;
sigma = .07;

N = 500;

rot_rate = mu + sigma*randn(1,N);

t = linspace(0,10*pi,100);
dt = diff(t);

Z = zeros(length(t),length(rot_rate));
for i = 1:length(rot_rate)
   qb = qBit([1;0]);
   for j = 1:length(t)-1
       qb.evolve([0,1;1,0]*rot_rate(i),dt(j));
       tmpPsi = qb.psi;
       Z(j,i) = tmpPsi'*sz*tmpPsi;
   end
end

figure(1); clf;
plot(t,mean(Z,2))