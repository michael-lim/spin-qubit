%%
%Demonstrations of qBit class
qb = qBit([1;0]) %Generates a qBit with psi=[1;0]
qb.rho %Displays the density matrix
isPure(qb) %Shows whether the qBit is pure

H = [0,-1i;1i,0]; %sy - this will rotate psi with respect to y-axis
evolve(qb,H,pi) %This will evolve qBit for pi duration
plotev(qb,H,pi,0.03) %This shows an animated plot of qBit evolution