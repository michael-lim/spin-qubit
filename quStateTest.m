%%
% Demonstration of quState class:
% First we check that the constructor with construct something according to
% the inputs:

qs1 = quState();
qs2 = quState([1;0]);
qs3 = quState([1,0;0,1]);
qs4 = quState([0;0;1;0]);
qs5 = quState([1,0,0,0;0,0,0,0;0,0,0,0;0,0,0,0]);
qs6 = quState([1,0;0,1],[0,1;1,0]);
qs7 = quState([1;0],[0;1]);
qs8 = quState([1;0;0]);
qs9 = quState([1;0],[1;0],[1;0]);

%%
% Checking Time evolution

clear;clc;
qs = quState([1;0],[0;-1]) %Generates an antisymmetric quState


sx=[0,1;1,0];
sy=[0,-1i;1i,0];
sz=[1,0;0,-1];
sx1 = kron([0 1;1 0],eye(2));
sy1 = kron([0 -1i; 1i 0],eye(2));
sz1 = kron([1 0;0 -1],eye(2));
sx2 = kron(eye(2),[0 1;1 0]);
sy2 = kron(eye(2),[0 -1i; 1i 0]);
sz2 = kron(eye(2),[1 0;0 -1]);

%%

H0 = sx1*sx2; %this will rotate psi 
H1 = sz1*sz2; %there is a random noise
t = linspace(0,10*pi,100);
steps = length(t)-1;
m = 1;
sg = .3;

N = 500;

H=qs.hamrunnoise(H0,H1,t,N,m,sg);

qs.nevolve(H,t,N,steps);
