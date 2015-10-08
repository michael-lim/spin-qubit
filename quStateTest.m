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
qs = quState([1;0],[0;1]) %Generates an antisymmetric quState


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

H0 = sx1+sx2; %this will rotate psi 
H1 = sx1+sx2; %there is a random noise
t = linspace(0,10*pi,100);   
steps = length(t)-1;
m = 1;
sg = 0.3;

N = 50;

H=qs.invpownoise(H0,H1,t,N,m,sg);
%H=qs.hamrunnoise(H0,H1,t,N,m,sg);

qs.nevolve(H,t,N,steps);

X = measureSiN12(qs,'x',N,1);
Y = measureSiN12(qs,'y',N,1);
Z = measureSiN12(qs,'z',N,1);
figure(1); clf;
hold on;
plot(t,mean(X,2),'r')
plot(t,mean(Y,2),'b')
plot(t,mean(Z,2),'k')
hold off;


%%


H0 = sx1+sx2; %this will rotate psi 
H1 = sx1+sx2; %there is a random noise
t = linspace(0,10*pi,100);   
steps = length(t)-1;
m = 1;
sg = 0.3;

N = 50;

%H=qs.invpownoise(H0,H1,t,N,m,sg);
H=qs.hamrunnoise(H0,H1,t,N,m,sg);

qs.nevolve(H,t,N,steps);

X = measureSiN12(qs,'x',N,1);
Y = measureSiN12(qs,'y',N,1);
Z = measureSiN12(qs,'z',N,1);
figure(2); clf;
hold on;
plot(t,mean(X,2),'r')
plot(t,mean(Y,2),'b')
plot(t,mean(Z,2),'k')
hold off;
%%

tl = 1e6;
invpow = zeros(1,tl);
invpow(1) = rand(1);
for i=1:(tl-1)
    invpow(i+1) = mod((invpow(i)+invpow(i)^(1.7)),1);
end
t=tic
pow = fft(invpow);
toc(t)