%%
% 
repeats = 1024;
matsize = 32;
tsize1 = zeros(repeats,1);
tsize2 = zeros(repeats,1);
tsize3 = zeros(repeats,1);
tavg1 = zeros(matsize,1);
tavg2 = zeros(matsize,1);
tavg3 = zeros(matsize,1);
tstd1 = zeros(matsize,1);
tstd2 = zeros(matsize,1);
tstd3 = zeros(matsize,1);
t = linspace(1,matsize,matsize)';

h = @(r) exp(1i*pi*rand(r,r));

%%
for j=1:matsize
    %expm1
    for i=1:repeats
        m = h(j);
        m = m'+m;
        t = tic;
        mm = expmdemo1(m);
        tsize1(i) = toc(t);
    end
    tavg1(j)=mean(tsize1);
    tstd1(j)=std(tsize1);
    %expm2
    for i=1:repeats
        m = h(j);
        m = m'+m;
        t = tic;
        mm = expmdemo2(m);
        tsize2(i) = toc(t);
    end
    tavg2(j)=mean(tsize2);
    tstd2(j)=std(tsize2);
    %expm3
    for i=1:repeats
        m = h(j);
        m = m'+m;
        t = tic;
        mm = expmdemo3(m);
        tsize3(i) = toc(t);
    end
    tavg3(j)=mean(tsize3);
    tstd3(j)=std(tsize3);
end
%%
t = linspace(1,matsize,matsize)';
clf;
figure(1);
hold on;
plot(tavg1, 'blue','DisplayName','expm');
%plot(tavg2, 'yellow','DisplayName','demo2');
plot(tavg3, 'green','DisplayName','demo3');
legend show
hold off;

figure(2);
hold on;
errorbar(t,tavg1,tstd1,'DisplayName','expm');
%errorbar(t,tavg2,tstd2,'DisplayName','demo2');
errorbar(t,tavg3,tstd3,'DisplayName','demo3');
legend show
hold off;

figure(3);
hold on;
errorbar(t,tavg1,tstd1,'DisplayName','expm');
errorbar(t,tavg3,tstd3,'DisplayName','demo3');
legend show
hold off;