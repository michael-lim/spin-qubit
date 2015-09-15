classdef qBit < handle
% qBit class represents a qubit which could be stored in global qstate
%     qBit Properties
%         psi - Spin state
%         rho - Density matrix
%         cpsi - Spin state time evolution
%         crho - Density matrix time evolution
%         dephase - not implemented
%         depolar - not implemented
%         parent - global parent quantum state, not implemented
%         epsilon - rounding error
%         bvec - Bloch vector of the qbit for pure states
%         purity - the purity of qubit
%     qBit Methods
%         qBit - constructor
%         isPure - returns boolean for the purity of qubit
%         evolve - evolves qBit according to the input Hamiltonian
%         plot - plots Bloch vector of the qBit
%         plotev - plots the evolution of the Bloch vector

    properties
        psi;
        rho;
        ipsi;
        irho;
        cpsi;
        crho;
        dephase; %Not implemented
        depolar; %Not implemented
        parent;
    end
    properties (SetAccess = immutable) %Rounding error
       epsilon = 1e-12; 
    end
    properties (Dependent = true, SetAccess = private)
        bvec;
        purity;
    end
    properties (Constant, Hidden)
        si = eye(2);
        sx = [0 1;1 0];
        sy = [0 -1i; 1i 0];
        sz = [1 0;0 -1];
    end
    
    methods
        function qb = qBit(psi) %Constructor
            if nargin <1 || isempty(psi)
               psi = [1;0]; 
            end
            [s1, s2]= size(psi);
            if s1==1 || s2 ==1 % input of a vector
               psit = psi;
               psit = psit/norm(psit);
               qb.psi = psit;
               qb.rho = psit*psit';
            elseif s1==s2 %input of a density matrix
                qb.rho = psi;
                qb.psi = [0;0];
            else
                error('input state must be vector or square matrix');
            end 
            qb.ipsi = qb.psi;
            qb.irho = qb.rho;
            qb.cpsi = struct();
            qb.crho = struct();
        end
        function bool = isPure(qb)
            if (abs(trace(qb.rho)-1)>qb.epsilon)
                warning('Trace rho not equal to 1')
            end
            bool = abs(qb.purity-1)<qb.epsilon;
        end
        function bvec = get.bvec(qb)
            bvec = real([trace(qb.rho*qb.sx), trace(qb.rho*qb.sy), trace(qb.rho*qb.sz)]);
        end
        function purity = get.purity(qb)
            purity = sqrt(trace(qb.rho^2));
        end
        function H = runnoise(qb,H0,H1,N,m,sg)
            rot_rate = sg * randn(1,N);
            mu = m * ones(1,N);
            for k = 1:N
                H(:,:,k) = mu(k)*H0 + rot_rate(k)*H1;
            end
        end
        function H = hamnoise(qb,H0,H1,t,m,sg)
            if length(t) == 1
                tl = t-1;
            else
                tl = length(t)-1;
            end
            rot_rate = sg * randn(1,tl);
            mu = m * ones(1,tl);
            for k = 1:tl
                H(:,:,k) = mu(k)*H0 + rot_rate(k)*H1;
            end
        end
        function H = hamrunnoise(qb,H0,H1,t,N,m,sg)
            if length(t) == 1
                tl = t-1;
            else
                tl = length(t)-1;
            end
            for z=1:N
                rot_rate = sg * randn(1,tl);
                mu = m * ones(1,tl);
                for k = 1:tl
                    H(:,:,k,z) = mu(k)*H0 + rot_rate(k)*H1;
                end
            end
        end
        function evolve(qb,H,t)
            qb.cpsi.runs(1).psi= qb.psi;
            qb.crho.runs(1).rho= qb.rho;
            if ndims(H) < 3
                if length(t) ==1
                    qb.psi = qb.stevolve(qb.psi,H,t);
                    qb.rho = qb.stevolve(qb.rho,H,t);
                    qb.cpsi.runs(1).psi= cat(3,qb.cpsi.runs(1).psi,qb.psi);
                    qb.crho.runs(1).rho= cat(3,qb.crho.runs(1).rho,qb.rho);
                else
                    dt = diff(t);
                    for i=1:length(dt)
                        qb.psi = qb.stevolve(qb.psi,H,dt(i));
                        qb.rho = qb.stevolve(qb.rho,H,dt(i));
                        qb.cpsi.runs(1).psi= cat(3,qb.cpsi.runs(1).psi,qb.psi);
                        qb.crho.runs(1).rho= cat(3,qb.crho.runs(1).rho,qb.rho);
                    end
                end
            else
                dt = diff(t);
                if (length(dt) ~= size(H,3))
                    warning('Time step does not equal the number of Hamiltonians.');
                end
                for i=1:length(dt)
                    qb.psi = qb.stevolve(qb.psi,H(:,:,i),dt(i));
                    qb.rho = qb.stevolve(qb.rho,H(:,:,i),dt(i));
                    qb.cpsi.runs(1).psi= cat(3,qb.cpsi.runs(1).psi,qb.psi);
                    qb.crho.runs(1).rho= cat(3,qb.crho.runs(1).rho,qb.rho);
                end
            end
        end

        function nevolve(qb,H,t,N,steps)
            for z = 1:N
                qb.psi = qb.ipsi;
                qb.rho = qb.irho;
                qb.cpsi.runs(z).psi= qb.psi;
                qb.crho.runs(z).rho= qb.rho;
                if length(t)==1
                    tl = linspace(1,t,steps+1);
                else
                    tl = t;
                end
                dt = diff(tl);
                switch ndims(H)
                    case 2
                        Hl = repmat(H,1,1,N,N);
                    case 3
                        Hl = repmat(H,1,1,1,N);
                    case 4
                        Hl = H;
                    otherwise
                        warning('Not a valid Hamiltonian.');
                end

                if length(dt) ~= steps || size(H,3) ~= steps
                    warning('Time step does not equal the increment number of Hamiltonians.');
                end

                if size(H,4) ~= N
                    warning('Number of experiments doest not equal the run number of Hamiltonians.');
                end

                for i=1:length(dt)
                    qb.psi = qb.stevolve(qb.psi,H(:,:,i,z),dt(i));
                    qb.rho = qb.stevolve(qb.rho,H(:,:,i,z),dt(i));
                    qb.cpsi.runs(z).psi= cat(3,qb.cpsi.runs(z).psi,qb.psi);
                    qb.crho.runs(z).rho= cat(3,qb.crho.runs(z).rho,qb.rho);
                end
            end
        end

        function s = measureSi(qb,Si)
            X = expm(i*sy*pi/4);
            Y = expm(i*sx*pi/4);
            if strcmpi(Si,'x') 
                s = (X*qb.psi)'*qb.sz*(X*qb.psi);
            elseif strcmpi(Si,'y') 
                s = (Y*qb.psi)'*qb.sz*(Y*qb.psi);
            elseif strcmpi(Si,'z') 
                s = qb.psi'*qb.sz*qb.psi;
            else
                warning('There is no such axis.');
            end
        end

        function h = plot(qb)
            clf;
            if (isPure(qb)<qb.epsilon)
                warning('State is not pure and will not plot Bloch vector');
            end
            c=[255/255 139/255 29/255];
            c1=[0 0 0];
            cmds={};
            cmds=[cmds struct('type','sphere','color',[c .1])];
            cmds=[cmds struct('type','equator','color',[c1 .2],'color2',c)];
            cmds=[cmds struct('type','spline','color',[c1 .2],'color2',c)];
            cmds=[cmds struct('type','label','val',[1.1 0 0],'label','|0+1>')];
            cmds=[cmds struct('type','label','val',[-1.1 0 0],'label','|0-1>')];
            cmds=[cmds struct('type','label','val',[0 0 1.1],'label','|0>')];
            cmds=[cmds struct('type','label','val',[0 0 -1.1],'label','|1>')];
            cmds=[cmds struct('type','vector','val',qb.bvec,'size',1,'color',[0 0 1])];
            h=plotBloch(cmds);
        end
        function h = plotev(qb,H,t)
            clf;
            if (isPure(qb)<qb.epsilon)
                warning('State is not pure and will not plot Bloch vector');
            end
            c=[255/255 139/255 29/255];
            c1=[0 0 0];
            cmds={};
            cmds=[cmds struct('type','sphere','color',[c .1])];
            cmds=[cmds struct('type','equator','color',[c1 .2],'color2',c)];
            cmds=[cmds struct('type','spline','color',[c1 .2],'color2',c)];
            cmds=[cmds struct('type','label','val',[1.1 0 0],'label','|0+1>')];
            cmds=[cmds struct('type','label','val',[-1.1 0 0],'label','|0-1>')];
            cmds=[cmds struct('type','label','val',[0 0 1.1],'label','|0>')];
            cmds=[cmds struct('type','label','val',[0 0 -1.1],'label','|1>')];
            h = plotBloch(cmds);
            dt = diff(t);
            for i = 1:length(dt)
                evolve(qb,H,dt(i));
                cmd={};
                cmd=[cmd struct('type','vector','val',qb.bvec,'size',1,'color',[0 0 1])];
                h = plotBloch(cmd);
                pause(0.001);
            end
        end
    end
    
    methods (Static)
        function rho = bloch2rho(v)
              x = v(1)*qb.sx;
              y = v(2)*qb.sy;
              z = v(3)*qb.sz;  
              rho = (qb.si+x+y+z)/2;
        end
        function y = stevolve(psii,H,t)
            psit = psii;
            Httemp = expm(1i*H*t);
            Httempm = expm(-1i*H*t);
            if size(psii,2)==2 
                psit=Httemp*psit*Httempm;
                y = psit;
            else
                psit=Httemp*psit;
                y = psit;
            end
        end
    end
end
