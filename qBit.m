classdef qBit < handle
% qBit class represents a qubit which could be stored in global qstate
%     qBit Properties
%         psi - spin state
%         rho - Density matrix
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
%         plot - plots Bloch vector of the qBit for pure states
%         plotev - plots the evolution of the Bloch vector

    properties
        psi;
        rho;
        dephase; %Not implemented
        depolar; %Not implemented
        parent;
        tpsi;
        trho;
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
            qb.trho = qb.rho;
            qb.tpsi = qb.psi;
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
        function evolve(qb,H,t)
            dt = diff(t);
            qb.trho = cat(3,qb.trho,qb.stevolve(qb.rho,H,dt));
            qb.tpsi = cat(3,qb.tpsi,qb.stevolve(qb.psi,H,dt));
            qb.rho = qb.trho(end);
            qb.psi = qb.tpsi(end);
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
        function y = stevolve(psii,H,dt)
            psit = psii;
            for i=1:length(dt)
                Httemp = expm(1i*H*dt(i));
                Httempm = expm(-1i*H*dt(i));
                if size(psii,2)==2
                    y=zeros(2,2,length(dt));
                    psit=Httemp*psit*Httempm;
                    y(:,:,i)=psit;
                else
                    y=zeros(2,1,length(dt));
                    psit=Httemp*psit;
                    y(:,:,i)=psit;
                end
            end
        end
    end
end
