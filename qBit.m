classdef qBit < handle
% qBit class represents a qubit which could be stored in global qstate
%     qBit Properties
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
        rho;
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
               qb.rho = psit*psit';
            elseif s1==s2 %input of a density matrix
                qb.rho = psi;
            else
                error('input state must be vector or square matrix');
            end 
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
            qb.rho = qb.stevolve(qb.rho,H,t);
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
        function h = plotev(qb,H,t,dt)
            clf;
            if nargin == 3
               dt = 0.1;
            end
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
            t = [0:dt:t];
            for i = 1:length(t)
                evolve(qb,H,dt);
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
