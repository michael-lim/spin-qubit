classdef quState < handle
% quState class represents a quantum state which could be stored in global quState
%     quState Properties
%         psi - Spin state
%         rho - Density matrix
%         cpsi - Spin state time evolution
%         crho - Density matrix time evolution
%         dephase - not implemented
%         depolar - not implemented
%         epsilon - rounding error
%         bvec - Bloch vector of the quantum state for pure states
%         purity - the purity of quantum state
%     quState Methods
%         quState - constructor
%         isPure - returns boolean for the purity of qubit
%         evolve - evolves the bits according to the input Hamiltonian
%         plot - plots Bloch vector of the qubits
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
        function qs = quState(psi1,psi2) %Constructor
            switch nargin
                case 2
                    if iscolumn(psi1)==1
                        ptemp1 = psi1;
                    else
                        ptemp1 = psi1';
                    end
                    if iscolumn(psi2)==1
                        ptemp2 = psi2;
                    else
                        ptemp2 = psi2';
                    end
                    if isequal(psi1,psi2)
                        s1 = size(psi1);
                        switch s1(1)
                            case 2
                                if s1(2) == 1
                                    qs.psi = kron(ptemp,ptemp);
                                    qs.rho = qs.psi*qs.psi';
                                    disp('Two 2 x 1 inputs detected; will initialize with the product vector');
                                elseif s1(2) == 2
                                    qs.psi = kron(ptemp,ptemp);
                                    qs.rho = qs.psi*qs.psi';
                                    disp('Two 2 x 2 inputs detected; will initialize with the product density matrix');
                            otherwise
                                warning('Not a valid input dimension');
                            end
                        end
                    else 
                        warning('Inputs have different dimensions');
                    end

                case 1
                    s = size(psi1);

                    if isequal(s,[2,2])
                        qs.psi = [0;0;0;0];
                        qs.rho = kron(psi1,psi1);
                        disp('One 2 x 2 input detected; will initialize with duplicate density matrices');
                    if isequal(s,[4,4])
                        qs.psi = [0;0;0;0];
                        qs.rho = psi1;
                        disp('One 4 x 4 input detected; will initialize with the density matrix');
                    else
                        if iscolumn(psi1)==1
                            ptemp = psi1;
                        else
                            ptemp = psi1';
                        end
                        
                        switch s(1)
                            case 2
                                qs.psi = kron(ptemp,ptemp);
                                qs.rho = qs.psi*qs.psi';
                                disp('One 2 x 1 input detected; will initialize with duplicate vectors');
                            case 4
                                qs.psi = ptemp;
                                qs.rho = qs.psi*qs.psi';
                                disp('One 4 x 1 input detected; will initialize with the product vector');
                            otherwise
                                warning('Not a valid input dimension');
                            end
                        end
                    end

                case 0
                    qs.psi = [0;0;0;0];
                    qs.rho = qs.psi*qs.psi';
                    disp('No input detected; will initialize with zero vectors');
                otherwise
                    warning('Too many inputs!');
            end

            qs.ipsi = qs.psi;
            qs.irho = qs.rho;
            qs.cpsi = struct();
            qs.crho = struct();
        end
        function bool = isPure(qs)
            if (abs(trace(qs.rho)-1)>qs.epsilon)
                warning('Trace rho not equal to 1')
            end
            bool = abs(qs.purity-1)<qs.epsilon;
        end
        function bvec = get.bvec(qs)
            bvec = real([trace(qs.rho*qs.sx), trace(qs.rho*qs.sy), trace(qs.rho*qs.sz)]);
        end
        function purity = get.purity(qs)
            purity = sqrt(trace(qs.rho^2));
        end
        function H = runnoise(qs,H0,H1,N,m,sg)
            rot_rate = sg * randn(1,N);
            mu = m * ones(1,N);
            for k = 1:N
                H(:,:,k) = mu(k)*H0 + rot_rate(k)*H1;
            end
        end
        function H = hamnoise(qs,H0,H1,t,m,sg)
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
        function H = hamrunnoise(qs,H0,H1,t,N,m,sg)
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
        function evolve(qs,H,t)
            qs.cpsi.runs(1).psi= qs.psi;
            qs.crho.runs(1).rho= qs.rho;
            if ndims(H) < 3
                if length(t) ==1
                    qs.psi = qs.stevolve(qs.psi,H,t);
                    qs.rho = qs.stevolve(qs.rho,H,t);
                    qs.cpsi.runs(1).psi= cat(3,qs.cpsi.runs(1).psi,qs.psi);
                    qs.crho.runs(1).rho= cat(3,qs.crho.runs(1).rho,qs.rho);
                else
                    dt = diff(t);
                    for i=1:length(dt)
                        qs.psi = qs.stevolve(qs.psi,H,dt(i));
                        qs.rho = qs.stevolve(qs.rho,H,dt(i));
                        qs.cpsi.runs(1).psi= cat(3,qs.cpsi.runs(1).psi,qs.psi);
                        qs.crho.runs(1).rho= cat(3,qs.crho.runs(1).rho,qs.rho);
                    end
                end
            else
                dt = diff(t);
                if (length(dt) ~= size(H,3))
                    warning('Time step does not equal the number of Hamiltonians.');
                end
                for i=1:length(dt)
                    qs.psi = qs.stevolve(qs.psi,H(:,:,i),dt(i));
                    qs.rho = qs.stevolve(qs.rho,H(:,:,i),dt(i));
                    qs.cpsi.runs(1).psi= cat(3,qs.cpsi.runs(1).psi,qs.psi);
                    qs.crho.runs(1).rho= cat(3,qs.crho.runs(1).rho,qs.rho);
                end
            end
        end

        function nevolve(qs,H,t,N,steps)
            for z = 1:N
                qs.psi = qs.ipsi;
                qs.rho = qs.irho;
                qs.cpsi.runs(z).psi= qs.psi;
                qs.crho.runs(z).rho= qs.rho;
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
                    warning('Number of experiments does not equal the run number of Hamiltonians.');
                end

                for i=1:length(dt)
                    qs.psi = qs.stevolve(qs.psi,H(:,:,i,z),dt(i));
                    qs.rho = qs.stevolve(qs.rho,H(:,:,i,z),dt(i));
                    qs.cpsi.runs(z).psi= cat(3,qs.cpsi.runs(z).psi,qs.psi);
                    qs.crho.runs(z).rho= cat(3,qs.crho.runs(z).rho,qs.rho);
                end
            end
        end

        function s = measureSi(qs,Si)
            X = expm(1i*qs.sy*pi/4);
            Y = expm(1i*qs.sx*pi/4);
            if strcmpi(Si,'x') 
                s = trace(qs.sz*X*qs.rho);
            elseif strcmpi(Si,'y') 
                s = trace(qs.sz*Y*qs.rho);
            elseif strcmpi(Si,'z') 
                s = trace(qs.sz*qs.rho);
            else
                warning('There is no such axis.');
            end
        end
        
        function sn = measureSiN(qs,Si,N)
            X = expm(1i*qs.sy*pi/4);
            Y = expm(1i*qs.sx*pi/4);
            sn = zeros(size(qs.crho.runs(1).rho,3),N);
            for k=1:N
            rtemp = qs.crho.runs(k).rho;
                for j= 1:size(rtemp,3)
                    rt(:,:) = rtemp(:,:,j);
                    if strcmpi(Si,'x') 
                        sn(j,k) = real(trace(qs.sz*X*rt));
                    elseif strcmpi(Si,'y') 
                        sn(j,k) = real(trace(qs.sz*Y*rt));
                    elseif strcmpi(Si,'z') 
                        sn(j,k) = real(trace(qs.sz*rt));
                    else
                        warning('There is no such axis.');
                    end
                end
            end
        end

        function beta = fitesin(qs,t,sn,beta0)
            expsine = @(b,x) (b(1).*exp(-b(2).*x).*sin(b(3).*x+b(4)))';
            if length(beta0) ~= 4
                b0 = [1,0.01,pi/2,0]; %Initial guess which works for most cases
            else
                b0 = beta0;
            end
            beta = nlinfit(t,sn,expsine,beta0);
        end

        function h = plot(qs)
            clf;
            if (isPure(qs)<qs.epsilon)
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
            cmds=[cmds struct('type','vector','val',qs.bvec,'size',1,'color',[0 0 1])];
            h=plotBloch(cmds);
        end
        function h = plotev(qs,H,t)
            clf;
            if (isPure(qs)<qs.epsilon)
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
                evolve(qs,H,dt(i));
                cmd={};
                cmd=[cmd struct('type','vector','val',qs.bvec,'size',1,'color',[0 0 1])];
                h = plotBloch(cmd);
                pause(0.001);
            end
        end
    end
    
    methods (Static)
        function rho = bloch2rho(v)
              x = v(1)*qs.sx;
              y = v(2)*qs.sy;
              z = v(3)*qs.sz;  
              rho = (qs.si+x+y+z)/2;
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