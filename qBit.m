classdef qBit < handle
    properties
        rho;
        dephase;
        
    end
    properties
    end
end

function s=def(s,def)
  f=fields(def);
  for i=1:length(f)
      l=f{i};
      if ~isfield(s,l) || isempty(s.(l))
          s.(l)=def.(l);
      end
  end
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