function M=AssemblyMassP1OptVS(Th)
% function M=AssemblyMassP1OptVS(Th)
%   Assembly of the Mass Matrix using P1-Lagrange finite elements
%   - OptVS version (see report).
%
% Parameters:
%  Th: mesh structure,
%
% Return values:
%  M: Global mass matrix, Th.nq-by-Th.nq sparse matrix.
%
% Example:
%    Th=HyperCube(2,10);
%    M=AssemblyMassP1OptVS(Th);
%
% Copyright (C) 2015  CJS (LAGA)
%   see README for details
  d=Th.d;Kg=1/((d+1)*(d+2))*Th.vols;
  M=spalloc(Th.nq,Th.nq,d*d*Th.nq);
  for il=1:d+1
    for jl=il+1:d+1
      M=M+sparse(Th.me(il,:),Th.me(jl,:),Kg,Th.nq,Th.nq);
    end
  end
  M=M+M';
  Kg=2*Kg;
  for il=1:d+1
    M=M+sparse(Th.me(il,:),Th.me(il,:),Kg,Th.nq,Th.nq);
  end
end
