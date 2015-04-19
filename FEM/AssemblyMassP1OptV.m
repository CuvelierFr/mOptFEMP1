function M=AssemblyMassP1OptV(Th)
% function M=AssemblyMassP1OptV(Th)
%   Assembly of the Mass Matrix using P1-Lagrange finite elements
%   - OptV version (see report).
%
% Parameters:
%  Th: mesh structure,
%
% Return values:
%  M: Global mass matrix, Th.nq-by-Th.nq sparse matrix.
%
% Example:
%    Th=HyperCube(2,10);
%    M=AssemblyMassP1OptV(Th);
%
% Copyright (C) 2015  CJS (LAGA)
%   see README for details
  d=Th.d;
  M=spalloc(Th.nq,Th.nq,d*d*Th.nq);
  for il=1:d+1
    for jl=1:d+1
      Kg=(1+(il==jl))/((d+1)*(d+2))*Th.vols;
      M=M+sparse(Th.me(il,:),Th.me(jl,:),Kg,Th.nq,Th.nq);
    end
  end
end
