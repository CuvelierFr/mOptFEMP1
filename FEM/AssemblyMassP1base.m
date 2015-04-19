function M=AssemblyMassP1base(Th)
% function M=AssemblyMassP1base(Th)
%   Assembly of the Mass Matrix using P1-Lagrange finite elements
%   - base version (see report).
%
% Parameters:
%  Th: mesh structure,
%
% Return values:
%  M: Global mass matrix, Th.nq-by-Th.nq sparse matrix.
%
% Example:
%    Th=HyperCube(2,10);
%    M=AssemblyMassP1base(Th);
%
% Copyright (C) 2015  CJS (LAGA)
%   see README for details
  d=Th.d;ndfe2=(d+1)*(d+1);
  M=spalloc(Th.nq,Th.nq,ndfe2*Th.nq);
  for k=1:Th.nme
    E=ElemMassP1(Th.d,Th.vols(k));
    for il=1:d+1
      i=Th.me(il,k);
      for jl=1:d+1
        j=Th.me(jl,k);
        M(i,j)=M(i,j)+E(il,jl);
      end
    end
  end
end
