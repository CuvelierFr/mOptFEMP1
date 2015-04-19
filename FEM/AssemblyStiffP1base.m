function M=AssemblyStiffP1base(Th)
% function M=AssemblyStiffP1base(Th)
%   Assembly of the Stiffness Matrix using P1-Lagrange finite elements
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
%    M=AssemblyStiffP1base(Th);
%
% Copyright (C) 2015  CJS (LAGA)
%   see README for details
  d=Th.d;ndfe2=(d+1)*(d+1);
  M=spalloc(Th.nq,Th.nq,ndfe2*Th.nq);
  for k=1:Th.nme
    G=Gradient(Th.q(:,Th.me(:,k)),Th.vols(k));
    E=ElemStiffP1(Th.d,G,Th.vols(k));
    for il=1:d+1
      i=Th.me(il,k);
      for jl=1:d+1
        j=Th.me(jl,k);
        M(i,j)=M(i,j)+E(il,jl);
      end
    end
  end
end
