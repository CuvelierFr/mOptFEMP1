function M=AssemblyMassWP1base(Th,W)
% function M=AssemblyMassWP1base(Th,W)
%   Assembly of the Weighted Mass Matrix using P1-Lagrange finite elements
%   - base version (see report).
%
% Parameters:
%  Th: mesh structure,
%  W : Array containing weigthed values at vertices,
%      1-by-Th.nq array.
%
% Return values:
%  M: Global weighted mass matrix, Th.nq-by-Th.nq sparse matrix.
%
% Example:
%    Th=HyperCube(2,10);w=@(x,y) cos(x+y);
%    W=w(Th.q(1,:),Th.q(2,:));
%    M=AssemblyMassWP1base(Th,W);
%
% Copyright (C) 2015  CJS (LAGA)
%   see README for details
  d=Th.d;ndfe2=(d+1)*(d+1);
  M=spalloc(Th.nq,Th.nq,ndfe2*Th.nq);
  for k=1:Th.nme
    E=ElemMassWP1(Th.d,Th.vols(k),W(Th.me(:,k)));
    for il=1:d+1
      i=Th.me(il,k);
      for jl=1:d+1
        j=Th.me(jl,k);
        M(i,j)=M(i,j)+E(il,jl);
      end
    end
  end
end
