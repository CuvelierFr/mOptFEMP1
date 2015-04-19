function M=AssemblyMassWP1OptV1(Th,W)
% function M=AssemblyMassWP1OptV1(Th,W)
%   Assembly of the Weighted Mass Matrix using P1-Lagrange finite elements
%   - OptV1 version (see report).
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
%    M=AssemblyMassWP1OptV1(Th,W);
%
% Copyright (C) 2015  CJS (LAGA)
%   see README for details
  d=Th.d;ndfe2=(d+1)*(d+1);
  Ig=zeros(Th.nme*ndfe2,1);
  Jg=zeros(Th.nme*ndfe2,1);
  Kg=zeros(Th.nme*ndfe2,1);
  l=1;
  for k=1:Th.nme
    E=ElemMassWP1(Th.d,Th.vols(k),W(Th.me(:,k)));
    for il=1:d+1
      for jl=1:d+1
        Ig(l)=Th.me(il,k);
        Jg(l)=Th.me(jl,k);
        Kg(l)=E(il,jl);
        l=l+1;
      end
    end
  end
  M=sparse(Ig(:),Jg(:),Kg(:),Th.nq,Th.nq);
end
