function M=AssemblyMassP1OptV2(Th)
% function M=AssemblyMassP1OptV2(Th)
%   Assembly of the Mass Matrix using P1-Lagrange finite elements
%   - OptV2 version (see report).
%
% Parameters:
%  Th: mesh structure,
%
% Return values:
%  M: Global mass matrix, Th.nq-by-Th.nq sparse matrix.
%
% Example:
%    Th=HyperCube(2,10);
%    M=AssemblyMassP1OptV2(Th);
%
% Copyright (C) 2015  CJS (LAGA)
%   see README for details
  d=Th.d;ndfe2=(d+1)*(d+1);
  Ig=zeros(Th.nme,ndfe2);
  Jg=zeros(Th.nme,ndfe2);
  Kg=zeros(Th.nme,ndfe2);
  l=1;
  for il=1:d+1
    for jl=1:d+1
      Ig(:,l)=Th.me(il,:);
      Jg(:,l)=Th.me(jl,:);
      Kg(:,l)=(1+(il==jl))/((d+1)*(d+2))*Th.vols;
      l=l+1;
    end
  end
  M=sparse(Ig(:),Jg(:),Kg(:),Th.nq,Th.nq);
end
