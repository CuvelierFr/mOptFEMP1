function M=AssemblyMassWP1OptV2(Th,W)
% function M=AssemblyMassWP1OptV2(Th,W)
%   Assembly of the Weighted Mass Matrix using P1-Lagrange finite elements
%   - OptV2 version (see report).
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
%    M=AssemblyMassWP1OptV2(Th,W);
%
% Copyright (C) 2015  CJS (LAGA)
%   see README for details
  d=Th.d;C=1/prod(d+1:d+3);ndfe2=(d+1)*(d+1);
  Wme=W(Th.me);Ws=sum(Wme,1);
  Ig=zeros(Th.nme,ndfe2);
  Jg=zeros(Th.nme,ndfe2);
  Kg=zeros(Th.nme,ndfe2);
  l=1;
  for il=1:d+1
    for jl=1:d+1
      Ig(:,l)=Th.me(il,:);
      Jg(:,l)=Th.me(jl,:);
      Kg(:,l)=C*(1+(il==jl))*Th.vols'.*(Ws+Wme(il,:)+Wme(jl,:));
      l=l+1;
    end
  end
  M=sparse(Ig(:),Jg(:),Kg(:),Th.nq,Th.nq);
end
