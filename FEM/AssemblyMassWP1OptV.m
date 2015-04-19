function M=AssemblyMassWP1OptV(Th,W)
% function M=AssemblyMassWP1OptV(Th,W)
%   Assembly of the Weighted Mass Matrix using P1-Lagrange finite elements
%   - OptV version (see report).
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
%    M=AssemblyMassWP1OptV(Th,W);
%
% Copyright (C) 2015  CJS (LAGA)
%   see README for details
  d=Th.d;C=1/prod(d+1:d+3)*Th.vols';
  Wme=W(Th.me);Ws=sum(Wme,1);
  M=spalloc(Th.nq,Th.nq,d*d*Th.nq);
  for il=1:d+1
    for jl=1:d+1
      Kg=(1+(il==jl))*C.*(Ws+Wme(il,:)+Wme(jl,:));
      M=M+sparse(Th.me(il,:),Th.me(jl,:),Kg,Th.nq,Th.nq);
    end
  end
end
