function M=AssemblyMassP1OptV1(Th)
% function M=AssemblyMassP1OptV1(Th)
%   Assembly of the Mass Matrix using P1-Lagrange finite elements
%   - OptV1 version (see report).
%
% Parameters:
%  Th: mesh structure,
%
% Return values:
%  M: Global mass matrix, Th.nq-by-Th.nq sparse matrix.
%
% Example:
%    Th=HyperCube(2,10);
%    M=AssemblyMassP1OptV1(Th);
%
% Copyright (C) 2015  CJS (LAGA)
%   see README for details
  d=Th.d;ndfe2=(d+1)*(d+1);
  Ig=zeros(Th.nme*ndfe2,1);
  Jg=zeros(Th.nme*ndfe2,1);
  Kg=zeros(Th.nme*ndfe2,1);
  l=1;
  for k=1:Th.nme
    E=ElemMassP1(Th.d,Th.vols(k));
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
