function M=AssemblyStiffP1OptV1(Th)
% function M=AssemblyStiffP1OptV1(Th)
%   Assembly of the Stiffness Matrix using P1-Lagrange finite elements
%   - OptV1 version (see report).
%
% Parameters:
%  Th: mesh structure,
%
% Return values:
%  M: Global Stiffness matrix, Th.nq-by-Th.nq sparse matrix.
%
% Example:
%    Th=HyperCube(2,10);
%    M=AssemblyStiffP1OptV1(Th);
%
% Copyright (C) 2015  CJS (LAGA)
%   see README for details
  d=Th.d;ndfe2=(d+1)*(d+1);
  Ig=zeros(Th.nme*ndfe2,1);
  Jg=zeros(Th.nme*ndfe2,1);
  Kg=zeros(Th.nme*ndfe2,1);
  l=1;
  for k=1:Th.nme
    G=Gradient(Th.q(:,Th.me(:,k)),Th.vols(k));
    E=ElemStiffP1(Th.d,G,Th.vols(k));
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
