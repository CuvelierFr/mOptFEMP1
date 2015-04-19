function M=AssemblyStiffP1OptV2(Th)
% function M=AssemblyStiffP1OptV2(Th)
%   Assembly of the Stiffness Matrix using P1-Lagrange finite elements
%   - OptV2 version (see report).
%
% Parameters:
%  Th: mesh structure,
%
% Return values:
%  M: Global Stiffness matrix, Th.nq-by-Th.nq sparse matrix.
%
% Example:
%    Th=HyperCube(2,10);
%    M=AssemblyStiffP1OptV2(Th);
%
% Copyright (C) 2015  CJS (LAGA)
%   see README for details
  d=Th.d;ndfe2=(d+1)*(d+1);
  G=gradientVecOpt(Th.q,Th.me,Th.vols);
  Ig=zeros(Th.nme,ndfe2);
  Jg=zeros(Th.nme,ndfe2);
  Kg=zeros(Th.nme,ndfe2);
  l=1;
  for il=1:d+1
    for jl=1:d+1
      Ig(:,l)=Th.me(il,:);
      Jg(:,l)=Th.me(jl,:);
      for i=1:d
        Kg(:,l)=Kg(:,l)+G(:,il,i).*G(:,jl,i);
      end
      Kg(:,l)=Th.vols.*Kg(:,l);
      l=l+1;
    end
  end
  M=sparse(Ig(:),Jg(:),Kg(:),Th.nq,Th.nq);
end
