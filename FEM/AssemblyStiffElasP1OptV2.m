function S=AssemblyStiffElasP1OptV2(Th,lam,mu)
% function M=AssemblyStiffElasP1OptV2(Th,lam,mu)
%   Assembly of the Elastic Stiffness Matrix using P1-Lagrange finite elements
%   - OptV2 version (see report).
%
% Parameters:
%  Th: mesh structure,
%  lam : Array containing Lamé's parameter of incompressibility "lambda" values at vertices,
%        1-by-Th.nq array.
%  mu  : Array containing Lamé's parameter of rigidity "mu" values at vertices,
%        1-by-Th.nq array.
% Return values:
%  M: Global Elastic Stiffness matrix, ndof-by-ndof sparse matrix (ndof=Th.d*Th.nq)
%
% Example:
%    Th=HyperCube(2,10);
%    mu=@(x,y) 7.2414e+05*(x<=0.5)+ 8.2031e+05*(x>0.5);
%    lam=@(x,y) 6.5172e+06*(x<=0.5)+ 1.0440e+06*(x>0.5);
%    M=AssemblyStiffElasP1OptV2(Th,setFdata(lam,Th),setFdata(mu,Th));
%
% Copyright (C) 2015  CJS (LAGA)
%   see README for details
  d=Th.d;N=(d+1)*(d+1)*d*d;
  assert((d==2)||(d==3))
  nme=Th.nme;
  [Q,S]=MatQS(d);
  lams=sum(lam(Th.me),1).*Th.vols'/(d+1);
  mus=sum(mu(Th.me),1).*Th.vols'/(d+1);
  G=gradientVecOpt(Th.q,Th.me,Th.vols);
  Kg=zeros(nme,N);
  Ig=zeros(nme,N);
  Jg=zeros(nme,N);
  l=1;
  for i=1:d
    for j=1:d
      for jl=1:d+1
        for il=1:d+1
          Kg(:,l)=lams.*dotMatVecG(Q{i,j},G,il,jl)+mus.*dotMatVecG(S{i,j},G,il,jl);
          Ig(:,l)=d*(Th.me(il,:)-1)+i;
          Jg(:,l)=d*(Th.me(jl,:)-1)+j;
          l=l+1;
        end
      end
    end
  end
  S=sparse(Ig(:),Jg(:),Kg(:),d*Th.nq,d*Th.nq);
end