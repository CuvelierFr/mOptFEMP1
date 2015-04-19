function M=AssemblyStiffElasP1OptV(Th,lam,mu)
% function M=AssemblyStiffElasP1OptV(Th,lam,mu)
%   Assembly of the Elastic Stiffness Matrix using P1-Lagrange finite elements
%   - OptV version (see report).
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
%    M=AssemblyStiffElasP1base(Th,setFdata(lam,Th),setFdata(mu,Th));
%
% Copyright (C) 2015  CJS (LAGA)
%   see README for details
  d=Th.d;N=(d+1)*d;
  assert((d==2)||(d==3))
  nme=Th.nme;
  [Q,S]=MatQS(d);
  lams=sum(lam(Th.me),1).*Th.vols'/(d+1);
  mus=sum(mu(Th.me),1).*Th.vols'/(d+1);
  G=gradientVecOpt(Th.q,Th.me,Th.vols);
  M=spalloc(d*Th.nq,d*Th.nq,N*N*Th.nq);
  for i=1:d
    for j=1:d
      for il=1:d+1
        ii=d*(il-1)+i;
        I=d*(Th.me(il,:)-1)+i;
        for jl=1:d+1
          jj=d*(jl-1)+j;
            M=M+sparse(I,d*(Th.me(jl,:)-1)+j, ...
                       lams.*dotMatVecG(Q{i,j},G,il,jl)+mus.*dotMatVecG(S{i,j},G,il,jl),d*Th.nq,d*Th.nq);
        end
      end
    end
  end
end
