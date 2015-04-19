function K=AssemblyStiffElasP1OptV1(Th,lam,mu)
% function M=AssemblyStiffElasP1OptV1(Th,lam,mu)
%   Assembly of the Elastic Stiffness Matrix using P1-Lagrange finite elements
%   - OptV1 version (see report).
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
%    M=AssemblyStiffElasP1OptV1(Th,setFdata(lam,Th),setFdata(mu,Th));
%
% Copyright (C) 2015  CJS (LAGA)
%   see README for details
  d=Th.d;ndfe2=(d+1)*(d+1)*d*d;
  [Q,S]=MatQS(d);
  Ig=zeros(Th.nme*ndfe2,1);
  Jg=zeros(Th.nme*ndfe2,1);
  Kg=zeros(Th.nme*ndfe2,1);
  p=1;
  for k=1:Th.nme
    mek=Th.me(:,k);
    G=Gradient(Th.q(:,mek),Th.vols(k));
    Cl=Th.vols(k)*sum(lam(mek))/(d+1);
    Cm=Th.vols(k)*sum(mu(mek))/(d+1);
    for l=1:d
      for n=1:d
	for il=1:d+1
	  r=d*(mek(il)-1)+l;
	  for jl=1:d+1
	    s=d*(mek(jl)-1)+n;
	    Kg(p,:)=Cl*G(jl,:)*Q{n,l}*G(il,:)'+Cm*G(jl,:)*S{n,l}*G(il,:)';
	    Ig(p,:)=r;
	    Jg(p,:)=s;
	    p=p+1;
	  end
	end
      end
    end
  end
  K=sparse(Ig(:),Jg(:),Kg(:),d*Th.nq,d*Th.nq);
end
