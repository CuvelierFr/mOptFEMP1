function M=AssemblyStiffP1OptVS(Th)
% function M=AssemblyStiffP1OptVS(Th)
%   Assembly of the Stiffness Matrix using P1-Lagrange finite elements
%   - OptVS version (see report).
%
% Parameters:
%  Th: mesh structure,
%
% Return values:
%  M: Global Stiffness matrix, Th.nq-by-Th.nq sparse matrix.
%
% Example:
%    Th=HyperCube(2,10);
%    M=AssemblyStiffP1OptVS(Th);
%
% Copyright (C) 2015  CJS (LAGA)
%   see README for details
  d=Th.d;
  G=gradientVecOpt(Th.q,Th.me,Th.vols);
  M=spalloc(Th.nq,Th.nq,d*d*Th.nq);
  for il=1:d+1
    for jl=il+1:d+1
      M=M+sparse(Th.me(il,:),Th.me(jl,:),Th.vols.*dotVecG(G,il,jl),Th.nq,Th.nq);
    end
  end
  M=M+M';
  for il=1:d+1
    M=M+sparse(Th.me(il,:),Th.me(il,:),Th.vols.*dotVecG(G,il,il),Th.nq,Th.nq);
  end
end
