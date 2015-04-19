function M=AssemblyStiffP1OptV(Th)
% function M=AssemblyStiffP1OptV(Th)
%   Assembly of the Stiffness Matrix using P1-Lagrange finite elements
%   - OptV version (see report).
%
% Parameters:
%  Th: mesh structure,
%
% Return values:
%  M: Global Stiffness matrix, Th.nq-by-Th.nq sparse matrix.
%
% Example:
%    Th=HyperCube(2,10);
%    M=AssemblyStiffP1OptV(Th);
%
% Copyright (C) 2015  CJS (LAGA)
%   see README for details
  d=Th.d;
  G=gradientVecOpt(Th.q,Th.me,Th.vols);
  M=spalloc(Th.nq,Th.nq,d*d*Th.nq);
  for il=1:d+1
    for jl=1:d+1
      M=M+sparse(Th.me(il,:),Th.me(jl,:),Th.vols.*dotVecG(G,il,jl),Th.nq,Th.nq);
    end
  end
end
