function fd=EvalFuncVFOnMesh(f,q)
  assert(iscell(f));
  m=length(f);
  nq=size(q,2);
  fd=zeros(m,nq);E=1:nq;
  I=1:nq;
  for i=1:m
    fd(i,I)=EvalFuncOnMesh(f{i},q);
  end
  fd=fd(:);
end