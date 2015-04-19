function X=dotVecG(G,il,jl)
  nme=size(G,1);
  d=size(G,3);
  X=zeros(nme,1);
  for i=1:d
    X=X+G(:,il,i).*G(:,jl,i);
  end
end