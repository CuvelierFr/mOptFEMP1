function X=dotMatVecG(A,G,il,jl)
  nme=size(G,1);
  d=size(A,1);
  X=zeros(nme,1);
  for i=1:d
    for j=1:d
      if (A(i,j)~=0)
        X=X+A(i,j)*(G(:,il,i).*G(:,jl,j));
      end
    end
  end
  X=X';
end