function vol=ComputeVolVec(d,q,me)
  % d-simplex in R^n (d<=n)
  nme=size(me,2);
  n=size(q,1);assert(d<=n && size(me,1)==d+1);
  X=zeros(d,n,nme);
  for i=1:d 
    X(i,:,:)=q(:,me(i+1,:))-q(:,me(1,:));
  end
  V=zeros(nme,d,d);
  for i=1:d
    V(:,i,i)=sum(reshape(X(i,:,:).*X(i,:,:),n,nme),1);
    for j=i+1:d
      V(:,i,j)=sum(reshape(X(i,:,:).*X(j,:,:),n,nme),1);
      V(:,j,i)=V(:,i,j);
    end
  end
  vol=sqrt(detVec(V)')/factorial(d);
end 

