function D=detVec(A)
% A : N-by-d-by-d
% D : N-by-1
% D(i) = det A(i,:,:)
n=size(A);
d=n(2);
if length(n)==3
  N=n(1);assert(n(2)==n(3))
else
  if (d==1)
    N=n(1);
  else
    N=1;assert(n(1)==n(2))
  end
end

D=zeros(N,1);
if d==1 
  D=A(:,1,1);   
else
  if d==2
   D=(A(:,1,1).*A(:,2,2)-A(:,1,2).*A(:,2,1));   
  else
    for ii=1:d
      D=D+(-1)^(ii+1)*(A(:,1,ii).*detVec(A(:,2:d,[1:ii-1,ii+1:d])));
    end
  end
end