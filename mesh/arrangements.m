function A=arrangements(V,d,varargin)
  p = inputParser; 
  p.addParamValue('equal',0,@isscalar);
  p.parse(varargin{:});
  k=p.Results.equal;
  if d==1
    A=V;
  else
    B=arrangements(V,d-1);
    N=size(B,1);
    nv=length(V);
    A=[];
    for i=1:nv
      P=V(i)*ones(N,1);
      A=[A;[B,P]];
    end
  end
  if k>=1
    I=find(sum(A,2)==k);
    A=A(I,:);
    
  end
end