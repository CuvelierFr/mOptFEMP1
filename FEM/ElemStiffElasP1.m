function E=ElemStiffElasP1(d,G,vol,lam,mu,Q,S)
  ndfe2=(d+1)*(d+1)*d*d;
  E=zeros(ndfe2,ndfe2);
  for 

end

% vectorial <A*u',v'>
% A d-by-d
% u,v 1-by-d
% x 1-by-nme
function x=dotMat(A,u,v)
  v*A*u'

end