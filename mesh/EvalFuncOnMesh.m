function fd=EvalFuncOnMesh(f,q,t)
  if nargin==2, t=[];end
  d=size(q,1);
  if isfhandle(f)
    if isempty(t)
      fd=eval(['f(q(1,:)',sprintf(',q(%d,:)',2:d),').'';']) ; 
    else
      fd=eval(['f(q(1,:)',sprintf(',q(%d,:)',2:d),',t).'';']) ;
    end
    if (length(fd)==1) % f == constant
      fd=fd*ones(size(q,2),1);
    end
    return
  end
  if isnumeric(f) && isscalar(f)
    fd=f*ones(size(q,2),1);
  end
    
end