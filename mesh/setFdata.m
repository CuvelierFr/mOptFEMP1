function Fh=setFdata(f,Th,varargin)
% 
  if isempty(f), Fh=[]; return; end
  if (nargin==3), Num=varargin{1};else Num=1; end
  Fh=[];
  if (isnumeric(f)&&isscalar(f)),  Fh=f*ones(Th.nq,1);           return, end
  if (isnumeric(f)&&(size(f,1)==Th.nq)); Fh=f; return, end
  if isfhandle(f), Fh=EvalFuncOnMesh(f,Th.q);return, end
  if iscell(f) % Vector Field case
    m=length(f);
    Fh=zeros(m*Th.nq,1);I=1:Th.nq;
    VFInd=getVFindices(0,m,Th.nq);
    for i=1:m
      Fh(VFInd(I,i))=setFdata(f{i},Th);
    end
    return
  end
  error('Unknown f type');
end

function VFInd=getVFindices(Num,m,nq)
  if (Num==1) 
    VFInd=@(I,i) I+(i-1)*nq;
  else 
    VFInd=@(I,i) m*I-(m-i);
  end
end