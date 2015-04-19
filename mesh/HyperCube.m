function Th=HyperCube(d,N,varargin)
% function Th=HyperCube(d,N,trans)
%   Build d-simplicial mesh of the hypercube [0,1]^d. Can be transform using 
%   optional parameter trans.
%
% Parameters:
%   d     : space dimension.
%   N     : numbers of discretisation in each space direction.
%           1-by-d array of integer or integer for same number in each space direction.
%   trans : (optional). function used to transfrom vertices of hypercube mesh.
%
% Return values:
%   Th    : d-simplicial mesh structure of fields (see report)
%      d     space dimension and dimension of simplices (mesh elements).
%      nq    number of vertices
%      nme   number of elements (d-simplices)
%      nbe   number of boundary elements ((d-1)-simplices)
%      q     d-by-nq array of vertices coordinates
%      me    (d+1)-by-nme connectivity array for mesh elements
%      be    d-by-nbe connectivity array for boundary elements
%      bel   1-by-nbe boundary elements label.
%      vols  1-by-nme array of mesh elements volumes
%      ...
% 
% Example:
%   Th=HyperCube(2,50);
%   Th=HyperCube(2,[100,10],@(q) [20*q(1,:);2*q(2,:)-1]);
%   Th=HyperCube(2,[100,10],@(q) [20*q(1,:);(2*q(2,:)-1)+cos(2*pi*q(1,:))]);
%   Th=HyperCube(3,[10,100,10],@(q) [(2*q(1,:)-1);20*q(2,:);(2*q(3,:)-1)+cos(2*pi*q(2,:))]);
  trans=[];
  if nargin==3, trans=varargin{:};  end
  if length(N)==1
    [q,me]=HypercubeKuhn(N*ones(1,d));
  else
    [q,me]=HypercubeKuhn(N);
  end

  Th.d=d;
  Th.q=q;
  Th.me=me;
  Th.nq=size(q,2);
  Th.nme=size(me,2);

  % Set Boundary
  V=nchoosek(1:d+1,d); % Get combinaisons
  BE=[];
  for i=1:size(V,1)    %trans=varargin{1}; % transformation function

    BE=[BE,Th.me(V(i,:),:)];
  end
  be1=sort(BE,1)';
  be=unique(be1,'rows')';
  nbe=size(be,2);
  bel=zeros(1,nbe);

  for i=1:d
    Qb{i}=Th.q(:,be(i,:));
  end

  tol=1e-12;

  label=1;

  for i=1:d % hypercube number of faces  2^d
    I=1:nbe;
    for j=1:d
      I=intersect(I,find(abs(Qb{j}(i,:)-0)<tol));
    end
    bel(I)=label;label=label+1;
    I=1:nbe;
    for j=1:d
      I=intersect(I,find(abs(Qb{j}(i,:)-1)<tol));
    end
    bel(I)=label;label=label+1;
  end
  I=find(bel==0);
  J=setdiff(1:nbe,I);

  [Bel1,ii]=sort(bel(J));
  Be1=be(:,J(ii));
  Th.be=Be1;
  Th.bel=Bel1;
  Th.nbe=length(Bel1);
  if isempty(trans)
    Th.vols=ones(Th.nme,1)/Th.nme;
  else
    Th.q=trans(Th.q);
    Th.vols=ComputeVolVec(d,Th.q,Th.me);
  end
  Th.h=GetMaxLengthEdges(q,me);
end

function [q,me]=HypercubeKuhn(N)
  d=length(N);
  BaseSimplex=[zeros(1,d);tril(ones(d))]';
  P=perms(1:d); 
  nmeBase=factorial(d);
  nt=nmeBase;
  meBase=zeros(d+1,nmeBase);
  for i=1:d
    T{i}=[0,1];
  end
  qBase=HypercubeVerticesV2(T{:});

  A=2.^(0:d-1);
  for i=1:nmeBase
    meBase(:,i)=A*BaseSimplex(P(i,:),:)+1;  
  end
  for k=1:nmeBase % Chech d-simplices orientations
    s=orientation(qBase(:,meBase(:,k)));
    if (s==-1)
      tmp=meBase(1,k);;
      meBase(1,k)=meBase(2,k);
      meBase(2,k)=tmp;
    end
  end
  % Maillage 
  for i=1:d
    T{i}=linspace(0,1,N(i));
  end
  q=HypercubeVerticesV2(T{:});

  Nhypercube=prod(N-1);
  C=[1,cumprod(N(1:(end-1)))];
  J=C*qBase;
  nme=Nhypercube*nmeBase; 
  I=[1:nmeBase];
  
  for i=1:d
    T{i}=[1:(N(i)-1)]-1;
  end
  qInd=HypercubeVerticesV2(T{:});
  ki=C*qInd+1;
  
  nJ=length(J); % 2^d number of points in d-hypercube
  Ind=zeros(nJ,Nhypercube);
  for l=1:nJ
    Ind(l,:)=ki+J(l);
  end
  K=1:nmeBase:nme;
  me=zeros(d+1,nme);
  for l=1:nmeBase
    me(:,K)=Ind(meBase(:,l),:);
    K=K+1;
  end
end 

function q=HypercubeVerticesV2(varargin)
  d=length(varargin);
  N = cellfun('length',varargin);
  c=cell(1,d);
  [c{1:d}]=ndgrid(varargin{:});
  q=zeros(d,prod(N));
  for i=1:d
    q(i,:)=c{i}(:);
  end
end

function s=orientation(ql)
  d=size(ql,1);
  D=[ones(1,d+1);ql];
  s=sign(det(D));
end
