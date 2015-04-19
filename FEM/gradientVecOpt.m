function G=gradientVecOpt(q,me,vol)
% G : nme x (d+1) x d
    d=size(me,1)-1;
nme=size(me,2);
switch d
  case 1
    G=zeros(nme,2,1);
    G(:,1,1)=-1./vol;
    G(:,2,1)=1./vol;
  case 2
    D12=q(:,me(1,:))-q(:,me(2,:));
    D13=q(:,me(1,:))-q(:,me(3,:));
    D23=q(:,me(2,:))-q(:,me(3,:));
  
    if size(vol,2)==1
      vol2=2*vol';
    else
      vol2=2*vol;
    end
  
    G=zeros(nme,3,2);
    G(:,1,1)=D23(2,:)./vol2;
    G(:,1,2)=-D23(1,:)./vol2;
    G(:,2,1)=-D13(2,:)./vol2;
    G(:,2,2)=D13(1,:)./vol2;
    G(:,3,1)=D12(2,:)./vol2;
    G(:,3,2)=-D12(1,:)./vol2;
  case 3
    
    D12=q(:,me(1,:))-q(:,me(2,:));
    D13=q(:,me(1,:))-q(:,me(3,:));
    D14=q(:,me(1,:))-q(:,me(4,:));
    D23=q(:,me(2,:))-q(:,me(3,:));
    D24=q(:,me(2,:))-q(:,me(4,:));
    
    if size(vol,2)==1
      C=1./(6*vol');
    else
      C=1./(6*vol);
    end
    
    G=zeros(nme,4,3);
    
    G(:,1,1)=(-D23(2,:).*D24(3,:) + D23(3,:).*D24(2,:)).*C;
    G(:,1,2)=( D23(1,:).*D24(3,:) - D23(3,:).*D24(1,:)).*C;
    G(:,1,3)=(-D23(1,:).*D24(2,:) + D23(2,:).*D24(1,:)).*C;

    G(:,2,1)=( D13(2,:).*D14(3,:) - D13(3,:).*D14(2,:)).*C;
    G(:,2,2)=(-D13(1,:).*D14(3,:) + D13(3,:).*D14(1,:)).*C;
    G(:,2,3)=( D13(1,:).*D14(2,:) - D13(2,:).*D14(1,:)).*C;

    G(:,3,1)=(-D12(2,:).*D14(3,:) + D12(3,:).*D14(2,:)).*C;
    G(:,3,2)=( D12(1,:).*D14(3,:) - D12(3,:).*D14(1,:)).*C;
    G(:,3,3)=(-D12(1,:).*D14(2,:) + D12(2,:).*D14(1,:)).*C;
    
    G(:,4,1)=( D12(2,:).*D13(3,:) - D12(3,:).*D13(2,:)).*C;
    G(:,4,2)=(-D12(1,:).*D13(3,:) + D12(3,:).*D13(1,:)).*C;
    G(:,4,3)=( D12(1,:).*D13(2,:) - D12(2,:).*D13(1,:)).*C;
  
  otherwise
    d=size(me,1)-1;
    nme=size(me,2);
    getQ=@(i,il) q(i,me(il,:));
    Grad=[-ones(d,1),eye(d)];
    
    K=zeros(d,d,nme);
    I=zeros(d,d,nme);
    J=zeros(d,d,nme);
    
    ii=d*[0:(nme-1)];
    for i=1:d % Composantes
      for j=1:d
	K(i,j,:)=getQ(i,j+1)-getQ(i,1);
	I(i,j,:)=ii+j;
	J(i,j,:)=ii+i;
      end
    end
    b=repmat(Grad,nme,1);
    spD=sparse(I(:),J(:),K(:),d*nme,d*nme);
    G=spD\b;
    G=reshape(G,[d,nme,d+1]);
    G=permute(G,[2,3,1]);
end