function [Q,S]=MatQS(d)
  if d==2
    C0=[1 1 0;1 1 0; 0 0 0];
    C1=[2 0 0; 0 2 0; 0 0 1];
  else
    C0=zeros(6,6);C1=zeros(6,6);
    C0(1:3,1:3)=ones(3,3);
    for i=1:6
      C1(i,i)=1+(i<4);
    end
  end
    B=MatB(d);
    Q=cell(d,d);
    S=cell(d,d);
    for n=1:d
      for l=1:d
        Q{n,l}=B{n}'*C0*B{l};
        S{n,l}=B{n}'*C1*B{l};
      end
    end
end

function B=MatB(d)
  B=cell(d,1);
  if d==2
    for i=1:d
      B{i}=[(i==1), 0 ; 0, (i==2); (i==2) , (i==1)];
    end
  else
    for i=1:d
      BB=zeros(6,3);
      BB(1,1)=(i==1);BB(2,2)=(i==2);BB(3,3)=(i==3);
      BB(4,1)=(i==2);BB(4,2)=(i==1);
      BB(5,2)=(i==3);BB(5,3)=(i==2);
      BB(6,1)=(i==3);BB(6,3)=(i==1);
      B{i}=BB;
    end
  end
end