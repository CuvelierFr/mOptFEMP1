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