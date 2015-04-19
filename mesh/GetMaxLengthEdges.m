function h=GetMaxLengthEdges(q,me)
  ne=size(me,1);
  h=0;
  for i=1:ne
    for j=i+1:ne
      h=max(h,max(sum((q(:,me(i,:))-q(:,me(j,:))).^2,1)));
    end
  end
  h=sqrt(h);
end