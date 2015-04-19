function E=ElemStiffP1(d,G,vol)
  ndfe=d+1;
  E=zeros(ndfe,ndfe);
  for il=1:ndfe
    E(il,il)=vol*G(il,:)*G(il,:)';
    for jl=il+1:ndfe
      E(il,jl)=vol*G(jl,:)*G(il,:)';
      E(jl,il)=E(il,jl);
    end
  end
end