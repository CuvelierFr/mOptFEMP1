function E=ElemMassWP1(d,vol,W)
  ndfe=d+1;
  ws=sum(W);
  C=vol/((d+1)*(d+2)*(d+3));
  E=zeros(ndfe,ndfe);
  for il=1:ndfe
    E(il,il)=2*C*(ws+2*W(il));
    for jl=il+1:ndfe
      E(il,jl)=C*(ws+W(il)+W(jl));
      E(jl,il)=E(il,jl);
    end
  end
end