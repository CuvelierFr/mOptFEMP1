function E=ElemMassP1(d,vol)
  ndfe=d+1;
  C=vol/((d+1)*(d+2));
  E=(vol/((d+1)*(d+2)))*ones(ndfe,ndfe);
  for il=1:ndfe
    E(il,il)=2*E(il,il);
  end
end