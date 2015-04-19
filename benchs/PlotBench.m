function PlotBench(Versions,Lndof,Tcpu,ctitle)
  nV=length(Versions);
  assert((size(Tcpu,2)==nV)&&(length(Lndof)==size(Tcpu,1)))
  Colors= select_colors(nV);
  figure(1)
  clf
  loglog(Lndof,Tcpu(:,1),'LineWidth',2,'Color',Colors(1,:));
  hold on
  for v=2:size(Tcpu,2);
    loglog(Lndof,Tcpu(:,v),'LineWidth',2,'Color',Colors(v,:));
  end
  legend(Versions{:})
  xlabel('n_{dof}')
  ylabel('cputime(s)')
  title(ctitle);
end