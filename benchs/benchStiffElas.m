function [Lndof,Tcpu]=benchStiffElas(d,Versions,LN,varargin)
%  benchMass : Benchmark for AssemblyStiffElasP1<version> functions using HyperCube(d,N) to generate meshes.
%
%  USAGE:
%    [Lndof,Tcpu]=benchStiffElas(d,Versions,LN)
%    [Lndof,Tcpu]=benchStiffElas(d,Versions,LN,'nbruns',3)
%    [Lndof,Tcpu]=benchStiffElas(d,Versions,LN,'isplot',true)
%    [Lndof,Tcpu]=benchStiffElas(d,Versions,LN,'lam',func_handle1,'mu',func_handle2)
%
%  INPUTS:
%    d        : space dimension (2 or 3)
%    Versions : List of Assembly version (cell array),  
%               for example {'OptVS','OptV','OptV2','OptV1'}
%    LN       : List of N values (1D array of integer values)
%    <Optional>
%    'nbruns' : Number of runs
%    'isplot' : if true, plot results 
%    'lam'      : Lamé's parameter of incompressibility lambda, function or constant value
%    'mu'       : Lamé's parameter of rigidity mu , function or constant value
%
%  OUTPUTS:
%    Lndof    : List of size matrices (#LN-by-1 array)
%    Tcpu     : #LN-by-#Versions array, cputimes for each Assembly versions and 
%               each HyperCube meshes. Tcpu(i,v), is the mean cputimes for
%               Assembly with version Versions{v} and HyperCube(d,LN[i]) mesh.
%
%  SAMPLES:
%    [Lndof,Tcpu]=benchStiffElas(2,{'OptVS','OptV','OptV2','OptV1'},[25:25:125]);
%    [Lndof,Tcpu]=benchStiffElas(3,{'OptVS','OptV','OptV2','OptV1'},[5:5:25],'nbruns',1);
%    [Lndof,Tcpu]=benchStiffElas(2,{'OptVS','OptV','OptV2'},[100:100:500],'nbruns',1,'isplot',true);
%
  assert( (d==2) || (d==3) )
  % Poisson and Lamé's parameters for rubber
  E = 21e5; nu = 0.45; %nu = 0.28; 
  mu= E/(2*(1+nu));
  lam = E*nu/((1+nu)*(1-2*nu));
  
  p = inputParser;
  if isOctave()
    %p=p.addParamValue('LN', [] , @isnumeric );
    p=p.addParamValue('nbruns', 5 , @isnumeric );
    p=p.addParamValue('isplot',false, @islogical);
    p=p.addParamValue('lam',lam);
    p=p.addParamValue('mu',mu);
    p=p.parse(varargin{:});
  else
    %p.addParamValue('LN', [] , @isnumeric );
    p.addParamValue('nbruns', 5 , @isnumeric );
    p.addParamValue('isplot',false, @islogical);
    p.addParamValue('lam',lam);
    p.addParamValue('mu',mu);
    p.parse(varargin{:});
end
  %LN=p.Results.LN;
  nbruns=p.Results.nbruns;
  isplot=p.Results.isplot;
  lam=p.Results.lam;
  mu=p.Results.mu;
  nLN=length(LN);
  nV=length(Versions);
  Tcpu=zeros(nLN,nV);
  Lndof=zeros(nLN,1);
  for i=1:nLN
    N=LN(i);
    fprintf(' -> (%d/%d) Creation of the mesh - HyperCube(%d,%d) ...\n',i,nLN,d,N);
    Th=HyperCube(d,N);
    fprintf('    %dd-mesh : nq=%d, nme=%d\n',Th.d,Th.nq,Th.nme);
    Lam=setFdata(lam,Th);Mu=setFdata(mu,Th);
    for v=1:nV
      fprintf('    -> (%d/%d) Run AssemblyStiffElasP1%s (%dd)\n',v,nV,Versions{v},d);
      Assembly=eval(sprintf('@(Th) AssemblyStiffElasP1%s(Th,Lam,Mu);',Versions{v}));
      for l=1:nbruns
	start=tic();
	M=Assembly(Th);
	t(l)=toc(start);
	fprintf('         run (%d/%d) : %.4f(s)\n',l,nbruns,t(l))
      end
      if v==1,M1=M;end
      Tcpu(i,v)=mean(t);
      fprintf('       AssemblyStiffElasP1%6s (%dd) : %d-by-%d matrix in %.4f(s)\n',Versions{v},d,size(M,1),size(M,2),Tcpu(i,v));
      if v>1, 
        fprintf('       Error with %6s matrix   : %.5e\n',Versions{1},full(max(max(abs(M-M1)))));
        fprintf('       %6s speedup             : x%.3f\n',Versions{1},Tcpu(i,1)/Tcpu(i,v));
      end
    end
    Lndof(i)=size(M1,1);
  end   
  if isplot
    PlotBench(Versions,Lndof,Tcpu,sprintf('Benchmark of AssemblyStiffElasP1 functions in %dD',d))
  end
end

