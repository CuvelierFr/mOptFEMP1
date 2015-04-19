function [Lndof,Tcpu]=benchMassW(d,Versions,LN,varargin)
%  benchMass : Benchmark for AssemblyMassWP1<version> functions using HyperCube(d,N) to generate meshes.
%
%  USAGE:
%    [Lndof,Tcpu]=benchMassW(d,Versions,LN)
%    [Lndof,Tcpu]=benchMassW(d,Versions,LN,'nbruns',3)
%    [Lndof,Tcpu]=benchMassW(d,Versions,LN,'isplot',true)
%    [Lndof,Tcpu]=benchMassW(d,Versions,LN,'w',func_handle)
%
%  INPUTS:
%    d        : space dimension
%    Versions : List of Assembly version (cell array),  
%               for example {'OptVS','OptV','OptV2','OptV1'}
%    LN       : List of N values (1D array of integer values)
%    <Optional>
%    'nbruns' : Number of runs
%    'isplot' : if true, plot results 
%    'w'      : weight function or constant value
%
%  OUTPUTS:
%    Lndof    : List of size matrices (#LN-by-1 array)
%    Tcpu     : #LN-by-#Versions array, cputimes for each Assembly versions and 
%               each HyperCube meshes. Tcpu(i,v), is the mean cputimes for
%               Assembly with version Versions{v} and HyperCube(d,LN[i]) mesh.
%
%  SAMPLES:
%    [Lndof,Tcpu]=benchMassW(2,{'OptVS','OptV','OptV2','OptV1'},[50:50:250]);
%    [Lndof,Tcpu]=benchMassW(3,{'OptVS','OptV','OptV2','OptV1'},[5:5:25],'nbruns',1);
%    [Lndof,Tcpu]=benchMassW(2,{'OptVS','OptV','OptV2'},[100:100:1000],'nbruns',1,'isplot',true,'w',pi);
%    [Lndof,Tcpu]=benchMassW(4,{'OptVS','OptV','OptV2'},[5:5:20],'nbruns',1,'isplot',true);
%    [Lndof,Tcpu]=benchMassW(5,{'OptVS','OptV','OptV2'},[2:2:8],'nbruns',1,'isplot',true);
%
  p = inputParser;
  if isOctave()
    %p=p.addParamValue('LN', [] , @isnumeric );
    p=p.addParamValue('nbruns', 5 , @isnumeric );
    p=p.addParamValue('isplot',false, @islogical);
    p=p.addParamValue('w',genericFunc(d,['cos(1',sprintf('+x%d',1:d),')']));
    p=p.parse(varargin{:});
  else
    %p.addParamValue('LN', [] , @isnumeric );
    p.addParamValue('nbruns', 5 , @isnumeric );
    p.addParamValue('isplot',false, @islogical);
    p.addParamValue('w',genericFunc(d,['cos(1',sprintf('+x%d',1:d),')']));
    p.parse(varargin{:});
end
  %LN=p.Results.LN;
  nbruns=p.Results.nbruns;
  isplot=p.Results.isplot;
  w=p.Results.w;
  nLN=length(LN);
  nV=length(Versions);
  Tcpu=zeros(nLN,nV);
  Lndof=zeros(nLN,1);
  for i=1:nLN
    N=LN(i);
    fprintf(' -> (%d/%d) Creation of the mesh - HyperCube(%d,%d) ...\n',i,nLN,d,N);
    Th=HyperCube(d,N);
    fprintf('    %dd-mesh : nq=%d, nme=%d\n',Th.d,Th.nq,Th.nme);
    W=setFdata(w,Th);
    for v=1:nV
      fprintf('    -> (%d/%d) Run AssemblyMassWP1%s (%dd)\n',v,nV,Versions{v},d);
      Assembly=eval(sprintf('@(Th) AssemblyMassWP1%s(Th,W);',Versions{v}));
      for l=1:nbruns
	start=tic();
	M=Assembly(Th);
	t(l)=toc(start);
	fprintf('         run (%d/%d) : %.4f(s)\n',l,nbruns,t(l))
      end
      if v==1,M1=M;end
      Tcpu(i,v)=mean(t);
      fprintf('       AssemblyMassWP1%6s (%dd) : %d-by-%d matrix in %.4f(s)\n',Versions{v},d,size(M,1),size(M,2),Tcpu(i,v));
      if v>1, 
        fprintf('       Error with %6s matrix   : %.5e\n',Versions{1},full(max(max(abs(M-M1)))));
        fprintf('       %6s speedup             : x%.3f\n',Versions{1},Tcpu(i,1)/Tcpu(i,v));
      end
    end
    Lndof(i)=size(M1,1);
  end   
  if isplot
    PlotBench(Versions,Lndof,Tcpu,sprintf('Benchmark of AssemblyMassWP1 functions in %dD',d))
  end
end

