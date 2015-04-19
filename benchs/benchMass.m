function [Lndof,Tcpu]=benchMass(d,Versions,LN,varargin)
%  benchMass : Benchmark for AssemblyMassP1<version> functions using HyperCube(d,N) to generate meshes.
%
%  USAGE:
%    [Lndof,Tcpu]=benchMass(d,Versions,LN)
%    [Lndof,Tcpu]=benchMass(d,Versions,LN,'nbruns',3)
%    [Lndof,Tcpu]=benchMass(d,Versions,LN,'isplot',true)
%
%  INPUTS:
%    d        : space dimension
%    Versions : List of Assembly version (cell array),  
%               for example {'OptVS','OptV','OptV2','OptV1'}
%    LN       : List of N values (1D array of integer values)
%    <Optional>
%    'nbruns' : Number of runs
%    'isplot' : if true, plot results 
%
%  OUTPUTS:
%    Lndof    : List of size matrices (#LN-by-1 array)
%    Tcpu     : #LN-by-#Versions array, cputimes for each Assembly versions and 
%               each HyperCube meshes. Tcpu(i,v), is the mean cputimes for
%               Assembly with version Versions{v} and HyperCube(d,LN[i]) mesh.
%
%  SAMPLES:
%    [Lndof,Tcpu]=benchMass(2,{'OptVS','OptV','OptV2','OptV1'},[50:50:250]);
%    [Lndof,Tcpu]=benchMass(3,{'OptVS','OptV','OptV2','OptV1'},[5:5:25]);
%    [Lndof,Tcpu]=benchMass(2,{'OptVS','OptV','OptV2'},[100:100:1000],'nbruns',1,'isplot',true);
%    [Lndof,Tcpu]=benchMass(4,{'OptVS','OptV','OptV2'},[5:5:20],'nbruns',1,'isplot',true);
%    [Lndof,Tcpu]=benchMass(5,{'OptVS','OptV','OptV2'},[2:2:8],'nbruns',1,'isplot',true);
%
  p = inputParser;
  if isOctave()
    p=p.addParamValue('nbruns', 5 , @isnumeric );
    p=p.addParamValue('isplot',false, @islogical);
    p=p.parse(varargin{:});
  else
    p.addParamValue('nbruns', 5 , @isnumeric );
    p.addParamValue('isplot',false, @islogical);
    p.parse(varargin{:});
end
  nbruns=p.Results.nbruns;
  isplot=p.Results.isplot;
  nLN=length(LN);
  nV=length(Versions);
  Tcpu=zeros(nLN,nV);
  Lndof=zeros(nLN,1);
  for i=1:nLN
    N=LN(i);
    fprintf(' -> (%d/%d) Creation of the mesh - HyperCube(%d,%d) ...\n',i,nLN,d,N);
    Th=HyperCube(d,N);
    fprintf('    %dd-mesh : nq=%d, nme=%d\n',Th.d,Th.nq,Th.nme);
    for v=1:nV
      fprintf('    -> (%d/%d) Run AssemblyMassP1%s (%dd)\n',v,nV,Versions{v},d);
      Assembly=eval(sprintf('@(Th) AssemblyMassP1%s(Th);',Versions{v}));
      for l=1:nbruns
	start=tic();
	M=Assembly(Th);
	t(l)=toc(start);
	fprintf('         run (%d/%d) : %.4f(s)\n',l,nbruns,t(l))
      end
      if v==1,M1=M;end
      Tcpu(i,v)=mean(t);
      fprintf('       AssemblyMassP1%6s (%dd) : %d-by-%d matrix in %.4f(s)\n',Versions{v},d,size(M,1),size(M,2),Tcpu(i,v));
      if v>1, 
        fprintf('       Error with %6s matrix   : %.5e\n',Versions{1},full(max(max(abs(M-M1)))));
        fprintf('       %6s speedup             : x%.3f\n',Versions{1},Tcpu(i,1)/Tcpu(i,v));
      end
    end
    Lndof(i)=size(M1,1);
  end   
  if isplot
    PlotBench(Versions,Lndof,Tcpu,sprintf('Benchmark of AssemblyMassP1 functions in %dD',d))
  end
end

