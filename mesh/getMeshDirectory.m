function meshdir=getMeshDirectory(d)
  name = getComputerName();
  switch lower(name)
    case {'fx-8350','hercule'}
      meshdir=sprintf('~/meshes/%dd',d);
    case {'gpuschwarz','gpucreos1'}
      meshdir=sprintf('/home/CJS/meshes/%dD',d);
    case {'box'}
      meshdir=sprintf('/home/scarella/meshes/%dd',d);
    otherwise
      meshdir=sprintf('/home/CJS/meshes/%dD',d);
  end