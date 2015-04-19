function InitmOptFEMP1(varargin)
  log=ver;
  bool=strcmp(log(1).Name,'Octave');
  if bool
    pkg load general;
  end
  p = inputParser;
  if bool
    p=p.addParamValue('pathm', pwd, @ischar );
    p=p.parse(varargin{:});
  else
    p.addParamValue('pathm', pwd, @ischar );
    p.parse(varargin{:});
  end

  addpath([p.Results.pathm,filesep,'common']);
  addpath([p.Results.pathm,filesep,'FEM']);
  addpath([p.Results.pathm,filesep,'mesh']);
  addpath([p.Results.pathm,filesep,'benchs']);
%  addpath([p.Results.pathm,filesep,'samples']);  
  addpath([p.Results.pathm,filesep,'graphic']);
  if isOctave()
    more off
    pkg load msh
    warning('off','Octave:fopen-file-in-path') % disable warning for mesh files in path
  end
end