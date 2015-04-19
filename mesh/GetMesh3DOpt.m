function Mesh=GetMesh3DOpt(cFileName,varargin)
% 'format' : 'freefem' (default) , 'gmsh' or 'medit'
  p = inputParser; 
  if isOctave()
    p=p.addParamValue('format', 'freefem', @ischar );
    p=p.addParamValue('order', 1, @isscalar );
    p=p.parse(varargin{:});
  else
    p.addParamValue('format', 'freefem', @ischar );
    p.addParamValue('order', 1, @isscalar );
    p.parse(varargin{:});
  end
  order=p.Results.order;
  Format=p.Results.format;
  % FC -> a ameliorer avec une liste des formats dispo et 
  %       utilisation + generique GetMesh%s...
  if (strcmp(Format,'freefem')) 
    assert(order==1);
    Mesh=GetFreefemMesh(cFileName);
  end

  if (strcmp(Format,'gmsh'))
    assert(order==1);
    Mesh=GetGmshMesh(cFileName);
  end
  if (strcmp(Format,'medit'))
    Mesh=GetMeditMesh(cFileName,order);
  end
end

% Read FreeFEM++ 3D meshes
function Mesh=GetFreefemMesh(cFileName)
[fid,message]=fopen(cFileName,'r');
  if ( fid == -1 )
      error([message,' : ',cFileName]);
  end
  for i=1:5
      tline = fgetl(fid);
  end
  if isOctave()
    nq=fscanf(fid,'%d',1);
    
    R=fscanf(fid,'%f %f %f %d',[4,nq]);
    q=R([1 2 3],:);
    ql=R(4,:);
    for i=1:3
        tline = fgetl(fid);
    end
    nme=fscanf(fid,'%d',1);
    R=fscanf(fid,'%d %d %d %d %d',[5,nme]);
    
    me=R([1:4],:);
    mel=R(5,:);
    for i=1:3
        tline = fgetl(fid);
    end
    nbf=fscanf(fid,'%d',1);
    R=fscanf(fid,'%d %d %d %d',[4,nbf]);
    
    bf=R([1 2 3],:);
    bfl=R(4,:);
  else % Matlab
    nq=fscanf(fid,'%d',1);
    
    R=textscan(fid,'%f %f %f %d',nq);
    q=[R{1},R{2},R{3}]';
    ql=R{4}';
    
    for i=1:3
        tline = fgetl(fid);
    end
    nme=fscanf(fid,'%ld',1);
    
    R=textscan(fid,'%d %d %d %d %d',nme);
    me=[R{1},R{2},R{3},R{4}]';
    mel=R{5}';
    for i=1:3
        tline = fgetl(fid);
    end
    nbf=fscanf(fid,'%d',1);
    
    R=textscan(fid,'%d %d %d %d',nbf);
    bf=[R{1},R{2},R{3}]';
    bfl=R{4}';
  end
  fclose(fid);
  Mesh=struct('d',3,'q',q,'me',double(me),'ql',ql,'mel',double(mel),'be',double(bf), ...
              'bel',double(bfl),'nq',nq,'nme',nme,'nbe',nbf, ...
              'vols',ComputeVolVec(3,q,me),'h',GetMaxLengthEdges(q,me));
end

% Read gmsh meshes
function Th=GetGmshMesh(cFileName)
  msh=load_gmsh(cFileName)
  Th.d=3;
  Th.nq=msh.nbNod;
  Th.q=msh.POS';
  Th.me=msh.TETS(1:msh.nbTets,1:4)';
  Th.nme=msh.nbTets
  Th.mel=msh.TETS(1:msh.nbTets,5)';
  Th.vols=ComputeVolVec(Th.d,Th.q,Th.me);
end

% Read medit meshes
function Mesh=GetMeditMeshOld(cFileName,order)
  [fid,message]=fopen(cFileName,'r');
   if ( fid == -1 )
     error([message,' : ',cFileName]);
   end
  if isOctave()
    % Read Vertices
    tline='';
    while ~strcmp(strtrim(tline),'Vertices')
        tline = fgetl(fid);
        if (tline ==-1)
          error('Error : Vertices not found');
        end
    end
    nq=fscanf(fid,'%d',1);
    R=fscanf(fid,'%f %f %f %d',[4,nq]);
    q=R([1 2 3],:);
    ql=R(4,:);
    
    % Read Triangles
    d=2;
    N=factorial(d+order)/(factorial(d)*factorial(order)); % Nb of Pk vertices on 2-simplex
    fclose(fid);
    [fid,message]=fopen(cFileName,'r');
    while ~strcmp(strtrim(tline),'Triangles')
        tline = fgetl(fid);
        if (tline ==-1)
          error('Error : Triangles not found');
        end
    end
    %for i=1:N+1, 
    nbe=fscanf(fid,'%d',1);
    R=fscanf(fid,'%d %d %d %d',[4,nbe]);   
    be=R([1 2 3],:);
    bel=R(4,:);
    
    % Read Tetrahedra
    fclose(fid);
    [fid,message]=fopen(cFileName,'r');
    tline='';
    while ~strcmp(strtrim(tline),'Tetrahedra')
        tline = fgetl(fid);
        if (tline ==-1)
          error('Error : Tetrahedra not found');
        end
    end
    nme=fscanf(fid,'%d',1);
    R=fscanf(fid,'%d %d %d %d %d',[5,nme]);  
    me=R([1:4],:);
    mel=R(5,:);
  else % Matlab
    tline='';
    while ~strcmp(strtrim(tline),'Vertices')
        tline = fgetl(fid);
    end
    nq=fscanf(fid,'%d',1);
    R=textscan(fid,'%f %f %f %d',nq);
    q=[R{1},R{2},R{3}]';
    ql=R{4}';
    fclose(fid);
    [fid,message]=fopen(cFileName,'r');
    tline='';
    % Triangles
    d=2;
    N=factorial(d+order)/(factorial(d)*factorial(order)); % Nb of Pk vertices on 2-simplex
    while ~strcmp(strtrim(tline),'Triangles')
        tline = fgetl(fid);
    end
    nbe=fscanf(fid,'%d',1);
    format=repmat('%d ',1,N+1);
    R=cell2mat(textscan(fid,format,nbe))';
    be=R(1:N,:);
    bel=R(N+1,:);
    fclose(fid);
    [fid,message]=fopen(cFileName,'r');
    % Tetrahedra
    d=3;
    N=factorial(d+order)/(factorial(d)*factorial(order)); % Nb of Pk vertices on 3-simplex
    tline='';
    while (~strcmp(strtrim(tline),'Tetrahedra'))
        tline = fgetl(fid);
    end
    % Apres Tetrahedra
    nme=fscanf(fid,'%ld',1);
    format=repmat('%d ',1,N+1);
    R=cell2mat(textscan(fid,format,nme))';
    %R=textscan(fid,'%d %d %d %d %d',nme);
    me=R(1:N,:);
    mel=R(N+1,:);
  end
  fclose(fid);
    
  Mesh=struct('d',3,'q',q,'me',double(me),'ql',ql,'mel',double(mel), ...
              'be',double(be),'bel',double(bel), ...
              'nq',nq, ...
              'nme',nme, ...
              'nbe',nbe, ...
              'vols',ComputeVolVec(3,q,me), ...
              'h',GetMaxLengthEdges(q,me));
end