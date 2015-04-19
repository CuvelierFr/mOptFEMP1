function Mesh=GetMeditMesh(cFileName,order,varargin)
  p = inputParser; 
  if isOctave()
    p=p.addParamValue('info', false, @islogical );
    p=p.addParamValue('verbose', false, @islogical );
    p=p.parse(varargin{:});
  else
    p.addParamValue('info', false, @islogical );
    p.addParamValue('verbose', false, @islogical );
    p.parse(varargin{:});
  end
  info=p.Results.info;verbose=p.Results.verbose;
  if verbose, fprintf('reading mesh file : %s\n',cFileName);end
  [nq,q,ql]=meditVertices(cFileName,info);
  [nme,me,mel]=meditElements(cFileName,order,'Tetrahedra',3,info);
  if nme==0 % 2d mesh
    [nme,me,mel]=meditElements(cFileName,order,'Triangles',2,info);
    if nme==0, error('Unable to read 2d mesh file %s\n Field "Triangles" not found !!!',cFileName); end
    [nbe,be,bel]=meditElements(cFileName,order,'Edges',1,info);
    if nbe==0, error('Unable to read 2d mesh file %s\n Field "Edges" not found !!!',cFileName); end
    d=2;
    N=factorial(d+order)/(factorial(d)*factorial(order));
    indve=[1,order+1,N]'; %index of triangles vertices
  else % 3d mesh
    [nbe,be,bel]=meditElements(cFileName,order,'Triangles',2,info);
    if nbe==0, error('Unable to read 3d mesh file %s\n Field "Triangles" not found !!!',cFileName); end
    d=3;
    N=factorial(d+order)/(factorial(d)*factorial(order));
    indve=[1,order+1,(order+1)*(order+2)/2,N]'; %index of tetrahedra vertices
  end
  if info
    Mesh=struct('order',order, ... % mesh order
              'd',d,'nq',nq,'nme',nme,'nbe',nbe,'indve',indve);
  else
    %Jb=changeOrderD(q(:,me(:,1)),order,d);
    %me=me(Jb,:);
    Mesh=struct('order',order, ... % mesh order
              'd',d,'q',q(1:d,:),'me',double(me),'ql',ql,'mel',double(mel), ...
              'be',double(be),'bel',double(bel), ...
              'nq',nq, ...
              'nme',nme, ...
              'nbe',nbe, ...
              'indve',indve, ...
              'vols',ComputeVolVec(d,q,me(indve,:))', ...
              'h',GetMaxLengthEdges(q,me(indve,:)));
  end
end

function Jb=changeOrderD(Q,order,d)
% Change l'ordre des noeuds locaux de 'medit' -> 'FC' 
%
% Au format 'medit' provenant de gmsh les 4 1er pts sont les
% sommets du tetra
% On peut faire mieux avec les distances pour trouver les 4 points
% 

  Qk=getPkPoints(Q(:,1:d+1),order); % A ameliorer

  [QQ1,Ik]=sortrows(Qk');
  [QQ,I]=sortrows(Q');
  %[L,IL]=sort(I);
  %Ja=Ik(IL);
  [L,IL]=sort(Ik);
  Jb=I(IL);
  %  max(max(abs(Qk(:,Ja)-Q)))
  %  max(max(abs(Q(:,Jb)-Qk)))
end

function [nq,q,ql]=meditVertices(cFileName,info)
  [fid,message]=fopen(cFileName,'r');
  if ( fid == -1 ), error([message,' : ',cFileName]);end
  tline='';
  while ~strcmp(strtrim(tline),'Vertices')
    tline = fgetl(fid);
    if (tline==-1), error('Unable to read mesh file %s\n Field "Vertices" not found !!!',cFileName); end
  end
  nq=fscanf(fid,'%d',1);
  if info, q=[];ql=[];fclose(fid);return;end
  if isOctave()
    R=fscanf(fid,'%f %f %f %d',[4,nq]);
    q=R([1 2 3],:);
    ql=R(4,:);
  else
    R=textscan(fid,'%f %f %f %d',nq);
    q=[R{1},R{2},R{3}]';
    ql=R{4}';
  end
  fclose(fid);
end

function [nme,me,mel]=meditElements(cFileName,order,ElementStr,d,info)
  [fid,message]=fopen(cFileName,'r');
  nme=0;me=[];mel=[];
  tline='';
  N=factorial(d+order)/(factorial(d)*factorial(order)); % Nb of Pk vertices on 2-simplex
  while ~strcmp(strtrim(tline),ElementStr)
    tline = fgetl(fid);
    if (tline==-1), return; end
  end
  nme=fscanf(fid,'%d',1);
  if info, me=[];mel=[];fclose(fid);return;end
  format=repmat('%d ',1,N+1);
  if isOctave()
    R=fscanf(fid,format,[N+1,nme]);
    me=R([1:N],:);
    mel=R(N+1,:);
  else
    R=cell2mat(textscan(fid,format,nme))';
    me=R(1:N,:);
    mel=R(N+1,:);
  end
  fclose(fid);
end
