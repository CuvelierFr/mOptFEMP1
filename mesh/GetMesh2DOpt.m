function Mesh=GetMesh2DOpt(cFileName)
  [fid,message]=fopen(cFileName,'r');
  if ( fid == -1 )
    error([message,' : ',cFileName]);
  end
  if isOctave()
    [n]=fscanf(fid,'%d %d %d',3);

    R=fscanf(fid,'%f %f %d',[3,n(1)]);
    q=R([1 2],:);
    ql=R(3,:);
    R=fscanf(fid,'%d %d %d %d',[4,n(2)]);

    me=R([1:3],:);
    mel=R(4,:);
    R=fscanf(fid,'%d %d %d',[3,n(3)]);
    
    be=R([1 2],:);
    bel=R(3,:);
  else % Matlab
    n=textscan(fid,'%d %d %d',1); % n(1) -> number of vertices
			% n(2) -> number of triangles
			% n(3) -> number of boundary edges
    
    R=textscan(fid,'%f %f %d',n{1});
    q=[R{1},R{2}]';
    ql=R{3}';
    R=textscan(fid,'%d %d %d %d',n{2});
    me=[R{1},R{2},R{3}]';
    mel=R{4}';

    R=textscan(fid,'%d %d %d',n{3});
    be=[R{1},R{2}]';
    bel=R{3}';
  end
  fclose(fid);

  Mesh=struct('d',2,'q',q,'me',double(me),'ql',ql,'mel',double(mel), ...
	      'be',double(be),'bel',double(bel), ...
	      'nq',size(q,2), ...
	      'nme',size(me,2), ...
	      'nbe',size(be,2), ...
	      'vols',ComputeVolVec(2,q,me)', ...
	      'h',GetMaxLengthEdges(q,me));
end