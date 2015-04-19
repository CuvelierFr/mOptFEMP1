function P=mOptFEMP1path()
  P = mfilename('fullpath');
  I=strfind(P,filesep);
  P=P(1:I(end-1)-1);
end