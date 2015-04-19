function bool=isfhandle(f)
  bool=strcmp(class(f),'function_handle');
end