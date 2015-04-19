function bool=isOctave()
% function bool=isOctave()
%   To determine whether Octave is used or not
%   
% Return values:
%  bool: if true Octave is used else Matlab is
%
log=ver;
bool=strcmp(log(1).Name,'Octave');
