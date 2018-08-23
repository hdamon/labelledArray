function varargout = subsref(obj,s)
% subsref method for labelledArray objects
%
%
% What should subreferencing an arrayDim look like?'
%
% 

switch s(1).type
  case '.'        
    varargout = {builtin('subsref',obj,s)};        
  case '()'
    newIdx = obj.findDimensions(s(1).subs);
    s(1).subs = {newIdx};
    varargout = {builtin('subsref',obj,s)};        
  case '{}'
    varargout = {builtin('subsref',obj,s)};        
  otherwise
    error('Not a valid indexing expression')
end
