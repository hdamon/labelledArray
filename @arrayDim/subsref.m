function varargout = subsref(obj,s)
% subsref method for labelledArray objects
%

switch s(1).type
  case '.'    
    for i = 1:nargout
     varargout(i) = {builtin('subsref',obj(i),s)};    
    end
  case '()'
    newIdx = obj.findDimensions(s(1).subs);
    s(1).subs = {newIdx};
    varargout = {builtin('subsref',obj,s)};        
  case '{}'
    varargout = {builtin('subsref',obj,s)};        
  otherwise
    error('Not a valid indexing expression')
end
