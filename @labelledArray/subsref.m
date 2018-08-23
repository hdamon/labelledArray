function varargout = subsref(obj,s)
% subsref method for labelledArray objects
%

switch s(1).type
  case '.'       
    if length(s)>1
      if ischar(s(1).subs)&&ismember(s(1).subs,methods(obj))
        % If we're calling a method, pass the arguments
        varargout = {builtin('subsref',obj,s)};
      else
        tmp = subsref(obj,s(1));
        
        objFields = {'array' 'dimNames' 'dimLabels' 'dimValues' 'dimUnits'};
        if ismember(s(1).subs,objFields)&&isequal(s(2).type,'()')         
          % Parenthetical referencing into a labelledArray
          s(2).subs = obj.getNumericIndex(s(2).subs{:});
          varargout = {subsref(tmp,s(2:end))};
        else
          varargout = {subsref(tmp,s(2:end))};
        end;

      end;
    else
      %disp(['Using builtin to access: ' s(1).subs]);
      varargout = {builtin('subsref',obj,s)};
    end
    
  case '()'
    if numel(obj) == 1
      %% Implement obj(indices)
      
      if (numel(s(1).subs)==1)&&(numel(s(1).subs{1})==1)&&(s(1).subs{1}==1)
        varargout = {builtin('subsref',obj,s)};
        return;
      end
      
      if numel(s)==1
        % Internal object indices can only be accessed individually        
        %disp(['Accessing single object: Single Reference']);
        %dimIdx = obj.getNumericIndex(s.subs{:});        
        varargout = {obj.subcopy(s.subs{:})};
      else
        % Get the right object, then reference into it.
        %disp(['Accessing single object: Multiple Reference']);
        tmp = subsref(obj,s(1));
        varargout = {subsref(tmp,s(2:end))};                
      end;
    else
        % Just use the builtin for arrays of objects
        
        if numel(s)==1
          varargout = {builtin('subsref',obj,s)};
        else          
          tmp = builtin('subsref',obj,s(1));
          varargout = {subsref(tmp,s(2:end))};
        end;
    end
    
  case '{}'
    s2.type = '()';
    s2.subs = s(1).subs;
    tmp = obj.subsref(s2);
    
    if length(s)==1
      varargout = {tmp};
    else
      varargout = {tmp.subsref(s(2:end))};
    end;

  otherwise
    error('Not a valid indexing expression')
end
