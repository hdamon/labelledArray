function varargout = subsref(obj,s)
% subsref method for labelledArray objects
%
%
% It's likely that this subsref() function is a bit more complicated than
% it needs to be.
%
% The main goal here is to enable referencing by dimension label and value,
% and to provide similar capability when referencing into the .array
% property.
%
% Part of the crlBundle
% Copyright 2009-2018
% Damon Hyde
%

%for i = 1:numel(s)

switch s(1).type
  case '.'       
    if length(s)>1
      
      hiddenMethods = {'findDimensions','isMutuallyConsistent','getConsistentDimensions'};
      if ischar(s(1).subs)&&(ismember(s(1).subs,methods(obj))||ismember(s(1).subs,hiddenMethods))
        % If we're calling a method, us the builtin
        [varargout{1:nargout}] = builtin('subsref',obj,s);
      else
        tmp = subsref(obj,s(1));
        
        objFields = {'array' 'dimNames' 'dimLabels' 'dimValues' 'dimUnits'};
        if ismember(s(1).subs,objFields)&&isequal(s(2).type,'()')         
          % Parenthetical referencing into a labelledArray
          s(2).subs = obj.getNumericIndex(s(2).subs{:});
          [varargout{1:nargout}] = subsref(tmp,s(2:end));
        else
          [varargout{1:nargout}] = builtin('subsref',tmp,s(2:end));%{subsref(tmp,s(2:end))};
        end;
      end;
    else
      %disp(['Using builtin to access: ' s(1).subs]);
      
      %varargout = {builtin('subsref',obj,s)};
      [varargout{1:nargout}] = builtin('subsref',obj,s);
    end
    
  case '()'
    if numel(obj) == 1
      %% Implement obj(indices) for numel(obj)==
      
      if (numel(s(1).subs)==1)&&isnumeric(s(1).subs{1})...
            &&(numel(s(1).subs{1})==1)&&(s(1).subs{1}==1)
        if numel(s)>1
        % varargout = {subsref(obj,s(2:end))};
         [varargout{1:nargout}] = subsref(obj,s(2:end));
        else
          %varargout = {builtin('subsref',obj,s)};
          [varargout{1:nargout}] = builtin('subsref',obj,s);
        end;
        return;
      end
      
      if numel(s)==1
        % Internal object indices can only be accessed individually        
        %disp(['Accessing single object: Single Reference']);
        %dimIdx = obj.getNumericIndex(s.subs{:});        
        %varargout = {obj.subcopy(s.subs{:})};        
        [varargout{1:nargout}] = obj.subcopy(s.subs{:});
      else
        % Get the right object, then reference into it.
        %disp(['Accessing single object: Multiple Reference']);
        tmp = subsref(obj,s(1));
        %varargout = {subsref(tmp,s(2:end))};                
        [varargout{1:nargout}] = subsref(tmp,s(2:end));
      end;
    else
        % obj(indices) for numel(obj)>1
                
        if numel(s)==1
          % Just use the builtin for arrays of objects
          %varargout = {builtin('subsref',obj,s)};
          [varargout{1:nargout}] = builtin('subsref',obj,s);
        else          
          % Otherwise, reference into the array using the builtin, then get
          % the appropriate values.
          tmp = builtin('subsref',obj,s(1));
          %varargout = {subsref(tmp,s(2:end))};
          [varargout{1:nargout}] = subsref(tmp,s(2:end));
        end;
    end
    
  case '{}'
    s2.type = '()';
    s2.subs = s(1).subs;
    tmp = obj.subsref(s2);
    
    if length(s)==1
      %varargout = {tmp};
      [varargout{1:nargout}] = tmp;
    else
      %varargout = {tmp.subsref(s(2:end))};
      [varargout{1:nargout}] = tmp.subsref(s(2:end));
    end;

  otherwise
    error('Not a valid indexing expression')
end
