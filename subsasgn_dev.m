function obj = subsasgn(obj,s,varargin)
   switch s(1).type
      case '.'
         if length(s) == 1
            % Implement obj.PropertyName = varargin{:};
            ...
         elseif length(s) == 2 && strcmp(s(2).type,'()')
            % Implement obj.PropertyName(indices) = varargin{:};
            ...
         else
            % Call built-in for any other case
            obj = builtin('subsasgn',obj,s,varargin);
         end
      case '()'
         if length(s) == 1
            % Implement obj(indices) = varargin{:};
         elseif length(s) == 2 && strcmp(s(2).type,'.')
            % Implement obj(indices).PropertyName = varargin{:};
            ...
         elseif length(s) == 3 && strcmp(s(2).type,'.') && strcmp(s(3).type,'()')
            % Implement obj(indices).PropertyName(indices) = varargin{:};
            ...
         else
            % Use built-in for any other expression
            obj = builtin('subsasgn',obj,s,varargin);
         end       
      case '{}'
         if length(s) == 1
            % Implement obj{indices} = varargin{:}
            ...
         elseif length(s) == 2 && strcmp(s(2).type,'.')
            % Implement obj{indices}.PropertyName = varargin{:}
            ...
            % Use built-in for any other expression
            obj = builtin('subsasgn',obj,s,varargin);
         end
      otherwise
         error('Not a valid indexing expression')
   end
end