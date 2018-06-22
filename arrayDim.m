classdef arrayDim
% Object class for labelledArray dimension information
%
%

  properties
    dimName
    dimSize
    dimLabels
    dimUnits
    dimValues
  end

  properties (Access=private)
    dimName_
    dimSize_
    dimLabels_
    dimUnits_
    dimValues_
  end
  
  
  methods


    function obj = arrayDim(varargin)
    % Class constructor
      if nargin>0
        if isa(varargin{1},'arrayDim')
          obj = varargin{1};
          return;
        end

        p = inputParser;
        p.addParameter('dimName',[],@ischar);
        p.addParameter('dimSize',[],@(x) isnumeric(x)&&isscalar(x));
        p.addParameter('dimLabels',[],@(x) ischar(x)||iscellstr(x));
        p.addParameter('dimUnits',[],@(x) ischar(x)||iscellstr(x));
        p.addParameter('dimValues',[],@(x) isnumeric(x)&&isvector(x));
        p.parse(varargin{:});

        obj.dimName   = p.Results.dimName;
        obj.dimSize   = p.Results.dimSize;
        obj.dimLabels = p.Results.dimLabels;
        obj.dimUnits  = p.Results.dimUnits;
        obj.dimValues = p.Results.dimValues;

      end
    end

    function set.dimName(obj,val)
    
    function isEqual = isequal(obj,a)
      if ~(isequal(obj.name,a.name))
       isEqual = false;
       warning('Dimension name mismatch');
       end 
       if ~(isequal(obj.labels,a.labels))
         isEqual = false;
         warning('Dimension label mismatch');
       end

       if ~(isequal(obj.units,a.units))
         isEqual = false;
         warning('Dimension unit mismatch');
       end

       if ~(isequal(obj.values,a.values))
         isEqual = false;
         warning('Dimension value mismatch');
       end
     
    end

    function isEmpty = isempty(obj)
      isEmpty = true;
      if numel(obj)==1
        % Empty only if all components are empty
        isEmpty = isEmpty&&isempty(obj.name);
        isEmpty = isEmpty&&isempty(obj.labels);
        isEmpty = isEmpty&&isempty(obj.units);
        isEmpty = isEmpty&&isempty(obj.values);
      end
    end

    end

    %% Private Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    
    methods (Access=private)
      function setProperty(obj,propName,val,castCharCell,isSingleOk,typeCheckFunc,typeCheckErr)
      % Protected utility function for internal setting of properties
      %
      % setProperty(obj,propName,val,castCharCell,isSingleOk,typeCheckFunc,typeCheckErr)
      %
      % NOTE: BECAUSE SUPERCLASS METHODS CANNOT ACCESS SUBCLASS PROTECTED
      % PROPERTIES, THIS FUNCTION MUST BE DUPLICATED IN ANY SUBCLASS!!!!
      %
      % Inputs
      % ------
      %            obj : labelledArray object 
      %       propName : String with propertyname to set
      %            val : Value to set to the property
      %   castCharCell : Flag to turn on typecasting of character string to
      %                   cellstrs.
      %     isSingleOK : Flag to permit singleton values
      %  typeCheckFunc : Type checking function handle
      %                   (Note: Empty is always valid)
      %  typeCheckErr  : Error to throw if type check fails
      %
      %
      
      
      if isempty(val)        
        % Empty cell for each dimension
        obj.([propName '_']) = cell(obj.ndims,1);
        return;
      end
      
      val = obj.validateCellInput(val);
      
      % Check Input
      for i = 1:size(val,1)
        currDim = val{i,1};
        currVal = val{i,2};
                     
        assert(isnumeric(currDim)&&isscalar(currDim),...
          'Dimension parameter must be a numeric scalar');
        
        % Convert char to cellstr if casting requested
        if castCharCell
          if ischar(currVal)
            currVal = {currVal};
            val{i,2} = strtrim(currVal);
          end
        end
       
        % Check Type
        assert(isempty(currVal)||typeCheckFunc(currVal),typeCheckErr);
        
        % Check Size
        assert(isempty(currVal)||...
               (numel(currVal)==size(obj,currDim))||...
               (isSingleOk&&(numel(currVal)==1)),...
          ['Must provide ' propName 'for the full dimension']);
      end
      
      % Assign dimLabels if All Checks Passed
      for i = 1:size(val,1)        
        if (numel(val{i,2})==1)&&castCharCell
          tmp = val{i,2};
          obj.([propName '_']){val{i,1}} = tmp{1};                
        else
          obj.([propName '_']){val{i,1}} = val{i,2};                
        end;
      end
      
    end   
    end
    
    
end
