classdef arrayDim
  % Object class for labelledArray dimension information
  %
  %
  
  properties (Dependent)
    dimName
    dimSize
    fixedSize
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
        
        %% Input Parsing
        p = inputParser;
        p.addParameter('dimName',  [],@(x) isempty(x)||iscell(x)||ischar(x));
        p.addParameter('dimSize',  [],@(x) isempty(x)||iscell(x)||(isnumeric(x)&&isscalar(x)));
        p.addParameter('dimLabels',[],@(x) isempty(x)||iscell(x)||ischar(x)||iscellstr(x));
        p.addParameter('dimUnits', [],@(x) isempty(x)||iscell(x)||ischar(x)||iscellstr(x));
        p.addParameter('dimValues',[],@(x) isempty(x)||iscell(x)||(isnumeric(x)&&isvector(x)));
        p.parse(varargin{:});
                        
        [parsed,nOutput] =  obj.parseInputsIntoCells(p.Results);
        
        %% Single output object
        if nOutput==1
          obj.dimSize   = parsed.dSize{1};
          obj.dimName   = parsed.dName{1};
          obj.dimLabels = parsed.dLabels{1};
          obj.dimUnits  = parsed.dUnits{1};
          obj.dimValues = parsed.dValues{1};
          return;
        end
        
        %% Array of output objects
        obj(nOutput) = arrayDim;
        for i = 1:nOutput
          obj(i) = arrayDim('dimName',parsed.dName{i},...
                            'dimSize',parsed.dSize{i},...
                            'dimLabels',parsed.dLabels{i},...
                            'dimUnits',parsed.dUnits{i},...
                            'dimValues',parsed.dValues{i});
        end;              
      end;
    end
    
    function name = get.dimName(obj)
      name = obj.dimName_;
    end;
    
    function obj = set.dimName(obj,val)
      assert(isempty(val)||ischar(val),'Dimension name must be a character string');
      obj.dimName_ = val;
      assert(obj.isInternallyConsistent);
    end
    
    %% Get/Set obj.dimSize
    %%%%%%%%%%%%%%%%%%%%%%
    function size = allDimSize(obj)
      % Return an array of all dimensions sizes
      size = cell2mat({obj.dimSize});
    end
    
    function isFixed = get.fixedSize(obj)
      if isempty(obj.dimSize_)
        isFixed = false;
      else
        isFixed = true;
      end;
    end
    
    function size = get.dimSize(obj)
      % Get the size of an individual dimension
      %
      % If not specifically defined, the dimension size is the length of
      % the longest assigned attribute.
      %
      % Otherwise empty dimensions have a size of 1.
      %
      if isempty(obj.dimSize_)
        if ischar(obj.dimLabels)
          nLabels = 1;
        else
          nLabels = numel(obj.dimLabels);
        end;
        if ischar(obj.dimUnits)
          nUnits = 1;
        else
         nUnits  = numel(obj.dimUnits);
        end;
        nValues = numel(obj.dimValues);
        size = max([nLabels nUnits nValues 1]);
      else
        size = obj.dimSize_;
      end;
    end
    
    function obj = set.dimSize(obj,val)
      assert(isempty(val)||(isnumeric(val)&&isscalar(val)),'Dimension size must be a numeric scalar');
      obj.dimSize_ = val;
    end
    
    %% Get/Set obj.dimLabels
    %%%%%%%%%%%%%%%%%%%%%%%%
    function labels = get.dimLabels(obj)                  
      labels = obj.dimLabels_;
    end;
    
    function obj = set.dimLabels(obj,val)
      obj.dimLabels_ = val;
      assert(obj.isInternallyConsistent);
    end;
    
    %% Get/Set obj.dimUnits
    %%%%%%%%%%%%%%%%%%%%%%%
    function units = get.dimUnits(obj)
      units = obj.dimUnits_;
    end
    
    function obj = set.dimUnits(obj,val)
      obj.dimUnits_ = val;
      assert(obj.isInternallyConsistent);
    end    
    
    %% Get/Set obj.dimValues
    %%%%%%%%%%%%%%%%%%%%%%%%
    function values = get.dimValues(obj)
      values = obj.dimValues_;
    end
    
    function obj = set.dimValues(obj,val)
      obj.dimValues_ = val;
      assert(obj.isInternallyConsistent);
    end;
                  
    %% Overloaded Functions
    %%%%%%%%%%%%%%%%%%%%%%%        
    function isEqual = isequal(obj,a)
      isEqual = true;
      if ~(isequal(obj.dimName,a.dimName))
        isEqual = false;
        warning('Dimension name mismatch');
      end
      if ~(isequal(obj.dimLabels,a.dimLabels))
        isEqual = false;
        warning('Dimension label mismatch');
      end
      
      if ~(isequal(obj.dimUnits,a.dimUnits))
        isEqual = false;
        warning('Dimension unit mismatch');
      end
      
      if ~(isequal(obj.dimValues,a.dimValues))
        isEqual = false;
        warning('Dimension value mismatch');
      end
    end
    
    function isEmpty = isempty(obj)
      isEmpty = false;
      if numel(obj)==1
        isEmpty = true;
        % Empty only if all components are empty
        isEmpty = isEmpty&&isempty(obj.dimName);
        isEmpty = isEmpty&&isempty(obj.dimLabels);
        isEmpty = isEmpty&&isempty(obj.dimUnits);
        isEmpty = isEmpty&&isempty(obj.dimValues);
      end
    end
    
    function isEmpty = isMostlyEmpty(obj)
      % Allows dimension name to have been set.
      isEmpty = false;
      if numel(obj)==1
        isEmpty = true;
        % Empty only if all components are empty        
        isEmpty = isEmpty&&isempty(obj.dimLabels);
        isEmpty = isEmpty&&isempty(obj.dimUnits);
        isEmpty = isEmpty&&isempty(obj.dimValues);
      end
    end    
    
    function obj = cat(dim,obj,a,varargin)
      
      assert(isMutuallyConsistent(obj,a,...
                                    'bsxValid',false,...
                                    'sizeMismatchValid',true,...
                                    'sizeMismatchDim',dim,...
                                    'nExtraDimsValid',0));
           
      % Copy info for dimensions that haven't been fully defined yet
      for i = 1:numel(obj)
        if isMostlyEmpty(obj(i))
          obj(i) = a(i);
        end
      end
                                  
      if ~isempty(obj(dim).dimSize_)&&~isempty(a(dim).dimSize_)
        obj(dim).dimSize_ = obj(dim).dimSize_ + a(dim).dimSize_;
      else
        % Clear if both aren't explicitly defined.
        obj(dim).dimSize_ = [];
      end;
            
      obj(dim).dimLabels_ = cat(2,checkType(obj(dim).dimLabels),checkType(a(dim).dimLabels));
      obj(dim).dimUnits_  = cat(2,checkType(obj(dim).dimUnits),checkType(a(dim).dimUnits));
      obj(dim).dimValues_ = cat(2,obj(dim).dimValues,a(dim).dimValues);
      assert(obj(dim).isInternallyConsistent);
      
      function out = checkType(in)
        if isempty(in)
          out = [];
        else
          out = cellstr(in);
        end
      end
      
    end
    
    function [obj,varargout] = subselectDimensions(obj,varargin)      
      idx = obj.getIndexIntoDimensions(varargin{:});
      
      for idxDim = 1:numel(obj)
        obj(idxDim) = obj(idxDim).subcopy(idx{idxDim});
      end;   
      
      if nargout>1
        varargout{1} = idx;
      end;      
    end
    
    
    function idxOut = getIndexIntoDimensions(obj,varargin)
      % Index into one or more dimensions
      %
      % idxOut = getIndexIntoDimensions(obj,varargin)
      % 
      % Inputs
      % ------
      %       obj : arrayDim object to index into
      %  varargin : Index into individual elements of the obj array.
      %                 numel(varargin)<=numel(obj) must be satisfied.      
      %
      % Output
      % ------
      %   idxOut : F
      %
      
      if numel(obj)>1
        idxOut = cell(1,numel(obj));
        for idxDim = 1:numel(obj)
          if idxDim<=numel(varargin)
            idxOut{idxDim} = obj(idxDim).getIndexIntoDimensions(varargin{idxDim});
          else
            idxOut{idxDim} = obj(idxDim).getIndexIntoDimensions(':');
          end;
        end    
        return;
      end
                                        
      if isempty(obj.dimLabels)
        cellIn = [];
        isStringValid = false;
      else
        cellIn = obj.dimLabels;
        isStringValid = true;
      end;
      
%       if isempty(obj.dimValues)
%         valIn= [];
%         isValueValid = false;
%       else
%         valIn = obj.dimValues;
%         isValueValid = false;
%       end;
      
       % Find the correct numeric indices
      if isequal(varargin{1},':')
        %% Requested Everything
        idxOut = varargin{1};
        return;
        
      elseif islogical(varargin{1})
        %% Logical Indexing
        assert(numel(varargin{1})==size(obj,dim),...
          'Logical indexing must match the size of the dimension');
        idxOut = find(varargin{1});
        return;
        
      elseif isnumeric(varargin{1})
        %% Numeric Reference
        if any(varargin{1}<1)||any(varargin{1}>obj.dimSize)
          error('Requested index outside of available range');
        end;
        idxOut = varargin{1};
        %idxOut(idxOut<1) = nan;
        %idxOut(idxOut>numel(cellIn)) = nan;
        return;
        
      elseif ischar(varargin{1})||iscellstr(varargin{1})
        %% String Reference
        if ~isStringValid
          error('String indexing unavailable for this dimension');
        end;
        
        if ischar(varargin{1}), varargin{1} = {varargin{1}}; end;
        cellIn = strtrim(cellIn);
        varargin{1} = strtrim(varargin{1});
        idxOut = zeros(1,numel(varargin{1}));
        for idx = 1:numel(idxOut)
          tmp = find(strcmp(varargin{1}{idx},cellIn));
          if isempty(tmp)
            error('Requested string does not appear in cell array');
          end;
          assert(numel(tmp)==1,'Multiple string matches in cellIn');
          idxOut(idx) = tmp;
        end
      elseif iscell(varargin{1})
        % Value Indexing
        error('Value indexing not yet implemented');
      else
        %% Otherwise, error.
        error('Incorrect reference type');
      end;
      
      
    end
    
    
    function idxOut = findDimensions(obj,idxIn)
      %% Get the numeric indices of dimensions by name or numeric reference
      %
      % idxOut = findDimensions(obj,idxin)
      %
      % Inputs
      % ------
      %    obj : labelledArray object
      %  idxIn : New dimension ordering for the object
      %           Can be:
      %             1) A numeric vector indexing the dimensions
      %             2) A cell array, with either numeric dimension indices
      %                 or character strings referencing index names.
      %
      % Outputs
      % -------
      %  idxOut : Numeric reference to each requested dimension.
      %
      if ischar(idxIn), idxIn = {idxIn}; end;
      if isnumeric(idxIn)
        idxOut = idxIn;
      elseif iscell(idxIn)
        idxOut = nan(1,numel(idxIn));
        for i = 1:numel(idxIn)
          if isnumeric(idxIn{i})&&isscalar(idxIn{i})
            idxOut(i) = idxIn{i};
          elseif ischar(idxIn{i})
            validNames = {obj.dimName};            
            validNames = validNames(~cellfun(@isempty,validNames));
            if ~isempty(validNames)
             inName = validatestring(idxIn{i},validNames);
            
             matched = cellfun(@(x) isequal(x,inName),{obj.dimName});
             idxOut(i) = find(matched);
            else
              error('Name indexing unavailable');
            end;
          else
            error('Invalid dimension reference');
          end;
        end
      else
        error('Invalid dimension referencing');
      end
    end;
    
    function [dimOut] = getConsistentDimensions(obj,a)
      % Get a new set of mutually consistent dimensions
      nDimsOut = max(numel(obj),numel(a));
      
      if nDimsOut==1
        assert(isMutuallyConsistent(obj,a),'Inconsistent decomposition sizes');
        
        if (obj.dimSize>1)||(a.dimSize==1)
          dimOut = obj;
        else
          dimOut = a;
        end;        
      else  
        dimOut(nDimsOut) = arrayDim;
        for i = 1:nDimsOut
          if i>numel(obj)
            dimOut(i) = a(i);
          elseif i>numel(a)
            dimOut(i) = obj(i);
          else
            dimOut(i) = getConsistentDimensions(obj(i),a(i));
          end;
        end
      end
                  
    end
    
    function isValid = isMutuallyConsistent(obj,a,varargin)
      % Check consistency between two sets of array dimensions
      %
      % function isValid = isMutuallyConsistent(obj,a,varargin)
      %
      % Inputs
      % ------
      %       obj : 
      %         a :
      %  bsxValid : 
      %
      % Output
      % ------
      %    isValid : Binary logical value
      %
      
      p = inputParser;
      p.addParameter('bsxValid',true,@islogical);
      p.addParameter('testStrict',true,@islogical);      
      p.addParameter('sizeMismatchValid',false,@islogical);
      p.addParameter('nExtraDimsValid',[],@(x) isnumeric(x)&&isscalar(x));
      p.addParameter('sizeMismatchDim',[],@(x) isnumeric(x)&&isscalar(x));
      p.parse(varargin{:});
                        
      assert(isa(a,'arrayDim'),'Second input must be an arrayDim object');
      
      if (numel(obj)==1)&&(numel(a)==1)
        % Directly compare two dimensions
                
        if p.Results.bsxValid
          % bsxfun can be applied if one of the dimension sizes is 1, or if
          % the dimensions otherwise match.
          isValid = xor((obj.dimSize==1),(a.dimSize==1));
          if isValid, return; end;
        end;
        
        if p.Results.sizeMismatchValid&&p.Results.testStrict
          % If mismatched sizes are valid, just want to check the
          % dimensions names under strict testing
          isValid = isEmptyOrEqual(obj.dimName,a.dimName);
          return;
        end
                                
        if isMostlyEmpty(obj)||isMostlyEmpty(a)
          % If one or the other is mostly empty (only a name assigned, if
          % anything), then they are compatible if they share the same
          % name.
          if ~isEmptyOrEqual(obj.dimName,a.dimName), return; end;
          isValid = true;
          return;
        end;
        
        % Size Must be Equal
        isValid = false;
        if ~isequal(obj.dimSize,a.dimSize), return; end;                
        % Strict Testing
        if p.Results.testStrict
          if ~isEmptyOrEqual(obj.dimName,  a.dimName),   return; end;
          if ~isEmptyOrEqual(obj.dimLabels,a.dimLabels), return; end;
          if ~isEmptyOrEqual(obj.dimUnits, a.dimUnits),  return; end;
          if ~isEmptyOrEqual(obj.dimValues,a.dimValues), return; end;
        end;
        isValid = true;
        return;
      else
        
        if ~isempty(p.Results.nExtraDimsValid)
          % Limit the difference in the number of extra dimensions
          isValid = abs(numel(obj)-numel(a))<=p.Results.nExtraDimsValid;
          if ~isValid, return; end;
        end;
        
        % Compare Multiple Dimensions        
        isValid = true;
        for idxDim = 1:min([numel(obj) numel(a)])          
          isValid = isValid && ...
                    isMutuallyConsistent(obj(idxDim),a(idxDim),...
                           'bsxValid',p.Results.bsxValid,...
                           'sizeMismatchValid',...
                               ismember(idxDim,p.Results.sizeMismatchDim));          
        end
      end
      
      function isValid = isEmptyOrEqual(A,B)
        isValid = isempty(A)||isempty(B);
        isValid = isValid||isequal(A,B);               
      end
           
    end
    
    
    
    function [obj,varargout] = subcopy(obj,varargin)
            
      idx = getIndexIntoDimensions(obj,varargin{:});
      
      if numel(obj)>1
        for idxDim = 1:numel(obj)
          obj(idxDim) = obj(idxDim).subcopy(idx{idxDim});
        end
        return;
      end
      
      if ~isempty(obj.dimSize_)&&~isequal(idx,':')
        obj.dimSize_ = numel(idx);
      end
      
      if ~isempty(obj.dimLabels)
       obj.dimLabels_ = obj.dimLabels(idx);
      end
      
      if ~isempty(obj.dimValues)
        obj.dimValues_ = obj.dimValues(idx);
      end
      
      if ~isempty(obj.dimUnits)&&~ischar(obj.dimUnits)
        obj.dimUnits_ = obj.dimUnits(idx);
      end
                    
      if nargout>=2
        varargout{1} = idx;
      end;
      
      assert(obj.isInternallyConsistent);
    end    
    
  end
  
  %% Private Methods
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
  
  methods 
    
    function isOk = isInternallyConsistent(obj)
      % Check the internal consistency of the object
      
      if numel(obj)>1
        isOk = true;
        for idxDim = 1:numel(obj)
          isOk = isOk&&obj(idxDim).isInternallyConsistent;
        end
        return;
      end
      
      assert(isempty(obj.dimName)||ischar(obj.dimName),...
              'Dimension name must be a character string');
            
      if ischar(obj.dimLabels)
        nLabels = 1;
      else
        nLabels = numel(obj.dimLabels);
      end
      if ischar(obj.dimUnits)
        nUnits = 1;
      else
        nUnits  = numel(obj.dimUnits);
      end
      nValues = numel(obj.dimValues);
         
      assert((nLabels==0)||(nLabels==obj.dimSize),...
                'Inconsistent number of labels');
      assert((nUnits==0)||(ischar(obj.dimUnits))||(nUnits==obj.dimSize),...
                'Inconsistent number of units');
      assert((nValues==0)||(nValues==obj.dimSize),...
                'Inconsistent number of values');
              
      assert(issorted(obj.dimValues)||...
             issorted(flip(obj.dimValues,2)),'Value array must be sorted');              
              
      isOk = true;      
    end
                   
  end
  
  methods (Static=true)
    
   function [validOut,nOutput] = parseInputsIntoCells(Results)
      % Input parser for the main object constructor. 
      %
      % Enables simultaneous construction of multiple objects
      %
      %
      [dName,nameSize]    = convertToCell(Results.dimName,@iscell);
      [dSize,sizeSize]    = convertToCell(Results.dimSize,@iscell);
      [dLabels,labelSize] = convertToCell(Results.dimLabels,...
        @(x) iscell(x)&&~iscellstr(x));
      [dUnits,unitSize]   = convertToCell(Results.dimUnits,...
        @(x) iscell(x)&&~iscellstr(x));
      [dValues,valueSize] = convertToCell(Results.dimValues, @iscell);
      
      nOutput = max([1 nameSize sizeSize labelSize unitSize valueSize]);
      
      validOut.dName   = validateCell(dName,nOutput);
      validOut.dSize   = validateCell(dSize,nOutput);
      validOut.dLabels = validateCell(dLabels,nOutput);
      validOut.dUnits  = validateCell(dUnits,nOutput);
      validOut.dValues = validateCell(dValues,nOutput);
      
      function [cellOut, size] = convertToCell(input,checkFHandle)
        % Convert input into correct cell form
        if checkFHandle(input)
          size = numel(input);
          cellOut = input;
        else
          size = 1;
          cellOut = {input};
        end;
      end
      
      function cellOut = validateCell(cellIn, nOutput)
        % Validate the size of each input cell.
        isEmpty = (numel(cellIn)==1)&&isempty(cellIn{1});
        sizeMatch = numel(cellIn)==nOutput;
        assert(isEmpty||sizeMatch,'Input size mismatch');
        if isEmpty
          cellOut = cell(nOutput,1);
        else
          cellOut = cellIn;
        end;
      end
      
    end    
    
    
  end
  
  
end
