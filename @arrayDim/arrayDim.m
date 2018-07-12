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
				isEmptyOrCell = @(x) isempty(x)||iscell(x);
				isNumericScalar = @(x) isnumeric(x)&&isscalar(x);
				isNumericVector = @(x) isnumeric(x)&&isscalar(x);
				isCharOrCellStr = @(x) isnumeric(x)||iscellstr(x);


        p = inputParser;
        p.addParameter('dimName',  [],@(x) isEmptyOrCell(x)||ischar(x));
        p.addParameter('dimSize',  [],@(x) isEmptyOrCell(x)||isNumericVector(x);
        p.addParameter('dimLabels',[],@(x) isEmptyOrCell(x)||isCharOrCellStr(x);
        p.addParameter('dimUnits', [],@(x) isEmptyOrCell(x)||isCharOrCellStr(x);
        p.addParameter('dimValues',[],@(x) isEmptyOrCell(x0||isNumericVector(x);
        p.parse(varargin{:});
                        
        [parsed,nOutput] =  obj.parseInputsIntoCells(p.Results);
        
        %% Multiple output object
        if nOutput>1
          %% Array of output objects
          obj(nOutput) = arrayDim;
          for i = 1:nOutput
            obj(i) = arrayDim('dimName',parsed.dName{i},...
                              'dimSize',parsed.dSize{i},...
                              'dimLabels',parsed.dLabels{i},...
                              'dimUnits',parsed.dUnits{i},...
                              'dimValues',parsed.dValues{i});
          end;              
					return;
				end;

          obj.dimSize   = parsed.dSize{1};
          obj.dimName   = parsed.dName{1};
          obj.dimLabels = parsed.dLabels{1};
          obj.dimUnits  = parsed.dUnits{1};
          obj.dimValues = parsed.dValues{1};
        
      end;
    end; % END Constructor
    
    function name = get.dimName(obj)
      name = obj.dimName_;
    end;
   
		%% Get/Set obj.dimName
		%%%%%%%%%%%%%%%%%%%%%% 
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
			% Determine if dimension has a fixed size
			% 
			% Not sure if this is entirely necessary.
      if isempty(obj.dimSize_)
        isFixed = false;
      else
        isFixed = true;
      end;
    end
    
    function sizeOut = get.dimSize(obj)
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
        sizeOut = max([nLabels nUnits nValues 1]);
      else
        sizeOut = obj.dimSize_;
      end;
    end
    
    function obj = set.dimSize(obj,val)
      assert(isempty(val)||(isnumeric(val)&&isscalar(val)),...
											'Dimension size must be a numeric scalar');
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
    
    function obj = set.dimLabels_(obj,val)
			%% Make sure dimLabels_ is always a row vector
      if (size(val,1)>1)
        val = val';
      end
			assert(isvector(val),'Dimension labels must be a row vector');
      obj.dimLabels_ = val;
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

    function obj = set.dimUnits_(obj,val)
      if (size(val,1)>1)
        val = val';
      end;
			assert(isempty(val)||isvector(val),...
								'Dimension units must be a row vector');
      obj.dimUnits_ = val;
    end;
    
    %% Get/Set obj.dimValues
    %%%%%%%%%%%%%%%%%%%%%%%%
    function values = get.dimValues(obj)
      values = obj.dimValues_;
    end
    
    function obj = set.dimValues(obj,val)
      obj.dimValues_ = val(:)';
      assert(obj.isInternallyConsistent);
    end;

    function obj = set.dimValues_(obj,val)
      if size(val,1)>1
          val = val';
      end;
			assert(isempty(val)||isvector(val),...
							'Dimension values must be a row vector');
      obj.dimValues_ = val;
    end
                  
    %% Overloaded Functions
    %%%%%%%%%%%%%%%%%%%%%%%        
    function isEqual = isequal(obj,a)

			if numel(obj)~=numel(a)
				isEqual = false;
				warning('Dimensionality mismatch');
				return;
			end;

			% Recurse for multiple dimensions
			isEqual = true;
			for i = 1:numel(obj)
  
  			% Check name equality
        isDimEqual = true;
        if ~(isequal(obj.dimName,a.dimName))
          isDimEqual = false;
          warning(['Dimension name mismatch in dimension:' num2str(i)]);
       end
  
  			% Check label equality
        if ~(isequal(obj.dimLabels,a.dimLabels))
          isDimEqual = false;
          warning(['Dimension label mismatch in dimension:' num2str(i)]);
        end
       
  			% Check Unit Equality 
        if ~(isequal(obj.dimUnits,a.dimUnits))
          isDimEqual = false;
          warning(['Dimension unit mismatch in dimension:' num2str(i)]);
        end
       
  			% Check Value Equality 
        if ~(isequal(obj.dimValues,a.dimValues))
          isDimEqual = false;
          warning(['Dimension value mismatch in dimension:' num2str(i)]);
        end;

				isEqual = isEqual&&isDimEqual;

				if ~isEqual
					break;
				end;
			end;
    end
    
    function isEmpty = isempty(obj)

			isEmpty = true;
			for i = 1:numel(obj)
        isDimEmpty = true;
        % Empty only if all components are empty
        isDimEmpty = isDimEmpty&&isempty(obj(i).dimName);
        isDimEmpty = isEmpty&&isempty(obj(i).dimLabels);
        isDimEmpty = isEmpty&&isempty(obj(i).dimUnits);
        isDimEmpty = isEmpty&&isempty(obj(i).dimValues);
	
				isEmpty = isEmpty&&isDimEmpty;
      end
    end
    
    function isEmpty = isMostlyEmpty(obj)
      % Allows dimension name to have been set.
      isEmpty = true;
		  for i = 1:numel(obj)
        isDimEmpty = true;
        % Empty only if all components are empty        
        isDimEmpty = isDimEmpty&&isempty(obj(i).dimLabels);
        isDimEmpty = isDimEmpty&&isempty(obj(i).dimUnits);
        isDimEmpty = isDimEmpty&&isempty(obj(i).dimValues);
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
          if (i~=dim)
            obj(i) = a(i);
          else
            % Copy the name for the concatenated dimension
            obj(i).dimName = a(i).dimName;
          end
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
      % Select a subset of samples along one or more dimensions
      %
      % Inputs:
      % ------
      %      obj : A single arrayDim object, or a 1D array of them
      % varargin : Indexing arguments
      %             The elements of varargin in are treated as 
      %             indexes into each of the components of the 
      %             input object. 
      %             If numel(varargin)<numel(obj),
      %               ':' indexing is used for all remaining dimensions.
      %             If numel(varargin)>numel(obj), extra indices are 
      %               ignored
      %
      % Outputs:
      % --------
      %   obj :  arrayDim object with subselection performed
      %
      % 
      % 
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
     
      % Recurse if more than one object input 
      if numel(obj)>1
        idxOut = cell(1,numel(obj));
        for idxDim = 1:numel(obj)
          if idxDim<=numel(varargin)
            idxOut{idxDim} = obj(idxDim).getIndexIntoDimensions(varargin{idxDim});
          else
            % Use ':' indexing if not all dimensions specified.
            idxOut{idxDim} = obj(idxDim).getIndexIntoDimensions(':');
          end;
        end    
        return;
      end
     
      % Should have numel(obj)==1 and numel(varargin)==1 below this line.
      assert(numel(obj)==1,'This shouldn''t happen');
      assert(numel(varargin)==1,'This shouldn''t happen');
     
      refIn = varargin{1};
      nOut = numel(refIn);
      
      % String referencing is valid if dimension labels are assigned
      isStringValid = true;
      if isempty(obj.dimLabels),  isStringValid = false;  end;
   
      % Value referencing is valid if values are assigned 
      isValueValid = false; % Turned off for now
      if isempty(obj.dimValues), isValueValid = false; end;
        
      %% Select According to Reference Type and Find Correct Numeric Indices
      if isequal(refIn,':')
        %% Requested Everything
        idxOut = refIn;
        return;
        
      elseif islogical(refIn)
        %% Logical Indexing
        assert(numel(refIn)==obj.dimSize,...
          'Logical indexing must match the size of the dimension');
        idxOut = find(refIn);
        return;
        
      elseif isnumeric(refIn)
        %% Numeric Indexing
        if any(refIn<1)||any(refIn>obj.dimSize)
          error('Requested index outside of available range');
        end;
        idxOut = refIn;
        return;
        
      elseif ischar(refIn)||iscellstr(refIn)
        %% String IndexingA
        assert(isStringValid,'String indexing unavailable for this dimension');
        
        if ischar(refIn), refIn = {refIn}; end;
        cellIn = strtrim(obj.dimLabels);
        refIn = strtrim(refIn);
        idxOut = zeros(1,numel(refIn));
        for idx = 1:numel(idxOut)
          tmp = find(strcmp(refIn{idx},cellIn));
          if isempty(tmp)
            error('Requested string does not appear in cell array');
          end;
          assert(numel(tmp)==1,'Multiple string matches in cellIn');
          idxOut(idx) = tmp;
        end
      elseif iscell(refIn)
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
				% Numeric Dimension Indexing
				assert(all(idxIn<numel(obj)), ...
									'Index exceeds object dimensions');
        idxOut = idxIn;
      elseif iscell(idxIn)
				% Cell referencing by dimension name or number
        idxOut = nan(1,numel(idxIn));
        for i = 1:numel(idxIn)
          if isnumeric(idxIn{i})&&isscalar(idxIn{i})
						% Numeric indexing within a cell array
						assert(idxIn{i}<=numel(obj),...
									'Index exceeds object dimensions');
            idxOut(i) = idxIn{i};
          elseif ischar(idxIn{i})
						% Referencing by dimension name
            validNames = {obj.dimName};            
            validNames = validNames(~cellfun(@isempty,validNames));
            if isempty(validNames)
							error('Name indexing unavailable');
						end;
            inName = validatestring(idxIn{i},validNames);
            
            matched = cellfun(@(x) isequal(x,inName),{obj.dimName});
            idxOut(i) = find(matched);
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
      %
      % [dimOut] = getConsistentDimensions(obj,a)
      %
      % Inputs
      % ------
      %  obj :
      %    a : arrayDim object
      %
      % Outputs
      % -------
      %  dimOut : New arrayDim object
      %

      nDimsOut = max(numel(obj),numel(a));
      
      if nDimsOut==1
        assert(isMutuallyConsistent(obj,a,...
                            'sizeMismatchValid',true,...
                            'nameMismatchValid',true),'Inconsistent decomposition sizes');
        
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
      p.addParameter('nameMismatchValid',false,@islogical);
      p.addParameter('nExtraDimsValid',[],@(x) isnumeric(x)&&isscalar(x));
      p.addParameter('sizeMismatchDim',[],@(x) isnumeric(x)&&isscalar(x));
      p.parse(varargin{:});
                        
      assert(isa(a,'arrayDim'),'Second input must be an arrayDim object');
      
      if (numel(obj)==1)&&(numel(a)==1)
        % Directly compare two dimensions

        if isMostlyEmpty(obj)||isMostlyEmpty(a)
          % If one or the other is mostly empty (only a name assigned, if
          % anything), then they are compatible.           
          isValid = true;
          return;
        end;                

        if ~p.Results.nameMismatchValid
          % By default, check that the names match
          isValid = isEmptyOrEqual(obj.dimName,a.dimName);
          if ~isValid, return; end;
        end;

        if p.Results.bsxValid
          % bsxfun can be applied if one of the dimension sizes is 1, or if
          % the dimensions otherwise match.
          isValid = xor((obj.dimSize==1),(a.dimSize==1));
          if isValid, return; end;
        end;
        
        if p.Results.sizeMismatchValid&&p.Results.testStrict
          % If mismatched sizes are valid, just want to check the
          % dimensions names under strict testing
          %isValid = isEmptyOrEqual(obj.dimName,a.dimName);
          isValid = true;
          return;
        end
                                       
        % Size Must be Equal
        isValid = false;
        if ~isequal(obj.dimSize,a.dimSize), return; end;                
        % Strict Testing
        if p.Results.testStrict          
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
                           'testStrict',p.Results.testStrict,...
                           'nameMismatchValid',p.Results.nameMismatchValid,...
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
