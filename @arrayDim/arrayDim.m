classdef arrayDim < handle & matlab.mixin.Copyable
  % Object class for labelledArray dimension information
  %
  %
  %
  % Properties
  % ----------
  %     dimName : String containing the dimension name
  %     dimSize : Number of elements along that dimension
  %   fixedSize : Flag to determine whether dimSize is allowed to change
  %   dimLabels : Cellstr of names for each element of the array
  %    dimUnits : Either:
  %                  String containing a uniform set of units for all
  %                  values
  %                  Cellstr containing units defined for each individiual
  %                  element (Typically used when naming elements)
  %   dimValues : Numeric value for each element. Must be sorted in
  %                  ascending order.
  %    dimRange : Minimum and maximum values of obj.dimValues
  %
  %
  
  
  %% Public property access is through the dependent properties
  properties (Dependent)
    dimName
    dimSize    
    dimLabels
    dimUnits
    dimValues
    dimRange
  end
  
  properties (Dependent, Hidden=true)
    fixedSize;
  end;
  
  
  %% Actual values are stored in private properties 
  properties (Access=private)
    dimName_
    dimSize_
    dimLabels_
    dimUnits_
    dimValues_
  end
  
  %% PUBLIC METHODS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  methods
    
    function obj = arrayDim(varargin)
      %% Class constructor
      if nargin>0
        if isa(varargin{1},'arrayDim')
          % Return a copy of the input object
          obj = varargin{1}.copy;
          return;
        end
        
        %% Input Parsing
        isEmptyOrCell = @(x) isempty(x)||iscell(x);
        isNumericScalar = @(x) isnumeric(x)&&isscalar(x);
        isNumericVector = @(x) isnumeric(x)&&isvector(x);
        isCharOrCellStr = @(x) ischar(x)||iscellstr(x);
        
        
        p = inputParser;        
        p.addParameter('dimName'  , [], @(x) isEmptyOrCell(x)||ischar(x));
        p.addParameter('dimSize'  , [], @(x) isEmptyOrCell(x)||isNumericScalar(x)||isNumericVector(x));
        p.addParameter('dimLabels', [], @(x) isEmptyOrCell(x)||isCharOrCellStr(x));
        p.addParameter('dimUnits' , [], @(x) isEmptyOrCell(x)||isCharOrCellStr(x));
        p.addParameter('dimValues', [], @(x) isEmptyOrCell(x)||isNumericVector(x));
        p.parse(varargin{:});
        
        [parsed,nOutput] =  obj.parseInputsIntoCells(p.Results);
        
        %% Multiple output object
        if nOutput>1
          %% Array of output objects
          obj(nOutput) = arrayDim;
          for i = 1:nOutput
            obj(i) = arrayDim('dimName'  , parsed.dName{i}   , ...
              'dimSize'  , parsed.dSize{i}   , ...
              'dimLabels', parsed.dLabels{i} , ...
              'dimUnits' , parsed.dUnits{i}  , ...
              'dimValues', parsed.dValues{i} );
          end;
          return;
        end;
        
        %% Construct Single Object
        obj.dimSize   = parsed.dSize{1};
        obj.dimName   = parsed.dName{1};
        obj.dimLabels = parsed.dLabels{1};
        obj.dimUnits  = parsed.dUnits{1};
        obj.dimValues = parsed.dValues{1};
        
      end;
    end; % END Constructor
    
    
    function trySet(obj,field,val)
      oldVal = obj.(field);
      obj.(field) = val;
      try
        assert(obj.isInternallyConsistent);
      catch
        obj.(field) = oldVal;
      end
    end
    
    %% Get/Set obj.dimName
    %%%%%%%%%%%%%%%%%%%%%%
    function name = get.dimName(obj)
      name = obj.dimName_;
    end;
    
    function set.dimName(obj,val)
      obj.trySet('dimName_',val);      
    end
    
    function set.dimName_(obj,val)
      assert(isempty(val)||ischar(val),'Dimension name must be a character string');
      obj.dimName_ = val;
    end;
    
    %% Get/Set obj.dimSize and related properties
    %%%%%%%%%%%%%%%%%%%%%%
    function size = allDimSize(obj)
      % Return an array of all dimensions sizes
      %
      % size = allDimSize(obj)
      %
      % For single arrayDim objects, this is the same as obj.dimSize. When
      % passed an array of a
      size = cell2mat({obj.dimSize});
    end
    
    function isFixed = get.fixedSize(obj)
      % Determine if dimension has a fixed size
      %
      % Not sure if this is entirely necessary.
      warning('arrayDim.fixedSize is deprecated');
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
    
    function set.dimSize(obj,val)
      obj.trySet('dimSize_',val);      
    end
    
    function set.dimSize_(obj,val)
      assert(isempty(val)||(isnumeric(val)&&isscalar(val)),...
        'Dimension size must be a numeric scalar');
      obj.dimSize_ = val;
    end
    
    %% Get/Set obj.dimLabels
    %%%%%%%%%%%%%%%%%%%%%%%%
    function labels = get.dimLabels(obj)
      labels = obj.dimLabels_;
    end;
    
    function set.dimLabels(obj,val)
      obj.trySet('dimLabels_',val);
      %obj.dimLabels_ = val;
      %assert(obj.isInternallyConsistent);
    end;
    
    function set.dimLabels_(obj,val)
      % dimLabels_ can be:
      %     1) Empty
      %     2) A cellstr with one cell per channel
      %
      if (size(val,1)>1)
        val = val';
      end
      assert(isempty(val)||(iscellstr(val)&&isvector(val)),...
        'Dimension labels must be a single dimensional cellstr');
      obj.dimLabels_ = val;
    end;
    
    %% Get/Set obj.dimUnits
    %%%%%%%%%%%%%%%%%%%%%%%
    function units = get.dimUnits(obj)
      units = obj.dimUnits_;
    end
    
    function set.dimUnits(obj,val)
      obj.trySet('dimUnits_',val);
      %obj.dimUnits_ = val;
      %assert(obj.isInternallyConsistent);
    end
    
    function set.dimUnits_(obj,val)
      % dimUnits_ can be:
      %     1) Empty
      %     2) A single string (same units for all channels)
      %     3) A cellstr, with one label per channel
      %
      if (size(val,1)>1)
        val = val';
      end;
      assert(isempty(val)||ischar(val)||(iscellstr(val)&&isvector(val)),...
        'Dimension units must be a row vector');
      obj.dimUnits_ = val;
    end;
    
    %% Get/Set obj.dimValues
    %%%%%%%%%%%%%%%%%%%%%%%%
    function values = get.dimValues(obj)
      values = obj.dimValues_;
    end
        
    function set.dimValues(obj,val)
      obj.trySet('dimValues_',val);
      %obj.dimValues_ = val(:)';
      %assert(obj.isInternallyConsistent);
    end;
    
    function set.dimValues_(obj,val)
      if size(val,1)>1
        val = val';
      end;
      assert(isempty(val)||isvector(val),...
        'Dimension values must be a row vector');
      obj.dimValues_ = val;
    end
    
    function out = get.dimRange(obj)
      if ~isempty(obj.dimValues)
        out = [obj.dimValues_(1) obj.dimValues_(end)];
      else
        out = [];
      end;
    end
    
    %% Overloaded Functions
    %%%%%%%%%%%%%%%%%%%%%%%
    function isEqual = isequal(obj,a)
      
      %% Check overall dimensionality
      if numel(obj)~=numel(a)
        isEqual = false;
        warning('Dimensionality mismatch');
        return;
      end;
      
      %% Recurse for multiple dimensions
      if numel(obj)>1
        isEqual = true;
        for i = 1:numel(obj)
          isEqual = isEqual && isequal(obj(i),a(i));
          if ~isEqual, break;	end;
        end
        return;
      end;
      
      %% Check an individual dimension
      
      % Check name equality
      isEqual = true;
      if ~(isequal(obj.dimName,a.dimName))
        isEqual = false;
        warning(['Dimension name mismatch in dimension:' num2str(i)]);
      end
      
      % Check label equality
      if ~(isequal(obj.dimLabels,a.dimLabels))
        isEqual = false;
        warning(['Dimension label mismatch in dimension:' num2str(i)]);
      end
      
      % Check Unit Equality
      if ~(isequal(obj.dimUnits,a.dimUnits))
        isEqual = false;
        warning(['Dimension unit mismatch in dimension:' num2str(i)]);
      end
      
      % Check Value Equality
      if ~(isequal(obj.dimValues,a.dimValues))
        isEqual = false;
        warning(['Dimension value mismatch in dimension:' num2str(i)]);
      end;
      
    end; % END isEqual()
    
    function isEmpty = isempty(obj)
      % Check if fully empty
      isEmpty = true;
      for i = 1:numel(obj)
        isDimEmpty = true;
        % Empty only if all components are empty
        isDimEmpty = isDimEmpty&&(obj(i).dimSize==1);
        isDimEmpty = isDimEmpty&&isempty(obj(i).dimName);
        isDimEmpty = isDimEmpty&&isempty(obj(i).dimLabels);
        isDimEmpty = isDimEmpty&&isempty(obj(i).dimUnits);
        isDimEmpty = isDimEmpty&&isempty(obj(i).dimValues);
        
        isEmpty = isEmpty&&isDimEmpty;
      end
    end
    
    function isEmpty = isMostlyEmpty(obj)
      % Check if empty except for dimension name
      isEmpty = true;
      for i = 1:numel(obj)        
        % Empty only if all components are empty
        isEmpty = isEmpty&&isempty(obj(i).dimLabels);
        isEmpty = isEmpty&&isempty(obj(i).dimUnits);
        isEmpty = isEmpty&&isempty(obj(i).dimValues);
      end
    end
    
    function isUniform = isUniformlySampled(obj)
      
      tol = 1e-3;
      
      isUniform = true;
      if numel(obj)>1
        for i = 1:numel(obj)
          isUniform = isUniform&&isUniformlySampled(obj(i));
        end
        return;
      end
            
      if ~isempty(obj.dimValues)
        delta = obj.dimValues(2:end)-obj.dimValues(1:end-1);
        isUniform = all(abs(delta-mean(delta))<tol);
      end
      
    end
    
    
    function objOut = cat(dim,objA,objB,varargin)
      % Concatenate dimensions
      
      %% Check for Mutual Consistency
      assert(isMutuallyConsistent(objA,objB,...
        'bsxValid',false,...
        'sizeMismatchValid',true,...
        'sizeMismatchDim',dim,...
        'nExtraDimsValid',1));
      
      objOut = objA.copy;
      
      % Copy info for dimensions that haven't been fully defined yet      
      for i = 1:numel(objA)
        if isMostlyEmpty(objA(i))
          if (i~=dim)
            objOut(i) = objB(i).copy;
          else
            % Copy the name for the concatenated dimension
            if i<=numel(objB)
             objOut(i).dimName = objB(i).dimName;
            end;
          end
        end
      end
            
      if isMostlyEmpty(objA(dim))
        if  dim <=numel(objB)
          objOut(dim) = objB(dim).copy;
        end;
      elseif isMostlyEmpty(objB(dim))
        % Do Nothing
      else      
        objOut(dim).dimSize_   = catSize;
        objOut(dim).dimLabels_ = catLabels;
        objOut(dim).dimUnits_  = catUnits;
        objOut(dim).dimValues_ = catValues;            
      end
      assert(objOut(dim).isInternallyConsistent);
      
      function out = checkType(in)
        if isempty(in)
          out = [];
        else
          out = cellstr(in);
        end
      end
      
      function newSize = catSize
        if ~isempty(objA(dim).dimSize_)&&~isempty(objB(dim).dimSize_)
          newSize = objA(dim).dimSize_ + objB(dim).dimSize_;
        else
          % Clear if both aren't explicitly defined.
          newSize = [];
        end;
      end
      
      function newLabels = catLabels
        % Concatenate Labels
        labelsA = checkType(objA(dim).dimLabels);
        labelsB = checkType(objB(dim).dimLabels);
        
        if isempty(labelsA)&&~isempty(labelsB)
          tmp(1:objA(dim).dimSize) = {''};
          labelsA = tmp;
        end;
        if isempty(labelsB)&&~isempty(labelsA)
          tmp(1:objB(dim).dimSize) = {''};
          labelsB = tmp;
        end;
        newLabels = cat(2,labelsA,labelsB);
      end
      
      function newUnits = catUnits
        % Concatenate Units
        unitsA = checkType(objA(dim).dimUnits);
        unitsB = checkType(objB(dim).dimUnits);
        
        
        if ischar(unitsA)&&ischar(unitsB)
          % Both have a single dimensional unit
          if isequal(unitsA,unitsB)
            newUnits = unitsA;
          else
            error('AAARGH');
          end                  
        elseif ischar(unitsA)&&iscellstr(unitsB)
          % Expand Dimensional Units
          tmp(1:objA(dim).dimSize) = {unitsA};
          unitsA = tmp;
          newUnits = cat(2,unitsA,unitsB);          
          
        elseif ischar(unitsB)&&iscellstr(unitsA)
          % Expand Dimensional Units
          tmp(1:objB(dim).dimSize) = {unitsB};
          unitsB = tmp;
          newUnits = cat(2,unitsA,unitsB);          
          
        elseif isempty(unitsA)&&~isempty(unitsB)
          tmp(1:objA(dim).dimSize) = {''};
          unitsA = tmp;
          newUnits = cat(2,unitsA,unitsB);
          
        elseif isempty(unitsB)&&~isempty(unitsA)
          tmp(1:objB(dim).dimSize) = {''};
          unitsB = tmp;
          newUnits = cat(2,unitsA,unitsB);
          
        else
          % Both are cellstr
          newUnits = cat(2,unitsA,unitsB);
        end
      end
      
      function newValues = catValues
        newValues = cat(2,objA(dim).dimValues_,objB(dim).dimValues_);
      end
      
    end
    
    function [objOut,varargout] = subselectDimensions(objIn,varargin)
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
      idx = objIn.getIndexIntoDimensions(varargin{:});
      if ~iscell(idx), idx = {idx}; end;
      
      objOut = objIn.copy;
      for idxDim = 1:numel(objOut)
        objOut(idxDim) = objOut(idxDim).subcopy(idx{idxDim});
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
            idxOut(idxDim) = {obj(idxDim).getIndexIntoDimensions(varargin{idxDim})};
          else
            % Use ':' indexing if not all dimensions specified.
            idxOut(idxDim) = {obj(idxDim).getIndexIntoDimensions(':')};
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
      isValueValid = true; % Turned off for now
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
        assert(isValueValid,'Value indexing unavailable for this dimension');
        
        % Value Indexing
        assert(numel(refIn{1})==2,'Value referencing must be done with range limits');
        
        range = refIn{1};
        [~,lowIdx] = min(abs(obj.dimValues-range(1)));
        [~,highIdx] = min(abs(obj.dimValues-range(2)));
        
        assert(lowIdx<highIdx,'Poorly defined value range');
        
        idxOut = lowIdx:highIdx;
                
        %error('Value indexing not yet implemented');
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
          dimOut = obj.copy;
        else
          dimOut = a.copy;
        end;
      else
        dimOut(nDimsOut) = arrayDim;
        for i = 1:nDimsOut
          if i>numel(obj)
            dimOut(i) = a(i).copy;
          elseif i>numel(a)
            dimOut(i) = obj(i).copy;
          else
            dimOut(i) = getConsistentDimensions(obj(i),a(i));
          end;
        end
      end
      
    end
    
    function isConsistent = isMutuallyConsistent(obj,a,varargin)
      % Check consistency between two sets of array dimensions
      %
      % function isConsistent = isMutuallyConsistent(obj,a,varargin)
      %
      % Inputs
      % ------
      %       obj :
      %         a :
      %
      % Param-Value Inputs
      % ------------------
      %   testStrict :
      %   sizeMismatchValid: 
      %   nameMismatchValid:
      %   nExtraDimsValid:
      %   sizeMismatchDim:
      %  bsxValid :
      %
      % Output
      % ------
      %    isConsistent : Binary logical value
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
      
      %% Recurse for Multiple Dimensions
      if (numel(obj)>1)||(numel(a)>1)
        if ~isempty(p.Results.nExtraDimsValid)
          % Limit the difference in the number of extra dimensions
          isConsistent = abs(numel(obj)-numel(a))<=p.Results.nExtraDimsValid;
          if ~isConsistent, return; end;
        end;
        
        % Compare Multiple Dimensions
        isConsistent = true;
        for idxDim = 1:min([numel(obj) numel(a)])
          isConsistent = isConsistent && ...
            isMutuallyConsistent(obj(idxDim),a(idxDim),...
            'bsxValid',p.Results.bsxValid,...
            'testStrict',p.Results.testStrict,...
            'nameMismatchValid',p.Results.nameMismatchValid,...
            'sizeMismatchValid',...
            ismember(idxDim,p.Results.sizeMismatchDim));
          if ~isConsistent, break; end;
        end;
        return;
      end;
      
      %% Compare Single Dimension Below This Line            
      if isMostlyEmpty(obj)||isMostlyEmpty(a)
        % If one or the other is mostly empty (only a name assigned, if
        % anything), then they are compatible.
        isConsistent = isEmptyOrEqual(obj.dimName,a.dimName);        
        return;
      end;
      
      if ~p.Results.nameMismatchValid
        % By default, check that the names match
        isConsistent = isEmptyOrEqual(obj.dimName,a.dimName);
        if ~isConsistent, return; end;
      end;
      
      if p.Results.bsxValid
        % bsxfun can be applied if one of the dimension sizes is 1, or if
        % the dimensions otherwise match.
        isConsistent = xor((obj.dimSize==1),(a.dimSize==1));
        if isConsistent, return; end;
      end;
      
      if p.Results.sizeMismatchValid&&p.Results.testStrict
        % If mismatched sizes are valid, just want to check the
        % dimensions names under strict testing
        isConsistent = isEmptyOrEqual(obj.dimName,a.dimName);
        %isConsistent = true;
        return;
      end
      
      % Size Must be Equal
      isConsistent = false;
      if ~isequal(obj.dimSize,a.dimSize), return; end;
      % Strict Testing
      if p.Results.testStrict
        if ~isEmptyOrEqual(obj.dimLabels,a.dimLabels), return; end;
        if ~isEmptyOrEqual(obj.dimUnits, a.dimUnits),  return; end;
        if ~isEmptyOrEqual(obj.dimValues,a.dimValues), return; end;
      end;
      isConsistent = true;
            
      function isValid = isEmptyOrEqual(A,B)
        % Require at least one empty value, or equality.
        isValid = isempty(A)||isempty(B);
        isValid = isValid||isequal(A,B);
      end
      
    end
    
    function [objOut,varargout] = subcopy(objIn,varargin)
      % Copies a subset of the dimension into a new arrayDim object
      %
      % Inputs
      
      %% Get Numeric Index
      idx = getIndexIntoDimensions(objIn,varargin{:});
      
      %% Recurse for Multiple Objects
      objOut = objIn.copy;
      if numel(objOut)>1
        for idxDim = 1:numel(objOut)
          objOut(idxDim) = objOut(idxDim).subcopy(idx{idxDim});
        end
        return;
      end
      
      %% Single Object Below this Point    
      if isequal(idx,':')
        % Do nothing if we're keeping the whole dimension
        return;
      end
      
      if ~isempty(objOut.dimSize_)
        objOut.dimSize_ = numel(idx);
      end
      
      if ~isempty(objOut.dimLabels)
        objOut.dimLabels_ = objOut.dimLabels(idx);
      end
      
      if ~isempty(objOut.dimValues)
        objOut.dimValues_ = objOut.dimValues(idx);
      end
      
      if ~isempty(objOut.dimUnits)&&~ischar(objOut.dimUnits)
        objOut.dimUnits_ = objOut.dimUnits(idx);
      end
      
      if nargout>=2
        varargout{1} = idx;
      end;
      
      assert(objOut.isInternallyConsistent);
    end
    
  end
  
  %% Private Methods
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
  
  methods
    
    function isOk = isInternallyConsistent(obj)
      % Check the internal consistency of the object
      %
      % Returns true if object is internally consistent. Errors otherwise.
      %
      % An arrayDim object is deemed internally consistent if:
      %
      % 1) obj.dimName is either empty or a character string
      % 2) obj.
      
      %% Recurse for multiple objects
      if numel(obj)>1
        isOk = true;
        for idxDim = 1:numel(obj)
          isOk = isOk&&obj(idxDim).isInternallyConsistent;
        end
        return;
      end
      
      %% Check Dimension Name
      assert(isempty(obj.dimName)||ischar(obj.dimName),...
        'Dimension name must be a character string');
      
      %% Check Dimension Labels
      if ischar(obj.dimLabels)
        nLabels = 1;
      else
        nLabels = numel(obj.dimLabels);
      end
      assert((nLabels==0)||(nLabels==obj.dimSize),...
        'Inconsistent number of labels');
      
      %% Check Dimension Units
      if ischar(obj.dimUnits)
        nUnits = 1;
      else
        nUnits  = numel(obj.dimUnits);
      end
      assert((nUnits==0)||(ischar(obj.dimUnits))||(nUnits==obj.dimSize),...
        'Inconsistent number of units');
      
      %% Check Dimension Values
      nValues = numel(obj.dimValues);
      assert((nValues==0)||(nValues==obj.dimSize),...
        'Inconsistent number of values');
      assert(issorted(obj.dimValues)||...
        issorted(flip(obj.dimValues,2)),'Value array must be sorted');
      
      %%
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
      [dSize,sizeSize]    = convertToCell(Results.dimSize,@iscell,true);
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
      
      function [cellOut, size] = convertToCell(input,checkFHandle,expand)
        % Convert input into correct cell form
        if ~exist('expand','var'),expand = false; end;
        
        if expand
          input = num2cell(input);
        end;
        
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
        isEmpty = isempty(cellIn)||(numel(cellIn)==1)&&isempty(cellIn{1});
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
