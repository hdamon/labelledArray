classdef labelledArray < handle & matlab.mixin.Copyable
  % Handle array array with labelled axes
  %
  % obj = labelledArray(array,varargin)
  %
  % Inputs
  % ------
  %   array : The array to put in the array
  %
  % Optional Param-Value Inputs
  % ---------------------------
  %  'dimNames'  : Cell string array of names for each dimension.
  %  'dimLabels' : A cell array defining labels for the elements along
  %                   one or more dimensions.
  %                 Labels can be assigned through the constructor in one
  %                 of two ways:
  %                   1) Provide a cell array of size {nDim x 1} or
  %                       {1 X nDim}, where nDim is the number of
  %                       dimensions. Each element in the cell must either
  %                       be empty, or be a cellstr with numel equal to the
  %                       size along that dimension.
  %                   2) Provide a cell array of size (N X 2), arranged as:
  %                        { <DIMENSION A> , <CellString of dimLabels> ;
  %                          <DIMENSION B> , <CellString of dimLabels> }
  %                       <DIMENSION> can be either a scalar numeric value,
  %                       or a character string matching one of the values
  %                       in obj.dimNames. Each cellstring of dimLabels
  %                       must have a number of elements equal to the
  %                       current size of the array along that dimension.
  %  'dimUnits'  : A cell array of units for one or more dimensions. This
  %                   be provided either as a single character string (same
  %                   units for all elements along that dimension), or as a
  %                   cell array of strings with length equal to the size
  %                   of that dimension (Different units for each element).
  %  'dimValues' : A cell array of dimension values along one or more dimensions. Cell
  %               array should be formatted as for 'dimLabels', except instead
  %               of the dimValues being cell strings, they must be numeric
  %               vectors.
  %
  % %%%%%%%%%%%%
  % NOTE: When referring to this set of additional properties below, the acronym
  % NLUV is used.
  % %%%%%%%%%%%%
  %
  % Overloaded Functions
  % --------------------
  %    size : Behaves in one of two ways:
  %             1) For an array of labelledarray objects, returns the size
  %                   of the array.
  %             2) For a single labelledarray object, returns the size of
  %                   the array.
  %
  %   ndims : Behavines in one of two ways:
  %             1) For an array of labelledarray objects, returns the size
  %                   of the array.
  %             2) For a single labelledarray object, returns the size of
  %                   the array.
  %
  % NOTE: The easiest way to determine which behavior will be elicited is
  %         to check numel(obj).  If numel(obj)>1, you get standard
  %         behavior.
  %
  % Referencing into labelledarray Objects
  % --------------------------------------
  %   Referencing for labelledarray objects behaves slightly differently
  %   than for normal Matlab arrays.
  %
  %  obj.<property> : Behaves normally.
  %
  %  Numeric Indexing:
  %
  %  obj(<indices>) : Behaves in one of two ways:
  %                     1) If obj is an array of labelledarray objects, this
  %                          references into the array
  %                     2) If obj is a single labelledarray object, this
  %                          returns a new object with the array all
  %                          associated information subselected according to the
  %                          indices.
  %
  %  obj{<indices>} : References into an array of labelledarray objects.
  %                     This works correctly in combination with (), so
  %                     that expressions such as:
  %                       obj{i}(x,y,z)       and
  %                       obj{i}.array(x,y,z)
  %                     work as expected.
  %
  % Label Indexing:
  %
  % For dimensions that have had dimLabels assigned, the array can also be
  % directly indexed using those labels:
  %
  %   IE:     obj(:,{'A' 'B' 'C'},:)
  %       or  obj(:,'A',:)
  %
  % Will both work for a three dimensional labelledArray where elements
  % along the second dimension have been assigned the labels 'A', 'B', and
  % 'C',
  %
  %
  properties (Hidden,Dependent)
    array
    arrayRange
    dimensions
  end
  
  properties (Access=protected)
    dimNames
    dimSize
    dimLabels
    dimValues
    dimUnits
    dimRange
  end
  
  properties (Access=protected) % Should possibly be private?
    array_      % The array array
    dimensions_ % Dimension information
    arrayRange_ %
  end
  
  events
    updatedOut
  end
  
  methods
    
    %%
    function obj = labelledArray(array,varargin)
      %% Main object constructor
      if nargin>0
        
        if isa(array,'labelledArray')
          obj = array.copy;
          return;
        end;
                  
        checkType = @(x) iscell(x) && ((size(x,2)==2)||(numel(x)>=ndims(array)));
        
        % Input parsing
        p = inputParser;
        p.KeepUnmatched = true;
        p.addRequired('array',@(x) (isnumeric(x)||isa(x,'labelledArray')));
        p.addParameter('dimNames' , [] , @(x) isempty(x)||checkType(x));
        %p.addParameter('dimLabels', [] , @(x) isempty(x)||checkType(x));
        %p.addParameter('dimUnits' , [] , @(x) isempty(x)||checkType(x));
        %p.addParameter('dimValues', [] , @(x) isempty(x)||checkType(x));
        p.parse(array,varargin{:});
        
        % Property Assignment
        obj.array   = p.Results.array;
        dims            = arrayDim('dimName',p.Results.dimNames,...
                                    p.Unmatched);
                                   %'dimLabels',p.Results.dimLabels,...
                                   %'dimUnits',p.Results.dimUnits,...          
                                   %'dimValues',p.Results.dimValues);  
                               
        if ~isempty(dims)
          if numel(dims)>=obj.ndims
            obj.dimensions_ = dims;
          else
            error('Insufficient number of dimensions provided');
          end                         
        end;
        
      end
      
    end
    
    %% Overloaded Functions
    %%%%%%%%%%%%%%%%%%%%%%%
    
    function out = abs(obj)
      out = obj.copy;
      out.array = abs(out.array);
    end
    
    %%
    function out = size(obj,dim)
      %% Return the size of a labelledArray object
      %
      % out = size(obj,dim)
      %
      % Behaves in one of two ways:
      %  1) If obj is a single labelledArray object, returns the
      %				size of obj.array_
      %	 2) If obj is an array of labelledArray objects, returns the
      %				size of the array.
      %
      % Size will have trailing singleton dimensions if NLUV have been
      % assigned for them.
      %
      if numel(obj)==1
        if ~exist('dim','var')
          out = ones(1,obj.ndims);
          tmp = size(obj.array);
          out(1:numel(tmp)) = tmp;
        else
          out = size(obj.array,dim);
        end;
      else
        out = builtin('size',obj);
      end
    end
    
    %%
    function out = ndims(obj)
      % Get dimensionality of array
      %
      % Ignores trailing singleton dimensions unless they have been
      % assigned NLUV
      %
      %
      
      if numel(obj)==1
        out = numel(obj.dimensions_);        
      else
        out = builtin('ndims',obj);
      end;
    end
    
    function ind = end(obj,k,n)
      ind = obj.size(k);
    end
    
    %%
    function out = permute(obj,order)
      % Permute the order of the array
      %
      % out = permute(obj,order)
      %
      % Inputs
      % ------
      %  obj : labelledArray object
      %  order : One of two forms:
      %						1) A numeric array of dimension indices
      %						2) A cell array combining numeric indices and
      %								dimension dimNames.
      %
      
      newOrder = obj.getDimensionOrdering(order);
      
      out = obj.copy;
      out.array_      = permute(obj.array_,newOrder);
      out.dimensions_ = obj.dimensions(newOrder);      
    end
    
    function out = cat(dim,obj,a,varargin)
      % Concatenate labelledArray objects
      %
      %
      
      assert(isa(a,class(obj)),'Can only concatenate like objects');
      
      % Build output object      
      out = obj.copy;
      
      if numel(out.dimensions_)<dim
        out.dimensions_(dim) = arrayDim;
      end;      
      
      out.dimensions_ = cat(dim,out.dimensions_,a.dimensions_);
      out.array_ = cat(dim,obj.array_,a.array_);
            
      if ~isempty(varargin)        
        if numel(varargin)<5
         out = cat(dim,out,varargin{:});
        else
         % Split recursion when concatenating multiple objects          
         nSplit = ceil(numel(varargin)/2);                
         blockA = cat(dim,out,varargin{1:nSplit});
         blockB = cat(dim,varargin{(nSplit+1):end});        
         out = cat(dim,blockA,blockB);
        end;
      end;
      
    end
    
    %% Functions that always operate on the data in place
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function out = power(obj,b)
      out = obj.copy;
      out.array_ = power(obj.array_,b);      
    end;
    
    function out = sqrt(obj)
      out = obj.copy;
      out.array_ = sqrt(obj.array_);
    end;
    
    %% Helper function for functions that collapse a dimension to a singleton
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function out = applyDimFunc(funcHandle,obj,idxDim,varargin)
      % Apply function handle that collapses a dimension to a singleton
      %
      % Inputs
      % ------
      %   funcHandle : Function handle
      %                  Will be called as:
      %                  funcHandle(obj.array,varargin{:});
      %         obj  : labelledArray object to operate on
      %       idxDim : Index into varargin of the dimension argument
      %     varargin : Additional arguments for funcHandle. One of these
      %                 must be a numeric dimension index.
      %
      % Outputs
      % -------
      %  out = labelledArray object after function application
      %
      %
      
      dim = varargin{idxDim};
      dim = obj.findDimensions(dim); % Enable name based referencing
      varargin{idxDim} = dim;
      out = obj.copy;
      out.array_ = funcHandle(obj.array_,varargin{:});      
      
      if ~isempty(obj.dimValues{dim})
        newVal = mean(obj.dimValues{dim});
      else
        newVal = [];
      end;
      
      newName = obj.dimNames{dim};
      if isempty(obj.dimUnits{dim})
        newUnits = [func2str(funcHandle) '(' obj.dimNames{dim} ')'];
      elseif ischar(obj.dimUnits{dim})
        newUnits = [func2str(funcHandle) '(' out.dimUnits{dim} ')'];
      else
        newUnits = [func2str(funcHandle) '(' out.dimUnits{dim}{1} ')'];
      end;
      newDim = arrayDim('dimName',newName,...
                        'dimUnits', newUnits,...
                        'dimValues', newVal);
      out.dimensions(dim) = newDim;
    end

    %% Functions that use applyDimFunc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function out = sum(obj,dim)
      if ~exist('dim','var'),dim = 1; end;
      out = applyDimFunc(@sum,obj,1,dim);
    end;
    
    function out = min(obj,~,dim)
      if ~exist('dim','var'), dim = 1; end;
      out = applyDimFunc(@min,obj,2,[],dim);
    end;
    
    function out = max(obj,~,dim)
      if ~exist('dim','var'), dim = 1; end;
      out = applyDimFunc(@max,obj,2,[],dim);
    end;
    
    function out = mean(obj,dim)
      if ~exist('dim','var'), dim = 1; end;
      out = applyDimFunc(@mean,obj,1,dim);
    end
    
    function out = std(obj,W,dim)
      if ~exist('dim','var'), dim = 1; end;
      if ~exist('W','var'), W = 0; end;
      out = applyDimFunc(@std,obj,2,W,dim);
    end
    
    function out = var(obj,W,dim)
      if ~exist('dim','var'), dim = 1; end;
      if ~exist('W','var'), W = 0; end;
      out = applyDimFunc(@var,obj,2,W,dim);
    end;
    
    %% Class-specific bsxfun function
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function objOut = bsxfun(funcHandle,obj,b)
      % Use bsxfun to apply funcHandle to obj.tfX
      %
      % objOut = bsxfun(fHandle,obj,b)
      %
      % Inputs
      % ------
      %        obj : labelledData object (or subclass)
      %          b : Value to apply
      % funcHandle : Function handle
      %
      % Outputs
      % -------
      %   objOut : Output object of class(obj), with array field operated
      %               on by bsxfun.
      %
      % Checks the consistency of the inputs, and then calls:
      %
      % tfX = bsxfun(funcHandle,obj,tfX,coeff)
      %
      % When b is:
      %   A timeFrequencyDecomposition obj:  coeff = b.tfX
      %                          OTHERWISE:  coeff = b;
      %
      %
      
      for idxObj = 1:numel(obj)
        if isa(b,'labelledArray')
          assert(isMutuallyConsistent(obj(idxObj),b,...
                                      'nameMismatchValid',true),...
                 'Inconsistent decomposition sizes');
          coeff = b.array;
          dimsOut = getConsistentDimensions(obj(idxObj),b);
        else
          coeff = b;
          % New dimensions need to be handled better here
          dimsOut = obj.dimensions;          
        end;
        
        objOut(idxObj) = copy(obj(idxObj));
        objOut(idxObj).array_ = bsxfun(funcHandle,obj(idxObj).array,coeff);
        objOut(idxObj).arrayRange_ = [];
        objOut(idxObj).dimensions_ = dimsOut;        
      end;
      
      %% Reshape if its an array of objects
      if numel(objOut)>1
        objOut = reshape(objOut,size(obj));
      end;
    end
    
    %% Functions that use bsxfun
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function out = plus(obj,b)
      out = bsxfun(@plus,obj,b);
    end
    
    function out = minus(obj,b)
      out = bsxfun(@minus,obj,b);
    end
    
    function out = rdivide(obj,b)
      out = bsxfun(@rdivide,obj,b);
    end
    
    function out = ldivide(obj,b)
      out = bsxfun(@ldivide,obj,b);
    end
    
    function out = times(obj,b)
      out = bsxfun(@times,obj,b);
    end;    
        
    %% dimensions get/set methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function out = get.dimensions(obj)
      out = obj.dimensions_;
    end
    
    function set.dimensions(obj,val)
      obj.setDimensions(val);
      %obj.dimensions_ = val;
      %assert(isInternallyConsistent(obj.dimensions_));
    end
    
    %% array Get/Set Methods
    %%%%%%%%%%%%%%%%%%%%%%%%
    function set.array(obj,val)
      % Done this way so it can be overloaded by subclasses.
      obj.setArray(val);
    end;
    
    function out = get.array(obj)
      out = obj.array_;
    end;
    
%     function set.array_(obj,val)
%       obj.array_ = val;
%       %% Adjust dimension sizes
%       for i = 1:ndims(val)
%         if ~isempty(val)
%          if i<=numel(obj.dimensions_)
%            obj.dimensions(i).dimSize = size(val,i);
%          else
%            % Add a new dimension
%            obj.dimensions(i) = arrayDim('dimSize',size(val,i));
%          end
%         else
%         end
%       end;      
%     end
    
    %% dimNames Get/Set Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    function set.dimNames(obj,val)
      if iscellstr(val)
        for idxDim = 1:numel(val)
          obj.dimensions_(idxDim).dimName = val{idxDim};
        end;
        return;
      end
      obj.dimensions_.dimName = val;      
    end
    
    function out = get.dimNames(obj)
      out = {obj.dimensions_.dimName};      
    end;
    
    function set.dimSize(obj,val)     
        assert(numel(val)==numel(obj.dimensions_),...
                  'Must define sizes for all dimensions simultaneously');
        for idxDim = 1:numel(val)
          obj.dimensions_(idxDim).dimSize = val(idxDim);
        end      
    end
    
    function out = get.dimSize(obj)
      out = [obj.dimensions_.dimSize];
    end;
    
    function out = get.dimRange(obj)
      out = {obj.dimensions_.dimRange};
    end;
    
    
    %% dimLabels Get/Set Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function set.dimLabels(obj,val)
      if iscell(val)&&~iscellstr(val)
        for idxDim = 1:numel(val)
          obj.dimensions_(idxDim).dimLabels = val{idxDim};
        end;
        return;
      end
      obj.dimensions_.dimLabels = val;      
    end
    
    %%
    function out = get.dimLabels(obj)
      out = {obj.dimensions_.dimLabels};      
    end;
    
    %% dimUnits Get/Set Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    function set.dimUnits(obj,val)
      if (iscell(val)&&~iscellstr(val))||(iscellstr(val)&&numel(val)==obj.ndims)
        for idxDim = 1:numel(val)
          obj.dimensions_(idxDim).dimUnits = val{idxDim};
        end
        return;
      end
      obj.dimensions_.dimUnits = val;
    end;
    
    function out = get.dimUnits(obj)
      out = {obj.dimensions_.dimUnits};      
    end;
    
    %% dimValue Get/Set Methods
    %%%%%%%%%%%%%%%%%%%%%%%%
    function set.dimValues(obj,val)
      if iscell(val)
        for idxDim = 1:numel(val)
          obj.dimensions_(idxDim).dimValues = val{idxDim};
        end
        return;
      end
      
      obj.dimensions_.dimValues = val;      
    end
    
    function out = get.dimValues(obj)
      out = {obj.dimensions_.dimValues};      
    end;
    
    %%  arrayRange 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Get/Set Methods for obj.arrayRange
    function rangeOut = get.arrayRange(obj)
      if isempty(obj.arrayRange_)
        % Update the array range if it hasn't been assigned yet.
        obj.updateArrayRange(true);
      end;
      rangeOut = obj.arrayRange_;           
    end;                  
    function set.arrayRange(obj,~)
      error('obj.arrayRange is derived from obj.array');
    end;
        
    function updateArrayRange(obj,force,restrictionList)
      % Update the stored array range
      %
      % Inputs
      % ------
      %              obj : labelledArray object
      %            force : Force computation of the range
      %  restrictionList : Restricted list of channels to consider
      %
      % Assigns a value of [min(obj.array_(:)) max(obj.array_(:))] to
      % obj.arrayRange_
      %
      % This should only get called when the range is actually requested.
      % Calling when otherwise unnecessary can increased computation time
      % when dealing with large arrays through the min() and max() calls.
      %
      % restrictionList is included here for support in subclasses.
      %
      
      if ~exist('force','var'),force = false; end;
      if isempty(obj.arrayRange_)&&~force  
        % Only update the range if we're actually requesting it, or it's
        % been set before.
        return;
      end;
      
      if ~exist('restrictionList','var'), restrictionList = {':'}; end;
      
      allIdx = obj.getNumericIndex(restrictionList{:});
      
      arrayMin = min(obj.array(allIdx{:}));
      while ~isscalar(arrayMin)
        arrayMin = min(arrayMin);
      end;
      
      arrayMax = max(obj.array(allIdx{:}));
      while ~isscalar(arrayMax)
        arrayMax = max(arrayMax);
      end;
      
      obj.arrayRange_ = [arrayMin arrayMax];      
    end    
    
  end
  
  %% Hidden Methods
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  methods (Hidden)
    
    function isValid = isMutuallyConsistent(obj,b,varargin)
      % Check consistency between timeFrequencyDecomposition objects
      %
      % Used to check if math operations can be applied between two
      % timefrequency decompositions
      %
      % Inputs
      % ------
      %       obj : labelledArray object
      %         b : labelledArray object
      %  checkBSX : Flag to allow dimension mismatch if it satisfies
      %               conditions for applying bsxfun.
      %
           
      isValid = isMutuallyConsistent(obj.dimensions_,b.dimensions_,varargin{:});      
    end
    
    function [dimsOut] = getConsistentDimensions(obj,b)
      % Get output dimensions for bsxfun
      %
      % Inputs
      % ------
      %   obj : labelledArray object
      %     b : labelledArray object
      %
      % Outputs
      % -------
      %  dimsOut : Array of arrayDim objects containing the consistent
      %               dimensions
      %
      
      dimsOut = getConsistentDimensions(obj.dimensions_,b.dimensions_);                
    end
    

    function idxOut = findDimensions(obj,dimRefs)
      %% Return the numeric index of the dimension, given a numeric or name based reference
      %
      % idxOut = findDimensions(obj,dimRef)
      %
      % Inputs
      % ------
      %		obj  : labelledArray object
      %	dimRef : Numeric vector of dimension indices, character string with
      %             a single name, or a cell array of numeric indices and
      %             character strings.
      %
      % Outputs
      % -------
      %  idxOut : Numeric index of the dimensions associated with
      %							the strings in dimNames
      %
      %
      
      idxOut = obj.dimensions_.findDimensions(dimRefs);      
    end
    
    
  end
  

  %% Protected Methods
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  methods (Access=protected)
    
    %%
    
    
    %%
    function idxOut = getNumericIndex(obj,varargin)
      %% Get numeric indexing into a single labelledArray object
      %
      % idxOut = getNumericIndex(obj,varargin)
      %
      % Inputs
      % ------
      %      obj : A labelledarray object
      % varargin : Cell array of indexes into each dimension.
      %              This must satisfy numel(varargin)==obj.ndims
      %              Provided indexing dimValues can be:
      %                       ':' : All values
      %                cellString : Reference by name
      %             numericVector : Reference by numeric index
      %
      % Outputs
      % -------
      %  idxOut : Cell array of numeric indices into each dimension.
      %
      % If the number of values in varargin is less than obj.ndims, the
      % remaining dimensions are referenced using ':'
      %
      %
      
      assert(numel(obj)==1,'Multiple objects passed. Not sure why we''re getting here');      
      idxOut = obj.dimensions_.getIndexIntoDimensions(varargin{:});            
    end
       
    %%
    function out = copyValuesFrom(obj,valObj)
      %% Individually copies dimValues from another object
      %
      % out = copyValuesFrom(obj,valObj)
      %
      % Inputs
      % ------
      %     obj : labelledArray object
      %  valObj : labelledArray object to copy dimValues from
      %
      % Output
      % ------
      %     out : New labelledArray object with dimValues copied from valObj
      %
      % This method is designed to allow copying of object dimValues from one
      % object to another without losing the class of the original passed
      % object. This is primarily used in the constructor when constructing
      % subclasses.
      %%%%
      
      assert(isa(valObj,'labelledArray'),...
        'Can only copy from a labelledArray object');
      
      out = obj.copy;
      out.array_      = valObj.array_;
      out.dimensions_ = valObj.dimensions_;
    end
    
    %%
    function [out,varargout] = subcopy(obj,varargin)
      %% Copy a subset of the object
      %
      % out = subcopy(obj,varargin);
      %
      % Inputs
      % ------
      %       obj : A labelledArray object
      %  varargin : Numeric indexes into each dimension.
      %               NOTE: Must satisfy numel(varargin)==obj.ndims
      %
      % Outputs
      % -------
      %  out : labelledArray object with array, dimLabels, and dimValues
      %           subselected based on the provided indices.
      %%%%%%%%%%%
            
      % Input assertions
      assert(numel(obj)==1,'Subcopy must be called with a single labelledArray object');      
   
      % Get New Dimensions and Indices
      [newDims,dimIdx] = obj.dimensions_.subselectDimensions(varargin{:});
      
      % Copy the object
      out = obj.copy;
      out.array_ = obj.array(dimIdx{:});
      out.dimensions_ = newDims;
      
      if nargout==2
        % Pass out the index, if requested.
        varargout{1} = dimIdx;
      end;
      
    end
        
    function setDimensions(obj,val)      
      obj.dimensions_ = val;
      assert(isInternallyConsistent(obj.dimensions_));
    end
    
    %%
    function setArray(obj,val)
      % Overloadable set function for obj.array
      %
      % setArray(obj,val)
      %
      % This gets called during the actual set.array() method. By writing
      % this protected method instead, this can be overloaded by
      % subclasses.
      %
      %
      
      if ~isempty(obj.array_)&&~isempty(val)
        % Check Overall Dimensionality
       
        dimsIn = builtin('ndims',val);
        maxDims = max([dimsIn, obj.ndims]);
                
        % Check Individual Dimension Sizes
        for idxDim = 1:maxDims
          
          if idxDim<=dimsIn
          
           dimSize = size(val,idxDim);
          
           if ( dimSize ~= size(obj.array_,idxDim) )
             error(['Input array does not match current size on dimension: ' num2str(idxDim)]);
           end
          else
            assert(obj.size(idxDim)==1,...
               ['Input array does not match current size on dimension: ' num2str(idxDim)]);
          end
        end
      end
      
      %resetDimensions = isempty(val)||isempty(obj.array_)||isempty(obj.dimensions_);
      
      if ~isequal(obj.array_,val)
        % Set Array, and Clear Range (Will be computed the next time its
        % requested)
        obj.array_ = val;
        obj.arrayRange_ = []; % Clear Array Range
        
        %% Adjust dimension sizes
        for i = 1:ndims(val)
          if ~isempty(val)
            if i<=numel(obj.dimensions_)
              obj.dimensions(i).dimSize = size(val,i);
            else
              % Add a new dimension
              if isempty(obj.dimensions_)
                obj.dimensions = arrayDim('dimSize',size(val,i));
              else
                obj.dimensions(i) = arrayDim('dimSize',size(val,i));
              end;
            end
          else
          end
        end;
        notify(obj,'updatedOut');
      end
    end
    
  end
  
  
  
end
