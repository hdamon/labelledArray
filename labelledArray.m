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
  %  'dimNames'  : dimNames for each dimension.
  %  'dimLabels' : An cell array defining dimLabels for one or more dimensions.
  %               This array must be of size (N X 2), arranged as:
  %                { <DIMENSION A> , <CellString of dimLabels> ;
  %                  <DIMENSION B> , <CellString of dimLabels> }
  %               Each dimension must be a scalar numeric value, with no
  %               repeats. Each cell string must have a number of elements
  %               equal to the current size of the array along that
  %               dimension.
  %  'dimValues' : A cell array of dimValues along one or more dimensions. Cell
  %               array should be formatted as for 'dimLabels', except instead
  %               of the dimValues beign cell strings, they must be numeric
  %               vectors.
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
  % -------------------------------------
  %   Referencing for labelledarray objects behaves slightly differently
  %   than for normal Matlab arrays.
  %
  %  obj.<property> : Behaves normally.
  %
  %  obj(<indices>) : Behaves in one of two ways:
  %                     1) If obj is an array of labelledarray objects, this
  %                          references into the array
  %                     2) If obj is a single labelledarray object, this
  %                          returns a new object with the array, dimValues,
  %                          and dimLabels subselected according to the
  %                          indices.
  %
  %  obj{<indices>} : References into an array of labelledarray objects.
  %                     This works correctly in combination with (), so
  %                     that expressions such as:
  %                       obj{i}(x,y,z)       and
  %                       obj{i}.array(x,y,z)
  %                     work as expected.
  %
  
  properties (Hidden,Dependent)
    array
    dimNames
    dimLabels
    dimValues
    dimUnits
  end
  
  properties (Access=protected) % Should possibly be private?
    array_   % The array array
    dimNames_  % dimNames for each dimension
    dimLabels_ % dimLabels for individual elements of each axis
    dimValues_ % dimValues for each axis (time,frequency,etc)
    dimUnits_ % Units for each axis
  end
  
  methods
    
    %%
    function obj = labelledArray(array,varargin)
      
      if nargin>0
        
        checkType = @(x) iscell(x) && ((size(x,2)==2)||(numel(x)==ndims(array)));
        
        % Input parsing
        p = inputParser;
        p.addRequired('array',@(x) (isnumeric(x)||isa(x,'labelledArray')));
        p.addParameter('dimNames',[],@(x) isempty(x)||checkType(x));
        p.addParameter('dimLabels',[],@(x) isempty(x)||checkType(x));
        p.addParameter('dimUnits',[],@(x) isempty(x)||checkType(x));
        p.addParameter('dimValues',[],@(x) isempty(x)||checkType(x));        
        p.parse(array,varargin{:});
        
        % Property Assignment
        obj.array   = p.Results.array;
        obj.dimNames  = p.Results.dimNames;
        obj.dimLabels = p.Results.dimLabels;
        obj.dimUnits  = p.Results.dimUnits;
        obj.dimValues = p.Results.dimValues;
      end
      
    end
    
    %% Overloaded Functions
    %%%%%%%%%%%%%%%%%%%%%%%
    
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
      if numel(obj)==1
        if ~exist('dim','var')
          out = size(obj.array);
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
      % Ignores trailing singleton dimensions.
      %
      
      if numel(obj)==1
%         if isempty(obj.array_)
%           out = 0;
%           return;
%         end
        out = builtin('ndims',obj.array_);
      else
        out = builtin('ndims',obj);
      end;
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
      out.array_   = permute(obj.array_,newOrder);
      out.dimNames_  = obj.dimNames_(newOrder);
      out.dimLabels_ = obj.dimLabels_(newOrder); 
      out.dimUnits_  = obj.dimUnits_(newOrder);
      out.dimValues_ = obj.dimValues_(newOrder);
    end
    
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
    
    function out = cat(dim,obj,a,varargin)
      
      assert(isa(a,class(obj)),'Can only concatenate like objects');
      
      % Concatenating with empty objects is quick.
      if isempty(obj)        
        out = a;
        return;
      end;
      
      if isempty(a)
        out = obj;
        return;
      end;
      
      newLabels = [];
      for idxDim = 1:obj.ndims
        assert(isequal(obj.dimNames{idxDim},a.dimNames{idxDim}),...
                  'Dimension names are inconsistent');
        if ( idxDim==dim )
          newLabels = cat(1,obj.dimLabels{idxDim},a.dimLabels{idxDim});
        else
          assert(isequal(obj.dimLabels{idxDim},a.dimLabels{idxDim}),...
                  'Dimension labels are inconsistent');
        end
      end            
      
      out = obj.copy;
      out.array_ = cat(dim,obj.array_,a.array_);
      out.dimLabels{dim} = newLabels;
                      
      if ~isempty(varargin)
        % Recurse when concatenating multiple objects
        out = cat(dim,out,varargin{:});
      end;
      
    end
    
    function out = mean(obj,dim)
      if ~exist('dim','var'), dim = 1; end;      
      out = applyDimFunc(@var,obj,1,dim);
    end
    
    function out = std(obj,W,dim)
      if ~exist('dim','var'), dim = 1; end;
      if ~exist('W','var'), W = 0; end;
      out = applyDimFunc(@var,obj,2,W,dim);      
    end
    
    function out = var(obj,W,dim)
      if ~exist('dim','var'), dim = 1; end;
      if ~exist('W','var'), W = 0; end;      
      out = applyDimFunc(@var,obj,2,W,dim);
    end;
    
    function out = applyDimFunc(funcHandle,obj,idxDim,varargin)      
      % For functions that collapse a dimension to a singleton
      %
      % Inputs
      % ------
      %   funcHandle : Function handle
      %                  Will be called as:
      %                  funcHandle(obj.array,varargin{:});
      %         obj
      %       idxDim : Index into varargin of the dimension argument
      %     varargin : Additional arguments for funcHandle. One of these
      %                 must be a numeric dimension index.
      %
      
      dim = varargin{idxDim};
      out = obj.copy;
      out.array_ = funcHandle(obj.array,varargin{:});
      out.dimLabels{dim} = func2str(funcHandle);
      out.dimUnits{dim} = [func2str(funcHandle) '(' out.dimUnits{dim} ')'];
      if ~isempty(obj.dimLabels{dim})
        out.dimLabels{dim} = [ out.dimLabels{dim} '(' strjoin(obj.dimLabels{dim}) ')'];            
      end;     
      if ~isempty(obj.dimValues{dim})
        % Return the average of the values.
        out.dimValues{dim} = mean(out.dimValues{dim});
      end
      
    end
    
    
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
          assert(isConsistent(obj(idxObj),b),'Inconsistent decomposition sizes');
          coeff = b.array;
          [names,labels] = getConsistentDimensions(obj(idxObj),b);
        else
          coeff = b;
          names = obj.dimNames;
          labels = obj.dimLabels;          
        end;
        
        objOut(idxObj) = copy(obj(idxObj));
        objOut(idxObj).array_ = bsxfun(funcHandle,obj.array,coeff);
        objOut(idxObj).dimNames = names;
        objOut(idxObj).dimLabels = labels;
        
      end;
      
      %% Reshape if its an array of objects
      if numel(objOut)>1
        objOut = reshape(objOut,size(obj));
      end;
    end      
    
    
    %% array Get/Set Methods
    function set.array(obj,val)
      % Done this way so it can be overloaded by subclasses.
      obj.setArray(val);
    end;    
    
    function out = get.array(obj)
      out = obj.array_;
    end;
    
    %% dimNames Get/Set Methods
    %%%%%%%%%%%%%%%%%%%%%%%
    function set.dimNames(obj,val)
      % Set Dimension Names
      %      
      typeCheckFunction = @(x) iscellstr(x);
      typeCheckErr = 'dimNames must be strings or cellstrings';
      
      obj.setProperty('dimNames',val,true,true,typeCheckFunction,typeCheckErr);      
    end
        
    function out = get.dimNames(obj)
      out = obj.dimNames_;
      [isChanged,out] = obj.fixDimLength(out);
      if isChanged
        obj.dimNames_ = out;
      end;
    end;
    
    %% dimLabels Get/Set Methods
    %%%%%%%%%%%%%%%%%%%%%%%%
    function set.dimLabels(obj,val)     
      typeCheckFunc = @(x) iscellstr(x);
      typeCheckErr = 'dimLabels must be strings or cellstrings';
      
      obj.setProperty('dimLabels',val,true,false,typeCheckFunc,typeCheckErr);        
    end
    
    %%
    function out = get.dimLabels(obj)
      out = obj.dimLabels_;
      [isChanged,out] = obj.fixDimLength(out);
      if isChanged
        obj.dimLabels_ = out;
      end;      
    end;
    
    %% dimUnits Get/Set Methods
    %%%%%%%%%%%%%%%%%%%%%%%%
    function set.dimUnits(obj,val)      
      typeCheckFunc = @(x) iscellstr(x);
      typeCheckErr = 'dimUnits must be strings or cellstrings';
      
      obj.setProperty('dimUnits',val,true,true,typeCheckFunc,typeCheckErr);      
    end
    
    function out = get.dimUnits(obj)
      out = obj.dimUnits_;
      
      [isChanged,out] = obj.fixDimLength(out);
      if isChanged
        obj.dimUnits_ = out;
      end;            
    end;    
   
    %% dimValue Get/Set Methods
    %%%%%%%%%%%%%%%%%%%%%%%%
    function set.dimValues(obj,val)      
      typeCheckFunc = @(x) isnumeric(x)&&isvector(x);
      typeCheckErr = 'Value inputs must be numeric vectors';
      
      obj.setProperty('dimValues',val,false,true,typeCheckFunc,typeCheckErr);   
    end
    
    function out = get.dimValues(obj)
      out = obj.dimValues_;
      [isChanged,out] = obj.fixDimLength(out);
      if isChanged
        obj.dimValues_ = out;
      end;
    end;    
    
  end
  
  %% Hidden Methods
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  methods (Hidden)
    
    function isValid = isConsistent(obj,b,checkBSX)
      % Check consistency between timeFrequencyDecomposition objects
      %
      % Used to check if math operations can be applied between two
      % timefrequency decompositions
      %
      % Inputs
      % ------
      %       obj : labelledArray object
      %         b : labelledArray object
      %  checkBSX : Flag to check compatibility 
      
      if ~exist('checkBSX','var'), checkBSX = true; end;
      
      isValid = false;
      
      if ~isa(b,'labelledArray'), return; end;
     % if ~(obj.ndims==b.ndims),  return; end; % Might not be needed?
      
      if checkBSX
        sizeValid = obj.validForBSXFUN(size(obj),size(b));
        if ~sizeValid, return; end;
      else
        if ~(obj.ndims==b.ndims)
          return;
        end;        
      end;      
      
      for idxDim = 1:obj.ndims
         % Dimension names must all be the same
         if ~(isequal(obj.dimNames{idxDim},b.dimNames{idxDim})), return; end;
         
         % Check Dimension Labels
         if checkBSX           
           if (size(obj,idxDim)~=1)&&(b.ndims<=idxDim)&&(size(b,idxDim)~=1)
             % Non-unitary dimensions need to be the same.
             if ~(isequal(obj.dimLabels{idxDim},b.dimLabels{idxDim})), 
               return; 
             end;
           end;
         else
           if ~(isequal(obj.dimLabels{idxDim},b.dimLabels{idxDim}))
             return;
           end;
         end
      end
                 
      % Passed all the tests, so its valid.
      isValid = true;      
    end      

   function [varargout] = getConsistentDimensions(obj,b)
     % Get output dimensions for bsxfun
     %
     % Inputs
     % ------
     %   obj : labelledArray object
     %     b : labelledArray object
     %
     % Outputs
     % -------
     %   names : dimNames for the 
     %  labels : dimLabels for the
     %
     
      assert(isa(b,'labelledArray'),...
                'Second input must be a labelledArray object');
              
      assert(isConsistent(obj,b),'Inconsistent decomposition sizes');  
      
      sizeEqual = crlEEG.util.validation.compareSizes(size(obj),size(b));
      
      nDimsOut = numel(sizeEqual);
      
      names  = cell(nDimsOut,1); 
      labels = cell(nDimsOut,1); 
      
      for idxDim = 1:nDimsOut        
        if sizeEqual(idxDim)||(size(b,idxDim)==1)
          names{idxDim}  = obj.dimNames{idxDim}; 
          labels{idxDim} = obj.dimLabels{idxDim};
          units{idxDim}  = obj.dimUnits{idxDim};
          values{idxDim} = obj.dimValues{idxDim};
        elseif size(obj,idxDim)==1
          names{idxDim}  = b.dimNames{idxDim};  
          labels{idxDim} = b.dimLabels{idxDim};
          units{idxDim}  = b.dimUnits{idxDim};
          values{idxDim} = b.dimValues{idxDim};          
        else
          error('Shouldn''t be getting here');
        end
      end   
      
      if nargout>=1, varargout{1} = names; end;
      if nargout>=2, varargout{2} = labels; end;
      if nargout>=3, varargout{3} = units; end;
      if nargout>=4, varargout{4} = labels; end;
      
   end      
    
   
  end  
  
  %% Hidden Static Methods
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  methods (Hidden,Static)
    function isValid = validForBSXFUN(sizeA,sizeB)
      % Check that matrix sizes are compatible for bsxfun operation
      
      maxlen = max([numel(sizeA) numel(sizeB)]);
      
      tmpA = ones(1,maxlen);
      tmpB = ones(1,maxlen);
      tmpA(1:numel(sizeA)) = sizeA;
      tmpB(1:numel(sizeB)) = sizeB;
      
      test = ( tmpA==tmpB );
      isValid = all((tmpA(~test)==1)|(tmpB(~test)==1));
      
    end

  end
  
  %% Protected Methods
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  methods (Access=protected)

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
        obj.([propName '_']){val{i,1}} = val{i,2};                
      end
      
    end    
    
    
    %%
    function newOrder = getDimensionOrdering(obj,order)
      %% Get the numeric indices of the new ordering.
      %
      % newOrder = getDimensionOrdering(obj,order)
      %
      % Inputs
      % ------
      %    obj : labelledArray object
      %  order : New dimension ordering for the object
      %           Can be:
      %             1) A numeric vector indexing the dimensions
      %             2) A cell array, with either numeric dimension indices
      %                 or character strings referencing index names.
      %
      % Outputs
      % -------
      %  newOrder : New order, with numeric dimension referencing
      %
      newOrder = nan(1,numel(order));
      if isnumeric(order)&&isvector(order)
        newOrder = order;
      elseif iscell(order)
        for i = 1:numel(order)
          if isnumeric(order{i})&&isscalar(order{i})
            newOrder(i) = order{i};
          elseif ischar(order{i})
            newOrder(i) = obj.getDimByName(order{i});
          else
            error('Invalid permutation argument');
          end;
        end
      else
        error('Invalid permutation argument');
      end
      
      assert(numel(unique(newOrder))==numel(newOrder),...
                'Each dimension can only be used once');
      
    end;
    
    %%
    function idxOut = getDimensionIndex(obj,dimRefs)
      %% Return the numeric index of the dimension, given a numeric or name based reference
      %
      % idxOut = getDimByName(obj,dimRef)
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
      
      assert((isnumeric(dimRefs)&&isvector(dimRefs))||...
              ischar(dimRefs)||iscell(dimRefs),...
        'See getDimensionIndex help for input format requirements');
      
      if ischar(dimRefs), dimRefs = {dimRefs}; end;
      
      if isnumeric(dimRefs)
        % Convert to cell for consistency
        dimRefs = num2cell(dimRefs);
      end
      
      % Get Valid Character Names
      validNames = obj.dimNames;
      validNames(cellfun(@isempty,validNames)) = [];
      
      % Get numeric indices from a cell array
      idxOut = nan(1,numel(dimRefs));
       
      for idxRef = 1:numel(dimRefs)
         currRef = dimRefs{idxRef};
         if ischar(currRef)     
           % Character Referencing                                   
           if ~isempty(validNames)
             matchedName = validatestring(currRef,validNames);
             idx{idxRef} = find(ismember(obj.dimNames,matchedName));
           end;

         elseif isnumeric(currRef)&&isscalar(currRef)
           % Numeric Indexing
           idxOut(idxRef) = currRef;
         else
           error('Unknown reference type');
         end           
       end
      
      
      assert(all(floor(idxOut)==idxOut)&&all(idxOut>=1)&&all(idxOut<=obj.ndims),...
                'Requested dimension is out of range');
      
    end
    
    
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
      %                       ':' : All dimValues
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
     
      % Get Indexing for each dimension
      idxOut = cell(obj.ndims,1);
      for idxDim = 1:obj.ndims
        if idxDim<=numel(varargin)          
          idxOut{idxDim} = obj.indexIntoDimension(idxDim,varargin{idxDim});
        else
          idxOut{idxDim} = ':';
        end;
      end
    end
    
    %%
    function idxOut = indexIntoDimension(obj,dim,index)
      % Indexing into a single dimension based on name or number
      %
      % idxOut = indexIntoDimension(obj,dim,index)
      %
      % Inputs
      % ------
      %   obj : labelledArray object
      %   dim : Selected dimension
      %  index: Index into the selected dimension
      %           This can be numeric (scalar or vector), logical, or character
      %           (string or cellstr). When provided with string/cellstr,
      %           the values in obj.dimLabels{dim} are used to index.
      %           Providing a string that does not match anything in the
      %           object will produce an error.      
      %
      % Outputs
      % -------
      %  idxOut : Numeric index into the selected dimension
      %
      
      if isempty(obj.dimLabels_{dim})
        cellIn = [];
        isStringValid = false;
      else
        cellIn = obj.dimLabels_{dim};
        isStringValid = true;
      end;
      
      if isequal(index,':')
        %% Requested Everything
        idxOut = index;
        return;
        
      elseif islogical(index)
        %% Logical Indexing
        assert(numel(index)==size(obj,dim),'FOOOO_-');
        idxOut = find(index);
        return;
        
      elseif isnumeric(index)
        %% Numeric Reference
        if any(index<1)||any(index>size(obj,dim))
          error('Requested index outside of available range');
        end;
        idxOut = index;
        %idxOut(idxOut<1) = nan;
        %idxOut(idxOut>numel(cellIn)) = nan;
        return;
                
      elseif ischar(index)||iscellstr(index)
        %% String Reference
        if ~isStringValid
          error('String indexing unavailable for this dimension');
        end;
        
        if ischar(index), index = {index}; end;
        cellIn = strtrim(cellIn);
        index = strtrim(index);
        idxOut = zeros(1,numel(index));
        for idx = 1:numel(idxOut)
          tmp = find(strcmp(index{idx},cellIn));
          if isempty(tmp)
            error('Requested string does not appear in cell array');
          end;
          assert(numel(tmp)==1,'Multiple string matches in cellIn');
          idxOut(idx) = tmp;
        end
        
      else
        %% Otherwise, error.
        error('Incorrect reference type');
      end;
      
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
      out.array_     = valObj.array_;
      out.dimNames_  = valObj.dimNames_;
      out.dimLabels_ = valObj.dimLabels_;
      out.dimValues_ = valObj.dimValues_;
      out.dimUnits_  = valObj.dimUnits_;
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
      
      % Get Numeric Indexing
      dimIdx = obj.getNumericIndex(varargin{:});
      
      % Input assertions
      assert(numel(obj)==1,'Passed multiple objects. Not sure why we''re getting here');      
      assert(numel(varargin)==obj.ndims,...
        'Must include referencing for all dimensions');
      
      % Allocate Memory
      tmparray = obj.array(varargin{:});
      tmpdimLabels = cell(obj.ndims,1);      
      tmpdimValues = cell(obj.ndims,1);
      tmpdimUnits  = cell(obj.ndims,1);
      
      % Get subsets
      for idxDim = 1:obj.ndims        
        if ~isempty(obj.dimLabels{idxDim})
          tmpdimLabels{idxDim} = obj.dimLabels{idxDim}(dimIdx{idxDim});
        end        

        currValues = obj.dimValues{idxDim};
        if ~isempty(currValues)
          tmpdimValues{idxDim} = currValues(dimIdx{idxDim});
        end;
        
        currUnits = obj.dimUnits{idxDim};
        if ~isempty(currUnits)
          if ischar(currUnits)
            tmpdimUnits{idxDim} = currUnits;
          else
            tmpdimUnits{idxDim} = currUnits(dimIdx{idxDim});
          end
        end;
        
      end;
      
      out = obj.copy;
      out.array_ = tmparray;
      out.dimLabels_ = tmpdimLabels;
      out.dimUnits_  = tmpdimUnits;
      out.dimValues_ = tmpdimValues;      
      
      if nargout==2
        % Pass out the index, if requested.
        varargout{1} = dimIdx;
      end;
                 
    end
    
    %%
    function validCell = validateCellInput(obj,val)
      %% Check cell input arrays, and convert style if necessary.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      assert(iscell(val)&&((size(val,2)==2)||...
        (numel(val)==obj.ndims)),'Incorrect input shape');
      
      % Check if input is already formatted correctly
      isValid = true;
      validCell = cell(size(val,1),2);
      
      nRows = size(val,1);
      nCols  = size(val,2);
      
      if ((nRows>1)&&(nCols==2))||...
          (isnumeric(val{1,1})&&(nCols==2))
        for i = 1:size(val,1)
          
          newDim = obj.getDimensionIndex(val{i,1});
          newVal = val{i,2};
          
          if isnan(newDim) %%||(~(size(obj,newDim)==numel(newVal)))
            isValid = false;
            break;
          end
          validCell{i,1} = newDim;
          validCell{i,2} = newVal;
        end
      else
        isValid = false;
      end;
      if isValid, return; end;
      
      % All dimensions are defined
      if numel(val)==obj.ndims
        validCell = cell(numel(val),2);
        for i = 1:obj.ndims
          validCell{i,1} = i;
          validCell{i,2} = val{i};
        end
        return;
      end
      
      error('Invalid input');
    end
    
    function [isChanged,newVal] = fixDimLength(obj,val)
      
      if numel(val)~=obj.ndims
        newVal = cell(obj.ndims,1);
        newVal(1:numel(val)) = val;
        isChanged=true;
      else
        newVal = val;
        isChanged = false;      
      end;        
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
      
      if ~isempty(obj.array_)
        % Check Overall Dimensionality
        nDims = builtin('ndims',val);
        if (nDims~=obj.ndims)&&(obj.ndims~=0)
          error('array dimensionality does not match existing size');
        end
        
        % Check Individual Dimension Sizes
        for idxDim = 1:nDims
          dimSize = size(val,idxDim);
          
          if ( dimSize ~= size(obj.array_,idxDim) )
            error(['Input array does not match current size on dimension: ' num2str(idxDim)]);
          end
        end
      end
                
      % Clear the labels and names if the array is getting cleared, or if
      % they haven't been set yet.
      resetLabels = isempty(val)||isempty(obj.array_)||isempty(obj.dimLabels_);
      resetNames  = isempty(val)||isempty(obj.array_)||isempty(obj.dimNames_);
      resetValues = isempty(val)||isempty(obj.array_)||isempty(obj.dimValues_);
      resetUnits  = isempty(val)||isempty(obj.array_)||isempty(obj.dimUnits_);
      
      obj.array_ = val;
      
      if resetLabels, obj.dimLabels = []; end;
      if resetNames,  obj.dimNames = []; end;             
      if resetValues, obj.dimValues = []; end;
      if resetUnits,  obj.dimUnits  = []; end;
    end    
    
  end
  
  
  
end
