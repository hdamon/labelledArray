classdef spatialGrid < handle & matlab.mixin.Copyable
  % classdef spatialGrid
  %
  % Group several arrayDims together to define a spatial grid in nD space
  %
  % obj = spatialGrid(dimensions,varargin)
  %
  % Inputs
  % ------
  %   dimensions : arrayDim object defining the dimensions
  %                   nDims = numel(dimensions)
  %
  % Param-Value Inputs
  % ------------------
  %        'origin' : (1 x nDims) Array Defining the Origin
  %    'directions' : (nDims x nDims) Array defining basis functions.
  %   'orientation' : Character string defining anatomic orientation of 
  %                     the grid. 
  %                     DEFAULT: 'Left-Posterior-Superior'
  %     'centering' : Defines whether grid points like in the center of
  %                     voxels ('cell'), or at the corners ('node')
  %                     DEFAULT: 'cell'
  %
  % Properties
  % ----------
  %   origin
  %  dimensions
  %  directions
  %  orientation
  %  centering
  %  voxSize
  %  boundingBox
  %  center
  %  aspect
  %
  % Methods
  % -------
  %
  
  properties (Dependent = true)
    origin
    dimensions
    directions
    orientation
    centering
    voxSize
    boundingBox
    center
    aspect
  end
  
  properties (Access=protected)
    % Defaults set in inputParser
    origin_
    directions_
    dimensions_
    orientation_
    centering_
  end
  
  
  methods
    
    function obj = spatialGrid(dimensions,varargin)
      % Main Object Constructor
      %
      if nargin>0
        %% Typecasting for subclassed objects
        if isa(dimensions,'spatialGrid')
          % Copy all object property values to the new object.
          obj.origin_      = dimensions.origin;
          obj.directions_  = dimensions.directions;
          obj.dimensions_  = dimensions.dimensions;
          obj.orientation_ = dimensions.orientation;
          obj.centering_   = dimensions.centering;
          return;
        end
        
        %% Input Parsing
        p = inputParser;
        p.KeepUnmatched = true;
        p.addRequired('dimensions',@(x) isa(x,'arrayDim'));
        p.addParameter('origin',zeros(1,numel(dimensions)),@(x) isnumeric(x)&&(numel(x)==numel(dimensions)));
        p.addParameter('directions',eye(numel(dimensions)));
        p.addParameter('orientation','Unknown');
        p.addParameter('centering','cell');
        p.parse(dimensions,varargin{:});
        
        %% Assign Property Values
        obj.dimensions = p.Results.dimensions;
        obj.origin = p.Results.origin;
        obj.directions = p.Results.directions;
        obj.orientation = p.Results.orientation;
        obj.centering = p.Results.centering;
      end
    end
    
    %% Set/Get Dimensions
    %%%%%%%%%%%%%%%%%%%%%
    function val = get.dimensions(obj)
      val = obj.dimensions_;
    end
    
    function set.dimensions(obj,val)
      if isnumeric(val)&&isvector(val)
        val = arrayDim('dimSize',num2cell(val));
      end
      
      assert(isa(val,'arrayDim'),'Must specify dimensions using an arraydim');
      
      %% Make sure values are populated
      for i = 1:numel(val)
        if isempty(val(i).dimValues)
          val(i).dimValues = 1:val(i).dimSize;
        end
      end
      
      obj.dimensions_ = val;
    end
    
    %% Set/Get Origin
    %%%%%%%%%%%%%%%%%
    function val = get.origin(obj)
      val = obj.origin_;
    end
    
    function set.origin(obj,val)
      assert(isnumeric(val)&&isvector(val)&&...
                (numel(val)==numel(obj.dimensions)),...
                'Origin size must match number of dimensions');      
      obj.origin_ = val;
    end
    
    function set.origin_(obj,val)      
      obj.origin_ = val(:)';
    end
    
    %% Set/Get Directions
    %%%%%%%%%%%%%%%%%%%%%
    function val = get.directions(obj)
      val = obj.directions_;
    end
    
    function set.directions(obj,val)
      isNumMat = isnumeric(val)&&ismatrix(val);
      isSizeCorrect = (size(val,1)==numel(obj.dimensions_))&&...
                      (size(val,2)==numel(obj.dimensions_));      

      assert(isNumMat&&isSizeCorrect,...
              'Mismatch in size of directions matrix');
                    
      obj.directions_ = val;
    end
    
    function set.directions_(obj,val)            
      obj.directions_ = val;
    end        
    
    %% Set/Get Orientation
    %%%%%%%%%%%%%%%%%%%%%%
    function val = get.orientation(obj)
      val = obj.orientation_;
    end
    
    function set.orientation(obj,val)
      obj.orientation_ = val;
    end
          
    
    %% Set/Get Centering
    %%%%%%%%%%%%%%%%%%%%
    function val = get.centering(obj)
      val = obj.centering_;
    end
    
    function set.centering(obj,val)
      obj.centering_ = val;
    end
    
    %% Get Aspect
    %%%%%%%%%%%%%
    function out = get.aspect(obj)
      out = sqrt(sum(obj.directions.^2,1));
    end    
    
    %% Get Grid Points
    %%%%%%%%%%%%%%%%%%    
    function pts = getGridPts(grid)
      % Get location in nD space of each point in the grid
      %
      % pts = getGridPts(grid);
      %
      % pts = obj.getGridPts;
      %
      % Inputs
      % ------
      %   grid : spatialGrid object
      %
      % Outputs
      % -------
      %  pts : Npts x 3 array of Euclidean coordinates for each point in
      %          the grid. Point ordering indexes across the contents of
      %          obj.dimensions in order
      %
      
      allVals = {grid.dimensions.dimValues};
      
      XYZ = cell(1,numel(allVals));
      [XYZ{:}] = ndgrid(allVals{:});
      
      pts = zeros(numel(XYZ{1}),numel(XYZ));
      for i = 1:numel(XYZ)
        pts(:,i) = XYZ{i}(:);
      end
      
      pts = pts*grid.directions' + grid.origin;
      
    end
    
    function isUniform = isUniformlySampled(grid)
      error('isUniformlySampled method not fully implemented');
      nDim = numel(grid.dimensions);        
    end
    
    %% Find Nearest Nodes
    %%%%%%%%%%%%%%%%%%%%%%
    function [idxOut] = getNearestNodes(grid,Positions,useNodes)
      % Find index of grid nodes nearest to specific locations
      %  
      % function [idxOut] = getNearestNodes(grid,Positions,useNodes)
      %
      % Inputs
      % ------
      %      grid : spatialGrid object
      % Positions : Nx3 Array of Positions to Project
      %  useNodes : Index array to restrict the nearest neighbor search
      %               to a limited set of voxels.
      %
      % Output
      % ------
      %  idxOut : Nx1 Array of indices into spatialGrid of the points
      %             closest to the input positions.
      %
      
      if ~exist('useNodes','var'),useNodes = ':'; end
      
      pts = grid.getGridPts;
      pts = pts(useNodes,:);
      
      % For each position to be shifted, find the closest node.
      idxOut = zeros(size(Positions,1),1);
      for i = 1:size(Positions,1)
        dist = Pts - repmat(Positions(i,:),size(Pts,1),1);
        dist = sqrt(sum(dist.^2,2));
        
        q = find(dist==min(dist));
        if numel(q)>1, % Pick a random node if there's more than one
          q = q(ceil(numel(q)*rand(1,1)));
        end
        idxOut(i) = UseNodes(q);
      end
      
    end
    
    %% Find Grid Bounding Box
    %%%%%%%%%%%%%%%%%%%%%%%%%
    function boxOut = getBoundingBox(grid)
      % Find the bounding box of a spatialGrid
      %
      % boxOut = getBoundingBox(grid)
      %
      % Inputs
      % ------
      %  grid : spatialGrid object
      %
      % Output
      % ------
      %  boxOut : 8x3 Array of bounding box corner coordinates
      %             Dimension ordering is the same as 
      %
      
      nDim = numel(grid.dimensions);
      boxOut = zeros(2^nDim,nDim);
      
      for idxDim = 1:nDim
        colIdx = getColIdx(nDim,idxDim);
        boxOut(:,idxDim) = grid.dimensions(idxDim).dimRange(colIdx);
      end
      
      boxOut = boxOut*grid.directions';
      boxOut = boxOut + grid.origin;
            
      if strcmpi(grid.centering,'cell')
        % There's a big assumption here that the grid is uniformly spaced.
        % Not sure how to deal with it otherwise.
        boxOut = boxOut - (grid.directions*[0.5 0.5 0.5]')';
      else        
      end
      
      function newCol = getColIdx(totDim,depth)        
        tmp = ones(2^(totDim-depth),1);
        tmp = [tmp ; 2*tmp];
        
        newCol = repmat(tmp,2^(depth-1),1);                
      end            
    end
    
  end
  
end

  
  
    
    