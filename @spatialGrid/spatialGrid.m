classdef spatialGrid < handle & matlab.mixin.Copyable
  % classdef spatialGrid
  %
  % Group several arrayDims together to define a spatial grid in nD space
  %
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
    origin_
    directions_
    dimensions_
    orientation_ = 'Left-Posterior-Superior';
    centering_ = 'cell';
  end
  
  
  methods
    
    function obj = spatialGrid(dimensions,varargin)
      
      if nargin>0
      if isa(dimensions,'spatialGrid')
        % Output copy of the input object
        obj = dimensions.copy;
        return;
      end
      
      obj.dimensions = dimensions;
      
      p = inputParser;   
      p.KeepUnmatched = true;
      p.addParameter('origin',[],@(x) isnumeric(x)&&(numel(x)==numel(dimensions)));
      p.addParameter('directions',[]);
      p.addParameter('orientation','Left-Posterior-Superior');
      p.addParameter('centering','cell');      
      p.parse(varargin{:});
                  
      if isempty(p.Results.origin)
        obj.origin_ = zeros(1,numel(obj.dimensions));
      else        
        obj.origin = p.Results.origin;
      end
      
      if isempty(p.Results.directions)
        obj.directions_ = eye(numel(obj.dimensions));
      else
        obj.directions = p.Results.directions;
      end
      
      obj.orientation = p.Results.orientation;
      obj.centering = p.Results.centering;
      end
    end
    
    %% Set/Get Dimensions
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
      nDim = numel(grid.dimensions);        
    end
    
    %% Find Nearest Nodes
    %%%%%%%%%%%%%%%%%%%%%%
    function [idxOut] = getNearestNodes(grid,Positions,useNodes)
      % Find index of grid nodes nearest to specific locations
      %  
      % function [idxOut] = getNearestNodes(grid,Positions,useNodes)
      %
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

  
  
    
    