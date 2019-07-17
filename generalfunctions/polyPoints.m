function [polyPnts] = polyPoints( fitter, varargin )
%POLYBDY Returns the  points inside a polygon or the
%sorted boundary points of a polygonal submesh
%   polyPnts(fitter, inPolyPnts, inBdyPnts)

    inPolyPnts = varargin{1};

    xBdy = fitter.fittedParam.submeshes{1,1}.u{1}(inPolyPnts,1);
    yBdy = fitter.fittedParam.submeshes{1,1}.u{1}(inPolyPnts,2);

    polyPnts = [xBdy, yBdy];

    if length(varargin) == 2
    
      inBdyPnts = varargin{2};
    
      xBdy = xBdy(inBdyPnts);
      yBdy = yBdy(inBdyPnts);
    
      polyPnts = [xBdy, yBdy];
      polyPnts = CounterClockWiseSort(polyPnts);
    
    end


end

