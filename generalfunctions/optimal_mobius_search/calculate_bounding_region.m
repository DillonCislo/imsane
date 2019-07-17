function [ rMin, rMax, thetaMin, thetaMax ] = calculate_bounding_region( cR, rR, Mj )
%CALCULATE_BOUNDING_REGION Calculates the bounding region
%   Input Parameters:
%       cR:     The center of the circumcircle of the input square region.
%               This is expected to be a complex variable.
%
%       rR:     The radius of the circular bounding region
%
%       Mj:     The input point being mapped

%==========================================================================
% Calculate the bounding radii
%==========================================================================
gamma = cR - rR^2 / conj( -1/conj(Mj) + cR );
cNR = ( gamma - Mj ) / ( 1 - conj(Mj) * gamma );
rNR = abs( cNR - ( cR + rR - Mj ) / ( 1 - conj(Mj) * (cR+rR) ) );

% Cap the radii should the exceed the appropriate limits
rMin = max(0, abs(cNR)-rNR);
rMax = min(1, abs(cNR)+rNR);

%==========================================================================
% Calculate the bounding angles
%==========================================================================
cAR = -cR + Mj; rAR = rR;
cBR = 1-conj(Mj)*cR; rBR = abs(Mj)*rR;

dAR = abs(cAR); angAR = angle(cAR);
dBR = abs(cBR); angBR = angle(cBR);

%--------------------------------------------------------------------------
% Check to see if either AR or BR contains the origin
%--------------------------------------------------------------------------
if (rAR > dAR) || (rBR > dBR)
    
    thetaMin = -pi;
    thetaMax = pi;
   
%--------------------------------------------------------------------------
% Otherwise calculate the angular shifts
%--------------------------------------------------------------------------
else
    
    dangAR = asin(rAR / dAR);
    dangBR = asin(rBR / dBR);
    
    
    thetaMin = wrapTo2Pi((angAR - dangAR)) + wrapTo2Pi((angBR - dangBR));
    thetaMax = wrapTo2Pi((angAR + dangAR)) + wrapTo2Pi((angBR + dangBR));
    
    thetaMin = wrapToPi(thetaMin);
    thetaMax = wrapToPi(thetaMax);

    
end
   
end

