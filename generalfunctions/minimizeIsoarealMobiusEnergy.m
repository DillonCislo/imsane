function z0 = minimizeIsoarealMobiusEnergy( face, R3D, R2D )
%MINIMIZEISOAREALMOBIUSENERGY Determines the Mobius automorphism of the
%unit disk that matches as closely the (re-scaled) areas of the mesh
%triangles in the plane to their counterparts in 3D (up to a rotation)
%   Input Parameters:
%       - face               #Fx3 connectivity list
%       - R3D                #Vx3 3D vertex coordinate list
%       - R2D                #Vx2 2D vertex coordinate list.  Assumed to be
%                            a triangulation of the unit disk
%   Output Parameters:
%       - z0                 The Mobius parameter of the best fit map
%
%--------------------------------------------------------------------------
% INPUT PROCESSING
%--------------------------------------------------------------------------

% Calculate the area of the 3D mesh ---------------------------------------
e1_3D = R3D( face(:,3), : ) - R3D( face(:,2), : );
e2_3D = R3D( face(:,1), : ) - R3D( face(:,3), : );

A3D = cross( e1_3D, e2_3D, 2 );
A3D = sum( sqrt( sum( A3D.^2, 2 ) ) ./ 2 );

AA = sqrt( A3D ./ pi );

% Calculate the target edge lengths ---------------------------------------
tr3D = triangulation( face, R3D );
eIDx = tr3D.edges;

L = R3D( eIDx(:,2), : ) - R3D( eIDx(:,1), : );
L = sqrt( sum( L.^2, 2 ) );

% Create complex representation of disk mapping ---------------------------
z = complex( R2D(:,1), R2D(:,2) );

%--------------------------------------------------------------------------
% MINIMIZATION PROCESSING
%-------------------------------------------------------------------------

fun = @isoarealMobiusEnergy;
x0 = [0 0];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon = @unitdisk;

options = optimoptions( 'fmincon', ...
    'Algorithm', 'interior-point', ...
    'Display', 'none', ...
    'ConstraintTolerance', 1e-6, ...
    'SpecifyObjectiveGradient', true, ...
    'SpecifyConstraintGradient', true, ...
    'CheckGradients', false );

z0 = fmincon( fun, x0, A, b, Aeq, beq, lb, ub, nonlcon, options );
    
%--------------------------------------------------------------------------
% ENERGY FUNCTION
%--------------------------------------------------------------------------
    function [E, EG] = isoarealMobiusEnergy( x )
        
        %------------------------------------------------------------------
        % Calculate the energy
        %------------------------------------------------------------------
        
        % Calculate the new vertex positions
        a = complex( x(1), x(2) );
        zeta = AA .* ( z - a ) ./ ( 1 - conj(a) .* z );
        
        % Calculate the new edge vectors
        zetaEdge = zeta(eIDx(:,2)) - zeta(eIDx(:,1));
        e2D = [ real( zetaEdge ), imag( zetaEdge ) ];
        
        % Calculate the new edge lengths
        l2D = sqrt( sum( e2D.^2, 2 ) );
        
        E = sum( ( l2D - L ).^2 );
        
        %------------------------------------------------------------------
        % Calculate the energy gradient
        %------------------------------------------------------------------
        if ( nargout > 1 )
            
            % Calculate the edge unit vectors
            e2DHat = e2D ./ l2D;
            
            % Calculate the derivatives
            dZdB = AA .* ( z.^2 - 2i .* x(2) .* z - 1 ) ./ ...
                ( 1 - conj(a) .* z ).^2;
            dZdC = -1i .* AA .* ( z.^2 - 2 .* x(1) .* z + 1 ) ./ ...
                ( 1 - conj(a) .* z ).^2;
            
            dZdCIJ = dZdC( eIDx(:,2) ) - dZdC( eIDx(:,1) );
            dZdBIJ = dZdB( eIDx(:,2) ) - dZdB( eIDx(:,1) );
            
            gradl = zeros( length(l2D), 2 );
            for i = 1:length(l2D)
                gradl(i,:) = e2DHat(i,:) * ...
                    [ real(dZdBIJ(i)), real(dZdCIJ(i)); ...
                    imag(dZdBIJ(i)), imag(dZdCIJ(i)) ];
            end
            
            EG = 2 .* sum( ( l2D - L ) .* gradl, 1 )';
            
        end
        
    end

%--------------------------------------------------------------------------
% CONSTRAINT FUNCTION
%--------------------------------------------------------------------------
    function [ c, ceq, DC, DCeq ] = unitdisk( x )
        
        c = x(1)^2 + x(2)^2 - 1;
        ceq = [];
        
        % Gradient of the constraints
        if nargout > 2
            
            DC = 2 .* [x(1); x(2)];
            DCeq = [];
            
        end
        
    end
        

end

