function [W, R, K] = beamAssembly(EA, EI, CNX, EQN, X, D, q)
% Function to compute and assemble the energy, nodal residual force vector,
% and stiffness matrix of an assembly of EB frame elements.
%
% (c) 2015 MAE M168
%
% Input parameters:
% EA: (Vector, NEL x 1) Stretching modulus for the elements
% EI: (Vector, NEL x 1) Bending modulus for the elements
% CNX: (Vector, 2 x NEL) Nodal connectivity array for the elements%
% EQN: (Matrix, 3 x NN) Global equation numbers
% X: (Vector, 2 x NN) Nodal positions
% D: (Vector, 2 x NN) Nodal displacements
%
%
% Output parameters
% W: (Scalar) Internal energy
% R: (Vector, N x 1) Internal nodal forces
% K: (Matrix, N x N) Stiffness matrix

% get number of elements
E = size(CNX,2);
N = sum(sum(EQN > 0));
elemDoF = 3;

% double check that input arrays have consistent sizes.
% assert(size(cnx,1)==E ,'size(cnx,1) = %d.  Should be %d.\n',size(cnx,1),E);
% assert(size(cnx,2)==2, 'size(cnx,2) = %d.  Should be 2.\n',size(cnx,2));

% initialize energy, residual, stiffness
W = 0.0;
R = zeros(N,1);
K = zeros(N,N);

% loop over elements

for e=1:E
    
    % Get both element nodes
    elementNodes = CNX(:,e);

    % Compute DoF locations in D, the global displacement array
    DoFLocations = [elementNodes(1)*elemDoF-2;
                    elementNodes(1)*elemDoF-1;
                    elementNodes(1)*elemDoF;
                    elementNodes(2)*elemDoF-2;
                    elementNodes(2)*elemDoF-1;
                    elementNodes(2)*elemDoF];

    % Copy global displacement values into element displacement array
    d = D(DoFLocations);
    
    x = X(DoFLocations);
    
    % Get element nodal positions, then call TrussBar

    [w,r,k] = beamElement(EA(e),EI(e), x, d, q(e));

    % Assemble element contributions to global energy, residual, stiffness
    W = W+w;

    % Get global DoF corresponding to both element nodes
    globalDoF = reshape(EQN(:,elementNodes),[6,1]);

    % Determine active DoF
    activeDoF = (globalDoF>0);

    % Remove inactive DoF so they aren't included in the assembly
    globalDoF(~activeDoF) = [];

    % Assembly
    R(globalDoF) = R(globalDoF) + r(activeDoF);
    K(globalDoF,globalDoF) = K(globalDoF,globalDoF) + k(activeDoF,activeDoF);
    
end

end


