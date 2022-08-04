function [stiffnessMatrix]=assemble1DStiffnessMatrix(elementArray,nodesPositionArray,structuralJointsArray,membersCrossSection,membersMaterial)
% 1D finite element matrix assembler
% 
% assembleStiffnessMatrix(elementType,elementArray.nodes,nodesPositionArray)
%
% stiffnessMatrix:      Assembled stiffness matrix
%
% elementArray.dof
% elementArray.nodes
% nodesPositionArray
% membersCrossSection
% membersMaterial
%

%% Definitions
nElements=size(elementArray.nodes,1);      %Number of elements
nTotalDof=max(max(elementArray.dof));    %Number of Dofs

%% Stiffness matrix assembly
stiffnessMatrix = zeros(nTotalDof);

for iElement = 1:nElements
    % Rotation
    V1 = nodesPositionArray(elementArray.nodes(iElement,2),:) - nodesPositionArray(elementArray.nodes(iElement,1),:);
    elementLength = norm(V1);
    V1 = V1/elementLength;
    V2 = structuralJointsArray(elementArray.auxiliarPoint(iElement),:) - nodesPositionArray(elementArray.nodes(iElement,1),:);
    V2 = V2/norm(V2);
    V3 = cross(V1,V2);
    V3 = V3/norm(V3);
    V2 = cross(V3,V1);
    V2 = V2/norm(V2);
    lambdaProjectionMatrix = [V1
                              V2
                              V3];
    projectionMatrix = blkdiag(lambdaProjectionMatrix,lambdaProjectionMatrix,lambdaProjectionMatrix,lambdaProjectionMatrix);
    
    % Coefficients
    X  = membersMaterial(elementArray.material(iElement),1)*membersCrossSection(elementArray.crossSection(iElement),1)/elementLength;
    S  = membersMaterial(elementArray.material(iElement),2)*membersCrossSection(elementArray.crossSection(iElement),4)/elementLength;
    Y1 = 12*membersMaterial(elementArray.material(iElement),1)*membersCrossSection(elementArray.crossSection(iElement),2)/elementLength^3;
    Y2 = 6*membersMaterial(elementArray.material(iElement),1)*membersCrossSection(elementArray.crossSection(iElement),2)/elementLength^2;
    Y3 = 4*membersMaterial(elementArray.material(iElement),1)*membersCrossSection(elementArray.crossSection(iElement),2)/elementLength^1;
    Y4 = 2*membersMaterial(elementArray.material(iElement),1)*membersCrossSection(elementArray.crossSection(iElement),2)/elementLength^1;
    Z1 = 12*membersMaterial(elementArray.material(iElement),1)*membersCrossSection(elementArray.crossSection(iElement),3)/elementLength^3;
    Z2 = 6*membersMaterial(elementArray.material(iElement),1)*membersCrossSection(elementArray.crossSection(iElement),3)/elementLength^2;
    Z3 = 4*membersMaterial(elementArray.material(iElement),1)*membersCrossSection(elementArray.crossSection(iElement),3)/elementLength^1;
    Z4 = 2*membersMaterial(elementArray.material(iElement),1)*membersCrossSection(elementArray.crossSection(iElement),3)/elementLength^1;
    
    elementalStiffnessMatrix=[  X   0   0  0   0   0 -X   0   0  0   0   0
                                0  Y1   0  0   0  Y2  0 -Y1   0  0   0  Y2
                                0   0  Z1  0 -Z2   0  0   0 -Z1  0 -Z2   0
                                0   0   0  S   0   0  0   0   0 -S   0   0
                                0   0 -Z2  0  Z3   0  0   0  Z2  0  Z4   0
                                0  Y2   0  0   0  Y3  0 -Y2   0  0   0  Y4
                               -X   0   0  0   0   0  X   0   0  0   0   0
                                0 -Y1   0  0   0 -Y2  0  Y1   0  0   0 -Y2
                                0   0 -Z1  0  Z2   0  0   0  Z1  0  Z2   0
                                0   0   0 -S   0   0  0   0   0  S   0   0
                                0   0 -Z2  0  Z4   0  0   0  Z2  0  Z3   0
                                0  Y2   0  0   0  Y4  0 -Y2   0  0   0  Y3];
    
    elementalStiffnessMatrix = projectionMatrix'*elementalStiffnessMatrix*projectionMatrix;
     
    % Matrix assembly
    
    stiffnessMatrix(elementArray.dof(iElement,:),elementArray.dof(iElement,:)) = stiffnessMatrix(elementArray.dof(iElement,:),elementArray.dof(iElement,:)) + elementalStiffnessMatrix;

end

