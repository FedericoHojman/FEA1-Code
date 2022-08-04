function [elementArray,nodesPositionArray]=trussFrameMeshGenerator(structuralMembersArray,structuralJointsArray);
% Mesh Generator for a Trusses / Frames
% 
% [elementArray,nodesPositionArray]=trussFrameMeshGenerator(elementType,structuralMembersArray,structuralJointsArray);
%
% elementArray.dof:             Element Dof conectivity matrix
% elementArray.nodes:           Element Nodes conectivity matrix
% elementArray.auxiliarPoint    Element auxiliar orientation point declaration. 
% elementArray.crossSection     Element cross section assignment.
% elementArray.material         Element material assignment.
% nodesPositionArray:           Nodal position in cartesian coordinates
%
% structuralMembersArray:       Structural information
% structuralJointsArray:        Structural joint locations
%


%% Sequential element construction
nStructuralMembers=size(structuralMembersArray.nodes,1);

% Multimember connection nodes lookup
nNodes=0;
for iStructuralJoints=1:max(max(structuralMembersArray.nodes(:,1:2))) %Only nodes are considered not auxiliar points
    nNodes=nNodes+1;
    nodesPositionArray(nNodes,:)=structuralJointsArray(iStructuralJoints,:);
end

% Structural members mesh generation
nElements=0;
for iMember=1:nStructuralMembers
    nElements=nElements+1;
    elementArray.dof(nElements,:)=[structuralMembersArray.nodes(iMember,1)*[6 6 6 6 6 6]-[5 4 3 2 1 0] structuralMembersArray.nodes(iMember,2)*[6 6 6 6 6 6]-[5 4 3 2 1 0]];
    elementArray.nodes(nElements,:)=[structuralMembersArray.nodes(iMember,1:2)];
    elementArray.auxiliarPoint(nElements)=[structuralMembersArray.nodes(iMember,3)];
    elementArray.crossSection(nElements)=structuralMembersArray.crossSection(iMember);
    elementArray.material(nElements)=structuralMembersArray.material(iMember);
end

return