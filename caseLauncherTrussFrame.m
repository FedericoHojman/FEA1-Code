% Case launcher
clear; close all; 

%% Structure definition
exercise = 3;
onlyBars = true;


switch exercise
    
    case 1
        
        structuralJointsArray=[ 0 0 0          % In [mm]
            2 0 0
            0 2 0
            3 3 0 ]*1000; 
        
        % Input of members connecting joints and which dof they connect
        % Begin | End | Cross Section Orientation
        structuralMembersArray.nodes=[1 2 4
            2 3 4
            3 1 4];
        planeStructure = true;
        
                                  
    case 2
        
        structuralJointsArray=   [0 0 0          % In [mm]
            2 0 0
            4 0 0
            0 2 0
            2 2 0
            4 2 0
            5 5 0]*1000; %For direction only
        
        % Input of members connecting joints and which dof they connect
        % Begin | End | Cross Section Orientation
        structuralMembersArray.nodes=[1 2 7
            2 3 7
            3 6 7
            5 6 7
            5 4 7
            6 2 7
            5 2 7
            4 2 7
            1 4 7
            ];
        planeStructure = true;
    case 3
        structuralJointsArray=[ 0 0 0          % In [mm]
                        1 0 0
                        0 1 0
                        1 1 0
                        0.5 0.5 1
                        1.5 0.5 1
                        0 0.5 0]*1000; %For direction only

% Input of members connecting joints and which dof they connect
% Begin | End | Cross Section Orientation
structuralMembersArray.nodes=[1 2 4
                              2 3 4
                              3 1 4
                              3 4 2
                              4 1 2
                              1 3 4
                              5 6 7
                              1 5 6
                              1 5 6
                              2 6 5
                              4 5 6
                              3 6 5
                              1 6 5
                              4 6 5];
                           planeStructure = false;

        
end


% Connected Dof                          
structuralMembersArray.dof = true(size(structuralMembersArray.nodes,1),12);

% Number of elements in member
structuralMembersArray.refinement = ones(size(structuralMembersArray.nodes,1));

% Member cross section number
structuralMembersArray.crossSection = ones(size(structuralMembersArray.nodes,1));

% Member material number
structuralMembersArray.material = ones(size(structuralMembersArray.nodes,1));
                    
% Cross sections definition
% Area | Inertia Moment in P123 plane | Inertia Moment orthogonal to P123 plane | Torsional Stiffness
membersCrossSection=[pi*20^2 0  0   0]; % Circular section
%membersCrossSection=[pi*20^2 pi*20^4/4 pi*20^4/4 pi*20^4/2]; % Circular section

% Material definition
% Young Modulus | Transverse Modulus | Density 
membersMaterial=[200000 76923 7800]; %MPa kg/m3 Steel

% Structure plot
linearMeshPlot(structuralMembersArray.nodes(:,1:2),structuralJointsArray,'b','Yes');
             
%% Preprocess
    
% Mesh generation
[elementArray,nodesPositionArray]=trussFrameMeshGenerator(structuralMembersArray,structuralJointsArray);

% Problem parameters
nElements=size(elementArray.nodes,1);    %Number of elements
nNodes=size(nodesPositionArray,1);       %Number of nodes
nTotalDof=max(max(elementArray.dof));    %Number of total dofs

switch exercise 
    case 1

        % Boundary conditions
        boundaryConditionsArray = false(nNodes,6);    % Boundary conditions array true=fixed
        boundaryConditionsArray(1,1) = true;
        boundaryConditionsArray(3,[1 2 ]) = true;
           % Rotations elimination

        % Load definition
        pointLoadsArray = zeros(nNodes,6);     % Point load nodal value for each direction
        pointLoadsArray(2,2) = -2000;

    case 2
        boundaryConditionsArray = false(nNodes,6);    % Boundary conditions array true=fixed
        boundaryConditionsArray(1,1) = true;
        boundaryConditionsArray(4,[1 2]) = true;
          % Rotations elimination

        % Load definition
        pointLoadsArray = zeros(nNodes,6);     % Point load nodal value for each direction
        pointLoadsArray([2 3],2) = -2000;
        
    case 3
         boundaryConditionsArray = false(nNodes,6);    % Boundary conditions array true=fixed
        boundaryConditionsArray(1,:) = true;
        boundaryConditionsArray(2,3) = true;
        boundaryConditionsArray(3,[1 3]) = true;
        boundaryConditionsArray(4,[ 2 3]) = true;
          % Rotations elimination

        % Load definition
        pointLoadsArray = zeros(nNodes,6);     % Point load nodal value for each direction
        pointLoadsArray(6,3) = -2000;
        
        
        
        
end

if onlyBars
	 boundaryConditionsArray(:,[ 4 5 6]) = true; 
end

if planeStructure
    boundaryConditionsArray(:,[3 6]) = true; 
end

%% Solver

% Stiffness calculation and assembly
[stiffnessMatrix]=assemble1DStiffnessMatrix(elementArray,nodesPositionArray,structuralJointsArray,membersCrossSection,membersMaterial);

% Matrix reduction
isFixed = reshape(boundaryConditionsArray',1,[])';
isFree = ~isFixed;

% Loads vector rearrangement
loadsVector = reshape(pointLoadsArray',1,[])';

% Equation solving
displacementsReducedVector = stiffnessMatrix(isFree,isFree)\loadsVector(isFree);

% Reconstruction
displacementsVector = zeros(nTotalDof,1);
displacementsVector(isFree) = displacementsVector(isFree) + displacementsReducedVector;

%% Postprocess
magnificationScale=10000;

% Nodal displacements rearranged
nodalDisplacements=reshape(displacementsVector,6,size(nodesPositionArray,1))';
nodalPositions=nodesPositionArray+nodalDisplacements(:,1:3)*magnificationScale;

% Deformed Structure plot
linearMeshPlot(elementArray.nodes(:,1:2),nodalPositions,'r','No');