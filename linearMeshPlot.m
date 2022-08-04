function linearMeshPlot(elementNodesArray,nodesPositionArray,color,numbering)
% Linear Mesh plotter
% 
% linearMeshPlot(elementNodesArray,nodesPositionArray,color,numbering)
% 
% Nodal order in elementNodesArray:
%       
%  1--3--4-...--2
%  
% elementNodesArray:    element conectivity matrix
% nodesPositionArray:   nodal position in cartesian coordinates
% color:                element borders color
% numbering:            'Yes' 'No'
%

%% Input control
if (nargin < 4) || isempty(numbering)
    numbering='No';
end

%% Definition
nElements = size(elementNodesArray,1);
nNodesPerElement = size(elementNodesArray,2);

nNodes = size(nodesPositionArray,1);
nDimensions = size(nodesPositionArray,2);

switch nDimensions
    case 1
        nodesPositionArray(:,[2 3])=0;
    case 2
        nodesPositionArray(:,[3])=0;
end

%% Element plotting and entity numbering
view ([0 0 1]); daspect([1 1 1]); hold on;

% Node location and numbering
for iNode = 1:nNodes;
    plot3(nodesPositionArray(iNode,1),nodesPositionArray(iNode,2),nodesPositionArray(iNode,3),'o','MarkerFaceColor','r','MarkerEdgeColor','b','MarkerSize',4);
    if strcmp(numbering,'Yes')
        text(nodesPositionArray(iNode,1),nodesPositionArray(iNode,2),nodesPositionArray(iNode,3),num2str(iNode),'VerticalAlignment','bottom','Color','r','FontSize',8);
    end
end

% Element ploting and numbering
for iElement = 1:nElements;
    plot3(nodesPositionArray(elementNodesArray(iElement,:),1),nodesPositionArray(elementNodesArray(iElement,:),2),nodesPositionArray(elementNodesArray(iElement,:),3),color);
    if strcmp(numbering,'Yes')
        xCoordinateElementCenter = mean(nodesPositionArray(elementNodesArray(iElement,:),1));
        yCoordinateElementCenter = mean(nodesPositionArray(elementNodesArray(iElement,:),2));
        zCoordinateElementCenter = mean(nodesPositionArray(elementNodesArray(iElement,:),3));
        text(xCoordinateElementCenter,yCoordinateElementCenter,zCoordinateElementCenter,num2str(iElement),'EdgeColor','k','Color','k','FontSize',8);
    end
end
