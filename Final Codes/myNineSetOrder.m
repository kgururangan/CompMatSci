function [orderVec] = myNineSetOrder(numnodes,Nx,Ny)


% global node numbers for boundary and interior points  
nodes_B = [2:Nx;
           2+numnodes:Nx+numnodes]; 
nodes_L = [Nx+2:Nx+1:numnodes-2*Nx-1;
           Nx+2+numnodes:Nx+1:numnodes-2*Nx-1+numnodes];
nodes_T = [numnodes-Nx+1:numnodes-1;
           numnodes-Nx+1+numnodes:numnodes-1+numnodes];
nodes_R = [2*(Nx+1):Nx+1:numnodes-Nx-1;
           2*(Nx+1)+numnodes:Nx+1:numnodes-Nx-1+numnodes];
nodes_LB = [1;
            1+numnodes];
nodes_LT = [numnodes-Nx;
            numnodes-Nx+numnodes];
nodes_RB = [Nx+1;
            Nx+1+numnodes];
nodes_RT = [numnodes;
            numnodes+numnodes];
NN = 1:numnodes;    
NN( [nodes_B(1,:), nodes_L(1,:), nodes_T(1,:), nodes_R(1,:), nodes_LB(1,:), nodes_LT(1,:), nodes_RB(1,:), nodes_RT(1,:)] ) = [];
nodes_int = [NN;
             NN + numnodes];
% order is chosen based on the ordering of the transformation matrix        
%orderVec = [nodes_L(:)', nodes_R(:)', nodes_B(:)', nodes_T(:)', nodes_LB(:)', nodes_RB(:)', nodes_LT(:)', nodes_RT(:)', nodes_int(:)'];

orderVec = [nodes_B(:)',nodes_L(:)',nodes_LB(:)',nodes_int(:)',nodes_T(:)',nodes_R(:)',nodes_RB(:)',nodes_LT(:)',nodes_RT(:)'];

end