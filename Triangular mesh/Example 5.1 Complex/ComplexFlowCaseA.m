function ComplexFlowCaseA % case (a) : vertical flow
% time: 2020.6.10 - 2020.6.10
% author : xuziyao
% Basis function is orthogonal polynomials(any degree) on the triangular element.
% Mixed boundary condition 
% Equation:
% s = -grad p
% (I+k_m*epsilon/k_n*delta*1*sigma*sigma)u = (k_m+epsilon*kt*delta*1*nu*nu)s
% q + div(u) = 0 (q=-f)
% =>
% A0*Sx=B1*P+D1+Cp11*Ux+Cp12*Uy, A0*Sy=B2*P+D2+Cp21*Ux+Cp22*Uy, 
% Au11*Ux+Au12*Uy=As11*Sx+As12*Sy  
% Au21*Ux+Au22*Uy=As21*Sx+As22*Sy
% A0*Q=C1*Ux+C2*Uy+Nv+Em*P+Ed
% =>
% FINISHED
clc,clear
format long
% set the geometry and triangulation manners:
%--------------------------------------------------------------------------
% mesh 1:(Distmesh)--------------------------------------------------------
%pv=[0,150; 400,100; 800,150; 1200,100; 1600,150; 1600,-1000; 0,-1000]; h0=5; geps=.0001*h0;
%[coordinates,elements3]=distmesh2d(@dpoly,@huniform,10,[0,-1000;1600,160],pv,pv);elements3 = elements3';
% mesh2:(existing mesh)----------------------------------------------------
coordinates = 0;
load('boxdfm_complex.mat')
%load('ccdfm_complex.mat');
elements3=delaunayn(coordinates); % List of triangles
elements3 = elements3';
bcol=[1,1,1];figure;simpplot(coordinates,elements3',[],bcol); hold on;
% mesh 3:(Matlab mesh, non-conforming)-----------------------------------------------------
%geodesc = [2 4 0 1 1 0 0 0 1 1 ]';
%h = 0.05;
%geometry = decsg(geodesc); % Create the geometry
%model = createpde;
%geometryFromEdges(model,geometry);
%generateMesh(model,'Hmax',h,'GeometricOrder','linear');
%coordinates = model.Mesh.Nodes'; elements3 =  model.Mesh.Elements;
%%pdeplot(coordinates',elements3,'EdgeColor','black','ElementLabels','off')
%bcol=[1,1,1];figure;simpplot(coordinates,elements3',[],bcol); hold on;
%--------------------------------------------------------------------------
% mesh 4:(Matlab mesh,conforming, domain decomposition)-----------------------------------------------------
%geodesc = zeros(10,10); % 10 subdomains
%geodesc(:,1) = [2;4;0;0.5;0.5;0;0.5;0.5;1;1];
%geodesc(:,2) = [2;4;0;0.5;0.5;0;0;0;0.5;0.5];
%geodesc(:,3) = [2;4;0.5;1;1;0.5;0;0;0.5;0.5];
%geodesc(:,4) = [2;4;0.5;0.75;0.75;0.5;0.75;0.75;1;1];
%geodesc(:,5) = [2;4;0.75;1;1;0.75;0.75;0.75;1;1];
%geodesc(:,6) = [2;4;0.75;1;1;0.75;0.5;0.5;0.75;0.75];
%geodesc(:,7) = [2;4;0.5;0.625;0.625;0.5;0.625;0.625;0.75;0.75];
%geodesc(:,8) = [2;4;0.625;0.75;0.75;0.625;0.625;0.625;0.75;0.75];
%geodesc(:,9) = [2;4;0.5;0.625;0.625;0.5;0.5;0.5;0.625;0.625];
%geodesc(:,10)= [2;4;0.625;0.75;0.75;0.625;0.5;0.5;0.625;0.625];
%ns = char('A','B','C','D','E','F','G','H','I','J'); ns = ns';
%sf = 'A+B+C+D+E+F+G+H+I+J';
%h = 0.2;
%geometry = decsg(geodesc,sf,ns); % Create the geometry
%model = createpde;
%geometryFromEdges(model,geometry);
%generateMesh(model,'Hmax',h,'GeometricOrder','linear');
%coordinates = model.Mesh.Nodes'; elements3 =  model.Mesh.Elements;
%pdeplot(coordinates',elements3,'EdgeColor','black','ElementLabels','off');
%bcol=[1,1,1];figure;simpplot(coordinates,elements3',[],bcol);hold on;
%--------------------------------------------------------------------------
EtoEmap = triangulation_neighbor_triangles ( elements3 ); 
[ Tarea , Elenth , Jacobimat ] = ComputeGeometryData( elements3 ,coordinates );
alpha = 5/min(Tarea)^1.5; % penalty parameter of pressure jump
beta = 5/min(Tarea)^1.5; % penalty parameter of flux jump
% Load the fracture and barrier data
% Nf xA yA xB yB
GeoData = [...
1 0.0500 0.4160 0.2200 0.0624;
2 0.0500 0.2750 0.2500 0.1350;
3 0.1500 0.6300 0.4500 0.0900;
4 0.1500 0.9167 0.4000 0.5000;
5 0.6500 0.8333 0.8500 0.1667;
6 0.7000 0.2350 0.8500 0.1675;
7 0.6000 0.3800 0.8500 0.2675;
8 0.3500 0.9714 0.8000 0.7143;
9 0.7500 0.9574 0.9500 0.8155;
10 0.1500 0.8363 0.4000 0.9727];
GeoData(:,1)=[]; GeoData = (GeoData-pi/6)*(1-1e-8)+pi/6;
% ParametersofFracture: x_c,y_c,length,theta,width,permeability_nu,permeability_sigma,xa,ya,xb,yb
% ,the # of elements passed by this fracture. 
NumberofFractures = 8; 
if NumberofFractures>0
ParametersofFracture = zeros(NumberofFractures, 12); 
ParametersofFracture(:,8:11) = GeoData([1:3,6:10],:);
ParametersofFracture(:,1) =  (GeoData([1:3,6:10],1)+GeoData([1:3,6:10],3))/2;% position of fractures
ParametersofFracture(:,2) =  (GeoData([1:3,6:10],2)+GeoData([1:3,6:10],4))/2;% position of fractures
ParametersofFracture(:,3) =  sqrt((GeoData([1:3,6:10],1)-GeoData([1:3,6:10],3)).^2+(GeoData([1:3,6:10],2)-GeoData([1:3,6:10],4)).^2);% position of fractures
ParametersofFracture(:,4) =  atan2(GeoData([1:3,6:10],4)-GeoData([1:3,6:10],2),GeoData([1:3,6:10],3)-GeoData([1:3,6:10],1));% position of fractures
ParametersofFracture(:,5) = 1e-4; % width of fractures
ParametersofFracture(:,6)= 1e4 ; % tangential permeability_nu of fractures
ParametersofFracture(:,7)= 1e4 ; % normal permeability_nu of fractures
end
% ParametersofBarrier: x_c,y_c,length,theta,width,permeability_nu,permeability_sigma,xa,ya,xb,yb
% ,the # of elements passed by this Barrier. 
NumberofBarriers = 2; 
if NumberofBarriers>0
ParametersofBarrier = zeros(NumberofBarriers, 12); 
ParametersofBarrier(:,8:11) = GeoData(4:5,:);
ParametersofBarrier(:,1) = (GeoData(4:5,1)+GeoData(4:5,3))/2;% position of fractures
ParametersofBarrier(:,2) = (GeoData(4:5,2)+GeoData(4:5,4))/2;% position of fractures
ParametersofBarrier(:,3) = sqrt((GeoData(4:5,1)-GeoData(4:5,3)).^2+(GeoData(4:5,2)-GeoData(4:5,4)).^2);% position of fractures
ParametersofBarrier(:,4) = atan2(GeoData(4:5,4)-GeoData(4:5,2),GeoData(4:5,3)-GeoData(4:5,1));% position of fractures
ParametersofBarrier(:,5) = 1e-4; % width of fractures
ParametersofBarrier(:,6)= 1e-4 ; % tangential permeability_nu of fractures
ParametersofBarrier(:,7)= 1e-4 ; % normal permeability_nu of fractures
end
% plot the triangulation and fractures and barriers
if NumberofFractures>0
plot(ParametersofFracture(:,[8,10])',ParametersofFracture(:,[9,11])',...
   'r-','LineWidth',1.5); 
end
if NumberofBarriers>0
plot(ParametersofBarrier(:,[8,10])',ParametersofBarrier(:,[9,11])',...
   'b-','LineWidth',1.5); 
end
% (2) compute basis data
HighOrderDegree = 1 ; % polynomial's order in each element 
HighOrderNp = (HighOrderDegree+1) * (HighOrderDegree+2) / 2; % degree of freedom in each element for high order polynomial
Gauss_num = 10 ; % the number of Gauss quadrature point in [-1,1]
Lobatto_num = 6 ; % the number of Lobatto quadrature point in [-1,1]
[Gauss_x,Gauss_w] = JacobiGQ(0,0,Gauss_num-1); 
Gauss_x = 0 + ( Gauss_x - (-1) ) / 2 ; % Gauss_x now is in [0,1].
Gauss_w = Gauss_w / 2 ;% sum(Gauss_w) = 1
[ quad_lamda1 , quad_lamda2 , quad_w ] = quad_rule1 ( 66 ) ; % positive quadrature rule
[basisP,GradrbasisP,GradsbasisP,boundaryP,GradrboundaryP,GradsboundaryP] = ...
    SetLegendreBasisData(HighOrderDegree,quad_lamda1,quad_lamda2,Gauss_x) ;
% compute fractures basic data
if NumberofFractures>0
[ParametersofFracture,FracturesPath,FracturesLength,basisPonFracture,GradrbasisPonFracture,GradsbasisPonFracture, ...
    FracturesStartRefCoord,FracturesEndRefCoord] = ...
    SetFractureBasisData(coordinates,elements3,ParametersofFracture,Gauss_x,HighOrderDegree);
end
% compute fractures and Barriers basic data
if NumberofBarriers>0
[ParametersofBarrier,BarriersPath,BarriersLength,basisPonBarrier,GradrbasisPonBarrier,GradsbasisPonBarrier, ...
    BarriersStartRefCoord,BarriersEndRefCoord] = ...
    SetFractureBasisData(coordinates,elements3,ParametersofBarrier,Gauss_x,HighOrderDegree);
end
FluxPenaltyEdges = JumpPenaltyEdges(EtoEmap,BarriersPath,BarriersStartRefCoord,BarriersEndRefCoord);
[HighOrderP,HighOrderUx,HighOrderUy] = EllipticLDGSolver;
% draw contour
ShowDG2DSolution_Ler(coordinates,elements3,HighOrderP);
colormap('jet');caxis([1,4]);cb=colorbar;cb.Position = cb.Position.*[1.07,1,1,1];
axis equal; axis off;ax=axis;axis(ax*1.001); view([0,90]);
hold on; plot3(GeoData(:,[1,3])',GeoData(:,[2,4])',100*ones(size(GeoData(:,[2,4])')),...
   'w-','LineWidth',0.05); hold off; 
drawnow
% profiles on slices
NumberofSlices = 1; 
ParametersofSlice = zeros(NumberofSlices, 12); 
ParametersofSlice(:,8:11)=([0.0,0.5,1.0,0.9]-pi/6)*(1-exp(-7))+pi/6;
ParametersofSlice(:,1) =  (ParametersofSlice(:,8)+ParametersofSlice(:,10))/2;% position of fractures
ParametersofSlice(:,2) =  (ParametersofSlice(:,9)+ParametersofSlice(:,11))/2;% position of fractures
ParametersofSlice(:,3) =  sqrt((ParametersofSlice(:,8)-ParametersofSlice(:,10)).^2+(ParametersofSlice(:,9)-ParametersofSlice(:,11)).^2);% position of fractures
ParametersofSlice(:,4) =  atan2(ParametersofSlice(:,11)-ParametersofSlice(:,9),ParametersofSlice(:,10)-ParametersofSlice(:,8));% position of fractures
[ParametersofSlice,SlicesPath,SlicesLength,basisPonSlice,GradrbasisPonSlice,GradsbasisPonSlice, ...
    SlicesStartRefCoord,SlicesEndRefCoord] = ...
    SetFractureBasisData(coordinates,elements3,ParametersofSlice,Gauss_x,HighOrderDegree);
for k = 1 : NumberofSlices
% slant line (0.0, 0.5) ---- (1.0, 0.9)
figure; hold on;
    for ell = 1 : ParametersofSlice(k,12)
        CurrentIndex = sum(ParametersofSlice(1:k-1,12)) + ell;
        jj = SlicesPath( CurrentIndex );
        SlicesStartPhysicCoord = SlicesStartRefCoord(CurrentIndex,:) * coordinates(elements3(:,jj),:); 
        SlicesEndPhysicCoord = SlicesEndRefCoord(CurrentIndex,:) * coordinates(elements3(:,jj),:); 
        Slice_XY = (1-[0;Gauss_x;1])*SlicesStartPhysicCoord+[0;Gauss_x;1]*SlicesEndPhysicCoord;
        Slice_Re = sqrt(sum((Slice_XY-ones(size([0;Gauss_x;1]))*ParametersofSlice(k,8:9)).^2,2));
        Slice_Ue = basisPonSlice((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),[end-1,1:end-2,end])'*HighOrderP(:,jj);
        handle1=plot(Slice_Re,Slice_Ue,'b','LineWidth',1.5);
    end    
    hold off;
% format the plot
xlabel('arc length', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
ylabel('pressure', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'FontSize', 10);
set(gcf, 'PaperPositionMode', 'auto')
fig_pos = get(gcf, 'PaperPosition');
set(gcf, 'PaperSize', [fig_pos(3) fig_pos(4)+1])
hold off;    
end

%---------------ELLIPTIC EQUATION LDG SOLVER------------------------------
function [HighOrderP,HighOrderUx,HighOrderUy] = EllipticLDGSolver
% define the variables
[HighOrderP,HighOrderUx,HighOrderUy] = deal( zeros( HighOrderNp , size(elements3,2) ) ); % unknowns 
[IMat,JMat] = ...
    deal(ones(HighOrderNp*size(elements3,2)*(size(EtoEmap,1)+1),HighOrderNp));
[VA0,VA0inv,VAu11,VAu12,VAu21,VAu22,VAs11,VAs12,VAs21,VAs22,VB1,VB2,VC1,VC2,VCp11,VCp12,VCp21,VCp22,VEm] = ...
    deal(zeros(HighOrderNp*size(elements3,2)*(size(EtoEmap,1)+1),HighOrderNp));
[D1,D2,Nv,Ed,A0Q] = deal(zeros(size( elements3,2 )*HighOrderNp,1));
% LDG solver:
for jj = 1 : size(elements3,2) % run element from first to end  
xy = quad_lamda1 * coordinates(elements3(1,jj),:)... 
     + quad_lamda2 * coordinates(elements3(2,jj),:)... 
     +(1-quad_lamda1-quad_lamda2)*coordinates(elements3(3,jj),:) ; 
I = (jj-1)*HighOrderNp+(1:HighOrderNp);
II = (jj-1)*(size(EtoEmap,1)+1)*HighOrderNp+(1:HighOrderNp);
[JMat(II,:),IMat(II,:)]   = meshgrid(I,I);
VA0(II,:) = VA0(II,:) + Tarea(jj)*basisP*diag(quad_w)*basisP';
VA0inv(II,:) = inv(Tarea(jj)*basisP*diag(quad_w)*basisP');
VAs11(II,:)  = VAs11(II,:) + Tarea(jj)*basisP*diag(k_m(xy).*quad_w)*basisP';
VAs22(II,:)  = VAs22(II,:) + Tarea(jj)*basisP*diag(k_m(xy).*quad_w)*basisP';
VAu11(II,:)  = VAu11(II,:) + Tarea(jj)*basisP*diag(quad_w)*basisP';
VAu22(II,:)  = VAu22(II,:) + Tarea(jj)*basisP*diag(quad_w)*basisP';
VB1(II,:) = VB1(II,:) + Tarea(jj)*(GradrbasisP*Jacobimat(1,1,jj)+GradsbasisP*Jacobimat(2,1,jj))*diag(quad_w)*basisP';
VB2(II,:) = VB2(II,:) + Tarea(jj)*(GradrbasisP*Jacobimat(1,2,jj)+GradsbasisP*Jacobimat(2,2,jj))*diag(quad_w)*basisP';
VC1(II,:) = VC1(II,:) + Tarea(jj)*(GradrbasisP*Jacobimat(1,1,jj)+GradsbasisP*Jacobimat(2,1,jj))*diag(quad_w)*basisP';
VC2(II,:) = VC2(II,:) + Tarea(jj)*(GradrbasisP*Jacobimat(1,2,jj)+GradsbasisP*Jacobimat(2,2,jj))*diag(quad_w)*basisP';
A0Q((jj-1)*HighOrderNp+(1:HighOrderNp)) = A0Q((jj-1)*HighOrderNp+(1:HighOrderNp)) + ...
    Tarea(jj)*basisP*(q(xy).*quad_w);
% Edge integral
    for ii = 1 : 3  % the ii-th edge of jj-th element
        vertex1 = mod(ii,3)+1 ; 
        vertex2 = mod(ii+1,3)+1 ;
        edge_vec =  coordinates(elements3(vertex2,jj),:)' - ...
            coordinates(elements3(vertex1,jj),:)' ;
        outer_vec = - [0,-1;1,0]*(edge_vec/norm(edge_vec));
        E = EtoEmap(ii,jj) ; % the index of neighborhood of element jj 
        Boundary_xy = (1-Gauss_x)*coordinates(elements3(vertex1,jj),:) + Gauss_x*coordinates(elements3(vertex2,jj),:);
        if E ~= -1 % the edge is not in boundary
            % find the index of the edge in neighborhood element:
            pp =  EtoEmap(:,E) == jj ;   
            J = (E-1)*HighOrderNp+1:E*HighOrderNp ;
            IJ = ((jj-1)*(size(EtoEmap,1)+1)+ii)*HighOrderNp+(1:HighOrderNp);
            [JMat(IJ,:),IMat(IJ,:)] = meshgrid(J,I);
            % central flux p^=(p1+p2)/2, u^=(u1+u2)/2
            VB1(II,:) = VB1(II,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(1)*1/2*diag(Gauss_w)*boundaryP(:,:,ii)';
            VB1(IJ,:) = VB1(IJ,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(1)*1/2*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)';
            VB2(II,:) = VB2(II,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(2)*1/2*diag(Gauss_w)*boundaryP(:,:,ii)';
            VB2(IJ,:) = VB2(IJ,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(2)*1/2*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)';
            VC1(II,:) = VC1(II,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(1)*1/2*diag(Gauss_w)*boundaryP(:,:,ii)';
            VC1(IJ,:) = VC1(IJ,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(1)*1/2*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)';
            VC2(II,:) = VC2(II,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(2)*1/2*diag(Gauss_w)*boundaryP(:,:,ii)';
            VC2(IJ,:) = VC2(IJ,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(2)*1/2*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)';
            if FluxPenaltyEdges(ii,jj) == 1 % penalty on jump of u:
            VCp11(II,:) = VCp11(II,:) - beta*Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(1)*diag(Gauss_w)*boundaryP(:,:,ii)'*(outer_vec(1));
            VCp11(IJ,:) = VCp11(IJ,:) - beta*Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(1)*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)'*(-outer_vec(1));
            VCp12(II,:) = VCp12(II,:) - beta*Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(1)*diag(Gauss_w)*boundaryP(:,:,ii)'*(outer_vec(2));
            VCp12(IJ,:) = VCp12(IJ,:) - beta*Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(1)*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)'*(-outer_vec(2));
            VCp21(II,:) = VCp21(II,:) - beta*Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(2)*diag(Gauss_w)*boundaryP(:,:,ii)'*(outer_vec(1));
            VCp21(IJ,:) = VCp21(IJ,:) - beta*Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(2)*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)'*(-outer_vec(1));
            VCp22(II,:) = VCp22(II,:) - beta*Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(2)*diag(Gauss_w)*boundaryP(:,:,ii)'*(outer_vec(2));
            VCp22(IJ,:) = VCp22(IJ,:) - beta*Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(2)*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)'*(-outer_vec(2));
            else % penalty on jump of p:
            VEm(II,:) = VEm(II,:) - alpha*Elenth(ii,jj)*boundaryP(:,:,ii)*diag(Gauss_w)*boundaryP(:,:,ii)';
            VEm(IJ,:) = VEm(IJ,:) + alpha*Elenth(ii,jj)*boundaryP(:,:,ii)*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)';
            end
        elseif abs(outer_vec(2))>1/2 % Dirichlet boundary.
            D1((jj-1)*HighOrderNp+(1:HighOrderNp)) = D1((jj-1)*HighOrderNp+(1:HighOrderNp)) - ...
                Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(1)*(Gauss_w.*p(Boundary_xy));
            D2((jj-1)*HighOrderNp+(1:HighOrderNp)) = D2((jj-1)*HighOrderNp+(1:HighOrderNp)) - ...
                Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(2)*(Gauss_w.*p(Boundary_xy));
            VC1(II,:) = VC1(II,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(1)*diag(Gauss_w)*boundaryP(:,:,ii)';
            VC2(II,:) = VC2(II,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(2)*diag(Gauss_w)*boundaryP(:,:,ii)';
            Ed((jj-1)*HighOrderNp+(1:HighOrderNp)) = Ed((jj-1)*HighOrderNp+(1:HighOrderNp)) + ...
                alpha*Elenth(ii,jj)*boundaryP(:,:,ii)*(Gauss_w.*p(Boundary_xy));
            VEm(II,:) = VEm(II,:) - alpha*Elenth(ii,jj)*boundaryP(:,:,ii)*diag(Gauss_w)*boundaryP(:,:,ii)';
        else % Neumann boundary
            VB1(II,:) = VB1(II,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(1)*diag(Gauss_w)*boundaryP(:,:,ii)';
            VB2(II,:) = VB2(II,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(2)*diag(Gauss_w)*boundaryP(:,:,ii)';
            Nv((jj-1)*HighOrderNp+(1:HighOrderNp)) = Nv((jj-1)*HighOrderNp+(1:HighOrderNp)) - ...
                Elenth(ii,jj)*boundaryP(:,:,ii)*(Gauss_w.*q_N(Boundary_xy,outer_vec));
        end
    end
end
% process the fractures
for k = 1 : NumberofFractures
    FractureNu = [cos(ParametersofFracture(k,4));sin(ParametersofFracture(k,4))];
    FractureSigma = [cos(pi/2+ParametersofFracture(k,4));sin(pi/2+ParametersofFracture(k,4))];
    for ell = 1 : ParametersofFracture(k,12)
    CurrentIndex = sum(ParametersofFracture(1:k-1,12)) + ell;
    jj = FracturesPath( CurrentIndex );
    if ~ismember(jj,BarriersPath)
    II = (jj-1)*(size(EtoEmap,1)+1)*HighOrderNp+(1:HighOrderNp);
    VAs11(II,:) = VAs11(II,:)  + ...
        ParametersofFracture(k,5)*ParametersofFracture(k,6)*FracturesLength( CurrentIndex )*...
        basisPonFracture((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)*diag(Gauss_w)*...
        basisPonFracture((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)'*FractureNu(1)*FractureNu(1);
    VAs12(II,:) = VAs12(II,:)  + ...
        ParametersofFracture(k,5)*ParametersofFracture(k,6)*FracturesLength( CurrentIndex )*...
        basisPonFracture((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)*diag(Gauss_w)*...
        basisPonFracture((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)'*FractureNu(1)*FractureNu(2);
    VAs21(II,:) = VAs21(II,:)  + ...
        ParametersofFracture(k,5)*ParametersofFracture(k,6)*FracturesLength( CurrentIndex )*...
        basisPonFracture((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)*diag(Gauss_w)*...
        basisPonFracture((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)'*FractureNu(2)*FractureNu(1);
    VAs22(II,:) = VAs22(II,:)  + ...
        ParametersofFracture(k,5)*ParametersofFracture(k,6)*FracturesLength( CurrentIndex )*...
        basisPonFracture((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)*diag(Gauss_w)*...
        basisPonFracture((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)'*FractureNu(2)*FractureNu(2);
    end
    end
end
% process the barrier
for k = 1 : NumberofBarriers
    BarrierNu = [cos(ParametersofBarrier(k,4));sin(ParametersofBarrier(k,4))];
    BarrierSigma = [cos(pi/2+ParametersofBarrier(k,4));sin(pi/2+ParametersofBarrier(k,4))];
    for ell = 1 : ParametersofBarrier(k,12)
    CurrentIndex = sum(ParametersofBarrier(1:k-1,12)) + ell;
    jj = BarriersPath( CurrentIndex );
    II = (jj-1)*(size(EtoEmap,1)+1)*HighOrderNp+(1:HighOrderNp);
    BarriersStartPhysicCoord = BarriersStartRefCoord(CurrentIndex,:) * coordinates(elements3(:,jj),:);
    BarriersEndPhysicCoord = BarriersEndRefCoord(CurrentIndex,:) * coordinates(elements3(:,jj),:);
    Barrier_xy = (1-Gauss_x)*BarriersStartPhysicCoord+Gauss_x*BarriersEndPhysicCoord;
    VAu11(II,:) = VAu11(II,:)  + ...
        ParametersofBarrier(k,5)/ParametersofBarrier(k,7)*BarriersLength( CurrentIndex )*...
        basisPonBarrier((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)*diag(Gauss_w.*k_m(Barrier_xy))*...
        basisPonBarrier((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)'*BarrierSigma(1)*BarrierSigma(1);
    VAu12(II,:) = VAu12(II,:)  + ...
        ParametersofBarrier(k,5)/ParametersofBarrier(k,7)*BarriersLength( CurrentIndex )*...
        basisPonBarrier((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)*diag(Gauss_w.*k_m(Barrier_xy))*...
        basisPonBarrier((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)'*BarrierSigma(1)*BarrierSigma(2);
    VAu21(II,:) = VAu21(II,:)  + ...
        ParametersofBarrier(k,5)/ParametersofBarrier(k,7)*BarriersLength( CurrentIndex )*...
        basisPonBarrier((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)*diag(Gauss_w.*k_m(Barrier_xy))*...
        basisPonBarrier((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)'*BarrierSigma(2)*BarrierSigma(1);
    VAu22(II,:) = VAu22(II,:)  + ...
        ParametersofBarrier(k,5)/ParametersofBarrier(k,7)*BarriersLength( CurrentIndex )*...
        basisPonBarrier((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)*diag(Gauss_w.*k_m(Barrier_xy))*...
        basisPonBarrier((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)'*BarrierSigma(2)*BarrierSigma(2);
    end
end
%A0 = sparse(IMat(:),JMat(:),VA0(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
A0inv = sparse(IMat(:),JMat(:),VA0inv(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
Au11 = sparse(IMat(:),JMat(:),VAu11(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
Au12 = sparse(IMat(:),JMat(:),VAu12(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
Au21 = sparse(IMat(:),JMat(:),VAu21(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
Au22 = sparse(IMat(:),JMat(:),VAu22(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
As11 = sparse(IMat(:),JMat(:),VAs11(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
As12 = sparse(IMat(:),JMat(:),VAs12(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
As21 = sparse(IMat(:),JMat(:),VAs21(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
As22 = sparse(IMat(:),JMat(:),VAs22(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
B1 = sparse(IMat(:),JMat(:),VB1(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
B2 = sparse(IMat(:),JMat(:),VB2(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
C1 = sparse(IMat(:),JMat(:),VC1(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
C2 = sparse(IMat(:),JMat(:),VC2(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
Cp11 = sparse(IMat(:),JMat(:),VCp11(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
Cp12 = sparse(IMat(:),JMat(:),VCp12(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
Cp21 = sparse(IMat(:),JMat(:),VCp21(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
Cp22 = sparse(IMat(:),JMat(:),VCp22(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
Em = sparse(IMat(:),JMat(:),VEm(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
mat_K11 = Au11-(As11*A0inv*Cp11+As12*A0inv*Cp21);
mat_K12 = Au12-(As11*A0inv*Cp12+As12*A0inv*Cp22);
mat_K21 = Au21-(As21*A0inv*Cp11+As22*A0inv*Cp21);
mat_K22 = Au22-(As21*A0inv*Cp12+As22*A0inv*Cp22);
mat_L1  = As11*A0inv*B1+As12*A0inv*B2;
mat_L2  = As21*A0inv*B1+As22*A0inv*B2;
vec_M1  = As11*A0inv*D1+As12*A0inv*D2;
vec_M2  = As21*A0inv*D1+As22*A0inv*D2;
HighOrderP(:) = ( [C1,C2]*([mat_K11,mat_K12;mat_K21,mat_K22]\[mat_L1;mat_L2])+Em )...
                    \ ( A0Q-Nv-Ed-[C1,C2]*([mat_K11,mat_K12;mat_K21,mat_K22]\[vec_M1;vec_M2]));
HighOrderU = ([mat_K11,mat_K12;mat_K21,mat_K22])\([mat_L1;mat_L2]*HighOrderP(:)+[vec_M1;vec_M2]);
HighOrderUx(:) = HighOrderU( 1 : HighOrderNp*size(elements3,2) );
HighOrderUy(:) = HighOrderU( HighOrderNp*size(elements3,2)+1 : end );
% flatten the solution on some particular cells :
%HighOrderP(2:end,BarriersPath)=0;
%HighOrderUx(2:end,FracturesPath)=0;
%HighOrderUy(2:end,FracturesPath)=0;
end
%--------------------------------------------------------------------------
function p = p(xy)
x = xy(:,1);
y = xy(:,2);
p = 1+3*y;
end
function f = f(xy)
x=xy(:,1);
y=xy(:,2);
f = zeros(size(x));
end
function k_m = k_m(xy)
x=xy(:,1);
y=xy(:,2);
k_m = ones(size(x));
end
function q_N = q_N(xy,n)
x=xy(:,1);
y=xy(:,2);
q_N = zeros(size(x));
end
function q = q(xy)
q = -f(xy);
end
%--------------------------------------------------------------------------
end