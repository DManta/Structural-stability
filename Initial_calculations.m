function [NEF, NDOF, DOF_rotacao, NGauss_r, NGauss_s, NGauss_t, V_n, V_1, V_2, X, J, J_inv, Jacobiano_0, T_bar_E, T_bar_S, J_tying] = Initial_calculations()

name = 'Nodes_C.txt'; %Geometry
table_Nodes = dlmread(name); clear name;

name = 'Link_EF_C.txt'; %Input EF e incidences
table_EF = dlmread(name); clear name;

NEF = size(table_EF,1)/2; %Nº of elements

Nodes = table_Nodes(1:3:end-2,1:4);
NNodes = size(Nodes,1);

%Pre-alocation
NDOF = 0; %NDOF total
DOF_Nodes = zeros(NNodes,1); %Nodal DOF
Incid = zeros(6,NNodes); %Incidences
k=0;

%For drawing purposes .......................................
name = 'Nodes_C_INDEF.txt';  %Geometry
table_Nodes_indef = dlmread(name); clear name;
Nodes_indef = table_Nodes_indef(1:3:end-2,1:4);

name = 'Link_EF_C_INDEF.txt'; 
table_EF_indef = dlmread(name); clear name; %Geometry
NEF_indef = size(table_EF_indef,1)/2;

%Incidences...................................................
NoIncid = table_EF(2:2:2*NEF,1:4);
NoIncid_indef = table_EF_indef(2:2:2*NEF_indef,1:4);

%Materials.................................................................
name = 'Material_C.txt';
table = dlmread(name); clear name;
h=table(2,4); %Thickness

% Numerical integration (r,s,t) with a grid of 2x2x2 points
% NGauss_r=2; NGauss_s=2; NGauss_t=2;

NGauss_r=2; NGauss_s=2; NGauss_t=5; 
Gauss_r=Gauss(NGauss_r); Gauss_s=Gauss(NGauss_s); Gauss_t=Gauss(NGauss_t);

[Psi_x, Psi_x_r, Psi_x_s] = deal(zeros(1,4,NGauss_r*NGauss_s));

cont = 0;
for i=1:NGauss_r
    r=Gauss_r(1,i);

    for ii=1:NGauss_s
        s=Gauss_s(1,ii); 
        cont = cont +1;

        Psi_x(:,:,cont)=1/4*[(1+r)*(1+s),(1-r)*(1+s),(1-r)*(1-s),(1+r)*(1-s)];
        Psi_x_r(:,:,cont)=1/4*[(s+1),-(s+1),(s-1),-(s-1)]; 
        Psi_x_s(:,:,cont)=1/4*[(r+1),-(r-1),(r-1),-(r+1)];         
    end
end

%..........................................................................

%Pre-alocation
uelement=zeros(24,NEF);
DOF_element=zeros(4,NEF);
[V_n, V_1, V_2, X]  =  deal(zeros(3,4,NEF));
[J, J_inv] = deal(zeros(3,3,NGauss_r*NGauss_s*NGauss_t,NEF)); 
Jacobian_0=zeros(NGauss_r*NGauss_s*NGauss_t,NEF);
T_bar_E = zeros(5,5,NGauss_r*NGauss_s*NGauss_t,NEF);
T_bar_S = zeros(5,5,NGauss_r*NGauss_s*NGauss_t,NEF);

%Calculations for each element
for el = 1:NEF
	
	%Nodes
    No1=NoIncid(el,1); No2=NoIncid(el,2); No3=NoIncid(el,3); No4=NoIncid(el,4);

    %Coordinates
    X1=Nodes(Nodes(:,1)==No1,2); X2=Nodes(Nodes(:,1)==No2,2); X3=Nodes(Nodes(:,1)==No3,2); X4=Nodes(Nodes(:,1)==No4,2);
    Y1=Nodes(Nodes(:,1)==No1,3); Y2=Nodes(Nodes(:,1)==No2,3); Y3=Nodes(Nodes(:,1)==No3,3); Y4=Nodes(Nodes(:,1)==No4,3);
    Z1=Nodes(Nodes(:,1)==No1,4); Z2=Nodes(Nodes(:,1)==No2,4); Z3=Nodes(Nodes(:,1)==No3,4); Z4=Nodes(Nodes(:,1)==No4,4);
    
    %Nodes
    X(:,:,el)=[X1, X2, X3, X4; Y1, Y2, Y3, Y4; Z1, Z2, Z3, Z4];
    
    %Normal vector
    V_n(:,1,el) = cross(X(:,2,el)-X(:,1,el),X(:,4,el)-X(:,1,el))/norm(cross(X(:,2,el)-X(:,1,el),X(:,4,el)-X(:,1,el)));
    V_n(:,2,el) = cross(X(:,3,el)-X(:,2,el),X(:,1,el)-X(:,2,el))/norm(cross(X(:,3,el)-X(:,2,el),X(:,1,el)-X(:,2,el)));
    V_n(:,3,el) = cross(X(:,4,el)-X(:,3,el),X(:,2,el)-X(:,3,el))/norm(cross(X(:,4,el)-X(:,3,el),X(:,2,el)-X(:,3,el)));
    V_n(:,4,el) = cross(X(:,1,el)-X(:,4,el),X(:,3,el)-X(:,4,el))/norm(cross(X(:,1,el)-X(:,4,el),X(:,3,el)-X(:,4,el)));
  
    %auxiliary
    cont=0; cont2 = 0; aux = [No1,No2,No3,No4];
    
    for i=1:NGauss_r
	
         for ii=1:NGauss_s
		 
            cont = cont +1;
 
            dx_dt = Psi_x(:,:,cont)*h/2*V_n(1,:,el)';
            dy_dt = Psi_x(:,:,cont)*h/2*V_n(2,:,el)';
            dz_dt = Psi_x(:,:,cont)*h/2*V_n(3,:,el)';
                
            for iii =1:NGauss_t
			
                cont2 = cont2 +1;
				
                t=Gauss_t(1,iii);
                
                %Derivatives
                dx_dr = Psi_x_r(:,:,cont)*(X(1,:,el)+t*h/2*V_n(1,:,el))';
                dy_dr = Psi_x_r(:,:,cont)*(X(2,:,el)+t*h/2*V_n(2,:,el))';
                dz_dr = Psi_x_r(:,:,cont)*(X(3,:,el)+t*h/2*V_n(3,:,el))';

                dx_ds = Psi_x_s(:,:,cont)*(X(1,:,el)+t*h/2*V_n(1,:,el))';
                dy_ds = Psi_x_s(:,:,cont)*(X(2,:,el)+t*h/2*V_n(2,:,el))';
                dz_ds = Psi_x_s(:,:,cont)*(X(3,:,el)+t*h/2*V_n(3,:,el))';
                                
                %Jacobian matrix
                J(:,:,cont2,el) = [dx_dr, dx_ds, dx_dt; dy_dr, dy_ds, dy_dt; dz_dr, dz_ds, dz_dt];
                
                %Jacobian = det(Jacobian)
                Jacobian_0(cont2,el) = det(J(:,:,cont2,el)); 
                
                %Inverse of Jacobian matrix           
                J_inv(:,:,cont2,el) = inv(J(:,:,cont2,el));
                
                %Covariant basis
                g_s_0 = J(:,2,cont2,el);
                g_t_0 = J(:,3,cont2,el);

                %Transformation matrices
                e_3_0 = g_t_0/norm(g_t_0);
                e_1_0 = cross(g_s_0,e_3_0)/norm(cross(g_s_0,e_3_0));
                e_2_0 = cross(e_3_0,e_1_0);
                              
                T_bar_E(:,:,cont2,el) = T_bar([e_1_0,e_2_0,e_3_0]'*J_inv(:,:,cont2,el)');
                T_bar_S(:,:,cont2,el) = T_bar(J(:,:,cont2,el)'*[e_1_0,e_2_0,e_3_0]);
            end
        end
    end
end
end

%Local function...........................................................

function T_b = T_bar(T)
    T_b=[T(1,1)^2,T(1,2)^2,T(1,1)*T(1,2),T(1,1)*T(1,3),T(1,2)*T(1,3); 
           T(2,1)^2,T(2,2)^2,T(2,1)*T(2,2),T(2,1)*T(2,3),T(2,2)*T(2,3);
           2*T(1,1)*T(2,1),2*T(1,2)*T(2,2),T(1,1)*T(2,2)+T(2,1)*T(1,2),T(1,1)*T(2,3)+T(2,1)*T(1,3),T(1,2)*T(2,3)+T(2,2)*T(1,3);  
           2*T(1,1)*T(3,1),2*T(1,2)*T(3,2),T(1,1)*T(3,2)+T(3,1)*T(1,2),T(1,1)*T(3,3)+T(3,1)*T(1,3),T(1,2)*T(3,3)+T(3,2)*T(1,3);
           2*T(2,1)*T(3,1),2*T(2,2)*T(3,2),T(2,1)*T(3,2)+T(3,1)*T(2,2),T(2,1)*T(3,3)+T(3,1)*T(2,3),T(2,2)*T(3,3)+T(3,2)*T(2,3)];
end