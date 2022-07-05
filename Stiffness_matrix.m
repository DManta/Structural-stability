function [X, V_n_t, V_1_t, V_2_t, E, S_trial, Se, Qt, Kt] = Stiffness_matrix(d, NEF, NGDL, NGauss_r, NGauss_s, NGauss_t, X_0, V_n_t, V_1_t, V_2_t, V_n_0, J_0_tying, Jacobian_0, J_0, T_bar_E, T_bar_S, S_0, E_0, S_trial)

%Stiffness_matrix
%This script perform the calculations of the elementary tangent stiffness matrix and internal force vector

global uelement h r_tying s_tying Psi_r_tying Psi_s_tying ...
       Psi_x_tying Psi_x_r_tying Psi_x_s_tying Psi_x_r Psi_x_s GDL_element Cct  

% Numerical integration (r,s,t) with a grid of 2x2x2 points
% NGauss_r=2; NGauss_s=2; NGauss_t=2;
Gauss_r=Gauss(NGauss_r); Gauss_s=Gauss(NGauss_s); Gauss_t=Gauss(NGauss_t);

%Pre-alocation 
X = zeros(3,4,NEF);
Kt = zeros(NGDL,NGDL);
Qt = zeros(NGDL,1);
E = zeros(5,NGauss_r*NGauss_s*NGauss_t, NEF);
Xi_D2_E = zeros(20,20,5);
Xi_D_E_tying = zeros(2,20,4);
Xi_D2_E_tying = zeros(20,20,2,4);
E_rt_tying = zeros(1,4); E_st_tying = zeros(1,4);

for el = 1:NEF

    %Restart the matrices for next element
    Ke = zeros(20,20);
    Qe = zeros(20,1);
	
    de=d(nonzeros(uelement(:,el)));
    cont=0;  %counter
    k=0; cont3=[]; %auxiliary
    
    for no = 1:4
	
        %Rotation matrix
		theta = de(cont+4)*V_1_t(:,no,el) + de(cont+5)*V_2_t(:,no,el); %rotation vector
		
        if norm(theta) < 1E-12 %To avoid ill conditioning cases     
            Lambda = eye(3,3);
        else
            Lambda = eye(3,3) + (sin(norm(theta))/norm(theta))*skew(theta) + ((1-cos(norm(theta)))/(norm(theta)^2))*skew(theta)*skew(theta);
        end
        
		%Update the normal vector
        V_1_t(:,no,el) = Lambda*V_1_t(:,no,el); 
        V_2_t(:,no,el) = Lambda*V_2_t(:,no,el);
        V_n_t(:,no,el) = Lambda*V_n_t(:,no,el);
    end
    
 
    %Nós do EF - Referencial Global 
    %Update da posição dos nós no plano médio (x_barra)
    X(:,:,el) = X_0(:,:,el) + [de(1:3),de(6:8),de(11:13),de(16:18)];%d_total
    
    %Variáveis na configuração anterior - tempo t..........................
    V_1_e = V_1_t(:,:,el); %V_1 no tempo t 
    V_2_e = V_2_t(:,:,el); %V_2 no tempo t 
    V_n_e = V_n_t(:,:,el); %V_n no tempo t 
   
    cont = 0;
    cont2 = 0;
    for i=1:NGauss_r
        r=Gauss_r(1,i);

        for ii=1:NGauss_s
            s=Gauss_s(1,ii);
            cont = cont +1;
                
            for iii =1:NGauss_t
                cont2 = cont2 +1;
                t=Gauss_t(1,iii);

                Psi_r = (1/4)*[(s+1)*eye(3,3), (t*h/2)*(s+1)*[-V_2_e(:,1),V_1_e(:,1)], -(s+1)*eye(3,3),(t*h/2)*-(s+1)*[-V_2_e(:,2),V_1_e(:,2)], (s-1)*eye(3,3),(t*h/2)*(s-1)*[-V_2_e(:,3),V_1_e(:,3)], -(s-1)*eye(3,3),(t*h/2)*-(s-1)*[-V_2_e(:,4),V_1_e(:,4)]];
                Psi_s = (1/4)*[(r+1)*eye(3,3), (t*h/2)*(r+1)*[-V_2_e(:,1),V_1_e(:,1)], -(r-1)*eye(3,3),(t*h/2)*-(r-1)*[-V_2_e(:,2),V_1_e(:,2)], (r-1)*eye(3,3),(t*h/2)*(r-1)*[-V_2_e(:,3),V_1_e(:,3)], -(r+1)*eye(3,3),(t*h/2)*-(r+1)*[-V_2_e(:,4),V_1_e(:,4)]];

                %Covariant basis -> time=0
                g_r_0 = J_0(:,1,cont2,el);
                g_s_0 = J_0(:,2,cont2,el);
         
                u_r = Psi_r*de + t*h/2 *[Psi_x_r(:,:,cont)*(V_n_t(1,:,el)-V_n_0(1,:,el))'; Psi_x_r(:,:,cont)*(V_n_t(2,:,el)-V_n_0(2,:,el))'; Psi_x_r(:,:,cont)*(V_n_t(3,:,el)-V_n_0(3,:,el))'];
                u_s = Psi_s*de + t*h/2 *[Psi_x_s(:,:,cont)*(V_n_t(1,:,el)-V_n_0(1,:,el))'; Psi_x_s(:,:,cont)*(V_n_t(2,:,el)-V_n_0(2,:,el))'; Psi_x_s(:,:,cont)*(V_n_t(3,:,el)-V_n_0(3,:,el))'];

                %Deformation
                Xi_D_E = [(g_r_0+u_r)'*Psi_r; (g_s_0+u_s)'*Psi_s; (g_r_0+u_r)'*Psi_s + (g_s_0+u_s)'*Psi_r; zeros(2,20)];
                                  
                Xi_D2_E(:,:,1) = PSI_r(s,t,h,V_n_e, g_r_0 + u_r) + Psi_r'*Psi_r;
                Xi_D2_E(:,:,2) = PSI_r(s,t,h,V_n_e, g_s_0 + u_s) + Psi_s'*Psi_s;
                Xi_D2_E(:,:,3) = PSI_s(r,t,h,V_n_e, g_r_0 + u_r) + PSI_r(s,t,h,V_n_e, g_s_0 + u_s) + Psi_r'*Psi_s + Psi_s'*Psi_r;
 
                %Derivatives
                dx_dr = Psi_x_r(:,:,cont)*(X(1,:,el)+t*h/2*(V_n_t(1,:,el)))';
                dy_dr = Psi_x_r(:,:,cont)*(X(2,:,el)+t*h/2*(V_n_t(2,:,el)))';
                dz_dr = Psi_x_r(:,:,cont)*(X(3,:,el)+t*h/2*(V_n_t(3,:,el)))';
				
                dx_ds = Psi_x_s(:,:,cont)*(X(1,:,el)+t*h/2*(V_n_t(1,:,el)))';
                dy_ds = Psi_x_s(:,:,cont)*(X(2,:,el)+t*h/2*(V_n_t(2,:,el)))';
                dz_ds = Psi_x_s(:,:,cont)*(X(3,:,el)+t*h/2*(V_n_t(3,:,el)))';
                                               
                %Covariant basis -> time=t
                g_r_t = [dx_dr; dy_dr; dz_dr];
                g_s_t = [dx_ds; dy_ds; dz_ds];
        
                %Strain tensor - covariant basis
                E_rr = g_r_t'*g_r_t - g_r_0'*g_r_0;
                E_ss = g_s_t'*g_s_t - g_s_0'*g_s_0;
                E_rs = g_r_t'*g_s_t - g_r_0'*g_s_0;

                E(:,cont2,el) = T_bar_E(:,:,cont2,el)*[E_rr;E_ss;2*E_rs]; %Deformação atual
                
                %Strain increment
                dE = E(:,cont2,el) - E_0(:,cont2,el);
                
                %Check the stress state at the integration point
                [S_trial(:,cont2,el),Cct,Se(cont2,el)] = Constitutive(S_0(:,cont2,el),dE);              

                %Numerical integration
                Qe = Qe + Gauss_r(2,i)*Gauss_s(2,ii)*Gauss_t(2,iii)*Jacobian_0(cont2,el)*(T_bar_E(:,:,cont2,el))'*S_trial(:,cont2,el);
                Ke = Ke + Gauss_r(2,i)*Gauss_s(2,ii)*Gauss_t(2,iii)*Jacobian_0(cont2,el)*((T_bar_E(:,:,cont2,el)*Xi_D_E)'*Cct*(T_bar_E(:,:,cont2,el)*Xi_D_E));               
            end
        end
    end
end
end