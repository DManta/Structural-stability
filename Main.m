clear; clc; close all; format long;

%Main program
%This script perform the initial calculations and call the necessaries sub-scripts to perform the analysis

%Perform the initial calculations, before the analysis begin
[NEF, NDOF, DOF_rotacao, NGauss_r, NGauss_s, NGauss_t, V_n_0, V_1_0, V_2_0, X_0, J_0, J_inv_0, Jacobian_0, T_bar_E, T_bar_S, J_0_tying] = Initial_calculations;

%Start from the initial load step -> t=0
V_n_t = V_n_0; %Normal vector, at each node
V_1_t = V_1_0; 
V_2_t = V_2_0; 
X_t = X_0;     %Coordinates of each node 
J_t = J_0;     %Jacobian_time t
J_t_tying = J_0_tying; %Jacobian_tying points_time t

%"time" = load increment factor
time=-1*[0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 , 110, 120, 130, 140, 150];

nincrements=length(time);

global Nodes Incid
DOF_support = sort(nonzeros(Incid(:, Nodes(:,2)==0))); %Cantilever - fixed at x=0
DOF_actives = setdiff(1:NDOF,DOF_support);

tol = 1E-4; %tolerance
iter_max = 20; %20 iterations max (check if more needed)

Results = zeros(nincrements,2); %U_z,time=F

%Pre-alocation
d = zeros(NDOF,1); %displacements 
dt = zeros(NDOF,1);
Lambda=0; %Load parameter
dLambda=0; 
S_0 = zeros(5,NGauss_r*NGauss_s*NGauss_t, NEF); %Initial stresses -> t =0
E_0 = zeros(5,NGauss_r*NGauss_s*NGauss_t, NEF); %Initial strains -> t =0
Gt_ext=zeros(NDOF+1,1); %Vetor of residuals 

DOF_mon = Incid(3,304); %Monitored Dof
one = zeros(NDOF,1);
one(DOF_mon) = 1;

%Start of the iterative procedure to calculate the non-linear equilibrium path -> Newton-Raphson scheme

for increment=2:nincrements
   
   %Start from zero
    iter=0; %iteration number

    S_trial = zeros(5,NGauss_r*NGauss_s*NGauss_t, NEF); %Stress trial
    Se = zeros(NGauss_r*NGauss_s*NGauss_t, NEF); %Effective stress

    
    %Update Kt, Qt.....................................................        
    [X_t, V_n_t, V_1_t, V_2_t, E, S_trial, Se, Qt, Kt] = Stiffness_matrix(d, NEF, NDOF, NGauss_r, NGauss_s, NGauss_t, X_0, V_n_t, V_1_t, V_2_t, V_n_0, J_0_tying, Jacobian_0, J_0, T_bar_E, T_bar_S, S_0, E_0, S_trial);
    
    error = tol + 1; %To start the iterative procedure

    while error > tol && iter < iter_max %iter_max, to avoid infinite loops
       
        iter = iter+1;
        
        tic %check the duration/iteration (performance)
               
        %Update Kt, Qt.....................................................        
        K_ext = Kt(DOF_actives,DOF_actives);
        
        Gt_ext = [Lambda*one(DOF_actives) - Qt(DOF_actives) ; time(increment)];
        
        d_ext= K_ext\Gt_ext; %Update/increment
        
        dt(DOF_actives) = d_ext(1:end-1); dLambda = d_ext(end); %Update/increment
        
        d(DOF_actives)=dt(DOF_actives); %Update/increment 
        Lambda=dLambda; %Update/increment

        %Update Kt, Qt.....................................................        
        [X_t, V_n_t, V_1_t, V_2_t, E, S_trial, Se, Qt, Kt] = Stiffness_matrix(d, NEF, NDOF, NGauss_r, NGauss_s, NGauss_t, X_0, V_n_t, V_1_t, V_2_t, V_n_0, J_0_tying, Jacobian_0, J_0, T_bar_E, T_bar_S, S_0, E_0, S_trial);
 
        error = norm(Gt_ext);
        
        if increment > 1 
            %Display the results - for control purposes
            disp(['inc: ',num2str(increment), ' | iter: ',num2str(iter), ' | displacement: ', num2str(time(increment)) , ' | Force: ', num2str(-Lambda) ,' | error: ', num2str(error)]);            
        end
        toc
    end
    
    %Convergence achieved
    S_0 = S_trial; %Save stresses
    E_0 = E; %Save strains

    %Save the results
    Results(increment,1) = -time(increment); %displacement
    save(['Results_', num2str(increment),'.mat'])
end