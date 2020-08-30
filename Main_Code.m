%Matrix of Elastic Constants
E = 1000; %psi
v = 0.3;  %poisson's ratio
rho=0.01; %density of the material
D = (1/((1+v)*(1-2*v)))*[E*(1-v), v*E; v*E, E*(1-v)]; %Matrix of Elastic Constants

% 4 pt gaussian integration data 
% weight = [0.652145154862546,0.652145154862546,0.347854845137453,0.347854845137453];
% intpt = [-0.339981043584856,0.339981043584856,-0.861136311594052,0.861136311594052];

%% Time based Parameters
t = 0.05;
t_end = 0.5; %end of time loop
% parameters for variation in delta-t
ct1 = 0.8;
ct2 = 1;
ct3 = 1.2;

%% Newmark Method Parameters
%parameters for central difference method
a1_beta = 0; a1_gamma = 0.5;  
%parameters for average acceleration method
a2_beta = 0.25; a2_gamma = 0.5;

%% Mesh Generation
Ri = 6; %inch
Ro = 9; %inch
n = 200;
nn = n+1;
np = linspace(Ri, Ro, n+1);
le = (Ro - Ri)/n; %length of each element
x1 = zeros(1,n); %coordinates of first node for each element
x2 = zeros(1,n); %coordinates of seconde node for each element
for i = 1:n
    x1(1,i) = x1(1,i) + 6 + (i-1)*le;
    x2(1,i) = x2(1,i) + (6+le) + (i-1)*le;
end
element_coord = [x1' x2'];

%% Stiffness Matrix, Force Vector, Mass Matrix generation
K_global = zeros(n+1,n+1);
F_global = zeros(n+1,1);
M_global = zeros(n+1, n+1);
eigen_global = zeros(n,1);
beta = [0 0.25];
gamma = [0.5 0.5];
for i = 1:n
    syms eps;
    N1 = 0.5*(1-eps);
    N2 = 0.5*(1+eps);
    r = N1*x1(1,i) + N2*x2(1,i);
    B = [-1/le  1/le; N1/r N2/r];
    L = @(eps) B'*D*B*r*(le/2);
    K_local = 2*pi*vpaintegral(L,eps,[-1 1]);
    M_local_1 = @(eps)[N1, N2]'*[N1 N2].*r.*((x2(1,i)-x1(1,i))/2);
    cons_mass_mat = 2*pi*rho*vpaintegral(M_local_1,eps,[-1 1]);
    for j = 1:2
        if abs(beta(1,j) - 0.25) <= 10^-10
            Mass_Matrix = cons_mass_mat;
        else
            M = diag(cons_mass_mat);
            M(1) = M(1) + cons_mass_mat(1,2);
            M(2) = M(2) + cons_mass_mat(2,1);
            Mass_Matrix = diag(M);
        end
    end 
    M_local = Mass_Matrix;
    lamda = eig(M_local\K_local);
    max_lamda = max(lamda);
    eigen_global(i,1) = max_lamda;
    id1 = i; id2 = i+1;
    K_global(id1:id2, id1:id2) = K_global(id1:id2,id1:id2) + K_local;
    M_global(id1:id2, id1:id2) = M_global(id1:id2,id1:id2) + M_local;
end

%% critical time step & evaluation of displacement, velocity and acceleration
big_lamda = max(eigen_global(:,1));

% for the given elements
% central difference calculation
[disp_t_n, vel_t_n, ac_t_n, dt_1,nt_1] = Eqsolver(nn,a1_beta,a1_gamma,K_global, F_global, M_global,big_lamda,Ri,t,t_end,ct1);
[disp_t_n_1, vel_t_n_1, ac_t_n_1, dt_2,nt_2] = Eqsolver(nn,a1_beta,a1_gamma,K_global, F_global, M_global,big_lamda,Ri,t,t_end,ct2);
[disp_t_n_2, vel_t_n_2, ac_t_n_2, dt_3, nt_3] = Eqsolver(nn,a1_beta,a1_gamma,K_global, F_global, M_global,big_lamda,Ri,t,t_end,ct3);

% average acceleration calculation
[disp_t_n_aa, vel_t_n_aa, ac_t_n_aa, dt_aa_1,nt_aa_1] = Eqsolver(nn,a2_beta,a2_gamma,K_global, F_global, M_global,big_lamda,Ri,t,t_end,ct1);
[disp_t_n_aa_1, vel_t_n_aa_1, ac_t_n_aa_1, dt_aa_2,nt_aa_2] = Eqsolver(nn,a2_beta,a2_gamma,K_global, F_global, M_global,big_lamda,Ri,t,t_end,ct2);
[disp_t_n_aa_2, vel_t_n_aa_2, ac_t_n_aa_2, dt_aa_3, nt_aa_3] = Eqsolver(nn,a2_beta,a2_gamma,K_global, F_global, M_global,big_lamda,Ri,t,t_end,ct3);

figure(1)
plot(disp_t_n_1(3,:))
hold on
plot(disp_t_n_aa_1(3,:))
xlabel('Time step')
ylabel('Displacement at node 3 in inches')
legend('Central Difference', 'Average Acceleration')
hold off

%% Radial and Tangential Stress Calculation
Be = [-1/le 1/le];
du_dr_cd = zeros(n,nt_1+1);
du_dt_cd = zeros(n,nt_1+1);
sigma_rr_cd = zeros(n,nt_1+1);
sigma_tt_cd = zeros(n,nt_1+1);
r1_cd = zeros(n,nt_1+1);
u_middle_cd = zeros(n,nt_1+1);

for i = 1: nt_1
    for j = 1:n
        u_e_cd = disp_t_n(j:j+1,i);
        du_dr_cd(j,i) = Be*u_e_cd;
        u_middle_cd(j,i) = 0.5*(disp_t_n(j) + disp_t_n(j+1));
        r1_cd(j,i) = x1(1,j)+ u_middle_cd(j,i) + le/2;
        du_dt_cd(j,i) = u_middle_cd(j,i)/r1_cd(j,i);
        sigma_rr_cd(j,i) = (E/((1+v)*(1-2*v)))*((1-v)*du_dr_cd(j,i) + v*du_dt_cd(j,i));
        sigma_tt_cd(j,i) = (E/((1+v)*(1-2*v)))*(v*du_dr_cd(j,i) + (1-v)*du_dt_cd(j,i));
    end
end


du_dr_aa = zeros(n,nt_aa_1+1);
du_dt_aa = zeros(n,nt_aa_1+1);
sigma_rr_aa = zeros(n,nt_aa_1+1);
sigma_tt_aa = zeros(n,nt_aa_1+1);
r1_aa = zeros(n,nt_aa_1+1);
u_middle_aa = zeros(n,nt_aa_1+1);

for i = 1: nt_aa_1
    for j = 1:n
        u_e_aa = disp_t_n_aa(j:j+1,i);
        du_dr_aa(j,i) = Be*u_e_aa;
        u_middle_aa(j,i) = 0.5*(disp_t_n_aa(j) + disp_t_n_aa(j+1));
        r1_aa(j,i) = x1(1,j)+ u_middle_aa(j,i) + le/2;
        du_dt_aa(j,i) = u_middle_aa(j,i)/r1_aa(j,i);
        sigma_rr_aa(j,i) = (E/((1+v)*(1-2*v)))*((1-v)*du_dr_aa(j,i) + v*du_dt_aa(j,i));
        sigma_tt_aa(j,i) = (E/((1+v)*(1-2*v)))*(v*du_dr_aa(j,i) + (1-v)*du_dt_aa(j,i));
    end
end


figure(2)
plot(sigma_rr_cd(3,:))
hold on
plot(sigma_rr_aa(3,:))
xlabel('Time step')
ylabel('Radial stress in psi at element 3')
legend('Central Difference Method','Average Acceleration Method')

figure(3)
% plot(sigma_tt_cd(3,:))
% hold on
plot(sigma_tt_aa(3,:))
xlabel('Time step')
ylabel('Tangential stress in psi at element 3')
legend('Central Difference Method','Average Acceleration Method')



    
    



