function[disp_t, vel_t, ac_t, dt, nt] = Eqsolver(nn,beta,gamma,K_global, F_global, M_global,big_lamda,Ri,t,t_end,ct)
dt = ct*2/sqrt(big_lamda); 

nt = ceil(t_end/dt);

def = zeros(nn,1);
vel = zeros(nn,1);
ac = zeros(nn,1);
disp_t = zeros(nn,nt+1);
vel_t = zeros(nn,nt+1);
ac_t = zeros(nn,nt+1);

%Start time loop
for j = 1:nt
    tn = j*dt;
    % predictor phase
    
        def1 = def + dt*vel+0.5*(dt^2)*(1-2*beta)*ac;
        def2(:,j)= def1;
        v1= vel+(1-gamma)*dt*ac;
        % solution phase
        A = (M_global + beta*(dt^2)*K_global);
        P = [2*pi*pressure(tn,t)*Ri];
        F_global(1) = P;
        F_star = F_global - K_global*def1;
        F_new(:,j) = F_global;
        ac = A\F_star;
        % correction phase
        def = def1 + beta*dt*dt*ac;
        vel = v1 + gamma*dt*ac;
        % final values of displacement, acceleration and velocity
        disp_t(:,j+1) = def;
        vel_t(:,j+1) = vel;
        ac_t(:,j+1) = ac;

end
end