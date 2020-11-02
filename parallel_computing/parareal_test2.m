%% parareal algorithm
% solving a system of ODEs, parallel in time
% y'(t) = A*y(t)
% y(t0) = y0

clear; close all;
%% problem definition

% define the system of ode in ode_operator
% f = @ode_operator;
% A = feval(f,t);

prob_n = 's2';

% global a_mat;
% 's0'
y0 = [1];
t_min = 0;
t_max = 14;
n_sub = 20;
n_coarse = 10;
scale_mesh = 10;

m = 1; %size(A,1);

tol = 1e-10;
error = 10*tol;
max_iter  = 1000;

%% coarse grid solution

dt_c = (t_max - t_min)/(n_sub*n_coarse);
t_c = t_min:dt_c:t_max;
n_c = length(t_c);

[y_c] = solve_ivp(m, t_c, y0);

t_sub = zeros(1,n_sub+1);
t_sub(1:n_sub) = t_c(1:n_coarse:n_c-1);
t_sub(n_sub+1) = t_c(n_sub*n_coarse+1);

s_mat = zeros(m,n_sub+1);
s_mat(m,1:n_sub) = y_c(m,1:n_coarse:n_c-1);
s_mat(m,n_sub+1) = y_c(m,n_sub*n_coarse+1);

s_mat_prev = zeros(m,n_sub+1);
s_mat_prev(m,1:n_sub) = y_c(m,1:n_coarse:n_c-1);
s_mat_prev(m,n_sub+1) = y_c(m,n_sub*n_coarse+1);

% s_mat_new = zeros(m,n_sub+1);
s_mat_new = s_mat_prev;

s_mat_all = zeros(max_iter+1,m,n_sub+1);
s_mat_fine_all = zeros(max_iter,m,n_sub+1);
s_mat_coarse_all = zeros(max_iter+1,m,n_sub+1);

s_mat_all(1,:,:) = s_mat;
s_mat_coarse_all(1,:,:) = s_mat_prev;

%% fine grid solution
n_fine = scale_mesh*n_coarse; %
n_f = n_sub*(n_fine) + 1;
dt = (t_c(n_coarse+1)-t_c(1))/(n_fine);
t = zeros(1,n_f);
for sub=1:n_sub
    t((sub-1)*n_fine+1:(sub)*n_fine+1) = t_c((sub-1)*n_coarse+1):dt:t_c((sub)*n_coarse+1);
end
y = zeros(m,n_f);

err_vec = zeros(max_iter,1);
y_all = zeros(max_iter,m,n_f);

iter = 1;
while (error > tol && iter < max_iter)
    
    % fine grid solution    
    for sub=1:n_sub
        t_sub_f = t((sub-1)*n_fine+1:(sub)*n_fine+1);
        [y_sub] = solve_ivp(m, ...
            t_sub_f, s_mat_prev(:,sub));
        y(:,(sub-1)*n_fine+1:(sub)*n_fine+1) = y_sub;
    end
    s_mat = y(:,1:n_fine:n_f);
    s_mat_fine_all(iter,:,:) = s_mat;

    % coarse grid solution
    for sub=1:n_sub
        t_sub_c = t_c((sub-1)*n_coarse+1:(sub)*n_coarse+1);
        [y_sub_c] = solve_ivp(m, ...
            t_sub_c, s_mat(:,sub));
        y_iter_c(:,(sub-1)*n_coarse+1:(sub)*n_coarse+1) = y_sub_c;
        
        s_mat_new(sub+1) = y_sub_c(n_coarse+1);
    end
    s_mat_coarse_all(iter+1,:,:) = s_mat_new;
    y_all(iter,:,:) = y;
    
    % update initial value
    s_mat(2:n_sub+1) = s_mat(2:n_sub+1) + s_mat_new(2:n_sub+1) ...
                        - s_mat_prev(2:n_sub+1);
    s_mat_prev = s_mat_new;
    s_mat_all(iter+1,:,:) = s_mat;
    
    error = norm(y(:,n_fine+1:n_fine:n_f) - s_mat(:,2:n_sub+1),2);
    err_vec(iter) = error;
    iter = iter + 1;
end

n_iters = iter -1;
err_vec(n_iters+1:max_iter) = [];
y_all(n_iters+1:max_iter,:,:) = [];

s_mat_all(n_iters+2:max_iter+1,:,:) = [];
s_mat_fine_all(n_iters+1:max_iter,:,:) = [];
s_mat_coarse_all(n_iters+2:max_iter+1,:,:) = [];

%% postprocessing

% solution
F(1:n_iters+1) = struct('cdata',[],'colormap',[]);
fig = figure(1);
grid on;
hold on;
for k=1:m
    p(k) = plot(t_sub,real(squeeze(s_mat_coarse_all(1,k,:))),'b','Marker','o');
    q(k) = plot(t_sub,real(squeeze(s_mat_coarse_all(1,k,:))),'r','Marker','+');
    r(k) = plot(t_sub,real(squeeze(s_mat_all(1,k,:))),'k','Marker','^');
    legend('Coarse', 'Fine', 'Corrected');
end
xlabel('$t$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
title('$y$','Interpreter','latex');
ylim([0,160]);
F(1) = getframe(gca);
pause(0.5);

% for k=1:m
%     set(h(k),'XData',t, 'YData',y_all(1,k,:));
% end
% F(2) = getframe(gca);
for a=1:n_iters
    for k=1:m
        set(p(k), 'YData', real(squeeze(s_mat_coarse_all(a+1,k,:))));
        set(q(k), 'YData', real(squeeze(s_mat_fine_all(a,k,:))));
        set(r(k), 'YData', real(squeeze(s_mat_all(a+1,k,:))));
        legend('Coarse', 'Fine', 'Corrected');
    end
    F(a+1) = getframe(gca);
    drawnow;
    pause(0.5);
end

filename = sprintf('%s_solution',prob_n);
print(fig,filename,'-dpng');

movie(F);
file_name = sprintf('%s_solution.avi',prob_n);
v = VideoWriter(file_name);
open(v);
writeVideo(v,F);
close(v);

% % progressive solution
% do it

% error plot
fig = figure(3);
plot(1:n_iters, err_vec);
grid on;
xlabel('Iterations \rightarrow');
ylabel('error');
title('global error at course grid');
filename = sprintf('%s_error',prob_n);
print(fig,filename,'-dpng');

%% helper functions

% ode operator
% function [A] = ode_operator(t)
% % evaluates A at t
% % user defined system of ode
% % t is one scalar
% % vectorize for t as a vector
% % n = 3;
% % A = zeros(n,n);
% 
% global a_mat;
% % 'test1'
% A = [a_mat];
% 
% end

% discrete ode solver

function [y] = solve_ivp(m, t, y0)
% explicit euler method
n_ = length(t);
dt = t(2) - t(1);
y = zeros(m,n_);
y(:,1) = y0;
for j=2:n_
%     A = feval(operator, t(j-1));
    %A = sin(t(k-1));
    y(:,j) = (1+ dt*sin(t(j-1)))*y(:,j-1) + dt*t(j-1);
    %% implicit method
%     y(:,k) = (1 - dt*A)\y(:,k-1);
    %% explicit method
%     y(:,j) = (1+ dt*A)*y(:,j-1);
    %% classical RK2
    k1 = sin(t(j-1))*y(:,j-1) + t(j-1);
    k2 = sin(t(j-1) + dt/2)*(y(:,j-1) + dt*k1/2) + t(j-1) + dt/2;
    k = (k1 + k2)/2;
    y(:,j) = y(:,j-1) + dt*k;
end
end
