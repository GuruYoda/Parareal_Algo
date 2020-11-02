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
% 's1'
% y0 = [2;1;-4];
% 's2'
% y0 = [4;11;-299];
% a_mat = -1i-0.8; %-0.9;
% y_exact = @(t)(y0(1)*exp(a_mat.*t));
t_min = 0;
t_max = 14;
scale_mesh = 10;

m = 1; %size(A,1);

%% coarse grid solution

n_sub = 8;
n_coarse = 4;
dt_c = (t_max - t_min)/(n_sub*n_coarse);
t_c = t_min:dt_c:t_max;
n_c = length(t_c);

[y_c] = solve_ivp(@ode_operator, m, t_c, y0);
s_mat = zeros(m,n_sub+1);
s_mat(m,1:n_sub) = y_c(m,1:n_coarse:n_c-1);
s_mat(m,n_sub+1) = y_c(m,n_sub*n_coarse+1);

s_mat_prev = zeros(m,n_sub+1);
s_mat_prev(m,1:n_sub) = y_c(m,1:n_coarse:n_c-1);
s_mat_prev(m,n_sub+1) = y_c(m,n_sub*n_coarse+1);

s_mat_new = zeros(m,n_sub+1);
s_mat_new(m,1:n_sub) = y_c(m,1:n_coarse:n_c-1);
s_mat_new(m,n_sub+1) = y_c(m,n_sub*n_coarse+1);

%% fine grid solution
n_fine = scale_mesh*n_coarse; %
n_f = n_sub*(n_fine) + 1;
dt = (t_c(n_coarse+1)-t_c(1))/(n_fine);
t = zeros(1,n_f);
for sub=1:n_sub
    t((sub-1)*n_fine+1:(sub)*n_fine+1) = t_c((sub-1)*n_coarse+1):dt:t_c((sub)*n_coarse+1);
end
y = zeros(m,n_f);

tol = 1e-5;
error = 10*tol;
iter = 1;
max_iter  = 1000;
err_vec = zeros(max_iter,1);
y_all = zeros(max_iter,m,n_f);

while (error > tol && iter < max_iter)
    
    
    for sub=1:n_sub
        % coarse grid solution
        t_sub_c = t_c((sub-1)*n_coarse+1:(sub)*n_coarse+1);
        [y_sub_c] = solve_ivp(@ode_operator, m, ...
            t_sub_c, s_mat(:,sub));
        %y(:,(sub-1)*n_fine+1:(sub)*n_fine+1) = y_sub;
        
        s_mat_new(sub+1) = y_sub_c(n_coarse+1);
        
        % fine grid solution
        t_sub_f = t((sub-1)*n_fine+1:(sub)*n_fine+1);
        [y_sub] = solve_ivp(@ode_operator, m, ...
            t_sub_f, s_mat(:,sub));
        y(:,(sub-1)*n_fine+1:(sub)*n_fine+1) = y_sub;
    end
    
    error = norm(y(:,n_fine+1:n_fine:n_f) - s_mat(:,2:n_sub+1),2);
    err_vec(iter) = error;
    y_all(iter,:,:) = y;
    
    % update initial value
    s_mat = y(:,1:n_fine:n_f);
    s_mat(2:n_sub+1) = s_mat(2:n_sub+1) + s_mat_new(2:n_sub+1) ...
                        - s_mat_prev(2:n_sub+1);
    s_mat_prev = s_mat_new;
    iter = iter + 1;
end

n_iters = iter -1;
err_vec(n_iters+1:max_iter) = [];
y_all(n_iters+1:max_iter,:,:) = [];

%% postprocessing

% yE = y_exact(t);
% for ploting only, we will interpolate the coarse grid solution
y_c_ = zeros(m,n_f);
for sub=1:n_sub
    y_c_((sub-1)*n_fine+1:(sub)*n_fine+1) = interp1(t_c((sub-1)*n_coarse+1:(sub)*n_coarse+1),...
                                    y_c((sub-1)*n_coarse+1:(sub)*n_coarse+1),t((sub-1)*n_fine+1:(sub)*n_fine+1));
end
%% plotting

% solution
F(1:n_iters+1) = struct('cdata',[],'colormap',[]);
fig = figure(1);
grid on;
hold on;
for k=1:m
    h(k) = plot(t,real(y_c_(k,:)),'r');
%     g(k) = plot(t,real(yE(k,:)),'b');
end
hold off;
xlabel('$t$','Interpreter','latex');
ylabel('$y$','Interpreter','latex');
title('$y$','Interpreter','latex');
% ylim([-2*abs(min(real(y_c_))),abs(y0)]);
F(1) = getframe(gcf);
pause(1);

% for k=1:m
%     set(h(k),'XData',t, 'YData',y_all(1,k,:));
% end
% F(2) = getframe(gca);
for a=1:n_iters
    hold on;
    for k=1:m
        set(h(k), 'YData', real(y_all(a,k,:)));
%         set(g(k), 'YData', real(yE(k,:)));
    end
    hold off;
    F(a+1) = getframe(gcf);
    drawnow;
    pause(1);
end

movie(F);

filename = sprintf('%s_solution',prob_n);
print(fig,filename,'-dpng');
file_name = sprintf('%s_solution.avi',prob_n);
v = VideoWriter(file_name);
open(v);
writeVideo(v,F);
close(v);

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
function [A] = ode_operator(t)
% evaluates A at t
% user defined system of ode
% t is one scalar
% vectorize for t as a vector
% n = 3;
% A = zeros(n,n);

% global a_mat;
% 's0'
A = [sin(t)];

% problem 2
% 's1'
% A = [0 1 0;
%     0 0 1;
%     6/t^3 -6/t^2 3/t];

% 's2'
% A = [0 1 0;
%     0 0 1;
%     100 -100 1];
end

% discrete ode solver

function [y] = solve_ivp(operator,m, t, y0)
% explicit euler method
n_ = length(t);
dt = t(2) - t(1);
y = zeros(m,n_);
y(:,1) = y0;
for k=2:n_
%     A = feval(operator, t(k-1));
%     y(:,k) = (1 - dt*A)\y(:,k-1);
    A = sin(t(k-1));
    y(:,k) = (1+ dt*A)*y(:,k-1) + dt*t(k-1);
end
end

