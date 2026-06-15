close all
clear
clc

%% Input
global a chi num_elements T cfl

a = 1;
chi = 0.05;
num_elements = 200;

T = 100;
cfl = 1.22;

flag_ini = 2; 
flag_tab = 3; % flag for the coefficients of the method, Table 2, Table 3

%% Core Part

g  = 9.81;
L = 500;
x   = linspace (0, L, num_elements+1) .';
x_cell = .5*(x(2:end) + x(1:end-1));
dx  = diff (x);

mass = [dx(1)/2;dx(2:end);dx(end)/2];


q0 = @(x) form_initial_cond(L, x, flag_ini);
q = q0(x);

exact_sol = @(x,t) q0(x-a*t).*exp(-chi*t);

figure()
plot(x, q, '--o')
ylim([0, max(q0(x))])

if (flag_ini==1)
    dt = .01;
else
    dt = time_step_selector(q, cfl, dx);
end

[a21, a32, a31, atilde32, atilde21, atilde31, gamma] = set_coefficients(flag_tab);

norm_inf = 0;

tic
figure()
time = 0;
plot(x, exact_sol(x,time), 'r', LineWidth=2)
hold on
time_vec = [time];
q_vec = [mean(q)];
norm_inf_vec = [];
error = 0;
while (time < T)
    % 1st step,
    F = flux(q);
    S = source(q);

    q_cell = .5*(q(2:end,:)+q(1:end-1,:));
    residual_2 = @(q_2) q_cell - a21.*dt.*(F(2:end,:) - F(1:end-1,:))./dx - atilde21.*dt/2*(S(2:end,:) + S(1:end-1,:)) - gamma*dt*source(q_2) - q_2;
    q_2 = fsolve(residual_2, q_cell);

    % 2nd step, low order sol,
    F_l = [flux(q(1,:)); F; flux(q(end,:))];
    F_m = [flux(q_2(1,:)); flux(q_2); flux(q_2(end,:))];
    F_ms = [source(q_2(1,:)); source(q_2); source(q_2(end,:))];
    residual_3 = @(q_3) q - a31/2*dt./mass .* (F_l(3:end,:)-F_l(1:end-2,:)) + a32*dt./mass.*(F_m(1:end-1,:)-F_m(2:end,:)) - atilde31.*dt.*S - atilde32*dt/2.*(F_ms(1:end-1,:)+F_ms(2:end,:)) - gamma*dt*source(q_3) - q_3;
    q_3 = fsolve(residual_3, q); 

    % update step,
    F_3 = [flux(q_3(1,:)); flux(q_3); flux(q_3(end,:))];
    if (flag_tab>1)
        q = q - atilde31*(dt*S + dt./mass./2.*(F_l(3:end,:)-F_l(1:end-2,:))) - atilde32*(dt/2*(F_ms(1:end-1,:)+F_ms(2:end,:)) - dt./mass.*(F_m(1:end-1,:)-F_m(2:end,:))) - gamma*(dt*source(q_3) + dt./mass./2.*(F_3(3:end,:)-F_3(1:end-2,:)));
    else
        q = q_3;
    end

    if (flag_ini==1)
        dt = .01;
    else
        dt = time_step_selector(q, cfl, dx);
    end
    cfl = compute_cfl(q, dt, dx);

    time = time + dt;
    disp("Simulation time: "+time+", cfl number "+cfl)

    time_vec = [time_vec; time];
    q_vec = [q_vec; mean(abs(q))];

    norm_inf = max(norm_inf, norm(q-exact_sol(x,time), inf));
    norm_inf_vec = [norm_inf_vec, norm(q-exact_sol(x,time), inf)];

    
    plot(x, q, 'ob','MarkerIndices', 1:12:length(x), LineWidth=1.3)
    hold on
    plot(x, exact_sol(x,time), 'r')
    ylim([-max(q0(x))*0, max(q0(x))])
    drawnow
    grid on
    xlabel("x (m)")
    ylabel("q")
    set(gca,"FontSize",16)

    error = max(error, norm(exact_sol(x,time)-q, inf));
end
legend('Exact','Numerical')
toc

figure
plot(x, exact_sol(x,time), 'r', LineWidth=2)
hold on
plot(x, q, '--b', LineWidth=1.3)
grid on
xlabel("x (m)")
ylabel("q")
axis equal
set(gca,"FontSize",16)

error

fin_time_err = norm(exact_sol(x,time)-q, inf)

tau = 1/chi

% plots

%% make plot for paper
figure()
plot(time_vec, exact_sol(a*time_vec,time_vec), 'r', LineWidth=2)
hold on
plot(time_vec, q_vec, '--ob','MarkerIndices', 1:5:length(x), LineWidth=1.2)
xlabel("t (s)")
ylabel("q")
xlim([0, T])
grid on
legend('Exact','Numerical')
set(gca,"FontSize",16)


return


%% Internal Functions


function [a21, a32, a31, atilde32, atilde21, atilde31, gamma] = set_coefficients(flag_tab)
if (flag_tab==0)
    % coefficients in https://doi.org/10.1016/j.jcp.2024.112798
    gamma = .5;
    atilde21 = 0;
    a31 = 0;
    atilde31 = gamma;
    atilde32 = 0;
    a21 = .5;
    a32 = 1;
elseif (flag_tab==1)
    % coefficients in https://doi.org/10.1016/j.camwa.2025.02.014
    gamma = (1.-sqrt(2.)*.5);
    atilde21 = (-1.+sqrt(2.))*.5;
    a31 = 0;
    atilde31 = gamma;
    atilde32 = (sqrt(2.)-1.);
    a21 = .5;
    a32 = 1;
elseif (flag_tab==2)
    % coefficients Table 2 in https://arxiv.org/abs/2509.09460
    a21 = 1 - sqrt(2)/2;
    a32 = 1;
    a31 = 0;
    atilde32 = sqrt(2)/2;
    atilde21 = 0;
    atilde31 = 0;
    gamma = 1 - sqrt(2)/2;
elseif (flag_tab==3)
    % coefficients Table 3 in https://arxiv.org/abs/2509.09460
    a21 = sqrt(2)/4;
    a32 = (sqrt(2)-1)/(3*sqrt(2)-4)/2;
    a31 = 0;
    atilde32 = sqrt(2)/2;
    atilde21 = 0;
    atilde31 = 0;
    gamma = 1 - sqrt(2)/2;
end

end

function dt = time_step_selector(q, cfl, dx)
dt = cfl*dx(1)/max_eigenvalues(q);
end

function cfl = compute_cfl(q, dt, dx)
cfl = dt/dx(1)*max_eigenvalues(q);
end

function q0 = form_initial_cond(L, x, flag)
if (flag==1)
    q0(:,1) = x.*0+1;
elseif (flag==2)
    q0(:,1) = 1+3.*exp(-5*( (x-L/4).^2 )/(0.1*L)^2);
elseif (flag==3)
    % Parameters for windowed checkerboard
    A     = 1.0;        % amplitude
    x0    = L/2;        % center of the Gaussian envelope
    sigma = 0.05*L;     % width of Gaussian envelope
    for j = 1:length(x)
        q0(j,1) = A * (-1)^j * exp(- ( (x(j)-x0)^2 ) / (2*sigma^2) );
    end
end
end

function F = flux(q)
% linear advection
global a
F = a*q;
end

function S = source(q)
% linear advection
global chi
S = chi*q;
end

function lam = max_eigenvalues(q)
% linear advection
global a
lam = a;
end
