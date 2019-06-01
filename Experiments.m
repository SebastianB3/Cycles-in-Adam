%% Sebastian Bock & Martin Weiß
% OTH Regensburg
% 02.06.2019
%% How To
% For the Experiments run firstly the chosen parameters and than the
% general Code
clc
close all
%% Parameters for Experiment 1
alpha_max = 0.001;
alpha_min = 10^(-4);
epsilon = 10^(-8);
beta_1 = 0.9;
beta_2 = 0.999;
alpha = 0.001;
eigenvalue = 10;
w_lsg = 0;
c = 10;
m_t_start =-0.1281144718e-4;
v_t_start =5.925207756*10^(-8);
w_t_start =0.00002434174964;

m = @(m_, w_) beta_1*m_+(1-beta_1)*c*w_;
v = @(v_, w_) beta_2*v_+(1-beta_2)*(c*w_)^2;
w = @(m_, v_, w_, t, alpha) w_-(alpha*m_/(sqrt(v_+epsilon)));

%% Parameters for Experiment 2

alpha_max = 0.01;
alpha_min = 10^(-4);
epsilon = 10^(-6);
beta_1 = 0.2;
beta_2 = 0.5;
alpha = 0.5;
c = 1;
m_t_start =0;
v_t_start =0;
w_t_start =eps;
eigenvalue = 1;
w_lsg = 0;

m = @(m_, w_) beta_1*m_+(1-beta_1)*c*w_;
v = @(v_, w_) beta_2*v_+(1-beta_2)*(c*w_)^2;
w = @(m_, v_, w_, t, alpha) w_-(alpha*m_/(sqrt(v_+epsilon)));

%% Parameters for Experiment 3

alpha_max = 1;
alpha_min = 10^(-4);
epsilon = 0.01;

beta_1 = 0.5;
beta_2 = 0.6;
alpha = 0.8;
c = 1;
m_t_start =0;
v_t_start =0;
w_t_start =eps;
eigenvalue = 1;
w_lsg = 0;

m = @(m_, w_) beta_1*m_+(1-beta_1)*c*w_;
v = @(v_, w_) beta_2*v_+(1-beta_2)*(c*w_)^2;
w = @(m_, v_, w_, t, alpha) w_-(alpha*m_/(sqrt(v_+epsilon)));

%% General Code

N = 10000;
Jacobian = @(m,v,w) [beta_1 0 (1-beta_1)*c; 0 beta_2 2*(1-beta_2)*c^2*w; -alpha*beta_1/sqrt(beta_2*v+(1-beta_2)*c^2*w^2+epsilon) 0.5*alpha*(beta_1*m+(1-beta_1)*c*w)*beta_2/(beta_2*v+(1-beta_2)*c^2*w^2+epsilon)^(2/3) 1-alpha*(1-beta_1)*c/sqrt(beta_2*v+(1-beta_2)*c^2*w^2+epsilon)+alpha*(beta_1*m+(1-beta_1)*c*w)*(1-beta_2)*c^2*w/(beta_2*v+(1-beta_2)*c^2*w^2+epsilon)^(2/3)];
Jacobian_bias = @(t,m,v,w) [beta_1 0 (1-beta_1)*c; 0 beta_2 2*(1-beta_2)*c^2*w; -alpha*sqrt(1-beta_2^t)*beta_1/((1-beta_1^t)*sqrt(beta_2*v+(1-beta_2)*c^2*w^2+epsilon)) 0.5*sqrt(1-beta_2^t)*alpha*(beta_1*m+(1-beta_1)*c*w)*beta_2/((1-beta_1^t)*(beta_2*v+(1-beta_2)*c^2*w^2+epsilon)^(2/3)) 1-alpha*sqrt(1-beta_2^t)*(1-beta_1)*c/((1-beta_1^t)*sqrt(beta_2*v+(1-beta_2)*c^2*w^2+epsilon))+alpha*sqrt(1-beta_2^t)*(beta_1*m+(1-beta_1)*c*w)*(1-beta_2)*c^2*w/((1-beta_1^t)*(beta_2*v+(1-beta_2)*c^2*w^2+epsilon)^(2/3))];

m_t = zeros(1,N);
v_t = zeros(1,N);
w_t = zeros(1,N);
m_t(1) = m_t_start;
v_t(1) = v_t_start; 
w_t(1) = w_t_start; 
eigenvalues = zeros(3,5);

for i=2:N
    m_t(i) = m(m_t(i-1),w_t(i-1));
    v_t(i) = v(v_t(i-1),w_t(i-1));
    w_t(i) = w(m_t(i), v_t(i), w_t(i-1),i, alpha);
    if i>N-5
        eigenvalues(1:3,i-N+5) = eig(Jacobian(m_t(i),v_t(i),w_t(i)));
    end
end

w_bias = @(m_, v_, w_, t) w_-(alpha*m_*sqrt(1-beta_2^t))/((1-beta_1^t)*(sqrt(v_+epsilon)));
m_t_bias = zeros(1,N);
v_t_bias = zeros(1,N);
w_t_bias = zeros(1,N);
m_t_bias(1) = m_t_start;
v_t_bias(1) = v_t_start;
w_t_bias(1) = w_t_start;
eigenvalues_bias = zeros(3,5);

for i=2:N
    m_t_bias(i) = m(m_t_bias(i-1),w_t_bias(i-1));
    v_t_bias(i) = v(v_t_bias(i-1),w_t_bias(i-1));
    w_t_bias(i) = w_bias(m_t_bias(i), v_t_bias(i), w_t_bias(i-1),i);
    if i>N-5
        eigenvalues_bias(1:3,i-N+5) = eig(Jacobian_bias(i, m_t_bias(i),v_t_bias(i),w_t_bias(i)));
    end
end

fontsize = 40;
fontsize_2D = 40;
markersize = 12;
figure(1)
hold on
plt1=plot3(w_t, m_t, v_t, 'g.', 'MarkerSize',markersize);
plt2=plot3(w_t_bias, m_t_bias, v_t_bias, 'r.', 'MarkerSize',markersize);
set(gca, 'FontSize', fontsize)
view(3)
xlabel('w_t', 'FontSize', fontsize)
ylabel('m_t', 'FontSize', fontsize)
zlabel('v_t', 'FontSize', fontsize)
legend([plt1(1),plt2(1)],{'without Bias-Correction', 'with Bias-Correction'}, 'FontSize', fontsize)
hold off

figure(2)
hold on
plot(w_t)
plot(w_t_bias)
set(gca, 'FontSize', fontsize_2D)
xlabel('iterations')
ylabel('w_t')
legend('without Bias-Correction', 'with Bias-Correction')
hold off

%% Plot alpha_t and w_t (Experiment 1, 2 and 3)
% First run the parameters

N = 100000;
N_alpha = 1000;
alphas = alpha_min:(alpha_max-alpha_min)/N_alpha:alpha_max;
fontsize_2D = 15;

m_t = zeros(N_alpha,N);
v_t = zeros(N_alpha,N);
w_t = zeros(N_alpha,N);
alpha_t = zeros(N_alpha,N);

for j=1:N_alpha
    alpha_t(j,1:N) = alphas(j);
    m_t(j,1) = m_t_start;
    v_t(j,1) = v_t_start; 
    w_t(j,1) = w_t_start; 
    for i=2:N
        m_t(j,i) = m(m_t(j,i-1),w_t(j,i-1));
        v_t(j,i) = v(v_t(j,i-1),w_t(j,i-1));
        w_t(j,i) = w(m_t(j,i), v_t(j,i), w_t(j,i-1),i, alpha_t(j,1));
    end
end

w_bias = @(m_, v_, w_, t, alpha) w_-(alpha*m_*sqrt(1-beta_2^t))/((1-beta_1^t)*(sqrt(v_+epsilon)));
m_t_bias = zeros(N_alpha,N);
v_t_bias = zeros(N_alpha,N);
w_t_bias = zeros(N_alpha,N);
alpha_t_bias = zeros(N_alpha,N);

figure(1)
hold on
plt1=plot(w_t(1:end, end-500:end), alpha_t(1:end, end-500:end),'g.');
plot(w_lsg,(sqrt(epsilon)*(2*beta_1+2))/(eigenvalue*(1-beta_1)),'rx')
set(gca, 'FontSize', fontsize_2D)
xlabel('w_t')
ylabel('\alpha_t')
hold off

%% Bifurcation in a multidimensional optimization (Figure 6)
clc
close all
epsilon = 10^(-8);
beta_1 = 0.1;
beta_2 = 0.5;
alpha = 0.5;
c = 1;

m_t_start = [-.1673553519 -.1673553519]';
v_t_start = [0.4183883298e-1 0.4183883298e-1]';
w_t_start = [0.2045454301 0.2045454301]';
w_lsg = 0;

phi= -0.2; Q = [cos(phi) -sin(phi); sin(phi) cos(phi)];
A = Q'*[1 0; 0 4]*Q;
eigenvalue = max(abs(eig(A)));
g = @(w) w'*A;

m = @(m_, w_) beta_1*m_+(1-beta_1)*g(w_)';
v = @(v_, w_) beta_2*v_+(1-beta_2)*(g(w_)').^2;
w = @(m_, v_, w_, t, alpha) w_-(alpha*m_./(sqrt(v_+epsilon)));

N = 100000;
N_alpha = 100;
alpha_max = 0.0005;
alpha_min = 10^(-5);
alphas = alpha_min:(alpha_max-alpha_min)/N_alpha:alpha_max;

m_t = zeros(N_alpha,2*N);
v_t = zeros(N_alpha,2*N);
w_t = zeros(N_alpha,2*N);
alpha_t = zeros(N_alpha,N);

for j=1:N_alpha
    alpha_t(j,1:N) = alphas(j);
    m_t(j,1:2) = m_t_start;
    v_t(j,1:2) = v_t_start; 
    w_t(j,1:2) = w_t_start; 
    for i=3:2:2*N-1
        m_t(j,i:i+1) = m(m_t(j,i-2:i-1)',w_t(j,i-2:i-1)');
        v_t(j,i:i+1) = v(v_t(j,i-2:i-1)',w_t(j,i-2:i-1)') ;  
        w_t(j,i:i+1) = w(m_t(j,i:i+1)', v_t(j,i:i+1)', w_t(j,i-2:i-1)',i, alpha_t(j,1));
    end
end

figure(1)
hold on
plt1=plot(w_t(1:end, end-500:2:end), alpha_t(1:end, end-250:end),'g.');
plt2=plot(w_t(1:end, end-501:2:end), alpha_t(1:end, end-250:end),'b.');
plt3 = plot(w_lsg,(sqrt(epsilon)*(2*beta_1+2))/(eigenvalue*(1-beta_1)),'rX');
legend([plt1(1), plt2(1),plt3(1)],{'First Component of m_t','Second Component of m_t','Our inequality'})
xlabel('w_t')
ylabel('\alpha_t')
hold off