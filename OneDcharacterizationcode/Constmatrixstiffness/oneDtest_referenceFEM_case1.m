clc;clear;

%% material properties

% Cater et al. 2018
E_matrix_avg = 3.67e9;%Unit: Pa
E_fiber = 19.5e9;%Unit: Pa
alpha_matrix_avg = 50*1e-6;%Unit: (mm/mm)/K
alpha_fiber = 5.6*1e-6;%Unit: (mm/mm)/K

%temperature load
dT = -100;

%mechanical load
epsilon_bar = 0.01;

%% element partitions
%total number of elements
elem_num = 30;
% material assignment, 1 stands for fiber, 0 for matrix
material_no = [1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 ...
                0 0 0 1 1 1 0 0 0 1 1 1 0 0 0];

% stiffness varies at each element using E_m = E_matrix_avg * E_ratio
E_ratio = [0 0 0 0.6 0.6 0.6 0 0 0 0.5 0.5 0.5 0 0 0 ...
                1.2 1.2 1.2 0 0 0 1.5 1.5 1.5 0 0 0 0.8 0.8 0.8];
            

% CTE varies at each element using alpha_m = alpha_matrix_avg * alpha_ratio
alpha_ratio = [0 0 0 0.6 0.6 0.6 0 0 0 0.5 0.5 0.5 0 0 0 ...
                1.2 1.2 1.2 0 0 0 1.5 1.5 1.5 0 0 0 0.8 0.8 0.8];

for i = 1 : elem_num
      E_matrix(i) = E_matrix_avg * E_ratio(i);
      alpha_matrix(i) = alpha_matrix_avg * alpha_ratio(i);
end

%total length of bar
total_length = 0.04*1e-3;%Unit: m

% element length
le = total_length/elem_num * (zeros(1,elem_num)+1);

% total number of nodes
node_num = elem_num + 1;

if length(material_no) ~= elem_num
    error('length of material assignment is incomplete');
end

if length(le) ~= elem_num
    error('length of element length assignment is incomplete');
end


%% global stiffness & force vector construction

% initialize dimensions of vector
vec_num = node_num;
[K,f] = globalstiffness(elem_num,node_num,le,dT,E_matrix,E_fiber,alpha_matrix,alpha_fiber,material_no,vec_num);

%% applying displacement Periodic Boundary condition

[K,f,vec_num] = dispPBC(K,f,vec_num,node_num,elem_num,epsilon_bar*total_length);

%% applying fixed disp. BC at left end

[K,f,vec_num] = fixedBC(K,f,vec_num,node_num,elem_num);

%% Computation of nodal displacement

u = inv(K)*f';

%% post processing

% get mid_fiber displacement
[l,u_fiber,l_fiber] = dispextraction(node_num,elem_num,material_no,u,le);

% plotting infomation for material properties distribution
[l,lm,m,alpha_matrix_plot,E_matrix_plot] = plotinfo(node_num,elem_num,material_no,alpha_matrix,E_matrix,le);
ff = K(1:node_num,1:node_num)*u(1:node_num);
reactionforce = ff(end);

figure(1)

plot(l*1e3,u(1:node_num)*1e3,'k-o');hold on
plot(l_fiber*1e3,u_fiber*1e3,'r-*');hold on
title('Nodal displacement')

set(gcf,'color','w');
xlabel('y/mm');
ylabel('u/mm');
xlim([0 total_length*1e3]);
% ylim([0 400]);
set(gca,'fontsize', 12);
get(gca,'fontname')
set(gca,'fontname','times')
% legend('Nodal displacement(ref)')

figure(2)

subplot(3,1,1)
plot(l*1e3,u(1:node_num)*1e3,'k--o');hold on
plot(l_fiber*1e3, u_fiber*1e3,'k--*');hold on
title('Nodal displacement')
% plot(l,ua,'-*');hold on
set(gcf,'color','w');
xlabel('y/mm');
ylabel('u/mm');
xlim([0 total_length*1e3]);
% ylim([0 400]);
set(gca,'fontsize', 12);
get(gca,'fontname')
set(gca,'fontname','times')
legend('Nodal displacement(ref)','Mid fiber displacement(ref)')

subplot(3,1,2)
plot(lm*1e3,(zeros(1,length(lm))+1)*E_matrix_avg/1e9,'k--');hold on
plot(lm*1e3,E_matrix_plot/1e9,'k-');hold on

title('Varying Matrix Stiffness')
set(gcf,'color','w');
% xlabel('y/mm');
ylabel('E_m/Gpa');
% ylabel('\alpha_m');
xlim([0 total_length*1e3]);
set(gca,'fontsize', 12);
get(gca,'fontname')
set(gca,'fontname','times')

subplot(3,1,3)
plot(lm*1e3,(zeros(1,length(lm))+1)*alpha_matrix_avg,'k--');hold on
plot(lm*1e3,alpha_matrix_plot,'k-');hold on

title('Varying Matrix CTE')
set(gcf,'color','w');
xlabel('y/mm');
% ylabel('E_m/Gpa');
ylabel('\alpha_m');
xlim([0 total_length*1e3]);
set(gca,'fontsize', 12);
get(gca,'fontname')
set(gca,'fontname','times')



