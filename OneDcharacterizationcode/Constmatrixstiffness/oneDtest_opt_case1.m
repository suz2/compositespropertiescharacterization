clc;clear

warning('off','all')
rmpath('folderthatisnotonpath')

%% element partitions
%total number of elements
elem_num = 30;
% material assignment, 1 stands for fiber, 0 for matrix
material_no = [1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 ...
                0 0 0 1 1 1 0 0 0 1 1 1 0 0 0];
            
E_matrix_avg = 3.67e9;
alpha_matrix_avg = 50*1e-6;

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

%%

% initialization for each matrix and fiber part
x0 = [2 2 2 2 2 20 20 20 20 20];
% x0 = [[0.6 0.5 1.2 1.5 0.8]*3.67,[0.6 0.5 1.2 1.5 0.8]*50];

% lower bound
lb = [0 0 0 0 0 0 0 0 0 0];

% upper bound
ub = [10 10 10 10 10 100 100 100 100 100];

A = [];
b = [];
Aeq = [];
beq = [];

nonlcon = @rfcon;
% nonlcon = [];

% selection of algorithm
% alg = 'interior-point';
alg = 'sqp';
% alg = 'active-set';

fname = ['E_alpha_objfunc_',alg,'.txt'];
fid = fopen(fname,'w+');% write data at the end of the file;
fclose(fid);
    
%% 
fname = ['E_alpha_objfunc_',alg,'_u.txt'];
fid = fopen(fname,'w+');
fclose(fid);

options = optimoptions('fmincon','Display','final-detailed','Algorithm',alg,...
    'FunctionTolerance',1e-16,'StepTolerance',1e-16,'ConstraintTolerance',1e-16,'MaxIterations',10000,...
    'MaxFunctionEvaluations',10000);
x = fmincon(@oneD,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

%% get optimized properties and displacement field

x0 = [2 2 2 2 2 20 20 20 20 20];
% alg = 'interior-point';
%     alg = 'sqp';
%     alg = 'active-set';
fname = ['E_alpha_objfunc_',alg,'.txt'];
data = importdata(fname);
x = data(end,:);

E_matrix = [0 0 0 x(1) x(1) x(1) 0 0 0 x(2) x(2) x(2) 0 0 0 ...
            x(3) x(3) x(3) 0 0 0 x(4) x(4) x(4) 0 0 0 x(5) x(5) x(5)]*1e9;

alpha_matrix = [0 0 0 x(6) x(6) x(6) 0 0 0 x(7) x(7) x(7) 0 0 0 ...
            x(8) x(8) x(8) 0 0 0 x(9) x(9) x(9) 0 0 0 x(10) x(10) x(10)]*1e-6;  
        
E_matrix_initial = [0 0 0 x0(1) x0(1) x0(1) 0 0 0 x0(2) x0(2) x0(2) 0 0 0 ...
    x0(3) x0(3) x0(3) 0 0 0 x0(4) x0(4) x0(4) 0 0 0 x0(5) x0(5) x0(5)]*1e9;

alpha_matrix_initial = [0 0 0 x0(6) x0(6) x0(6) 0 0 0 x0(7) x0(7) x0(7) 0 0 0 ...
            x0(8) x0(8) x0(8) 0 0 0 x0(9) x0(9) x0(9) 0 0 0 x0(10) x0(10) x0(10)]*1e-6;  


        
fname = ['E_alpha_objfunc_',alg,'_u.txt'];
data_u = importdata(fname);

u = data_u(end,1:(end-1));

% extract mid-fiber displacement
[l,u_fiber,l_fiber] = dispextraction(node_num,elem_num,material_no,u,le);
[l,lm,m,alpha_matrix_plot,E_matrix_plot] = plotinfo(node_num,elem_num,material_no,alpha_matrix,E_matrix,le);

[l,lm,m,alpha_matrix_initial_plot,E_matrix_initial_plot] = plotinfo(node_num,elem_num,material_no,alpha_matrix_initial,E_matrix_initial,le);


%% plot iteration step vs. x/objfunc

fname = ['E_alpha_objfunc_',alg,'.txt'];
data = importdata(fname);

figure(1)
plot(l*1e3,u(1:node_num)*1e3,'-o');hold on
plot(l_fiber*1e3, u_fiber*1e3,'-*');hold on
% plot(l,ua,'-*');hold on
set(gcf,'color','w');
xlabel('y/mm');
ylabel('u/mm');
xlim([0 total_length*1e3]);
% ylim([0 400]);
set(gca,'fontsize', 12);
get(gca,'fontname')
set(gca,'fontname','times')
% legend('Ref nodal disp.','Ref Mid fiber disp.','Nodal disp.','Mid fiber disp.')
legend('Nodal disp.(ref)','Fiber centroid disp.(ref)','Nodal disp.(opt)','Fiber centroid disp.(opt)')


%% 

figure(2)

subplot(3,1,1)
plot(l*1e3,u(1:node_num)*1e3,'-o');hold on
plot(l_fiber*1e3, u_fiber*1e3,'-*');hold on
% plot(l,ua,'-*');hold on
set(gcf,'color','w');
xlabel('y/mm');
ylabel('u/mm');
xlim([0 total_length*1e3]);
% ylim([0 400]);
set(gca,'fontsize', 12);
get(gca,'fontname')
set(gca,'fontname','times')
% legend('Ref nodal disp.','Ref Mid fiber disp.','Nodal disp.','Mid fiber disp.')
legend('Nodal disp.(ref)','Fiber centroid disp.(ref)','Nodal disp.(opt)','Fiber centroid disp.(opt)')

subplot(3,1,2)
plot(lm*1e3,E_matrix_initial_plot/1e9,'b--');hold on
plot(lm*1e3,E_matrix_plot/1e9,'b-');hold on

title('Stiffness')
set(gcf,'color','w');
xlabel('y/mm');
ylabel('E_m/Gpa');
xlim([0 total_length*1e3]);
set(gca,'fontsize', 12);
get(gca,'fontname')
set(gca,'fontname','times')
legend('Nominal','Ref.','Initial','Predicted')

subplot(3,1,3)
plot(lm*1e3,alpha_matrix_initial_plot,'b--');hold on
plot(lm*1e3,alpha_matrix_plot,'b-');hold on
title('CTE')
set(gcf,'color','w');
xlabel('y/mm');
ylabel('\alpha_m');
xlim([0 total_length*1e3]);
set(gca,'fontsize', 12);
get(gca,'fontname')
set(gca,'fontname','times')
legend('Nominal','Ref.','Initial','Predicted')


function objfunc = oneD(x)

%% element partitions 
%element and material information have to be redefined in fmincon function, 
%because x/objfunc is the only variables that can be transpass between main func and fmincon

    %total number of elements
    elem_num = 30;
    % material assignment, 1 stands for fiber, 0 for matrix
    material_no = [1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 ...
                    0 0 0 1 1 1 0 0 0 1 1 1 0 0 0];

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

    %% material properties
    
    E_matrix = [0 0 0 x(1) x(1) x(1) 0 0 0 x(2) x(2) x(2) 0 0 0 ...
                x(3) x(3) x(3) 0 0 0 x(4) x(4) x(4) 0 0 0 x(5) x(5) x(5)]*1e9;
    E_fiber = 19.5e9;%Unit: Pa
    alpha_matrix = [0 0 0 x(6) x(6) x(6) 0 0 0 x(7) x(7) x(7) 0 0 0 ...
                x(8) x(8) x(8) 0 0 0 x(9) x(9) x(9) 0 0 0 x(10) x(10) x(10)]*1e-6; 
    
    alpha_fiber = 5.6*1e-6;%Unit: (mm/mm)/K

    %% Loading information

    dT = -100;
    epsilon_bar = 0.01;

    
    %% global stiffness & force vector construction

    % initialize dimensions of vector
    vec_num = node_num;

    [K,f] = globalstiffness(elem_num,node_num,le,dT,E_matrix,E_fiber,alpha_matrix,alpha_fiber,material_no,vec_num);

    %% applying displacement Periodic Boundary condition

    [K,f,vec_num] = dispPBC(K,f,vec_num,node_num,elem_num, epsilon_bar*total_length);

    %% applying fixed disp. BC at left end

    [K,f,vec_num] = fixedBC(K,f,vec_num,node_num,elem_num);

    %% Computation of nodal displacement

    u = inv(K)*f';

    %% post processing

    % extract mid-fiber displacement
    [l,u_fiber,l_fiber] = dispextraction(node_num,elem_num,material_no,u,le);
    
    % u_ref is obtained from oneDtest_referenceFEM.m
    u_ref= [5.31576929778756e-09,1.17932322704304e-07,2.55345879073010e-07,2.98969925074055e-07,3.25195469594006e-07];
    reactionforce_ref = [51004750.6534239];
    
    ff = K(1:node_num,1:node_num)*u(1:node_num);
    reactionforce = ff(end);

    objfunc = sum(abs((u_fiber - u_ref)./u_ref))/length(u_ref);
    
%     alg = 'interior-point';
    alg = 'sqp';
%     alg = 'active-set';

%% print x and objective function at each iteration

    fname = ['E_alpha_objfunc_',alg,'.txt'];
    fid = fopen(fname,'at');% write data at the end of the file;
    fprintf(fid,'%d     %d     %d       %d     %d     %d        %d     %d      %d      %d      %d      %d\n',...
        x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),objfunc,abs(reactionforce - reactionforce_ref));
    fclose(fid);
    
%% print displacement field at each iteration

    fname = ['E_alpha_objfunc_',alg,'_u.txt'];
    fid = fopen(fname,'at');
    for i = 1: length(u)
        if i == length(u)
            fprintf(fid,'%d     %d\n',u(i),reactionforce);
        else
            fprintf(fid,'%d     ',u(i));
        end
    end
    
    fclose(fid);

end

function [c,ceq] = rfcon(x)
    %% element partitions 
%element and material information have to be redefined in fmincon function, 
%because x/objfunc is the only variables that can be transpass between main func and fmincon

    %total number of elements
    elem_num = 30;
    % material assignment, 1 stands for fiber, 0 for matrix
    material_no = [1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 ...
                    0 0 0 1 1 1 0 0 0 1 1 1 0 0 0];

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

    %% material properties
    
    E_matrix = [0 0 0 x(1) x(1) x(1) 0 0 0 x(2) x(2) x(2) 0 0 0 ...
                x(3) x(3) x(3) 0 0 0 x(4) x(4) x(4) 0 0 0 x(5) x(5) x(5)]*1e9;
    E_fiber = 19.5e9;%Unit: Pa
    alpha_matrix = [0 0 0 x(6) x(6) x(6) 0 0 0 x(7) x(7) x(7) 0 0 0 ...
                x(8) x(8) x(8) 0 0 0 x(9) x(9) x(9) 0 0 0 x(10) x(10) x(10)]*1e-6; 
    
    alpha_fiber = 5.6*1e-6;%Unit: (mm/mm)/K

    %% Loading information
    %Thermal load
    dT = -100;
    
    %mechanical load
    epsilon_bar = 0.01;
    
    %% global stiffness & force vector construction

    % initialize dimensions of vector
    vec_num = node_num;

    [K,f] = globalstiffness(elem_num,node_num,le,dT,E_matrix,E_fiber,alpha_matrix,alpha_fiber,material_no,vec_num);

    %% applying displacement Periodic Boundary condition

    [K,f,vec_num] = dispPBC(K,f,vec_num,node_num,elem_num, epsilon_bar*total_length);

    %% applying fixed disp. BC at left end

    [K,f,vec_num] = fixedBC(K,f,vec_num,node_num,elem_num);

    %% Computation of nodal displacement

    u = inv(K)*f';

    %% post processing

    % extract mid-fiber displacement
    [l,u_fiber,l_fiber] = dispextraction(node_num,elem_num,material_no,u,le);
    %case 1
    u_ref= [5.31576929778756e-09,1.17932322704304e-07,2.55345879073010e-07,2.98969925074055e-07,3.25195469594006e-07];
    reactionforce_ref = [51004750.6534239];
    
    ff = K(1:node_num,1:node_num)*u(1:node_num);
    reactionforce = ff(end);
    c = abs(reactionforce - reactionforce_ref)/reactionforce_ref-0.01;
    ceq = [];
end