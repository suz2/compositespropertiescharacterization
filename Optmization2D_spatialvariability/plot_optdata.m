clc;clear;
totallevel = 4;
totallevel = totallevel - 1;
ii = 1;

errorbnd = [1 2 4 8];
global directory 
directory = './';

x_struc = {};
objfunc_struc = {};
iter_struc = {};

x_vec = [];
Ei_vec_mean = [];
a_vec_mean = [];
num_vec_mean = [];

Ei_vec_std = [];
a_vec_std = [];
num_vec_std = [];

Ei_vec_minobj = [];
a_vec_minobj = [];
num_vec_minobj = [];

e_ind_vec = [];
e_ind_vec2 = [];
objmin_ind = [];
for err_ind = 1 : length(errorbnd)
    e_ind = errorbnd(err_ind);

    string = ['e',int2str(e_ind)];
    [x_struc,objfunc_struc,iter_struc] = databuilding(x_struc,objfunc_struc,iter_struc,string,ii,e_ind,totallevel);
    x_vec = [x_vec;x_struc(ii).value];
    e_ind_vec = [e_ind_vec; e_ind*ones(length(x_struc(ii).value(:,1)),1)];

    Ei_vec_mean = [Ei_vec_mean;mean(x_struc(ii).value(:,1))];
    a_vec_mean = [a_vec_mean;mean(x_struc(ii).value(:,2))];
    num_vec_mean = [num_vec_mean;mean(x_struc(ii).value(:,3))];

    Ei_vec_std = [Ei_vec_std;std(x_struc(ii).value(:,1))];
    a_vec_std = [a_vec_std;std(x_struc(ii).value(:,2))];
    num_vec_std = [num_vec_std;std(x_struc(ii).value(:,3))];

    ind_objmin = find(objfunc_struc(ii).value(:,3) == min(objfunc_struc(ii).value(:,3)));

    Ei_vec_minobj = [Ei_vec_minobj;x_struc(ii).value(ind_objmin,1)];
    a_vec_minobj = [a_vec_minobj;x_struc(ii).value(ind_objmin,2)];
    num_vec_minobj = [num_vec_minobj;x_struc(ii).value(ind_objmin,3)];
    objmin_ind = [objmin_ind;ind_objmin];

    e_ind_vec2 = [e_ind_vec2; e_ind];
        
    ii = ii + 1;

end

for i = 1:length(Ei_vec_minobj)
    maxdiff(i) = modulusdiff(Ei_vec_minobj(i),a_vec_minobj(i));
end

%%
figure(1)
plot([0,8],[1,1]*7.5426,'k--');hold on
plot(e_ind_vec2,Ei_vec_minobj,'-x');hold on
set(gcf,'color','w');
set(gca,'fontsize', 18);
get(gca,'fontname');
set(gca,'fontname','times');
xlabel('Error bound')
ylabel('E_{inter}')
% legend('Ref','1 sample','10 samples','20 samples','50 samples','100 samples')
% lgd.FontSize = 12;

%%
figure(2)
plot([0,8],[1,1]*0.23465,'k--');hold on
plot(e_ind_vec2,a_vec_minobj,'-x')
set(gcf,'color','w');
set(gca,'fontsize', 18);
get(gca,'fontname');
set(gca,'fontname','times');
xlabel('Error bound')
ylabel('\alpha')
% legend('Ref','1 sample','10 samples','20 samples','50 samples','100 samples')
% lgd.FontSize = 12;

%%
figure(3)
plot([0,8],[1,1]*0.34,'k--');hold on
plot(e_ind_vec2,num_vec_minobj,'-x')
set(gcf,'color','w');
set(gca,'fontsize', 18);
get(gca,'fontname');
set(gca,'fontname','times');
xlabel('Error bound')
ylabel('\nu_m')
% legend('Ref','1 sample','10 samples','20 samples','50 samples','100 samples')


%%
figure(4)
plot(e_ind_vec2,(Ei_vec_minobj - 7.5426)./7.5426,'-o');hold on
ylabel('E_{inter} error')
set(gcf,'color','w');
set(gca,'fontsize', 18);
get(gca,'fontname');
set(gca,'fontname','times');
xlabel('Error bound')
% legend('1 sample','10 samples','20 samples','50 samples','100 samples')
% lgd.FontSize = 12;

%%
figure(5)
plot(e_ind_vec2,(a_vec_minobj - 0.23465)./0.23465,'-o');hold on
ylabel('\alpha error')
set(gcf,'color','w');
set(gca,'fontsize', 18);
get(gca,'fontname');
set(gca,'fontname','times');
xlabel('Error bound')
% legend('1 sample','10 samples','20 samples','50 samples','100 samples')
% lgd.FontSize = 12;

%%
figure(6)
plot(e_ind_vec2,(num_vec_minobj - 0.34)./0.34,'-o');hold on
ylabel('\nu_m error')
set(gcf,'color','w');
set(gca,'fontsize', 18);
get(gca,'fontname');
set(gca,'fontname','times');
xlabel('Error bound')
% legend('1 sample','10 samples','20 samples','50 samples','100 samples')
% lgd.FontSize = 12;

%%
figure(7)
plot(errorbnd,maxdiff,'-x');hold on
set(gcf,'color','w');
set(gca,'fontsize', 18);
get(gca,'fontname');
set(gca,'fontname','times');
xlabel('Error bound')
ylabel('max error of stiffness distribution')
% legend('1 sample','10 samples','20 samples','50 samples','100 samples')
% lgd.FontSize = 12;

%%
function [x_struc,objfunc_struc,iter_struc] = databuilding(x_struc,objfunc_struc,iter_struc,string,ii,e_ind,totallevel)
    x_struc(ii).name = string;
    objfunc_struc(ii).name = string;
    iter_struc(ii).name = string;
    j = 1;
    for Ei_ind = 0 : totallevel
       for a_ind = 0 : totallevel
           for nu = 0 : totallevel
               [x,objfunc,iteration] = extractdata(e_ind,Ei_ind,a_ind,nu);
               x_struc(ii).value(j,:) = x;
               objfunc_struc(ii).value(j,:) = objfunc;
               iter_struc(ii).value(j) = iteration;
               j = j + 1;
           end
       end
    end
end

function [x,objfunc,iteration] = extractdata(e_ind,Ei_ind,a_ind,nu)
    global directory
    string = [directory,'objfunc_e',int2str(e_ind),...
                'Ei',int2str(Ei_ind),'a',int2str(a_ind),'num',int2str(nu),'.dat'];
    data = dlmread(string);
    iteration = data(end,1);
    x = data(end-1,2:end);
    objfunc = data(end,2:end);
    x(1) = x(1) * 10;
end

function maxdiff = modulusdiff(Einter,alpha)

    fibercentroid = importdata('../fibercoord2.dat');
    r = 5;
    xcoord = linspace(0,100,201);
    ycoord = linspace(0,100,201);
    for i = 1 : length(xcoord)
        for j = 1 : length(ycoord)
            point((i - 1) * length(xcoord) + j,:) = [xcoord(i),ycoord(j)];
        end
    end

    for i = 1 : length(point(:,1))
        flag = 0;
        j = 0;
        while flag == 0 && j < length(fibercentroid(:,1))
            j = j + 1;
            if sqrt((point(i,1) - fibercentroid(j,1))^2 + (point(i,2) - fibercentroid(j,2))^2) <= r
                flag = 1;
            end
        end

        if flag == 0
           mindist(i) = sqrt(min((point(i,1) - fibercentroid(:,1)).*(point(i,1) - fibercentroid(:,1)) + ...
                            (point(i,2) - fibercentroid(:,2)).*(point(i,2) - fibercentroid(:,2)) )) - r;
        else
            mindist(i) = NaN;
        end
    end
    aref = 2.4826;
    bref = -0.23465;
    cref = 5.06;
    modulus1 = aref * exp(bref * mindist) + cref;
    modulus1 = reshape(modulus1,[length(xcoord), length(ycoord)]);

    a = Einter - 5.06;
    b = -alpha;
    c = 5.06;
    modulus2 = a * exp(b * mindist) + c;
    modulus2 = reshape(modulus2,[length(xcoord), length(ycoord)]);
    % contourf(xcoord,ycoord,modulus)
    maxdiff = max(max((modulus1-modulus2)./modulus1));

end