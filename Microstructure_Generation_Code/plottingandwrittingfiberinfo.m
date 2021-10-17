for i = 1: nparticles
    x = circles(i,1);
    y = circles(i,2);
    r = 5;
    plotcircle(x,y,r);
    for j = 1:nparticles
        dist(i,j) = sqrt((circles(i,1) - circles(j,1))^2 + (circles(i,2) - circles(j,2))^2);
        if i == j;dist(i,j) = image_size*100;end
    end
    
    if circles(i,1) < r
        x = circles(i,1) + image_size;
        y = circles(i,2);
        plotcircle(x,y,r);
    end
    
    if image_size - circles(i,1) < r
        x = - image_size + circles(i,1);
        y = circles(i,2);
        plotcircle(x,y,r);
    end
    
    if circles(i,2) < r
        x = circles(i,1);
        y = circles(i,2) + image_size;
        plotcircle(x,y,r);
    end
    
    if image_size - circles(i,2) < r
        x = circles(i,1);
        y = - image_size + circles(i,2);
        plotcircle(x,y,r);
    end

    if (circles(i,1) < r) && (circles(i,2) < r)
        x = circles(i,1) + image_size;
        y = circles(i,2) + image_size;
        plotcircle(x,y,r);
    end
    
    if (image_size - circles(i,1) < r) && (circles(i,2) < r)
        x = - image_size + circles(i,1);
        y = circles(i,2) + image_size;
        plotcircle(x,y,r);
    end    
    
    if (image_size - circles(i,1) < r) && (image_size - circles(i,2) < r)
        x = - image_size + circles(i,1);
        y = - image_size + circles(i,2);
        plotcircle(x,y,r);
    end
    
    if (circles(i,1) < r) && (image_size - circles(i,2) < r)
        x = circles(i,1) + image_size;
        y = - image_size + circles(i,2);
        plotcircle(x,y,r);
    end   

end
xlim([0 image_size]);    
ylim([0 image_size]);    
set(gcf,'color','w');
set(gca,'fontsize', 18);
get(gca,'fontname')
set(gca,'fontname','times')

mindist = min(min(dist))-5*2;


fid = fopen('fibercentroid_r5vf55domain50.txt','w+');
for i = 1: nparticles  

    fprintf(fid,'%2d,%2d,%2.2f,%2d\n',circles(i,1)-25, circles(i,2)-25,circles(i,1)+circles(i,3)-25, circles(i,2)-25);

end
fclose(fid);


function plotcircle(x,y,r)

    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    plot(xunit, yunit,'k-');hold on
end