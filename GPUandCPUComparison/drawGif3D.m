function drawGif3D( timeIterator,deltaT, filename,  start, h )
    %color_range=[min(min(background )),max(max(background))];
    fig = h;
    
    title([filename, ' t = ',num2str( timeIterator * deltaT),'s'],'FontName','Times New Roman','FontSize', 16);
    %xlabel('X Point','FontName','Times New Roman','FontSize',14);
   %ylabel('Y Point','FontName','Times New Roman','FontSize',14); 
   %set(gca,'YDir','normal');

    fimage = getframe( fig );  
    imind=frame2im(fimage);  
    [imind,cm] = rgb2ind(imind,256);

    if ( timeIterator / start )== 1    
       imwrite(imind,cm, [ strtrim( filename) '.gif'],'gif', 'Loopcount',inf,'DelayTime',0.05);
    else
       imwrite(imind,cm, [ strtrim( filename) '.gif'],'gif','WriteMode','append','DelayTime',0.15);
    end