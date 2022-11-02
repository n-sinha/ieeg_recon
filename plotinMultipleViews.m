function plotinMultipleViews(plotObjs)

viewangles = [0,90;0,-90;-90,0;90,0;180,0;0,0];

for j = 1: size(viewangles,1)
    
    subplot(2,3,j)
    copyobj(plotObjs,gca);
    view(viewangles(j,:));
    axis square
    grid on
    axis off
    
end


end
