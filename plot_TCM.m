function plot_TCM(ntoud,ntolr,szt,szmv,cover,nelx,nely,nod_coor,loop,...
    xPhys,x_t,ele_nod,nele_t,filename1)

figure(1)
clf
hold on
axis equal; 
axis([0 (ntoud-1)*szt+szmv 0 (ntolr-1)*szt+cover*2]);
colormap(gray); imagesc(0.5,0.5,1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
% % plot vertical lines showing elements
% for ii = 1:nelx-1
%     h = plot( [ ii ii], [0 nely], 'k');
% end
% % plot horizontal lines showing elements
% for ii = 1:nely-1
%     h = plot( [0 nelx],[ii ii],  'k');
% end
% plot truss elements
for e=1:nele_t
    if x_t(e)>0.01
        xplot=[nod_coor(ele_nod(e,1),1) nod_coor(ele_nod(e,2),1)];
        yplot=[nod_coor(ele_nod(e,1),2) nod_coor(ele_nod(e,2),2)];
        h = plot(xplot,yplot,'Color',"#4DBEEE");
        set(h,{'LineWidth'},{x_t(e)*10})
    end
end
hold off

if loop == 0
    frame1 = getframe(1);
      im = frame2im(frame1);
      [imind,cm] = rgb2ind(im,256);
      imwrite(imind,cm,filename1,'gif', 'Loopcount',inf);
else
    frame1 = getframe(1);
      im = frame2im(frame1);
      [imind,cm] = rgb2ind(im,256);
      imwrite(imind,cm,filename1,'gif','WriteMode','append');
end

pause(0.1)