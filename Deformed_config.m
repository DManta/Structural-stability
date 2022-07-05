function Deformaed_config(NEF, X_0_plot, X_0, X_t, Escala, Psi_plast, NGauss_r, NGauss_s, NGauss_t, Plast, Se)

%Plot the deformed configuration

%Indeformed...............................................................

%Coordinates of nodes
F_indef=zeros(NEF_indef,4);
V_indef=zeros(4*NEF_indef,3);

for el = 1:NEF_indef  
    V_indef(1+(el-1)*4:4+(el-1)*4,:) = X_0_plot(:,:,el)';
    F_indef(el,:) = 1+(el-1)*4:4+(el-1)*4;    
end

%Indeformed
axis equal;
patch('Vertices',V_indef,'Faces',F_indef,'FaceColor','none','EdgeColor',[0 0 1],'LineWidth',1.0);
set(gcf,'color','w');
set(gca,'visible','off');
set(gcf, 'Renderer','painters') %renders

%Deformed.................................................................

%Material properties
name = 'Material_C.txt';
table = dlmread(name); clear name;
Sy=table(2,5); %Yield stress

%Patches coordinates
F_gauss=zeros(NEF*NGauss_r*NGauss_s,4);
V_gauss=zeros(4*NEF*NGauss_r*NGauss_s,3);
F_color=zeros(NEF*NGauss_r*NGauss_s,3);


Colors = jet;

for el = 1:NEF

    cont = 0;   %Counter

    for i=1:NGauss_r
         for ii=1:NGauss_s       
            cont = cont +1; %contador em (r,s)
            cont2 = (el-1)*NGauss_r*NGauss_s + (i-1)*NGauss_s + ii;

            cor = round(size(Colors,1)*Se((cont-1)*NGauss_t + Gauss_t,el)/Sy);
          
            F_color(cont2,:) = Colors(cor,:);
            
            V_gauss(1+(cont2-1)*4:4+(cont2-1)*4,:) = Psi_plast(:,:,cont)*(X + Escala*[d(uelement(1:3))'; d(uelement(7:9))'; d(uelement(13:15))'; d(uelement(19:21))']);
            F_gauss(cont2,:) = 1+(cont2-1)*4:4+(cont2-1)*4; 
         end
    end
end

%Deformed
patch('Vertices',V_gauss,'Faces',F_gauss,'FaceVertexCData',F_color,'FaceColor','flat','EdgeColor',[0 0 0]);
set(gcf,'color','w');
set(gcf, 'Renderer','painters')
set(gca,'visible','off');

% To plot with the efective stresses
colormap jet
caxis([0, Sy*10^-3]);
colorbar

hold off

end