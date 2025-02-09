clear all;
format long;
digits(64);
tic;

% load data from the origininal file

position = xlsread('bacteria.xlsx');                 %node position information
ele = xlsread('bacteria_element.xlsx');   %element information

%Define the bond strength of the particle.(harmonic bond)
bondcoeff = 1.0;

dia = 10;  %diameter
Dx=62; % box size  old value for 27 cells is 100
Dy=62;
Dz=62;

%recenter 
position(:,1)=position(:,1) - mean(position(:,1)) ;
position(:,2)=position(:,2) - mean(position(:,2)) ;
position(:,3)=position(:,3) - mean(position(:,3)) ;

% Define the box size 
box=zeros(3,1);
box(1,1)=0;
box(1,2)=Dx;
box(2,1)=0;
box(2,2)=Dy;
box(3,1)=0;
box(3,2)=Dz;

num_tot = 2;   %total number of the molecule old value = 27
num_atom = size(position,1);
atom_xyz = zeros(num_tot*num_atom,6);
num_ele = size(ele,1);
num_points = zeros(num_tot*num_ele,3);

x_tot=2; % old value = 3
y_tot=1;
z_tot=1;
L_x=Dx/x_tot;
L_y=Dy/y_tot;
L_z=Dz/z_tot;
NP_center = zeros(num_tot,3);


%without random translation

for x_n=1:x_tot              %distribute position of centers
    for y_n=1:y_tot
        for z_n=1:z_tot
            num = x_n+x_tot*(y_n-1)+x_tot*y_tot*(z_n-1);
for n_atom=1:num_atom
atom_xyz((num-1)*num_atom+n_atom,1)=n_atom+(num-1)*num_atom;  %node number
atom_xyz((num-1)*num_atom+n_atom,2)=num;                             %molecule type
atom_xyz((num-1)*num_atom+n_atom,3)=1;                             %atom type
atom_xyz((num-1)*num_atom+n_atom,4)=position(n_atom,1)+(x_n-1)*L_x+L_x/2.0;     %node position, 4->6: x,y,z
atom_xyz((num-1)*num_atom+n_atom,5)=position(n_atom,2)+(y_n-1)*L_y+L_y/2.0;
atom_xyz((num-1)*num_atom+n_atom,6)=position(n_atom,3)+(z_n-1)*L_z+L_z/2.0;
end

        end
    end
end

   
for num=1:num_tot
    %node transformation
   
   %element transformation
   for n_ele=1:num_ele
       num_points((num-1)*num_ele+n_ele,1)=ele(n_ele,1)+num_atom*(num-1);
       num_points((num-1)*num_ele+n_ele,2)=ele(n_ele,2)+num_atom*(num-1);
       num_points((num-1)*num_ele+n_ele,3)=ele(n_ele,3)+num_atom*(num-1);
   end
   % compute the center of the molecule
  NP_center(num,:)=[mean(atom_xyz((num-1)*num_atom+1:num*num_atom,4)),mean(atom_xyz((num-1)*num_atom+1:num*num_atom,5)),mean(atom_xyz((num-1)*num_atom+1:num*num_atom,6))];   %center of the model
   
end

toc
tic;

%bond infromation
NN=zeros(1,4);
count=0;
for kk=1:num_ele
 temp=num_points(kk,:);
 temp=sort(temp);
 count=count+1;
 NN(count,3)=temp(1);
 NN(count,4)=temp(2);
 count=count+1;
 NN(count,3)=temp(2);
 NN(count,4)=temp(3);
 count=count+1;
 NN(count,3)=temp(1);
 NN(count,4)=temp(3);
end
 [Bond, ia, ic] = unique(NN,'rows');
 Bond(:,1) = (1:size(Bond,1))';
 Bond(:,2) = Bond(:,1);
 
% store the bond length for stress free model

Len_bond=zeros(size(Bond,1),1);
for kk =1:size(Bond,1)
index1=find(atom_xyz(:,1)==Bond(kk,3));
index2=find(atom_xyz(:,1)==Bond(kk,4));
p1=[atom_xyz(index1,4), atom_xyz(index1,5), atom_xyz(index1,6)];
p2=[atom_xyz(index2,4), atom_xyz(index2,5), atom_xyz(index2,6)];
vector=[p1(1)-p2(1),p1(2)-p2(2),p1(3)-p2(3)];
Len_bond(kk,1)=norm(vector,2);
end
toc
tic;

 Bond_single = Bond;
 Bond_temp = Bond;
 for kk=1:num_tot-1
     Bond_temp(:,1:2) = Bond_single(:,1:2) + kk*size(Bond_single,1);
     Bond_temp(:,3:4) = Bond_single(:,3:4) + kk*num_atom;
     Bond = [Bond;Bond_temp]; 
 end
 
n_bond = size(Bond_single,1);
for i=1:num_tot
    Bond_len_temp((i-1)*n_bond+1:i*n_bond,1)=Len_bond;
end
Len_bond=Bond_len_temp;

toc
tic;

%add the angle part
Angle_single = zeros(size(ele,1),5);
Angle=zeros(size(num_points,1),5);
count=0;
for i=1:size(ele,1)
   count=count+1;
   Angle_single(count,3)=ele(i,1); 
   Angle_single(count,4)=ele(i,2);
   Angle_single(count,5)=ele(i,3);
   
   %ensure that the points in the anglelist is anticlocwise orders
   index1=find(atom_xyz(:,1)==ele(i,1));
   index2=find(atom_xyz(:,1)==ele(i,2));
   index3=find(atom_xyz(:,1)==ele(i,3));
   
   p1=[atom_xyz(index1,4), atom_xyz(index1,5), atom_xyz(index1,6)];
   p2=[atom_xyz(index2,4), atom_xyz(index2,5), atom_xyz(index2,6)];
   p3=[atom_xyz(index3,4), atom_xyz(index3,5), atom_xyz(index3,6)];
   tri_center = (p1+p2+p3)/3.0-NP_center(atom_xyz(index1,2),:);
   tri_center(1)=0;
   tri_center(2)=0;
   a21=p2-p1;
   a31=p3-p1;
   v_nor=cross(a21,a31);
   vvv=dot(tri_center,v_nor);
   if vvv<0 % angle potentil in Fedosov 
       Angle_single(count,4)=ele(i,3);
       Angle_single(count,5)=ele(i,2);
   end     
end
Angle_single(:,1)=(1:count)';
Angle_single(:,2)=Angle_single(:,1);
toc
tic;
%combination angle
 Angle = Angle_single;
 Angle_temp = Angle_single;
 for kk=1:num_tot-1
     Angle_temp(:,1:2) = Angle_single(:,1:2) + kk*size(Angle_single,1);
     Angle_temp(:,3:5) = Angle_single(:,3:5) + kk*num_atom;
     Angle = [Angle;Angle_temp]; 
 end
toc
tic;
%store the Angle area and volume for stress free
Area_tri = zeros(size(Angle_single,1),1);
Vol_tri = zeros(size(Angle_single,1),1);

T_area = zeros(num_tot,1);
Volume = zeros(num_tot,1);
T_area_temp = zeros(num_tot,1);
Volume_temp = zeros(num_tot,1);
for kk =1:size(Angle_single,1)
       
index1 = find(atom_xyz(:,1)==Angle_single(kk,3));
index2 = find(atom_xyz(:,1)==Angle_single(kk,4));
index3 = find(atom_xyz(:,1)==Angle_single(kk,5));
   
   imol = atom_xyz(index1,2);

   p1=[atom_xyz(index1,4), atom_xyz(index1,5), atom_xyz(index1,6)];
   p2=[atom_xyz(index2,4), atom_xyz(index2,5), atom_xyz(index2,6)];
   p3=[atom_xyz(index3,4), atom_xyz(index3,5), atom_xyz(index3,6)];
   a21=p2-p1;
   a31=p3-p1;
   epsilon= cross(a21,a31);
   Area_tri(kk,1) = 0.5*norm(epsilon,2);
   tri_center = (p1+p2+p3)/3.0-NP_center(atom_xyz(index1,2),:);
   Vol_tri(kk,1) = abs(dot(tri_center,epsilon)/6);
   
   T_area_temp(imol,1) = T_area_temp(imol,1) + Area_tri(kk,1);
   Volume_temp(imol,1) =  Volume_temp(imol,1) + Vol_tri(kk,1);
end
for i=1:num_tot
T_area(i) =   T_area_temp(1,1);
Volume(i) =   Volume_temp(1,1);
end

n_Angle = size(Angle_single,1);
for i=1:num_tot
    Area_tri_temp((i-1)*n_Angle+1:i*n_Angle,1)=Area_tri;
end
Area_tri=Area_tri_temp;

T_area
Volume
toc
tic;
Dihe = zeros(1,6);
count=0;

% for mm=1:num_tot*num_atom
%     
% ID1 =find(num_points(:,1)==mm);
% ID2 =find(num_points(:,2)==mm);
% ID3 =find(num_points(:,3)==mm);
% temp=[num_points(ID1,:);num_points(ID2,:);num_points(ID3,:)]; 
% 
% for kk=1:size(temp,1)
%       id1=temp(kk,:);    
%     for k=(kk+1):size(temp,1)
%       id2=temp(k,:);
%       atoms=intersect(id1,id2);
%        % if the two triangle is connected, they will have only four unique points
%       if length(atoms)==2 && atoms(1)>=mm && atoms(2)>=mm; 
%           
%        count=count+1;
%        
%        Dihe(count,4)=atoms(1);
%        Dihe(count,5)=atoms(2);
%        atom1=setxor(id1,atoms);
%        atom4=setxor(id2,atoms);
%        Dihe(count,3)=atom1;
%        Dihe(count,6)=atom4;
%        
%        % direction jugment 
%        index1=find(atom_xyz(:,1)==atom1);
%        index2=find(atom_xyz(:,1)==atoms(1));
%        index3=find(atom_xyz(:,1)==atoms(2));
%        p1=[atom_xyz(index1,4), atom_xyz(index1,5), atom_xyz(index1,6)];
%        p2=[atom_xyz(index2,4), atom_xyz(index2,5), atom_xyz(index2,6)];
%        p3=[atom_xyz(index3,4), atom_xyz(index3,5), atom_xyz(index3,6)];
%        tri_center = (p1+p2+p3)/3.0-NP_center(atom_xyz(index1,2),:);
%        tri_center(1)=0;
%        tri_center(2)=0;
%        a21=p2-p1;
%        a31=p3-p1;
%        v_nor=cross(a21,a31);
%        vvv=dot(tri_center,v_nor); 
%        if vvv<0    %switch the position of dihedral 
%        Dihe(count,4)=atoms(2);
%        Dihe(count,5)=atoms(1);
%        end 
%       else 
%       continue;
%       end      
%     end
% end
% end 

for mm=1:num_atom
    
ID1 =find(num_points(:,1)==mm);
ID2 =find(num_points(:,2)==mm);
ID3 =find(num_points(:,3)==mm);
temp=[num_points(ID1,:);num_points(ID2,:);num_points(ID3,:)]; 

for kk=1:size(temp,1)
      id1=temp(kk,:);    
    for k=(kk+1):size(temp,1)
      id2=temp(k,:);
      atoms=intersect(id1,id2);
       % if the two triangle is connected, they will have only four unique points
      if length(atoms)==2 && atoms(1)>=mm && atoms(2)>=mm; 
          
       count=count+1;
       
       Dihe(count,4)=atoms(1);
       Dihe(count,5)=atoms(2);
       atom1=setxor(id1,atoms);
       atom4=setxor(id2,atoms);
       Dihe(count,3)=atom1;
       Dihe(count,6)=atom4;
       
       % direction jugment 
       index1=find(atom_xyz(:,1)==atom1);
       index2=find(atom_xyz(:,1)==atoms(1));
       index3=find(atom_xyz(:,1)==atoms(2));
       p1=[atom_xyz(index1,4), atom_xyz(index1,5), atom_xyz(index1,6)];
       p2=[atom_xyz(index2,4), atom_xyz(index2,5), atom_xyz(index2,6)];
       p3=[atom_xyz(index3,4), atom_xyz(index3,5), atom_xyz(index3,6)];
       tri_center = (p1+p2+p3)/3.0-NP_center(atom_xyz(index1,2),:);
       tri_center(1)=0;
       tri_center(2)=0;
       a21=p2-p1;
       a31=p3-p1;
       v_nor=cross(a21,a31);
       vvv=dot(tri_center,v_nor); 
       if vvv<0    %switch the position of dihedral 
       Dihe(count,4)=atoms(2);
       Dihe(count,5)=atoms(1);
       end 
      else 
      continue;
      end      
    end
end
end 
toc

Dihe(:,1)=(1:count)';
Dihe(:,2)=Dihe(:,1);

tic;
% store the dihedral information 
Dihe_angle = zeros(size(Dihe,1),1);
for kk =1:size(Dihe,1)
index1 = find(atom_xyz(:,1)==Dihe(kk,3));
index2 = find(atom_xyz(:,1)==Dihe(kk,4));
index3 = find(atom_xyz(:,1)==Dihe(kk,5));
index4 = find(atom_xyz(:,1)==Dihe(kk,6));

p1 =[atom_xyz(index1,4),atom_xyz(index1,5),atom_xyz(index1,6)];
p2 =[atom_xyz(index2,4),atom_xyz(index2,5),atom_xyz(index2,6)];
p3 =[atom_xyz(index3,4),atom_xyz(index3,5),atom_xyz(index3,6)];
p4 =[atom_xyz(index4,4),atom_xyz(index4,5),atom_xyz(index4,6)];
a21=p2-p1;
a31=p3-p1;
a34=p3-p4;
a24=p2-p4;

tc1=(p1+p2+p3)/3-NP_center(atom_xyz(index1,2));
tc2=(p4+p2+p3)/3-NP_center(atom_xyz(index1,2));

xi=cross(a21,a31);
xi=xi/norm(xi,2);
zeta=cross(a34,a24);
zeta=zeta/norm(zeta,2);

c=dot(xi,zeta);
if c>1.0
    c;
end
if c<-1.0
    c;
end
    
d=dot((xi-zeta),(tc1-tc2));

if d<0
    theta0 = -acos(c);
else
    theta0 = acos(c);
end
Dihe_angle(kk)=theta0*180/pi;
end

Dihe_single = Dihe;
 Dihe_temp = Dihe;
 for kk=1:num_tot-1
     Dihe_temp(:,1:2) = Dihe_single(:,1:2) + kk*size(Dihe_single,1);
     Dihe_temp(:,3:6) = Dihe_single(:,3:6) + kk*num_atom;
     Dihe = [Dihe;Dihe_temp]; 
 end
 
n_dihe = size(Dihe_angle,1);
for i=1:num_tot
    Dihe_angle_temp((i-1)*n_dihe+1:i*n_dihe,1)=Dihe_angle;
end
Dihe_angle=Dihe_angle_temp;
toc
tic;

% %  rotation of the molecules using given random angle from -2*pi to 2*pi
% rotate_sita = -pi + (pi+pi).*rand(num_tot,1);
% for num=1:num_tot
%     rotatex = [1,0,0;0,cos(rotate_sita(num)),-sin(rotate_sita(num));0,sin(rotate_sita(num)),cos(rotate_sita(num))];
%     rotatey = [cos(rotate_sita(num)),0,sin(rotate_sita(num));0,1,0;-sin(rotate_sita(num)),0,cos(rotate_sita(num))];
%     rotatez = [cos(rotate_sita(num)),-sin(rotate_sita(num)),0;sin(rotate_sita(num)),cos(rotate_sita(num)),0;0,0,1];
%     
%     for n_atom=1:num_atom
%     atom_xyz((num-1)*num_atom+n_atom,4:6) = NP_center(num,:).'+rotatex*rotatey*rotatez*((atom_xyz((num-1)*num_atom+n_atom,4:6)-NP_center(num,:)).');
%     end
% end



% %  rotation of the molecules using given random angle from 0 to pi/12; % a + (b-a).*rand(N,1).
% rotate_sita = -pi + (pi+pi).*rand(num_tot,1);
% for num=1:num_tot
%     rotatex = [1,0,0;0,cos(rotate_sita(num)),-sin(rotate_sita(num));0,sin(rotate_sita(num)),cos(rotate_sita(num))];
%     rotatey = [cos(rotate_sita(num)),0,sin(rotate_sita(num));0,1,0;-sin(rotate_sita(num)),0,cos(rotate_sita(num))];
%     rotatez = [cos(rotate_sita(num)),-sin(rotate_sita(num)),0;sin(rotate_sita(num)),cos(rotate_sita(num)),0;0,0,1];
%     
%     for n_atom=1:num_atom
%     atom_xyz((num-1)*num_atom+n_atom,4:6) = NP_center(num,:).'+rotatex*rotatey*rotatez*((atom_xyz((num-1)*num_atom+n_atom,4:6)-NP_center(num,:)).');
%     end
% end

%  rotation of the molecules using given random angle from 0 to pi/12; % a + (b-a).*rand(N,1).
rotate_sita = -pi + (pi+pi).*rand(num_tot,1);
for num=1:num_tot
    rotatex = [1,0,0;0,cos(rotate_sita(num)),-sin(rotate_sita(num));0,sin(rotate_sita(num)),cos(rotate_sita(num))];
    rotatey = [cos(rotate_sita(num)),0,sin(rotate_sita(num));0,1,0;-sin(rotate_sita(num)),0,cos(rotate_sita(num))];
    rotatez = [cos(rotate_sita(num)),-sin(rotate_sita(num)),0;sin(rotate_sita(num)),cos(rotate_sita(num)),0;0,0,1];
    
    for n_atom=1:num_atom
    atom_xyz((num-1)*num_atom+n_atom,4:6) = NP_center(num,:).'+rotatex*rotatey*rotatez*((atom_xyz((num-1)*num_atom+n_atom,4:6)-NP_center(num,:)).');
    end
end

toc
tic;
print=1;
if print==1
name1=num2str(num_tot);
name2='bac.data';
name=[name1,name2];
S_angle=180*acos((sqrt(3)*(size(atom_xyz,1)-2)-5*pi)/(sqrt(3)*(size(atom_xyz,1)-2)-3*pi))/pi;
% write the input file
fid = fopen(name,'w');
fprintf(fid,'Generation of multi rbc models.\n');
fprintf(fid,'\n');
fprintf(fid,'%d atoms\n', size(atom_xyz,1));
fprintf(fid,'%d bonds\n', size(Bond,1));
fprintf(fid,'%d angles\n',size(Angle,1));
fprintf(fid,'%d dihedrals\n', size(Dihe,1));
fprintf(fid,'%d impropers\n', 0);
fprintf(fid,'\n');
fprintf(fid,'%d atom types\n', 1);
fprintf(fid,'%d bond types\n', size(Bond,1));
fprintf(fid,'%d angle types\n', size(Angle,1));
fprintf(fid,'%d dihedral types\n',size(Dihe,1));
fprintf(fid,'%d improper types\n', 0);
fprintf(fid,'\n');
fprintf(fid,'%4.3f %4.3f xlo xhi\n', box(1,1), box(1,2));
fprintf(fid,'%4.3f %4.3f ylo yhi\n', box(2,1), box(2,2));
fprintf(fid,'%4.3f %4.3f zlo zhi\n', box(3,1), box(3,2));
fprintf(fid,'\n');
fprintf(fid,'Masses\n');
fprintf(fid,'\n');
for mas_point=1:1
fprintf(fid,'%d %f \n',mas_point, 1.0);
end


% fprintf(fid,'\n');
% fprintf(fid,'Bond Coeffs #wlc\n');
% fprintf(fid,'\n');
% scale=2.0;
% for i=1:size(Bond,1)
%     %normLength=Len_bond(i,1)/1.2;
%     lmax = scale*Len_bond(i,1);
%     %temp(i,1:8)=[i, 1.1e-4,1/scale,lmax,0.01,2.0,30.0,90.0];
%     fprintf(fid,'%d %d %10.8f %10.8f %d %10.8f %10.8f %10.8f\n',Bond(i,2), 6.57e-4,1/scale,lmax,0.1,2.0,30.0,90.0);
% end


fprintf(fid,'\n');
fprintf(fid,'Bond Coeffs #harmonic\n');
fprintf(fid,'\n');
for i=1:size(Bond,1)    
fprintf(fid,'%d %10.8f %10.8f\n',i, bondcoeff, Len_bond(i,1)); 
end


fprintf(fid,'\n');
% fprintf(fid,'Angle Coeffs  #rbc\n');
% fprintf(fid,'\n');
% for i =1:size(Angle,1)
% %index1 = find(atom_xyz(:,1)==Angle(i,3));
% %imol =   atom_xyz(index1,2);
% fprintf(fid,'%d %10.8f %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f\n', Angle(i,2), 0.0, 1, 0.075, T_area(1,1),0.966,Volume(1,1), 3.67, Area_tri(i,1));
% end
% 
% fprintf(fid,'\n');
% fprintf(fid,'Dihedral Coeffs  #bend\n');
% fprintf(fid,'\n');
% for i =1:size(Dihe,1)
% fprintf(fid,'%d %d %10.8f\n', Dihe(i,2),0.133,Dihe_angle(i));
% end

fprintf(fid,'\n');
fprintf(fid,'Angle Coeffs  #rbc\n');
fprintf(fid,'\n');
for i =1:size(Angle,1)
%index1 = find(atom_xyz(:,1)==Angle(i,3));
%imol =   atom_xyz(index1,2);
fprintf(fid,'%d %10.8f %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f\n', Angle(i,2), 0.0, 1, 0.075, T_area(1,1),0.395,Volume(1,1), 3.67, Area_tri(i,1));
end
fprintf(fid,'\n');
fprintf(fid,'Dihedral Coeffs  #bend\n');
fprintf(fid,'\n');
for i =1:size(Dihe,1)
fprintf(fid,'%d %d %10.8f\n', Dihe(i,2),0.7937,Dihe_angle(i));
end


fprintf(fid,'\n');
fprintf(fid,'Atoms\n');
fprintf(fid,'\n');
fprintf(fid,'%d %d %d %10.8f %10.8f %10.8f\n',atom_xyz(:,1:6).');  
fprintf(fid,'\n');
fprintf(fid,'Bonds\n');
fprintf(fid,'\n');
fprintf(fid,'%d %d %d %d\n',Bond(:,1:end).');
fprintf(fid,'\n');
fprintf(fid,'Angles\n');
fprintf(fid,'\n');
fprintf(fid,'%d %d %d %d %d\n',Angle(:,1:end).');
fprintf(fid,'\n');
fprintf(fid,'Dihedrals\n');
fprintf(fid,'\n');
fprintf(fid,'%d %d %d %d %d %d\n',Dihe(:,1:end).');
fclose(fid);
end
toc