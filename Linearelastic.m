clc
clearvars
%%%%%% finding global K %%%%%%%%
%there is no boundary or source force
NdCoord=[0 0;2 0.5; 2 1; 0 1];
zeta=[-1/sqrt(3), 1/sqrt(3)];
etta=[-1/sqrt(3), 1/sqrt(3)];
E=3e11;
niu=0.3;
x=E/(1-niu*niu);
D= x*[1  niu    0;
     niu  1     0;
      0   0 (1-niu)/2];
K=zeros(8);
for i=1:2
    for j=1:2
        derivative_shape=0.25*[etta(j)-1 1-etta(j) 1+etta(j) -etta(j)-1;
                               zeta(i)-1 -zeta(i)-1 1+zeta(i) 1-zeta(i)];
        J=derivative_shape*NdCoord;

        B_sigma= inv(J)*derivative_shape;

        B=[B_sigma(1,1)     0       B_sigma(1,2)     0      B_sigma(1,3)      0         B_sigma(1,4)      0;
                   0         B_sigma(2,1)     0        B_sigma(2,2)    0         B_sigma(2,3)     0        B_sigma(2,4);
                B_sigma(2,1) B_sigma(1,1) B_sigma(2,2) B_sigma(1,2) B_sigma(2,3) B_sigma(1,3) B_sigma(2,4) B_sigma(1,4)];

        K=K+det(J)*B'*D*B;
    end
end
global_K=K
%{
%swaping the column 3 and 4 for transforming the all essential boundary at
%-the top of the matrix
global_K(:,[3,4])=global_K(:,[4,3]); % it will change the essential boundary displacement at the top
global_K([3,4],:)=global_K([4,3],:); % it will change the unknown reaction force at the top
modified_K= global_K


%Partitioning the stiffness matrix
KE= modified_K(1:3,1:3);
KEF= modified_K(1:3,4:8);
%here swaping change the symmetry of the matrix So there is no  Transpose
%of KEF.This will be written as KM
KM=modified_K(4:8,1:3);
KF= modified_K(4:8,4:8);
rF=[500;-1500;-1500;-500;-500];
dE=[0;0;0];
dF=inv(KF)*(rF-KM*dE)

%reaction force calculation
BC_force=[1500; 1500; 500] %external boundary force
%rE=[1500+r1x; 1500+r1y; 500+r2y];
rE=KE*dE+KEF*dF %total force 
Reaction_force=rE-BC_force


%Checking for the displacement is correct or not
d=[0;0;-.1230;0;-.1849;.0638;-.1256;.1124];
F=global_K*d;
%lets check what will be the actual displacement value of since we know all
%8 forces
f=[1500; 1500;500; 500;-1500;-1500;-500;-500];
du=inv(global_K)*f;

%}




