%clear all;
format LONG;

%J=4 model
%number of levels (2J+1) = 9
levels = 9;

%Choose coordinate frame Should be same as HCEF
a = [0,1,0];
b = [0,0,1];
c = [1,0,0];

%Define Hcef
B = [0.2485418 0.7430229 0.1039840 -0.0027785 -0.0143670 0.0861763 -0.0011329 -0.0011329 0.0000086 -0.0002414 0.0002160 0.0001899 -0.0001599 0.0000046];
Hcef = B(1,1)*OM(2,0)+B(1,2)*OM(2,2)+B(1,3)*OM(2,-2)+B(1,4)*OM(4,0)+B(1,5)*OM(4,2)+B(1,6)*OM(4,-2)+B(1,7)*OM(4,4)+B(1,8)*OM(4,-4)+B(1,9)*OM(6,0)+B(1,10)*OM(6,2)+B(1,11)*OM(6,-2)+B(1,12)*OM(6,4)+B(1,13)*OM(6,-4)+B(1,14)*OM(6,6)+B(1,15)*OM(6,-6);


%rotation matrices : see 'Magnetism of PrPtAl' Kinoshita.et.al
U = @(phi) expm(-1i*phi*JOp(3));
Uinv = @(phi) inv(U(phi));
Hion = @(phi) Uinv(phi)*Hcef*U(phi);

%Define Jx, Jy and Jz from function JOp()
Jx = JOp(1);
Jy = JOp(2);
Jz = JOp(3);
J = [Jx;Jy;Jz];

%define Pr3+ position coordinates from H Kitazawa et al.
y = 1/4;
x = 0.0254;
z = 0.318;

%define atomic positions (4 atoms in the unit cell) Kitazawa et al.
atomA = [z,x,1/4];
atomAp = [-z,-x,-1/4];
atomBp = [z-1/2,1/2-x,-1/4];
atomB = [-z+1/2,x-1/2,1/4];

%define lattice parameters Kitazawa et al.
latA = 7.114;
latB = 4.46;
latC = 7.785;

%calculate reciprocal lattice vectors
rlu = 2*pi*[1/latC;1/latA;1/latB];

%define positions of nearest neighbours in neighbouring unit cells
r1 = c + [-2*z, -2*x, 0];
r2 = [-2*z+1/2,0,0];
r3 = -c +[2*z,-2*x,0];
r4 = [0,1/2-2*x,0];


%Calculate energies via DMD formalism Check if its working first
energiesCalc = Energies_J4(32*pi/180, -0.21, 0.19, 0.1, 0.092, 0, 0.52, 0.05, 0.15, 0.22, 0, 11.6/7, [2; 0; 0], [0; 0; 0]);

%select range over which to calculate
h = [-2.5:0.05:2.5];

%do the calculation and save results to matrix dispersion
dispersion = zeros(length(energiesCalc),2,length(h));
for i = 1:length(h)
    dispersion(:,:,i) = Energies_J4(32*pi/180, -0.21, 0.19, 0.1, 0.092, 0, 0.52, 0.05, 0.15, 0.22, 0, 11.6/5.8, [2; h(1,i); 0], [0; 0; 0]);
end

