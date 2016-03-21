function EnergiesOut = Energies_J4( phi, JAA1, JAA2, JAA3, JAB1, JAB2, JAAp1,JAAp2,JAAp3,JABp1,JABp2, beta, q ,H )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

levels = evalin('base','levels');
take = 0;
trans = zeros(levels*(levels-1),2);

for i = 1:levels
    for j = 1:levels
        if i == j
            take = take+1;
        elseif i ~= j
            numb = levels*(i-1)+j-take;
            trans(numb,1) = i;
            trans(numb,2) = j;
        end
    end 
end


%Evaluate from base
a = evalin('base','a');
b = evalin('base','b');
c = evalin('base','c');

%evalueate from base
r1 = evalin('base','r1');
r2 = evalin('base','r2');
r3 = evalin('base','r3');
r4 = evalin('base','r4');

%evaluate from base
rlu = evalin('base','rlu');

%evalueate from base
atomA = evalin('base','atomA');
atomB = evalin('base','atomB');
atomAp = evalin('base','atomAp');
atomBp = evalin('base','atomBp');


%define Jex Matrix within this function Energies()
%The functions JAA to JApBp were calculated by Professor Andrew Huxley
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
JAA = @(Q) JAA1*cos(2*pi*dot(Q,b)) + JAA2*cos(2*pi*dot(Q,a)) + JAA3*cos(2*pi*dot(Q,c));
JAB = @(Q) JAB1*cos(2*pi*dot(Q,a)/2)*exp(2*pi*1i*dot(Q,r1)) +JAB2*cos(2*pi*dot(Q,a))*exp(2*pi*1i*dot(Q,(c+r2)));
JAAp = @(Q) JAAp1*cos(2*pi*dot(Q,b)/2)*exp(2*pi*1i*dot(Q,r1)) + JAAp2*cos(2*pi*dot(Q,b)/2)*exp(2*pi*1i*dot(Q,(r1-c))) +JAAp3*cos(2*pi*dot(Q,b)/2)*exp(2*pi*1i*dot(Q,(r1+a)));
JABp = @(Q) JABp1*exp(2*pi*1i*dot(Q,r4))*cos(2*pi*dot(Q,b)/2)*cos(2*pi*dot(Q,c)/2)+JABp2*exp(2*pi*1i*dot(Q,(r4-a)))*cos(2*pi*dot(Q,b)/2)*cos(2*pi*dot(Q,c)/2);
JBAp = @(Q) JABp(Q);
JBBp = @(Q) JAAp1*cos(2*pi*dot(Q,b)/2)*exp(2*pi*1i*dot(Q,r3))+JAAp2*cos(2*pi*dot(Q,b)/2)*exp(2*pi*1i*dot(Q,(c+r3)))+JAAp3*exp(2*pi*1i*dot(Q,(a+r3)));
JApBp = @(Q) conj(JAB(Q));


%Construct Jex for use in A()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jex1 = [JAA(q), JAB(q), JAAp(q), JABp(q)];
Jex2 = [conj(JAB(q)), JAA(q), JABp(q), JBBp(q)];
Jex3 = [conj(JAAp(q)), conj(JBAp(q)), JAA(q), JApBp(q)];
Jex4 = [conj(JABp(q)), conj(JBBp(q)), conj(JApBp(q)), JAA(q)];
Jex = [Jex1;Jex2;Jex3;Jex4];

%construct nL and nK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nL = JAA([0;0;0]) + JAAp([0;0;0]);
nK = JAB([0;0;0]) + JABp([0;0;0]);


%calculate moments for use in function Gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculating Moments');
momCalc = moments( H, beta, phi, nL, nK );



%calculate eigensystem for use in functions Gamma ans U
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculating H energy eigensystem');
ssystemACalc = ssystemA(momCalc(1,:)',momCalc(2,:)',H,phi,nL,nK);
ssystemBCalc = ssystemB(momCalc(1,:)',momCalc(2,:)',H,phi,nL,nK);


%Calculate Boltamann weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculating Occupation Weights');
nA1 = zeros(levels,1);
nA = zeros(levels,1);
nAZ = 0;
for i = 1:levels
    nA1(i,1) = exp(-beta*ssystemACalc(i,1));
    nAZ = nAZ + nA1(i,1);
end
nA = nA1/nAZ; %note: nA(end,1) is probability of occupation of ground state


%Calculate Lambda for sign of transition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculating LAMBDA');
Lambda = zeros(length(trans),length(trans),4,4);

for i = 1:4  
    for j = 1:4      
        for k = 1:length(trans)    
            for l = 1:length(trans)
                if trans(k,1)<=trans(k,2) && i == j && l == k  Lambda(l,k,j,i) = 1;
                elseif  trans(k,1)>trans(k,2) && i == j && l == k  Lambda(l,k,j,i) = -1;
                end
            end
        end
    end
end

L1 = [Lambda(:,:,1,1),Lambda(:,:,1,2),Lambda(:,:,1,3),Lambda(:,:,1,4)];
L2 = [Lambda(:,:,2,1),Lambda(:,:,2,2),Lambda(:,:,2,3),Lambda(:,:,2,4)];
L3 = [Lambda(:,:,3,1),Lambda(:,:,3,2),Lambda(:,:,3,3),Lambda(:,:,3,4)];
L4 = [Lambda(:,:,4,1),Lambda(:,:,4,2),Lambda(:,:,4,3),Lambda(:,:,4,4)];
LAMBDA = [L1;L2;L3;L4];



%construct vU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculating vUD');
vUD = zeros(length(trans),3,4);
for i = 1:4
    for k = 1:length(trans)
        vUD(k,:,i) = U(trans(k,1),trans(k,2),i);
    end   
end

vU1 = vUD(:,:,1);
vU2 = vUD(:,:,2);
vU3 = vUD(:,:,3);
vU4 = vUD(:,:,4);
vU = [vU1;vU2;vU3;vU4];


%construct kspace atom positions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
atoms = [exp(2*pi*1i*dot(q,atomA));exp(2*pi*1i*dot(q,atomB));exp(2*pi*1i*dot(q,atomAp));exp(2*pi*1i*dot(q,atomBp))];


%construct Generalized matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculating A Matrix');

AmD = zeros(length(trans),length(trans), 4, 4);
for i = 1:4
    for j = 1:4
        for k = 1:length(trans)
            for l = 1:length(trans)
                
                AmD(l,k,j,i) = A(trans(l,1),trans(l,2),trans(k,1),trans(k,2),j,i);
                
            end
        end
    end
end


Am1 = [AmD(:,:,1,1),AmD(:,:,1,2),AmD(:,:,1,3),AmD(:,:,1,4)];
Am2 = [AmD(:,:,2,1),AmD(:,:,2,2),AmD(:,:,2,3),AmD(:,:,2,4)];
Am3 = [AmD(:,:,3,1),AmD(:,:,3,2),AmD(:,:,3,3),AmD(:,:,3,4)];
Am4 = [AmD(:,:,4,1),AmD(:,:,4,2),AmD(:,:,4,3),AmD(:,:,4,4)];
Am = [Am1;Am2;Am3;Am4]


%Calculate Generlaized matrix for diagonalisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LAm = LAMBDA*Am;


%construct vGamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculating vGamma');
vGammaD = zeros(length(trans),1,4);
for i = 1:4
    for k = 1:length(trans)
        vGammaD(k,:,i) = Gamma(trans(k,1),trans(k,2),i);
    end   
end

vGamma = [vGammaD(:,:,1);vGammaD(:,:,2);vGammaD(:,:,3);vGammaD(:,:,4)];


%calculate Eigensystem using generalised eignevalue problem function **which
%does not normailise eigenvectors**
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('Solving Generalzed Eigenvalue Problem');
[VV,EE] = eig(Am,LAMBDA);


Eigenvals = real(diag(EE));
Emin = min(Eigenvals);
EEigenvals = Eigenvals - Emin;
eigenSystem = [transpose(EEigenvals);transpose(-EEigenvals+Eigenvals);VV];
eigenSystem = sortrows(transpose(eigenSystem),1);
eigenSystemOut = [eigenSystem(:,1)+eigenSystem(:,2),eigenSystem(:,3:end)];


%sorted eigenvectors in rows
eigenVectors = eigenSystemOut(:,2:end);

%Column of sorted energy eigenvalues
eigenValues = eigenSystemOut(:,1);


%sTm = eigenvectors in columns and column of energy eigenvalues
sTm = transpose(eigenVectors);
nsTm = sTm;
energies = eigenValues;



%Normalise eigenvectors
for i = 1:length(sTm)
    nsTm(:,i) = sTm(:,i)/(sqrt(sTm(:,i)'*Am*sTm(:,i)));
end


%Calculate qA = q*rlu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qA = q.*rlu;

%calculate selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
selection = eye(3) - kron(transpose(qA),qA)/(transpose(qA)*qA);


%Caclculate Intensity
disp('Calculating Intensity');
intensity = zeros(4*length(trans),1);

for level = 1:4*length(trans)

    int = [nsTm(:,level).*vU(:,1).*conj(sqrt(vGamma)),nsTm(:,level).*vU(:,2).*conj(sqrt(vGamma)),nsTm(:,level).*vU(:,3).*conj(sqrt(vGamma))];
    intSum = sum(int,1);
    intHold = selection.*kron(intSum,intSum');
    intHoldFlatten = [intHold(:,1);intHold(:,2);intHold(:,3)];
    occup = energies(level,1)/(1-exp(-beta*energies(level,1)));
    intensity(level,1) = Chop(sum(intHoldFlatten)*occup);
    
end

EnergiesOut = zeros(length(trans),2);
num = 1;

for i = 2*length(trans)+1:4*length(trans)
EnergiesOut(num,1) = energies(i,1);
EnergiesOut(num,2) = intensity(i,1);
num = num+1;
end

end

