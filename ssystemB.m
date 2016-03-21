function eigenSystemOut = ssystemB( jA,jB,H,phi,IAA,IAB )
%ssystem: solve eigensystem and sort 

[VB,EB] = eig(HB(jA,jB,H,phi,IAA,IAB));


%convert to same as Mathematica for debugging
%{
Mconvert = [1,0,0;0,1,0;0,0,-1];
VB = Mconvert*transpose(VB);
VB = transpose(VB);

%}
Eigenvals = diag(real(EB));
Emin = min(Eigenvals);
EEigenvals = Eigenvals - Emin;
eigenSystem = [transpose(EEigenvals);transpose(-EEigenvals+Eigenvals);VB];
eigenSystem = sortrows(transpose(eigenSystem),1);
eigenSystemOut = [eigenSystem(:,1)+eigenSystem(:,2),eigenSystem(:,3:end)];


end

