function eigenSystemOut = ssystemA( jA,jB,H,phi,IAA,IAB )
%ssystem: solve eigensystem and sort 

%returns sorted eigenvalues and eigenvectors
%[eValue1,eVec1;eValue2,eVec2;......;eValueN, eVecN]

%columns of VA are right eigenvectros

[VA,EA] = eig(HA(jA,jB,H,phi,IAA,IAB));

%convert to same as Mathematica for debugging
%{
Mconvert = [-1,0,0;0,1,0;0,0,1];
VA = Mconvert*transpose(VA);
VA = transpose(VA);
%}

Eigenvals = real(diag(EA));
Emin = min(Eigenvals);
EEigenvals = Eigenvals - Emin;
eigenSystem = [EEigenvals';(-EEigenvals+Eigenvals)';VA];
eigenSystemSort = sortrows(transpose(eigenSystem),1);
eigenSystemOut = [eigenSystemSort(:,1)+eigenSystemSort(:,2),eigenSystemSort(:,3:end)];

end

