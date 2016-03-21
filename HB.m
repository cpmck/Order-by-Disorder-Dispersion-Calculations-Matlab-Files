function returnHB = HB( jA, jB, H, phi, IAA, IAB )
%HB Calculates mean field Hamiltonian B
%   HB() calculates the mean field hamiltonian at A-position ion

%import from base script
Hion = evalin('base','Hion');
Jx = evalin('base','Jx');
Jy = evalin('base','Jy');
Jz = evalin('base','Jz');

%calculate exchange interaction between A-A sites and A-B sites
jAB = H + IAA*jB + IAB*jA;

%construct mean field Hamiltonian
HHB = zeros(length(Jx),length(Jx));
for i = 1:length(Jx)
    for j = 1:length(Jx)   
        HHB(i,j) = jAB(1,1)*Jx(i,j) + jAB(2,1)*Jy(i,j) + jAB(3,1)*Jz(i,j);  
    end
end

%return mean field hamiltonian
returnHB = conj(Hion(phi)) - HHB;
end
