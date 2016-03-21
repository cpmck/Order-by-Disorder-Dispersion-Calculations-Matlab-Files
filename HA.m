function returnHA = HA( jA, jB, H, phi, IAA, IAB )
%HA Calculates mean field Hamiltonian A
%   HA() calculates the mean field hamiltonian at A-position ion

%import Hion and angular momentum operators from base script
Hion = evalin('base','Hion');
Jx = evalin('base','Jx');
Jy = evalin('base','Jy');
Jz = evalin('base','Jz');

%calculae exchange interaction between A-A sites and A-B sites
%jA are the moments at A and B, H is the magnetic field
jAB = H + IAA*jA + IAB*jB;

%construct the mean field hamiltonian
HHA = zeros(length(Jx),length(Jx));
for i = 1:length(Jx)
    for j = 1:length(Jx)   
        HHA(i,j) = jAB(1,1)*Jx(i,j) + jAB(2,1)*Jy(i,j) + jAB(3,1)*Jz(i,j);
    end
end

%return the mean field Hamiltonian for A site
returnHA = Hion(phi) - HHA;

end

