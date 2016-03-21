
function O = OM(l,m)
format long;
% OM() returns the Stevens equivalent operator O_l^m
% source:  E Bauer and M Rotter. Magnetism of complex metallic alloys: crystalline electric field effects (2007)
J = 4;

%row vector of m values
setVals = [-4:1:4];

%calculate values of raising and lowering operators
Jp = zeros(9,9);
Jm = zeros(9,9);
for i = 1:9 
    for j = 1:9
        if j-i == 1
            Jp(i,j) = sqrt(J*(J+1) - (i-5)*(i-4));
        elseif i-j == 1
            Jm(i,j) = sqrt(J*(J+1) - (i-5)*(i-6));
        end
    end
end


%create Jz operator with m vals along diagonal
%row vector of m values
setVals = [-4:1:4];
Jz = diag(setVals);


Jp2 = Jp*Jp;
Jm2 = Jm*Jm;
Jz2 = Jz*Jz;

X = J*(J+1) * eye(2*J+1);

I = eye(2*J+1);

Jm3 = Jm2*Jm;
Jp3 = Jp2*Jp;
Jz3 = Jz2*Jz;

if l == 2 && m == 2
O = (Jp*Jp + Jm*Jm)/2;
elseif l == 2 && m == 0
O = 3*Jz*Jz - X;
elseif l == 2 && m == 1
O = (1/4) * (Jz*(Jp+Jm) + (Jp+Jm)*Jz);
elseif l == 2 && m == -1
O = (1/(4*1i)) * (Jz*(Jp-Jm) + (Jp-Jm)*Jz);
elseif l == 2 && m == -2
O = (Jp*Jp - Jm*Jm)/(2*1i);
elseif l == 4 && m == 4
O = (Jp2*Jp2 + Jm2*Jm2)/2;
elseif l == 4 && m == 2
O = ((7*Jz2 - X - 5*I)*(Jp2 + Jm2) + (Jp2 + Jm2) * (7*Jz2 - X - 5*I))/4;
elseif l == 4 && m == -2
O = ((7*Jz2 - X - 5*I)*(Jp2-Jm2) + (Jp2 - Jm2) * (7*Jz2 - X - 5*I))/(4*1i);
elseif l == 4 && m == 0
O = 35*Jz2*Jz2 - (30*X - 25*I)*Jz2 + 3*X*X - 6*X;
elseif l == 4 && m == -4
O = (Jp2*Jp2 -Jm2*Jm2)/(2*1i);
elseif l == 4 && m == 1
O = (1/4)*((Jp+Jm)*(7*Jz3 - (3*X + I)*Jz) + (7*Jz3 - (3*X + I)*Jz)*(Jp+Jm));
elseif l == 4 && m == 3
O = (1/4)*((Jp3+Jm3)*Jz + Jz*(Jp3+Jm3));
elseif l == 4 && m == -1
O = (1/(4*1i))*((Jp-Jm)*(7*Jz3 - (3*X + I)*Jz) + (7*Jz3 - (3*X + I)*Jz)*(Jp-Jm));
elseif l == 4 && m == -3
O = (1/(4*1i))*((Jp3 - Jm3)*Jz + Jz*(Jp3 -Jm3));
elseif l == 6 && m == 0
O = 231*Jz2*(Jz2*Jz2) - (315*X - 735*I)*(Jz2*Jz2) + (105*X*X - 525*X + 294*I)*Jz2 - 5*X*(X*X) + 40*X*X - 60*X;
elseif l == 6 && m == 2
O = (1/4)*((33* Jz2*Jz2 - 18*Jz2*X - 123*Jz2 + X*X + 10*X + 102*I)*(Jp2+Jm2) +(Jp2+Jm2)*(33*Jz2*Jz2 - 18*Jz2*X - 123*Jz2 + X*X + 10*X + 102*I));
elseif l == 6 && m == -2
O = (1/(4*1i))*( (33*Jz2*Jz2 - 18*Jz2*X - 123*Jz2 + X*X + 10*X + 102*I)*(Jp2-Jm2)+(Jp2-Jm2)*(33*Jz2*Jz2 - 18*Jz2*X - 123*Jz2 + X*X + 10*X +102*I) );
elseif l == 6 && m == 4
O = (1/4)*((11*Jz2 - X - 38*I)*(Jp2*Jp2+Jm2*Jm2) + (Jp2*Jp2+Jm2*Jm2)*(11*Jz2 -X - 38*I));
elseif l == 6 && m == -4
O = (1/(4*1i))*((11*Jz2-X-38*I)*(Jp2*Jp2-Jm2*Jm2)+(Jp2*Jp2-Jm2*Jm2)*(11*Jz2-X-38*I));
elseif l == 6 && m == 6
O = (1/2)*((Jp2*Jp2)*Jp2 + (Jm2*Jm2)*Jm2);
elseif l == 6 && m == -6
O = (1/(2*1i))*((Jp2*Jp2)*Jp2 - (Jm2*Jm2)*Jm2);
elseif l == 6 && m == 1
O = (1/4)*((Jp+Jm)*(33*Jz2*Jz3 - (30*X - 15*I)*Jz3 + (5*X*X - 10*X + 12*I)*Jz) + (33*Jz2*Jz3 - (30*X - 15*I)*Jz3 + (5*X*X - 10*X + 12*I)*Jz)*(Jp+Jm));
elseif l == 6 && m == 3
O = (1/4)*((Jp3+Jm3)*(11*Jz3 - (3*X + 59*I)*Jz) + ((11*Jz3 - (3*X +59*I)*Jz))*(Jp3+Jm3));
elseif l == 6 && m == 5
O = (1/4)*((Jp3*Jp2 + Jm3*Jm2)*Jz + Jz*(Jp3*Jp2 +Jm3*Jm2));
elseif l == 6 && m == -1
O = (1/(4*1i))*((Jp-Jm)*(33*Jz2*Jz3 - (30*X - 15*I)*Jz3 + (5*X*X -10*X +12*I)*Jz) +(33*Jz2*Jz3 - (30*X-15*I)*Jz3 + (5*X*X - 10*X + 12*I)*Jz)*(Jp-Jm));
elseif l == 6 && m == -3
O = (1/(4*1i))*((Jp3 -Jm3)*(11*Jz3 - (3*X + 59*I)*Jz) + ((11*Jz3 - (3*X + 59*I)*Jz))*(Jp3-Jm3));
elseif l == 6 && m == -5
O = (1/(4*1i))*((Jp3*Jp2-Jm3*Jm2)*Jz + Jz*(Jp3*Jp2 - Jm3*Jm2));
else
    disp('no Stevens Operator')
end    


end

