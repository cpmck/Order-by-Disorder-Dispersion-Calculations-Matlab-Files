function newj = newj( jA,jB,H,beta,phi,IAA,IAB )
%newj calculate moment <J> from mean field.
%Mean fields are defiend by HA and HB 
%HA is the mean feild at jA and HB is the mean field at jB
%Expectation values are calculated using Trace[] (see density matrices)

%evaluate from base script
Jx = evalin('base','Jx');
Jy = evalin('base','Jy');
Jz = evalin('base','Jz');

%define function to calculate numerators of expeactation values
numerHA = @(J) trace(expm( -beta*HA(jA,jB,H,phi,IAA,IAB) )*J );
numerHB = @(J) trace(expm( -beta*HB(jA,jB,H,phi,IAA,IAB) )*J );

%calculate the partition functions
ZHA = trace(expm( -beta*HA(jA,jB,H,phi,IAA,IAB) ) );
ZHB = trace(expm( -beta*HB(jA,jB,H,phi,IAA,IAB) ) );

%calculate thermal expeactation values of moments
newjA = [numerHA(Jx),numerHA(Jy),numerHA(Jz)]/ZHA;
newjB = [numerHB(Jx),numerHB(Jy),numerHB(Jz)]/ZHB;

%return moments
newj = [newjA',newjB'];

end

