%Author: Georgios Stergianakos
%Date: 26 May 2021
%Description: This function performs inner room hourly temperature 
%calculations, according to the ISO EN 13790 - simple hourly method. 
%All variables below are coherent with the respective ones from 
%ISO EN 13790 methodology 

function [Tm,Ts,Tair,Top] = Tair_calc(I,D,Hve,His,Hw,Hms,Hem,Fm,Fst,Fia,the,thsup,thmt,FHCnd,Cm)
    H1 = 1/((1/Hve)+(1/His)); H2 = H1+Hw; H3 = 1/((1/H2)+(1/Hms));
    Fmtot = Fm(I,D)+Hem*the(I,D)+H3*(Fst(I,D)+Hw*the(I,D)+H1*(((Fia(I,D)+FHCnd(I,D))/Hve)+thsup(I,D)))/H2;
    thmt2 = (thmt*((Cm/3600)-0.5*(H3+Hem))+Fmtot)/((Cm/3600)+0.5*(H3+Hem));
    Tm = (thmt+thmt2)/2;
    Ts = (Hms*Tm+Fst(I,D)+Hw*the(I,D)+H1*(thsup(I,D)+(Fia(I,D)+FHCnd(I,D))/Hve))/(Hms+Hw+H1);
    Tair = (His*Ts+Hve*thsup(I,D)+Fia(I,D)+FHCnd(I,D))/(His+Hve);
    Top = 0.3*Tair+0.7*Ts;
end
