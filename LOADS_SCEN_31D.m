%Author: Georgios Stergianakos
%Date: 26 May 2021
%Description: This function performs analytical required hourly thermal power calculations,
%according to the ISO EN 13790 - simple hourly method. All variables below 
%are coherent with the respective ones from ISO EN 13790 methodology 
%and can be looked into.
function [Tairac, FHCndac, FHCndun,Tm2,Vinf] = LOADS_SCEN_31D(Af,Apar,Aw,Awpar,Awpor,Asalot,Uoppar,Uw,Uf1,Uf2,Ta,Isol,TintH,TintC,FHmax,FCmax,Fshob,Fgl,asc,ggl,Uop,Uor,Udap,Ror,Radap,Raop,Cm,Am,Vinf)
    %% Calculation of the final transmission coefficients H = U*A (W/K)
    Hoppar = Uoppar*(Apar-Aw);
    Hw = Uw*45.6; Hf1 = Af*Uf1; Hf2 = Af*Uf2;
    Hop = Hoppar + Hf1 + Hf2; 
    Hve = Vinf/3600*1200;
    At = 4.5*2*Af; Hms = 9.1*Am; Hem = 1/((1/Hop)-(1/Hms)); His=3.45*At;
    %% Solar heat gain calculation
    Asolw = Fgl*ggl*(1-0.2)*Aw;
    Asolop = asc*Ror*Uor*Af+asc*Radap*Udap*Af+asc*Raop*Uop*(Apar-Aw);
    for i = 1:24
        for j=1:28
            Fsolop(i,j) = Fshob*Asolop*Isol(i,j)*0.0036;
            Fsolw(i,j) = Fshob*Asolw*Isol(i,j)*0.0036;
            Fsol(i,j) = Fsolop(i,j)+Fsolw(i,j);
        end
    end
    %% Internal heat gain calculation
    ss = 0;
    for j = 1:28
        ss = ss+1;
        for i = 1:24
            if  ss ~= 6 && ss ~= 7
                if (i >= 8 && i <= 17)
                    Qint(i,j) = (8*Asalot+1*(1-Asalot))*2*Af;
                elseif i>18 && i<24
                    Qint(i,j) = (20*Asalot+1*(1-Asalot))*2*Af;
                else
                    Qint(i,j) = (2*Asalot+6*(1-Asalot))*2*Af;
                end
            else
                if (i >= 8 && i <= 17)
                    Qint(i,j) = (8*Asalot+2*(1-Asalot))*2*Af;
                elseif i>18 && i<24
                    Qint(i,j) = (20*Asalot+4*(1-Asalot))*2*Af;
                else
                    Qint(i,j) = (2*Asalot+6*(1-Asalot))*2*Af;
                end
            end
        end
        if ss == 7
            ss=0;
        end
    end
    Fint = Qint*0.0036;
    %% Total inner rrom gains: Internal & solar gains
    Fia = 0.5*Fint;
    Fm = Am/At*(0.5*Fint+Fsol);
    Fst = (1-(Am/At)-Hw/(9.1*At))*(0.5*Fint+Fsol);
    %% Air and operative temperature - Thermal power calculation (Q)
    thmt1=10; Tm = thmt1;
    for j = 1:28
        for i = 1:24
            thmt1 = Tm;
            FHCnd(i,j) = 0;
            [Tm,Ts,Tair,Top] = Tair_calc(i,j,Hve,His,Hw,Hms,Hem,Fm,Fst,Fia,Ta,Ta,thmt1,FHCnd,Cm);
            Tair0(i,j) = Tair;
            Tm2(i,j) = Tm;
            if Tair0(i,j)<=TintC(i,j) && Tair0(i,j)>=TintH(i,j) 
                FHCndac(i,j) = 0;
                FHCndun(i,j) = 0;
                Tairac(i,j) = Tair0(i,j);
            elseif Tair0(i,j)>TintC(i,j)
                Tairset = TintC(i,j);
                FHCnd10(i,j) = -10*2*Af;
                FHCnd(i,j) = FHCnd10(i,j);
                [Tm,Ts,Tair,Top] = Tair_calc(i,j,Hve,His,Hw,Hms,Hem,Fm,Fst,Fia,Ta,Ta,thmt1,FHCnd,Cm);
                Tm2(i,j) = Tm;
                Tair10(i,j) = Tair;
                FHCndun(i,j) = FHCnd10(i,j).*((Tairset-Tair0(i,j))./(Tair10(i,j)-Tair0(i,j)));
                if FHCndun(i,j)>=FCmax
                    FHCndac(i,j) = FHCndun(i,j);
                    [Tm,Ts,Tair,Top] = Tair_calc(i,j,Hve,His,Hw,Hms,Hem,Fm,Fst,Fia,Ta,Ta,thmt1,FHCndac,Cm); %%%%
                    Tm2(i,j) = Tm;
                    Tairac(i,j) = Tairset;
                else
                    FHCndac(i,j) = FCmax;
                    [Tm,Ts,Tair,Top] = Tair_calc(i,j,Hve,His,Hw,Hms,Hem,Fm,Fst,Fia,Ta,Ta,thmt1,FHCndac,Cm);
                    Tm2(i,j) = Tm;
                    Tairac(i,j)=Tair;
                end
            elseif Tair0(i,j)<TintH(i,j)
                Tairset = TintH(i,j);
                FHCnd10(i,j) = 10*2*Af;
                FHCnd(i,j) = FHCnd10(i,j);
                [Tm,Ts,Tair,Top] = Tair_calc(i,j,Hve,His,Hw,Hms,Hem,Fm,Fst,Fia,Ta,Ta,thmt1,FHCnd,Cm);
                Tm2(i,j) = Tm;
                Tair10(i,j) = Tair;
                FHCndun(i,j) = FHCnd10(i,j).*((Tairset-Tair0(i,j))./(Tair10(i,j)-Tair0(i,j)));
                if FHCndun(i,j)<=FHmax
                    FHCndac(i,j) = FHCndun(i,j);
                    [Tm,Ts,Tair,Top] = Tair_calc(i,j,Hve,His,Hw,Hms,Hem,Fm,Fst,Fia,Ta,Ta,thmt1,FHCndac,Cm); %%%%
                    Tm2(i,j) = Tm;
                    Tairac(i,j) = Tairset;
                else
                    FHCndac(i,j) = FHmax;
                    [Tm,Ts,Tair,Top] = Tair_calc(i,j,Hve,His,Hw,Hms,Hem,Fm,Fst,Fia,Ta,Ta,thmt1,FHCndac,Cm);
                    Tm2(i,j) = Tm;
                    Tairac(i,j) = Tair;
                end
            end
        end
    end