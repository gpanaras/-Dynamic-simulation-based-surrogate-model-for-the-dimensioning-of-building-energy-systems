%Author: Georgios Stergianakos
%Date: 26 May 2021
%Description: This program calculates the hourly thermal loads, according
%to the ISO EN 13790 - simple hourly method. For a certain building, it
%generates random parameter values for Um, Cm and Cinf. All variables below 
%are coherent with the respective ones from ISO EN 13790 methodology 
%and can be looked into.
clear; clc; format long
close all
%% INPUT
%Opaque block side surface area 
Af = 119.4; A1 = 58.8425; A2 = 60.58; Apar = 346.785; Aw = 45.6;
%Windows/doors are (in contact with outdoor air)
Awpar = 4*0.5*0.8+2*0.3*0.8+2*0.25*0.8; Awpor = 45.6-Awpar;
%Kitchen and living room floor areaa percentage 
Asalot = (17.4*4+7*4)/(2*Af);
mys = 1; 
%Thermal capacity
Cm2 = [80000,110000,165000,260000,370000]*2*Af; 
Am2 = [2.5,2.5,2.5,3,3.5]*2*Af; 
%Air infiltration rate
Vinf = [100,400,700];
%% Calculation loop initiates
for xx = 1:40
    Uoppar = 0.3+(2.5-0.3)*rand(1);  Uw = 2+(6-2)*rand(1); Uf1 = Uoppar;
    Uf2 = 0.3+(2.5-0.3)*rand(1);
    Um(xx) = (Uoppar*(Apar-Aw) + Uw*Aw + Uf1*Af + Uf2*Af)/(Apar + 2*Af);
    for yy = 1:length(Am2)
        for ww = 1:length(Vinf) 
            Ta = xlsread('ATHIhour45.xls','I2:I745'); %Gets hourly outdoor 
            %temperature data from excel spreadsheet (use an .xls file of
            %your choice)
            Ta = reshape(Ta,[24,31]);
            Isol = xlsread('ATHIhour45.xls','F8018:F8761');% Gets hourly 
            %solar radiation data from excel spreadsheet (use an .xls file
            %of your choice)
            Isol = reshape(Isol,[24,31]);
            FHmax = 20000; FCmax = -5000; %max heating and cooling power
            % heat transmission coefficient - resistances of opaque blocks
            % for calculating building solar gains
            Fshob = 1; Fgl = 1; asc = 0.6; ggl = 0.68;
            Ror = 0.04; Radap = 0.17; Raop = 0.04;
            % Calculation of the inner temperature settings, according to 
            %the occupancy schedule
            Tcomf = 18:1:20; %Inner room setpoint for heating
            for xxx = 1:length(Tcomf)
                ss = 0;
                for j = 1:31
                    ss = ss+1;
                    for i = 1:24
                        if  ss ~= 6 && ss ~= 7 %Weekday occupancy schedule
                            if (i >= 8 && i < 16)
                                TintH(i,j) = -inf;
                                TintC(i,j) = inf;
                            else
                                TintH(i,j) = Tcomf(xxx);
                                TintC(i,j) = 26;
                            end
                        else       %Weekend occupancy schedule
                            if (i >= 8 && i < 16) 
                                TintH(i,j) = Tcomf(xxx);
                                TintC(i,j) = 26; 
                            else
                                TintH(i,j) = Tcomf(xxx);
                                TintC(i,j) = 26;
                            end
                        end
                    end
                    if ss == 7
                        ss=0;
                    end
                end
                %Call LOAD_SCEN_31D
                FHmax = 10000:500:40000;
                mysum = 1;       
                for zz = 1:length(FHmax)  %FHMAX
                    [Tairac_diak_HA, FHCndac_diak_HA, FHCndun_diak_HA] = LOADS_SCEN_31D(Af,Apar,Aw,Awpar,Awpor,Asalot,Uoppar,Uw,Uf1,Uf2,Ta,Isol,TintH,TintC,FHmax(zz),FCmax,Fshob,Fgl,asc,ggl,Uoppar,Uf1,Uf2,Ror,Radap,Raop,Cm2(yy),Am2(yy),Vinf(ww));
                    % ENERGY
                    En_diak_HA_monthly(xxx,zz,ww,yy,xx) = sum(abs(trapz(1:24,FHCndac_diak_HA)*10^(-3)));     
                    En_diak_HA_monthly_un(xxx,zz,ww,yy,xx) = sum(abs(trapz(1:24,FHCndun_diak_HA)*10^(-3)));  
                    [r,c] = size(Tairac_diak_HA);
                    sss = 0; ss = 0;
                    for jjj = 1:c   % 31 days
                        sss = sss + 1;
                        for iii = 1:r   % 24 hours
                            if sss ~= 6 && sss ~= 7 % Monday - Friday
                                if (iii <= 7 || iii >= 16)
                                    ss = ss + abs(Tcomf(xxx)-Tairac_diak_HA(iii,jjj));
                                end
                            else % Saturday - Sunday 
                                ss = ss + abs(Tcomf(xxx)-Tairac_diak_HA(iii,jjj));
                            end
                        end
                        if sss == 7  % when Sunday
                            sss = 0;
                        end
                    end
                    HDH(xxx,zz,ww,yy,xx) = ss/(2*(31*16+8*2*4));
                end
            end
        end
    end
end
%% Final generated data to be exported are organized into arrays
[r,c,d,e,f] = size(En_diak_HA_monthly);
for p = 1:f
    for t = 1:e
        for k = 1:d
            for i = 1:r
                for j = 1:c
                    if j<=c-1
                        if (abs(En_diak_HA_monthly(i,j,k,t,p)-En_diak_HA_monthly(i,j+1,k,t,p)))/(abs(FHmax(j+1)-FHmax(j))) <= 0.002
                            En_diak_HA_monthly_N(i,j,k,t,p) = -123;
                        else
                            En_diak_HA_monthly_N(i,j,k,t,p) = En_diak_HA_monthly(i,j,k,t,p);
                        end
                    elseif j==c
                        if (abs(En_diak_HA_monthly(i,j,k,t,p)-En_diak_HA_monthly(i,j-1,k,t,p)))/(abs(FHmax(j-1)-FHmax(j))) <= 0.002
                            En_diak_HA_monthly_N(i,j,k,t,p) = -123;
                        else
                            En_diak_HA_monthly_N(i,j,k,t,p) = En_diak_HA_monthly(i,j,k,t,p);
                        end
                    end
                end
            end
        end
    end
end
%---------------
[FHmax_m, Tcomf_m, Vinf_m, Cm_m, Um_m] = ndgrid(FHmax,Tcomf,Vinf,Cm2,Um);
FHmax_m = permute(FHmax_m,[2,1,3,4,5]);
Tcomf_m = permute(Tcomf_m,[2,1,3,4,5]);
Vinf_m = permute(Vinf_m,[2,1,3,4,5]);
Um_m = permute(Um_m,[2,1,3,4,5]);
Cm_m = permute(Cm_m,[2,1,3,4,5]);
xData = reshape(FHmax_m,[numel(FHmax_m),1]);
yData = reshape(Tcomf_m,[numel(Tcomf_m),1]);
zData = reshape(Um_m,[numel(Um_m),1]);
tData = reshape(Cm_m,[numel(Cm_m),1]);
pData = reshape(Vinf_m,[numel(Vinf_m),1]);
wData = reshape(En_diak_HA_monthly_N,[numel(En_diak_HA_monthly_N),1]);
wData_N = wData(find(wData~=-123));
zData_N = zData(find(wData~=-123));
xData_N = xData(find(wData~=-123));
yData_N = yData(find(wData~=-123));
tData_N = tData(find(wData~=-123));
pData_N = pData(find(wData~=-123));
xlswrite('C:\Users\leozo\Desktop\6D_EN_LCR_new_more_ENERGY.xlsx', [xData_N,yData_N,zData_N,tData_N,pData_N,wData_N])
%---------------
%Final generated data to be used
[r,c,d,e,f] = size(HDH); 
for p = 1:f
    for t = 1:e
        for k = 1:d
            for i = 1:r
                for j = 1:c
                    if HDH(i,j,k,t,p) <= 0.005 ||  HDH(i,j,k,t,p)>=0.02
                        HDH_N(i,j,k,t,p) = -123;
                    else
                        HDH_N(i,j,k,t,p) = HDH(i,j,k,t,p);
                    end
                end
            end
        end
    end
end
%-----------------
[FHmax_m2, Tcomf_m2, Vinf_m2, Cm_m2, Um_m2] = ndgrid(FHmax,Tcomf,Vinf,Cm2,Um);
FHmax_m2 = permute(FHmax_m2,[2,1,3,4,5]);
Tcomf_m2 = permute(Tcomf_m2,[2,1,3,4,5]);
Vinf_m2 = permute(Vinf_m2,[2,1,3,4,5]);
Um_m2 = permute(Um_m2,[2,1,3,4,5]);
Cm_m2 = permute(Cm_m2,[2,1,3,4,5]);
xData = reshape(FHmax_m2,[numel(FHmax_m2),1]);
yData = reshape(Tcomf_m2,[numel(Tcomf_m2),1]);
pData = reshape(Vinf_m2,[numel(Vinf_m2),1]);
zData = reshape(Um_m2,[numel(Um_m2),1]);
tData = reshape(Cm_m2,[numel(Cm_m2),1]);
wData = reshape(HDH_N,[numel(HDH_N),1]);
wData_N = wData(find(wData~=-123));
zData_N = zData(find(wData~=-123));
xData_N = xData(find(wData~=-123));
yData_N = yData(find(wData~=-123));
tData_N = tData(find(wData~=-123));
pData_N = pData(find(wData~=-123));
