clc
TypeLidf=2;
LIDFa	=	70;
LIDFb	=	0;
% Cab		= 40;		% chlorophyll content
Car		= 8;		% carotenoid content
Ant		= 0;		% carotenoid content 
Cbrown	= 0.0;	    % brown pigment content (arbitrary units)
% Cw		= 0.01;   	% EWT (cm)
% Cm		= 0.009;	% LMA (g.cm-2)
% N		= 2.275;	    % structure coefficient
data    = dataSpec_PDB_GSV;

Rsoil1  = data(:,11);
Rsoil2=data(:,12);
psoil	= 0; % soil factor (psoil=0: wet soil / psoil=1: dry soil)
% rsoil0  = psoil*Rsoil1+(1-psoil)*Rsoil2;
GSV1  = data(:,13);
GSV2  = data(:,14);
GSV3  = data(:,15);
GSV_SM  = data(:,16);
i=0
gsvsoil=zeros(2101,141);
for c1=0.25:0.2:0.75;
    for c2=-0.2:0.1:0.1;
        for c3=0:0.05:0.1;
            for csm=-0.4:0.1:0;
                for rho = c1*GSV1 + c2*GSV2 + c3*GSV3+ csm*GSV_SM;
                    if size(find(rho<0),1)==0
                        i=i+1;
                        gsvsoil(:,i)=rho;
                        %                     plot(band,rho)
                    end
                    hold on
                end
            end
        end
    end
end
meansoil=mean(gsvsoil,2);

hspot	=	0.1;       % hot spot
N=1.5;
tto		=	55.;		% observer zenith angle (?
sateazimuth=170;       %observer azimuth angle (?
%
SRF=textread('.\SRF\Himawari400_700.txt'); %Spectral response function of AHI
SRF=SRF(1:475,2);
Es  = data(:,9);  % coefficient props(um-1)
Ed  = data(:,10);

temp2=0;
%%
aLAIlist=0.1:0.2:7;
aCablist=20:10:60;
aCwlist=[0.004,0.0075,0.04];
aCmlist=[0.0019,0.0058,0.0169];
apsoillist=0.1:0.2:1.6;
attslist=25:5:80;%SOZ
apsilist=0:10:180;%RAA 
askyllist=0.01:0.1:0.7;
aLIDFalist=[50,60,70,80];
%%

for tts=attslist
    tts
    if tts<10
        fid=fopen(['HongheLUTSOZ0',num2str(tts*10),'.txt'],'w');%дȫ΄¼þ·¾¶
    else
        fid=fopen(['HongheLUTSOZ',num2str(tts*10),'.txt'],'w');%дȫ΄¼þ·¾¶
    end
	for psi=apsilist

        for LAI=aLAIlist
          
            for Cab=aCablist
                for Cw=aCwlist
                    for Cm=aCmlist
                        for psoil=apsoillist
                            rsoil  = psoil*meansoil;
                            for LIDFa=aLIDFalist
                                [rdot,rsot,rddt,rsdt,tdd,tss,tsd]=PRO4SAIL_FPAR(N,Cab,Car,Ant,Cbrown,Cw,Cm,LIDFa,LIDFb,TypeLidf,LAI,hspot,tts,tto,psi,rsoil);
                                rd  = pi/180;
                                tempi=0;
                                for skyl=askyllist
                                    % temp
                                    % rdot: hemispherical-directional reflectance factor in viewing direction    
                                    % rsot: bi-directional reflectance factor
                                    % rsdt: directional-hemispherical reflectance factor for solar incident flux
                                    % rddt: bi-hemispherical reflectance factor
                                    PARdiro	= (1-skyl)*Es;
                                    PARdifo	= (skyl)*Ed;

                                    hdresv	= (rdot.*PARdifo+rsot.*PARdiro)./(PARdiro+PARdifo);  
                                    %%
                                    %band3
                                    hdband3=sum((hdresv(202:279).*SRF(202:279)))/sum(SRF(202:279));
                                    %%
                                    %band4
                                    hdband4=sum((hdresv(441:474).*SRF(441:474)))/sum(SRF(441:474));                            
                                    %%
                                    abs_dif=1-rddt-(1-rsoil+tdd.*rsoil).*tdd./(1-rddt.*rsoil); %Ad
                                    abs_dir=1-rdot-(1-rsoil+tdd.*rsoil).*(tss+tsd)./(1-rddt.*rsoil); %As
                                    fPARdif=abs_dif(1:301,1);
                                    fPARdir=abs_dir(1:301,1);
                                    fPARdif1=mean(abs_dif(1:301,1));
                                    fPARdir1=mean(abs_dir(1:301,1));
                                    fPAR1 = fPARdir1*(1-skyl) + fPARdif1*skyl;

                                    %%
                                        if fPARdif1>0 && fPARdir1>0 && fPAR1>0
                                            fprintf(fid,'%d %d %f %f %f %f %f %f %f %f %f %f %f\r',...
                                                [tts,psi,...
                                                LAI,Cab,Cw,Cm,psoil,skyl,...
                                                fPARdif1,fPARdir1,fPAR1,...
                                                hdband3,hdband4]);
                                        end
                                end
                                %                                 fclose(fid);
                                % dlmwrite('LUT.txt',[sunazimuth,psi,band1,band2,band3,LAI,fPARdif,fPARdir,fPAR],'delimiter','\t','precision',5)
                            end
                        end
                    end
                end
            end
        end
    end
    fclose(fid);
end
