
% Algorithm to make climatology plots of satellite observations of
% disturbances in ionosphere  -- Dev Joshi 


%% Make the plots with new Data from CNOFS7pre2 ( folder name - saved by algorithmn CNOFS7pre2(a)) 

Cka         = clock;
dataDir     = 'C:\Users\joshide\Desktop\CNOFSDATAn1\Final\CNOFSpre2\N2\';
dataDir12   = 'C:\Users\joshide\Desktop\CNOFSDATAn1\Final\CNOFSpre2\Leap2012\2012a\';

folders   = dir(dataDir);
isub      = [folders(:).isdir]; %# returns logical vector
nameFolds = {folders(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = []; %

OM     = zeros(12,18); % Orbit Matrix
EM     = zeros(12,18); % Event Matrix
DOM    = zeros(12,18); % Orbit Matrix
DEM    = zeros(12,18); % Event Matrix
DEM1   = zeros(12,18); % Event Matrix

Sigmall1 = [];
Sigmall  = [];
Apexall  = [];
Lonall   = [];
Dayall   = [];
Monall   = [];
Yearall  = [];
Magall   = [];
Locall   = [];
Latall   = [];
Dennall  = [];
Denmall  = [];
UTCall   = [];
Delall   = [];


for k       = 4:7 % only for test : years 2011 - 2014
 
    

   nF            = nameFolds{k};
   
   if k ~= 5

   dataDir1      = strcat(dataDir,nF,'\'); % subfolder name : 2008, 2009, etc.
   
   elseif k == 5

   dataDir1      = dataDir12; % subfolder name : 2008, 2009, etc.
      
   end
   
 
  
  Year           = str2double(nF(1:4));
  
  load([dataDir1, 'Long1.mat']);
  
  Dn1a(k) = length(Densityn2long);

   
  
 % pause
  % 11.5 year
   if k == 4
        

     
     MT2                  = find(MON2medianlong > 6); %LT  = find(Test >= 23 | Test <=1)
    
        
   Sigmalong              = Sigmalong(MT2);
   DELNlong               = DELNlong(MT2); 
   DELN1long              = DELN1long(MT2);
   
   Apexalt2medianlong     = Apexalt2medianlong(MT2);
   Alt2medianlong         = Alt2medianlong(MT2);
   Magneticlat2medianlong = Magneticlat2medianlong(MT2);
   LocalTime2medianlong   = LocalTime2medianlong(MT2);
   Lon2medianlong         = Lon2medianlong(MT2); 
   Orbit2medianlong       = Orbit2medianlong(MT2);
   
   DAY2medianlong         = DAY2medianlong(MT2);
   MON2medianlong         = MON2medianlong(MT2);
   Lat2medianlong         = Lat2medianlong(MT2);
        
    end
     
   
   %% Filter - Local Time

   LT2                    = find(LocalTime2medianlong >= 20  &  LocalTime2medianlong <= 24); %LT  = find(Test >= 23 | Test <=1)
   

   Sigmalong              = Sigmalong(LT2);
   Sigmamedlong           = Sigmamedlong1(LT2); % Sigmamedlong1(LT2);
   DELNlong               = DELNlong(LT2); 
   DELN1long              = DELN1long(LT2);
   
   Apexalt2medianlong     = Apexalt2medianlong(LT2);
   Alt2medianlong         = Alt2medianlong(LT2);
   Magneticlat2medianlong = Magneticlat2medianlong(LT2);
   LocalTime2medianlong   = LocalTime2medianlong(LT2);
   Lon2medianlong         = Lon2medianlong(LT2); 
   Orbit2medianlong       = Orbit2medianlong(LT2);
   
   DAY2medianlong         = DAY2medianlong(LT2);
   MON2medianlong         = MON2medianlong(LT2);
   Lat2medianlong         = Lat2medianlong(LT2);
   UTC2medianlong         = UTCsecs2medianlong(LT2);
   
   
   
   Densityn2long          = Densityn2long1(LT2);
   Densitymn2long         = Densitymn2long1(LT2);
   delDen2meanlong        = delDen2meanlong(LT2);
   
   
   Dn(k) = length(Densityn2long);
   On(k) = length(Orbit2medianlong);
   
   %%
%    %
   ML2                    = find(Magneticlat2medianlong >= -10  &  Magneticlat2medianlong <= +10);
   
   Sigmalong              = Sigmalong(ML2);
   Sigmamedlong           = Sigmamedlong(ML2);
   DELNlong               = DELNlong(ML2); 
   DELN1long              = DELN1long(ML2);
   
   Apexalt2medianlong     = Apexalt2medianlong(ML2);
   Alt2medianlong         = Alt2medianlong(ML2);
   Magneticlat2medianlong = Magneticlat2medianlong(ML2);
   LocalTime2medianlong   = LocalTime2medianlong(ML2);
   Lon2medianlong         = Lon2medianlong(ML2); 
   Orbit2medianlong       = Orbit2medianlong(ML2);
   
   DAY2medianlong         = DAY2medianlong(ML2);
   MON2medianlong         = MON2medianlong(ML2);
   Lat2medianlong         = Lat2medianlong(ML2);
   UTC2medianlong         = UTC2medianlong(ML2);
   
   
   Densityn2long          = Densityn2long(ML2);
   Densitymn2long         = Densitymn2long(ML2);
   delDen2meanlong        = delDen2meanlong(ML2);

      Dn1(k) = length(Densityn2long);
   %%
%    
   AL2                    = find(Apexalt2medianlong >=  400  &  Apexalt2medianlong <= 500);
   % AL2                    = find(Apexalt2medianlong >=  400);
    
   Sigmalong              = Sigmalong(AL2);
   Sigmamedlong           = Sigmamedlong(AL2);
   DELNlong               = DELNlong(AL2); 
   DELN1long              = DELN1long(AL2);
   
   Apexalt2medianlong     = Apexalt2medianlong(AL2);
   Alt2medianlong         = Alt2medianlong(AL2);
   Magneticlat2medianlong = Magneticlat2medianlong(AL2);
   LocalTime2medianlong   = LocalTime2medianlong(AL2);
   Lon2medianlong         = Lon2medianlong(AL2); 
   Orbit2medianlong       = Orbit2medianlong(AL2);
   
   DAY2medianlong         = DAY2medianlong(AL2);
   MON2medianlong         = MON2medianlong(AL2);
   Lat2medianlong         = Lat2medianlong(AL2);
   UTC2medianlong         = UTC2medianlong(AL2);

   Year2medianlong        = Year*ones(length(DAY2medianlong),1);
   
   Densityn2long          = Densityn2long(AL2);
   Densitymn2long         = Densitymn2long(AL2);
   delDen2meanlong        = delDen2meanlong(AL2);

      Dn2(k) = length(Densityn2long);
    

  
%%


   Dnl1(k) = length(Densityn2long);
   Onl1(k) = length(Orbit2medianlong);
 
   Lon2med                =  Lon2medianlong/20;
   Lon2meds               =  ceil(Lon2med); % one longitude number per bin
   
   Apexall                = [Apexall Apexalt2medianlong];
   Sigmall1               = [Sigmall1 Sigmalong];
   Sigmall                = [Sigmall Sigmamedlong];

   Lonall                 = [Lonall  Lon2medianlong];
   Dayall                 = [Dayall  DAY2medianlong];
   Monall                 = [Monall  MON2medianlong];
   Yearall                = [Yearall Year2medianlong'];
   Magall                 = [Magall  Magneticlat2medianlong];
   Locall                 = [Locall  LocalTime2medianlong];
   Latall                 = [Latall  Lat2medianlong];
   UTCall                 = [UTCall  UTC2medianlong];
   
   Dennall                = [Dennall   Densityn2long];
   Denmall                = [Denmall   Densitymn2long];
   Delall                 = [Delall   delDen2meanlong];

%% Sigma -- excluding enhancements and smooth density changes

 for sel      = 1:18
  
      LonOrb      = find(Lon2meds == sel); % Position of the bins being used
         
   if isempty(LonOrb) == 0  % Only if the Longitude bin 'sel' isn't empty
       
       OrbLon                  = Orbit2medianlong(LonOrb); % Orbit selected for that longitude bin

       SigmaLon                = Sigmamedlong(LonOrb);        % Sigma selected for that longitude bin
       MonLon                  = MON2medianlong(LonOrb);   % Month selected for that longitude bin
       
       Denn2Lon                = Densityn2long(LonOrb);
       Denmn2Lon               = Densitymn2long(LonOrb);
       delDen2Lon              = delDen2meanlong(LonOrb);
   
       [MonLon1,ML,ML1]        = unique(MonLon);
      
        
      for mon1                 = 1:length(MonLon1) % no. of unique months for a longitude bin -- ideally should be all months 
      
          mon                  = MonLon1(mon1); % This is what goes in the matrices OM and EM. It is critically important when mon ? mon1 - which is the case when 
                                                % data for selected longitude sector(bin) aren't available for all months.   
                                                % which month
          olm                  = find(ML1 == mon1); % check - MonLon1(mon1) = unique(MonLon(olm)) ; 
          OrbLonMon            = OrbLon(olm);       % Orbit selected for that longitude-month bin
          SigmaLonMon          = SigmaLon(olm);     % Sigma selected for that longitude-month bin
          
          Denn2LonMon          = Denn2Lon(olm);
          Denmn2LonMon         = Denmn2Lon(olm);
          delDen2LonMon        = delDen2Lon(olm);
           
          [OrbNum1,ON,ON1]     = unique(OrbLonMon);
           OM(mon,sel)         = OM(mon,sel) + length(OrbNum1);  %% how many orbits
          
          SigmaCount  = 0;
          Sigmaloma1  = [];
          for jel              = 1:length(OrbNum1)
              OrbNum           = OrbNum1(jel);
              
              lom              = find(ON1 == jel);
              Sigmalom         = SigmaLonMon(lom);  %% value of Sigmas in that orbit
              
              Denn2lom          = Denn2LonMon(lom);
              Denmn2lom         = Denmn2LonMon(lom);
              delDen2lom        = delDen2LonMon(lom);
             
              Den3              = Denn2lom - Denmn2lom;
              
              %
              Sigmalomm         = Sigmalom(Den3 <= 0);  % excluding ( enhancements  Den3 > 0 )
              delDen2lomm       = delDen2lom(Den3 <= 0);
              
              Sigmalomn         = Sigmalomm(delDen2lomm >= 0.02); % excluding smooth density changes
               %
               
              Sigmaloma        = sum(double((Sigmalomn > 1.2))); 
              Sigmaloma1       = [Sigmaloma1 Sigmaloma];
              
                     if    Sigmaloma >= 1
                                SigmaCount = SigmaCount + 1;
                     end
          end

         EM(mon,sel)          = EM(mon,sel) + SigmaCount;  %% how many orbits have at least one instance of the criteria being fulfilled ?

       end 
       
    
     
   end
   
 end

 end


%Sigma
DM   = (EM./OM)*100;
DM2  = [DM(:,10:18) DM(:,1:9)]; % this is p-colored
EM2  = [EM(:,10:18) EM(:,1:9)]; % events(occurrence events)
OM2  = [OM(:,10:18) OM(:,1:9)]; % orbits

Ckb = clock;
Ckc = Ckb - Cka


%% Conversion of RGB

cm1       = cm;
cm1(19,:) = cm1(20,:);
cm1(17,:) = cm1(18,:);
cm1(15,:) = cm1(16,:);
cm1(1,:)  = [0.65 0.65 0.65];

figure;p1 = sanePColor(DM2);cm2 = colormap(cm1); caxis([0 100]);cb1 = colorbar ; shading interp

ax            = gca;
ax.YTick      = 1:12;
ax.XTick      = 2:2:18;
ax.XTickLabel = -150:40:180; % figure this out  150 or - 140
ax.YTickLabel = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
ax.TickDir = 'out';
ylabel('Month','FontSize',12);xlabel('Longitude','FontSize',12);
title({'Sigma > 1.2 %, 20 - 24 LT','400 - 500 km, 2011 - 2014'}, 'FontSize',12) % change this for every command


