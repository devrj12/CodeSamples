
% Algorithm to compare ground observations of radio wave scintillation (
% degradation ) with in-situ satellite observations of ionospheric
% irregularties    -- Dev Joshi

%% new
  dataDir2   = 'C:\Users\joshide\Desktop\CNOFSnew\';
  s          = what(dataDir2);
  S = s.mat;
  Num = [];
  for sj = 1:size(S,1)
  S1 = cell2mat(S(sj));
  Nu   = str2double(regexp(S1,'\d*','match'));
  Num  = [Num; Nu];
  end
  NUM = sort(Num);
  
  %%
  dirnew   = 'E:\CNOFSLocal\New\Localdays\';
  dirnew1  = 'E:\CNOFSLocal\New\LocalTimedays\';
  dirSig   = 'E:\CNOFSLocal\New\Sigmaa\';    %  (1) This is changed with r.s.t. Groundspace5ab 
  dirdel   = 'E:\CNOFSLocal\New\Denn\';
  
  load ASIchrisnew1.mat  % check two ground conditions
  load latlonipp1.mat % lat-lon vectors for latipp and lonipp created in Window.m ( latlonipp1.mat -- extends to -13 to + 13)

  NUMday   = []; 
  NUMdayy  = [];
  NUMday1  = [];
  Passnum  = [];
  Passnumm = [];
  
  NUMdaya  = []; %  days which have passes
  Passnuma = [];   %  Counts in each day with passes
  NUMdayb  = []; %  days which have passes
  Passnumb = [];   %  Counts in each day with passes
  NUMdayc  = []; %  days which have passes
  Passnumc = [];   %  Counts in each day with passes
  NUMdayd  = []; %  days which have passes
  Passnumd = [];   %  Counts in each day with passes
  
  
for nm  = 2:length(NUM)-1  % for a given year pick the indices of NUM

     st1       = strcat('Local_',num2str(NUM(nm))) % Days
     load([dirnew,   st1]);
 
     st        = num2str(NUM(nm));
      
     st3       = strcat('LT_',num2str(NUM(nm))); % LocalTime
     load([dirnew1, st3]);
     
     st5      = strcat('Denmn_',num2str(NUM(nm)));
     load([dirdel, st5]);
     

 % Assign quantities
Density      = lrz1(:,3);
Lat1         = oa1(:,11);
Lon1         = oa1(:,12);
% LocalTime1   = oa1(:,27);
Alt1         = oa1(:,13);
Apexalt1     = oa1(:,35);
Orbit1       = oa1(:,1);
Maglat1      = oa1(:,38);
UTC          = 1:length(Density);

LocalTime1   = LT;
minta        = @intersect;


 % find the window

latipp1(1)   = min(Lat1);   % relaxing latitude window
latipp1(end) = max(Lat1);


% changing lon-window according to the magnetic field line from the ipp : created at window.m
LATLON = [];
for lc = 2:length(latipp1)
   
Lafind4      = find(latipp1(lc - 1) < Lat1 & Lat1 <= latipp1(lc));  % use latipp1a to skip every 2nd latipp1 point -- fewer intervals


Lonfind4    = find(Lon1 >= (lonipp1(lc-1) - 0.75) & Lon1 <= (lonipp1(lc-1) + 0.75)) ;

Latlon      = minta(Lafind4,Lonfind4);
LATLON      = [LATLON Latlon' ];

end

Locfind3    = find(LocalTime1 >= 20  &  LocalTime1 <= 24);               % 3.   tbcl

minta1      = minta(LATLON,Locfind3);
UTC1        = minta1;
      

% find Sigmamed > 1.2 places in the window/ count the passes

st2 = strcat('LocSigMed_',num2str(NUM(nm)));  % this SigmaMed is expanded
load([dirSig, st2])
Sigma2 = Sigmamed(UTC1);
UTC2   = UTC1(Sigma2 > 1.2);


% pass counting

UTCdh  = find((diff(UTC1))> 1);
UTCdh1 = [0; UTCdh ; length(UTC1)];


A2 = [];
for pp  = 1:length(UTCdh)+1 % length(UTCdh) + 1 = total number of passes
    A   = length(UTCdh1(pp)+ 1: UTCdh1(pp+1));
    A1  = pp*ones([1, A]) ; 
    A2  = [A2 A1]; % pass number in UTC1
end

A2un1 = unique(A2);  % unique passes

A2un = unique(A2(Sigma2 > 1.2)); % unique passes with irregualarities
% TA2  = length(A2un); % pass number with Sigma  > 1 or passes with scintillation
TA2  = length(A2un1); % only the count of the passes irrespective of irregularities ( Sigma > 1);

% record days with no-pass at all
 if TA2 == 0    
     NUMday1 = [NUMday1 NUM(nm)]   
     

 end


% make ground-space comparison &&  calculate and save the quantities for each pass
TP = [];
for tp = 1:TA2 % TA2 = unique number of passes in a day 
    
  NUMday = [NUMday NUM(nm)]; %  days which have passes
  Passnum = [Passnum   TA2];   %  Counts in each day with passes

  
    tpv    = find(A2 == tp);
    UTC3   = UTC1(tpv(1):tpv(end)); % positions for a pass
    Sigma3 = Sigmamed(UTC3);
    
    % here delden and enhancement filters
    Dennloc3   =  Dennloc(UTC3);
    Demnloc3   =  Demnloc(UTC3);
    delDenloc3 =  delDenloc(UTC3);
    
    nanfin     = find(isnan(Sigma3) == 1);
    Dennloc3(nanfin)  =  [];
    Demnloc3(nanfin)  =  [];
    delDenloc3(nanfin)=  [];
    Sigma3(nanfin)    =  [];
    UTC3(nanfin)      =  [];
    
    
    delfin  = find(delDenloc3 < 0.02);
    
    Dennloc3(delfin) = [];
    Demnloc3(delfin) = [];
    Sigma3(delfin)   = [];
    UTC3(delfin)     = [];
    
    denfin  = find(Dennloc3 > Demnloc3);
    Sigma3(denfin)   = [];
    UTC3(denfin)     = [];
    
  if isempty(Sigma3) == 0
      
  NUMdayy  = [NUMdayy NUM(nm)]; %  days which have passes
  Passnumm = [Passnumm   TA2];   %  Counts in each day with passes
   

    UTC3a  = UTC3(Sigma3 == max(Sigma3));
    
    
    UTC4   = UTC3(Sigma3 > 1.2); % position of the pass being ( Sigma > 1)
    Sig5   = Sigmamed(UTC4);
    UTC5   = UTC4(Sig5 == max(Sig5));
 

%     Choose ground data
      GD       = find(ASIall(:,1) == NUM(nm)); % change for every year
      GD1      = ASIall(GD,:); % change for every year
      GD1(:,3) = GD1(:,3) - 0.96; % convert the UT to LT
      GDLT     = GD1(:,3); % Local Time column
      GDLT1    = find(GDLT >= 20  &  GDLT <= 24); % positions for the Local-time filter
      GD2      = GD1(GDLT1,:); % final ground data for the day for the filters
  

    if isempty(UTC4) == 0 %  if there is irregularity in space (pass is already assumed) : Check1
          
       
        TP   = [TP tp]; % check -- shall be same as A2un
        LTav = mean(LocalTime1(UTC5));
       
         % pick up date and local time filter the data

         Xtest = abs(LTav - GD2(:,3));
         Gs    = fix(mean(find(Xtest == min(Xtest))));  % LTav closest in local-time with the ground data
         
         Gsin1 = Gs - 1:Gs + 1;
%          Gsin1 = Gs - 2:Gs + 2;
         Gsin2 = 1:length(GD2(:,3)); % any column of GD2 - taken 3 here for instance
         Gsin  = intersect(Gsin1,Gsin2);
         
         if isempty(GD2) ~= 1 % if GD2 isn't empty
         GS = mean(GD2(Gsin,4)); % mean of the 15-min S4
         end
         

         if GS > 0.6
             
              NUMdaya = [NUMdaya NUM(nm)]; %  days which have passes
              Passnuma = [Passnuma   tp];   %  Counts in each day with passes
             
              PAltpos       = mean(Alt1(UTC5));  % the ground has scintillated
              PApexpos      = mean(Apexalt1(UTC5));
              PLatpos       = mean(Lat1(UTC5));
              PLonpos       = mean(Lon1(UTC5));
              PMaglatpos    = mean(Maglat1(UTC5));
              PLocTimepos   = mean(LocalTime1(UTC5));
              PSigpos       = mean(Sigmamed(UTC5));
              
              pp            = [PAltpos PApexpos PLatpos PLonpos PMaglatpos PLocTimepos PSigpos  GS];
              
              PP.(['D' num2str(tp) '_' num2str(NUM(nm))]) = pp;
              

         else
   
              NUMdayb       = [NUMdayb NUM(nm)]; %  days which have passes
              Passnumb      = [Passnumb   tp];   %  Counts in each day with passes
             
              NAltpos       = mean(Alt1(UTC5));    % the ground hasn't scintillated
              NApexpos      = mean(Apexalt1(UTC5));
              NLatpos       = mean(Lat1(UTC5));
              NLonpos       = mean(Lon1(UTC5));
              NMaglatpos    = mean(Maglat1(UTC5)); 
              NLocTimepos   = mean(LocalTime1(UTC5));
              NSigpos       = mean(Sigmamed(UTC5));
              
              np            = [NAltpos NApexpos NLatpos NLonpos NMaglatpos NLocTimepos NSigpos GS];
              
              NP.(['D' num2str(tp) '_' num2str(NUM(nm))]) = np;


         end
         
        
        
    elseif isempty(UTC4) == 1 % if no irregularity in space but passes
          

        
         LTav1 = mean(LocalTime1(UTC3a));
       
         % pick up date and local time filter the data

         Xtest1 = abs(LTav1 - GD2(:,3)); % 2nd column of GD2 is local-time at the IPP
         Gs1 = fix(mean(find(Xtest1 == min(Xtest1)))); 
         
         Gsin1a = Gs1 - 1:Gs1 + 1;
%          Gsin1a = Gs1 - 2:Gs1 + 2;
         Gsin2a = 1:length(GD2(:,3)); % any column of GD2 - taken 3 here for instance
         Gsina  = intersect(Gsin1a,Gsin2a);
         
         if isempty(GD2) ~= 1
         GS1 = mean(GD2(Gsina,4));
         end
         
         if GS1 > 0.6
              
              NUMdayc       = [NUMdayc NUM(nm)]; %  days which have passes
              Passnumc      = [Passnumc   tp];   %  Counts in each day with passes
              
             GS = GS1;
             
              PAltneg       = mean(Alt1(UTC3a));  % the ground has scintillated
              PApexneg      = mean(Apexalt1(UTC3a));
              PLatneg       = mean(Lat1(UTC3a));
              PLonneg       = mean(Lon1(UTC3a));
              PMaglatneg    = mean(Maglat1(UTC3a));
              PLocTimeneg   = mean(LocalTime1(UTC3a)); 
              PSigneg       = mean(Sigmamed(UTC3a));
              
              pn            = [PAltneg PApexneg PLatneg PLonneg PMaglatneg PLocTimeneg PSigneg GS];
              
              PN.(['D' num2str(tp) '_' num2str(NUM(nm))]) = pn;

         else
             
              NUMdayd       = [NUMdayd NUM(nm)]; %  days which have passes
              Passnumd      = [Passnumd   tp];   %  Counts in each day with passes
             
             GS  = GS1;
   
              NAltneg       = mean(Alt1(UTC3a));    % the ground hasn't scintillated
              NApexneg      = mean(Apexalt1(UTC3a));
              NLatneg       = mean(Lat1(UTC3a));
              NLonneg       = mean(Lon1(UTC3a));
              NMaglatneg    = mean(Maglat1(UTC3a)); 
              NLocTimeneg   = mean(LocalTime1(UTC3a));
              NSigneg       = mean(Sigmamed(UTC3a));
              
              nn            = [NAltneg NApexneg NLatneg NLonneg NMaglatneg NLocTimeneg NSigneg GS];
              
              NN.(['D' num2str(tp) '_' num2str(NUM(nm))]) = nn;

         end
 
    end
    
   % pause
  end
  
end

end

%% Unique number of days with number of passes in each day

Passnum2 = [];
NUMday2  = [];
for j = 10260:11110   % 10260:11110
    
   ak =  find(NUMday ==j);
   
   if isempty(ak)~= 1
   
   Passnum2 = [Passnum2 Passnum(ak(1))];
   NUMday2 = [NUMday2 NUMday(ak(1))];
   
   end
end


Passnum22 = [];
NUMday22 = [];
for jk = 10260:11110   % 
    
   akk =  find(NUMdayy ==jk);
   
   if isempty(akk)~= 1
   
   Passnum22 = [Passnum22 Passnumm(akk(1))];
   NUMday22 = [NUMday22 NUMdayy(akk(1))];
   
   end
end



  