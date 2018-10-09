
% Modification of 'PHaRLAP' ray-tracing MATLAB package to simulate the
% high-frequency radio wave propagation in the METAL OXIDE SPACE CLOUD
% (MOSC) experiment conducted by NASA and Air Force Research Laboratory in
% 2013 

%                                                      -- Dev Joshi

% PHaRLAP toolbox developed by Defence Science and Technology Organisation (DSTO),
% Australia. The toolbox is available by request from its author, Manuel Cervera.


% PHaRLAP code used in this algorithm : 
%   raytrace_3d
%
% Purpose :
%   3D magneto-ionic numerical raytrace for a multihop ray. The Hamiltonian
%   for the ray is that set down by Haselgrove and Hasegrove (1960), and
%   Hasegrove (1963). Geomagnetic field effects are considered for O and X 
%   polarized radio waves. No-field case is also considered if required. WGS84
%   coordinate system is assumed.




%%

CKa = clock; 


tic
load x5.mat %scaling factor
y = [1.75,1.5344,1.3959,1.2812,1.1881,1.1146,1.1146,1.1146,1.1224,1.1146,1.1146,1.1146];

load Delaylaunch1RW0741.mat
feqlist = freqRW1;
flist = [];
timed = [];
BK = [];
FH = [];
TM = [];
R1 = [];
PL = [];
la = [];
lo = [];
ht = [];
nu = [];
Dac = [];

m = -2:0.05:2;
elevs= 45:0.05:52; 



for ki = 1:length(feqlist)
    
    load('RWcLAUNCH1.mat');
    freq = feqlist(ki);
    
    
clearvars K K5;
for j = 1:115
ionopf1 = squeeze(iono_pf_grid(j,:,:));
Ca = bsxfun(@times,ionopf1,(y(ki)*x5b));
K(j,:,:) = Ca;
end

for j = 1:115
ionopf2 = squeeze(iono_pf_grid_5(j,:,:));
Cb = bsxfun(@times,ionopf2,(y(ki)*x5b));
K5(j,:,:) = Cb;
end

iono_pf_grid = K;
iono_pf_grid_5 = K5;
 


% convert plasma frequency grid to  electron density in electrons/cm^3
iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;


load Seconds5.txt –ascii; % This is ALTAIR profile
   
P = Seconds5(:,1); % height
Q =(10).^(Seconds5(:,2));  % number density : change logarithm into numbers /(cm)^3
 
C = [P,Q];
h = C(:,1);
n = C(:,2); 
[p, q] = min(abs(h-156.8126));                         
[e, g] = min(abs(h-184.8497));                        

B = C(q:g,:);                % MOSC cloud
size(B);

% A = iono_en_grid(276:294,323:341,102:120); % 1.5510 km :  LWnewcloud.mat
A = iono_en_grid(21:39,146:164,102:120); % 1.5510 km :  1st release 
S =  156.8127:1.5510:184.7307;
K = A(1,1,:);
L = squeeze(K);
A1 = [S',L]; % testnew1

% Interpolate 
x   = B(:,1);   %height in radar data (corresponding to MOSC)
v   = B(:,2);   %electron density in radar data (corresponding to MOSC): these are 136 in number.
xq  = A1(:,1); %height in iono_en_grid
vq1 = interp1(x,v,xq);  % electron density to be added in the heights of  iono_en_grid according to electron densities in radar data

vq2 = flip(vq1);
vq  = (vq1+vq2)/2;
vq  = vq(1:10);
vq  = flip(vq);  % flipped as I go farther : electron density decreases. It's highest at the center.


%Add Cloud
G = zeros(size(iono_en_grid));

for i = -9:1:9
   for j = -9:1:9
   for k = -9:1:9

          R = (13.9590/9)*((i)^(2) + (j)^(2)+(k)^2)^(1/2); 
    if R <= 13.9590
%       if R <= 10

  G(31+i,155+j,111+k) = vq(round(9*R/13.9590)+1);  %1st Release


    end
    end
    end
end

Gpf = ((G)* 80.6164e-6);
Gpfgrid = (Gpf.^(1/2));


% pause
iono_en_grid = iono_en_grid + G ;% adds the cloud
iono_en_grid_5 = iono_en_grid_5 + G;

% Update iono_pf_grid
iono_pf = ((iono_en_grid)* 80.6164e-6);
iono_pf_grid = (iono_pf.^(1/2));  
iono_pfa = ((iono_en_grid_5)* 80.6164e-6);
iono_pf_grid_5 = (iono_pfa.^(1/2));

nhops = 1;           % number of hops
num_elevs = length(elevs);

ray_O = [];


for mi = 1:length(m)

ray_bear =  144.9246 + m(mi);

% initial bearing of radar ray : rongelap-MOSC

figure(3)
clf
plot3(10.0042,167.6568,170,'r*') %Mid
hold on;grid on
plot3(10.1740,166.0046,0,'b*') %wotho
plot3(11.1523,166.8378,0,'b*') %rongelap


[X,Y,Z] = sphere(18);
surf(X*0.1269 + 10.0042,Y*0.1269 + 167.6580, Z*13.9590 + 170.7);
shading interp ;

xlabel('latitude (deg)')
ylabel('longitude (deg)')
zlabel('Height (km)')


tol = [1e-7 0.01 25];

for x = 1:1
  OX_mode = 1;
  [ray_data, ray_path_data, nhops_done, ray_label, ray_state_vec] = ...
      raytrace_3d(origin_lat, origin_long, origin_ht, elevs(1), ray_bear, ...
                  freq, OX_mode, nhops, tol, iono_en_grid, iono_en_grid_5, ...
	          collision_freq, iono_grid_parms, Bx, By, Bz, ...
	          geomag_grid_parms);
          
      rayl = ray_label;
      
      B1 = [ray_path_data(1,end) ray_path_data(2,end) ray_path_data(3,end)]; 

end


% Generate the O mode rays
tol = [1e-7 0.01 25];
OX_mode = 1;
rl = [];
BO = [];

 for elev_idx= 1:num_elevs
  elev = elevs(elev_idx);
  [ray_data, ray_path_data, nhops_done, ray_label, ray_state_vec] = ...
      raytrace_3d(origin_lat, origin_long, origin_ht, elev, ray_bear, freq, ...
                  OX_mode, nhops, tol);
  la1 = [];
  lo1 = [];
  ht1 = [];
 	 
  initial_elev = elev;	   
  final_elev = ray_data(7);
  
  ray_O(elev_idx).lat = ray_path_data(1, :);
  ray_O(elev_idx).lon = ray_path_data(2, :);
  num = length(ray_path_data(1, :));
  
  nu = [nu num];

  ray_O(elev_idx).gndrng = zeros(1, num);
  ray_O(elev_idx).gndrng(1) = 0;
  for ii = 2:num
    lat = ray_path_data(1, ii); % defining latitude along the path for all points for given elevation
    lon = ray_path_data(2, ii);
    ground_range = latlon2raz(lat, lon, origin_lat, origin_long, ...
	'wgs84') / 1000.0;
    ray_O(elev_idx).gndrng(ii) = ground_range;
  end
  ray_O(elev_idx).height = ray_path_data(3, :);
  ray_path_data(3, :);
  Aab = size(ray_path_data(3, :));

  Bab = size(ray_path_data(1, :));
  Cab = size(ray_path_data(2,:));
  Dab = ray_path_data(3,end);
  Dac = [Dac ; Dab];



figure(3)
plot3(ray_O(elev_idx).lat, ray_O(elev_idx).lon, ray_O(elev_idx).height, 'g', 'markersize', 5)


  la  = [ la ray_path_data(1, :)];
  lo  = [ lo ray_path_data(2, :)];
  ht =  [ ht  ray_path_data(3, :)];
  
  la1 = [ la1 ray_path_data(1, :)];
  lo1 = [ lo1 ray_path_data(2, :)];
  ht1 = [ ht1  ray_path_data(3, :)];
  
  Powerloss = ray_data(11,:);
  PL = [Powerloss PL];
  rayl = ray_label;
  rl = [rl rayl];
  R1 = [R1;rayl];
  C1 = [ray_path_data(1,end) ray_path_data(2,end) ray_path_data(3,end)]; 
  



  BO = [BO;C1] ;% coordinate of last points for this ray bearing( all angles)
  BK = [BK;C1] ; % all last points
 end
 

FI = find(rl==1);% find if it is within the region
BI = [];
for uv = 1:length(FI)
BJ = [BO(FI(uv),:)];
BI = [BI;BJ];  % last points within the region
end

BI;

if ~isempty(BI)
 D = [10.1740 166.0046 0]; % Wotho site
 F = [];
 
 for k = 1:size(BI,1) %(BI,1) - rows ; k = 1:11
     
F1 = abs(BI(k,:) - D);
F = [F; F1];

 end
 
 F2 = F*110;    % difference in km from the target site to Wotho - height multiplied with 110 presuming it would hit 0.
 FH = [FH;F2];  % there are all F2s like all BOs in BK.
 F3 = F2(:,1);
 F4 = F2(:,2);
 F5 = F2(:,3);

TN =  find(F3 < 10 & F4 < 10 & F5 == 0); % Target number : finds position of angles at which target is hit

le = length(TN);

if ~isempty(TN)
    
    clear ray
    clear ray_O 
    

     for in = 1:length(TN)
          
          [ray_data, ray_path_data, nhops_done, ray_label, ray_state_vec] = ...
      raytrace_3d(origin_lat, origin_long, origin_ht, elevs(FI(TN(in))), ray_bear, freq, ...
                  OX_mode, nhops, tol);  % elevs within the target site within the region  = FI(TN(in)
         
     ray(in).gndrng = ray_path_data(1, :);
     ray(in).height = ray_path_data(2, :);
  

     grouprange(in) = ray_data(4,:);
     geometricaldistance = ray_data(15,:); 
     timedelay(in)  =  grouprange(in)/(3*10^2);
 
     end


 timeadd = timedelay;
 fadd = freq*ones(1,length(TN));
 flist = [flist fadd];
 timed = [timed timeadd];
 timedelay = [];
 

end

end
clear ray
end

end


figure
plot(flist,timed,'r*')
grid on

toc
CKb = clock;
CKc = CKb - CKa;





