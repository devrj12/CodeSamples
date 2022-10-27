
% Algorithm to Optimize the difference between anchor profile and IRI profile by changing IRI variables

%%

function [x, fval, history,exitflag,output] =  Naldermead(x0) 

 history = [];

 options = optimset('OutputFcn', @outfun,'Display','iter');
 [x,fval,exitflag,output] = fminsearch(@(x)ban2(x),x0,options);
 
 function stop = outfun(x,~,state)
        stop = false;

        if isequal(state,'iter')
          history = [history; x];
        end
 end

function ban1 = ban2(x)
 
    
      xtest = x;
      assignin('base','xtest',xtest);
 
      A =       []; % variables in avea1.mat .. need to be initialized while using nested function
      B =       [];
      D6xa =    []; % ALTAIR profile
      alt3 =    [];
      yfinal1 = [];
  
load avea1.mat
lat = 10.0042;
lon  = 167.6568;
R12 = 190;
UT = [2013 5 01 07 41];
ht_start = 0.1617;

ht_step = 4 ;
num_hts = 175;



alt = 0.1617:4:174*4+0.1617;
 
assignin('base','alt',alt);
assignin('base','lat',lat);
 
iono_layer_parms = [x(1), x(2), -1, -1, -1, -1]; 
[iono, ~] = iri2012(lat, lon, R12, UT, ht_start, ht_step,num_hts,iono_layer_parms);

assignin('base','iono_layer_parms',iono_layer_parms)
ionofreq = 0.008978663597663*sqrt((10.^(-6))*(iono(1,18:end))); % before 18th column, it's -1.

assignin('base','iono',iono);

% pause
D6xb = interp1(alt3(101:697),D6xa(101:697),alt(26:175));

assignin('base','ionofreq',ionofreq);
D6xc = 0.008978663597663*sqrt(D6xb);

C = 0.008978663597663*sqrt((10.^(-6)*(B(9:end))));
assignin('base','C',C);
assignin('base','A',A);
assignin('base','D6xc',D6xc);


ban = ionofreq(9:158) - D6xc;
% pause

assignin('base','ban',ban)
ban1 = sum((abs(ban)));

end
end
