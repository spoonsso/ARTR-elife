%% Run simulations and plot number of resources gathered over swim bout number
clear all;close all;
%load('13Nov14_posHMM_50_expi_137_newHMMprobs.mat'); %load in simulated trajectories
%load('22Oct14_possave_normalized_diffusion_bias_fixed_expi17.mat');
%load 7Jan14_400sims_50_75_expi17_difflocked
%load 8Nov2014_randomfish_uncorrected
%load 7Nov14_posHMM_50_swimdist_16026_diff_locked
%load 13Nov14_posHMM_74_50_swimdist_174_newHMMprobs
%load 14Jan14_95_HMMdiffusion;

%Use nump = 1520 for angle normalized diffusion (the 1e-5 var/mean local
%threshold for the mean resources after 40 swims for HMM75 fish

%Use nump = 1304 for dist normalized diffusion (the 1e-5 var/mean local
%threshold for the mean resources after 40 swims for HMM50 / distance
%corrected fish

load 15Jan14_posHMM_50_and_real_basedon3_23_1_newangledistribution_expi_118
%load 15Jan14_posHMM_50_swimdist_15081_good_angle_hist

pos_HMM_50 = single(pos_HMM_50); %reduce vector file size
pos_HMM_75 = single(pos_HMM_real);

rng(2);
%Extend traces via replication
% pos_HMM_50(41:80,:,:) = bsxfun(@plus, pos_HMM_50,pos_HMM_50(end,:,:));
% pos_HMM_75(41:80,:,:) = bsxfun(@plus, pos_HMM_75,pos_HMM_75(end,:,:));

%We want the boundaries of the resource field to remain constant even as we
%change the density and size of resources. In the original simulations and
%plots, resources extended from ~-50.9 to 50.9, i.e. resR 7.27/2 and xRes =
%-resR*7*2:2*resR:resR*7*2

%initialize structure that will hold the diameter, density, and processed
%resource fficiencies over bout number. The first two fields were used for
%a screen where resD and resDensity were varied iteratively.
eff_struct.resD = 0;
eff_struct.resDensity = 0;
eff_struct.eff_50 = 0;
eff_struct.eff_75 = 0;
eff_struct.totsum = 0;
eff_cnt = 1;

% try a random distribution but force a certain density, no overlap
resR = 13/2; %For all initial simulations in jan. 2014, resR = 13/2;
%Use 16 resources

%For random angle starts: -------
randangs = rand(1,size(pos_HMM_50,3))*2*pi;
randangs = repmat(randangs,[size(pos_HMM_50,1) 1]);
pos_(:,1,:) = squeeze(pos_HMM_50(:,1,:)).*cos(randangs) - squeeze(pos_HMM_50(:,2,:)).*sin(randangs);
pos_(:,2,:) = squeeze(pos_HMM_50(:,1,:)).*sin(randangs) + squeeze(pos_HMM_50(:,2,:)).*cos(randangs);
pos_HMM_50 = single(pos_);

pos_(:,1,:) = squeeze(pos_HMM_75(:,1,:)).*cos(randangs) - squeeze(pos_HMM_75(:,2,:)).*sin(randangs);
pos_(:,2,:) = squeeze(pos_HMM_75(:,1,:)).*sin(randangs) + squeeze(pos_HMM_75(:,2,:)).*cos(randangs);
pos_HMM_75 = single(pos_);
%-------------

%Perform 100 resource distributions/gatherings for all 1e6 trajectories
for j = 1%:100%:100
overlaps = 1;
cnttest = [];
    xRes = 50.9*2*rand(64,1)-50.9; %
    yRes = 50.9*2*rand(64,1)-50.9;
while overlaps

    

    % For ensuring that no 2 resource areas overlap (i.e. each center is at
    % least 2*resR away from another center)
    %if length(find(sqrt(bsxfun(@minus,xRes,xRes').^2+bsxfun(@minus,yRes,yRes').^2)<=resR*2)) == 16 & ...
    
    %Prevent resources from being placed in center
    if isempty(find(sqrt(xRes.^2 + yRes.^2) <= resR*2))
        overlaps = 0;
    else
        badz = sqrt(xRes.^2 + yRes.^2) <= resR*2;
        xRes(badz) = 50.9*2*rand(length(find(badz)),1)-50.9;
        yRes(badz) = 50.9*2*rand(length(find(badz)),1)-50.9;
    end

    %figure(9);plot(xRes,yRes,'o');
%     cnttest = [cnttest length(find(sqrt(bsxfun(@minus,xRes,xRes').^2+bsxfun(@minus,xRes,xRes').^2)<resR*2))];
%     figure(9);plot(cnttest);drawnow;
end

load 9Jan16_resourceloc;


for resD = 13 %was 3:2:20 for initial screen
% resR = 7.27/2;
for resDensity = 2 %was 2:5 for initial screen
    resDensity
resR = resD/2

%More symmetric distribution ---------
% xRes =resDensity*resR/2:resDensity*resR:50.9;%*2;%*3;
% xRes = [-1*xRes(end:-1:1) xRes];
% yRes =resDensity*resR/2:resDensity*resR:50.9;%*2;%*3;
% yRes = [yRes(end:-1:1) -1*yRes(1:end)];
%More symmetric distribution ----------

% Original grid distribution -----------
% xRes =-50.9:resDensity*resR:50.9;
% yRes =50.9:-resDensity*resR:-50.9;
%----------

% Code for iterative efficiency screen -------
% xRes =-resR*7*2:2*resR:resR*7*2;%0:2*resR:resR*7*2;%
% yRes =resR*7*2:-2*resR:-resR*7*2;
%----------

%Needed for grid distributions -----------
% xRes = repmat(xRes,[length(yRes) 1]);
% 
% yRes = repmat(yRes',[1 size(xRes,2)]);
% yRes = reshape(yRes,[1 size(yRes,1)*size(yRes,2)]);
% 
% xRes = reshape(xRes,[ size(xRes,1)*size(xRes,2),1]);

%-----------------


pos = pos_HMM_50; %use the pos variable flexibly

% pos = pos(:,:,1:1e5);

% efficiency = zeros([size(pos,1),size(pos,3),length(xRes)],'uint8');
efficiency = false([size(pos,1),size(pos,3),length(xRes)]); %initalize efficiency matrix

        
        for i=1:length(xRes) %count how many times simulated trajectory passed within a resource radius
            foundRes = bsxfun(@minus,squeeze(pos(:,1,:)),xRes(i)).^2 + bsxfun(@minus,squeeze(pos(:,2,:)),yRes(i)).^2  <= resR^2;
            
            [a b] = max(foundRes,[],1);
     
           onez = sub2ind(size(foundRes),b(a~=0),find(a));
           foundRes = false(size(foundRes));
           foundRes(onez)  = 1;
            efficiency(:,:,i) = foundRes;
            
        end
        
                % Find place where the mean resources gathered after 40 swims
        % stabilizes; use this to make errorbars
        foundRes = squeeze(sum(efficiency,3));
meanp = zeros(1,1e4);for i=1:1e4;meanp(i) = mean(sum(foundRes(:,round(1:i)+1),1));end
figure(71);hold on;plot(meanp,'b');
xxx = zeros(1,1e4);for i=201:1e4-200;xxx(i) = var(meanp(i-200:i+200))/mean(meanp(i-200:i+200));end
xxx(1:200) = 100;
r = find(xxx<=1e-5);
% nump = r(1);
nump = 1520;%1304; %hard-coded, see above

eff_50 = cumsum(sum(sum(efficiency,3),2))/size(pos,3);
figure(10);hold on;plot(eff_50,'b')
% eff_50 = squeeze(sum(efficiency,3));
% [a b] = max(eff_50,[],1);
% eff_50 = b(a~=0);


sub50 = mean(cumsum(sum(efficiency(:,1:nump,:),3),1),2);
sub50std = std(cumsum(sum(efficiency(:,1:nump,:),3),1),0,2)/sqrt(nump);
figure(69);hold on;errorbar(sub50,sub50std);

hold50 = cumsum(sum(efficiency(:,1:nump,:),3),1);


% clear pos pos_HMM_50 efficiency
efficiency_50 = efficiency;
clear efficiency pos;

pos = pos_HMM_75;
% pos = pos(:,:,1:1e5);

% efficiency = zeros([size(pos,1),size(pos,3),length(xRes)],'uint8');
efficiency = false([size(pos,1),size(pos,3),length(xRes)]);

        for i=1:length(xRes) %count how many times simulated trajectory passed within a resource radius
            foundRes = bsxfun(@minus,squeeze(pos(:,1,:)),xRes(i)).^2 + bsxfun(@minus,squeeze(pos(:,2,:)),yRes(i)).^2  <= resR^2;
            [a b] = max(foundRes,[],1);
     

            

            
           onez = sub2ind(size(foundRes),b(a~=0),find(a));
           foundRes = false(size(foundRes));
           foundRes(onez)  = 1;
            efficiency(:,:,i) = foundRes;
            
        end

        % Find place where the mean resources gathered after 40 swims
        % stabilizes
        foundRes = squeeze(sum(efficiency,3));
meanp = zeros(1,1e4);for i=1:1e4;meanp(i) = mean(sum(foundRes(:,round(1:i)+1),1));end
figure(71);hold on;plot(meanp,'r');
xxx = zeros(1,1e4);for i=201:1e4-200;xxx(i) = var(meanp(i-200:i+200))/mean(meanp(i-200:i+200));end
xxx(1:200) = 100;
r = find(xxx<=1e-5);
% nump = r(1)
nump = 1520;%1304; %hard-coded, see above

%         eff_75 = squeeze(sum(efficiency,3));
% [a b] = max(eff_75,[],1);
% eff_75 = b(a~=0);

        eff_75 = cumsum(sum(sum(efficiency,3),2))/size(pos,3);
        
figure(10);hold on;plot(eff_75,'r')

figure(11);hold on;plot(eff_75./eff_50-1,'g');drawnow
figure(12);hold on;plot(eff_75-eff_50,'m');drawnow


hold75 = cumsum(sum(efficiency(:,1:nump,:),3),1);

    sub75 = mean(cumsum(sum(efficiency(:,1:nump,:),3),1),2);
sub75std = std(cumsum(sum(efficiency(:,1:nump,:),3),1),0,2)/sqrt(nump);
figure(69);hold on;errorbar(sub75,sub75std,'r');

eff_struct(eff_cnt).resD = resD;
eff_struct(eff_cnt).resDensity = resDensity;
eff_struct(eff_cnt).eff_50 = eff_50;
eff_struct(eff_cnt).eff_75 = eff_75;
eff_struct(eff_cnt).totsum = nansum(eff_75./eff_50-1);
eff_cnt = eff_cnt + 1;

end
end
end

%% For revision, plot fraction of total sims (constrained by stable N) in which efficiency curves are larger for HMM75 than HMM50
e50 = cumsum(sum(efficiency_50(:,1:nump,:),3),1);
e75 = cumsum(sum(efficiency(:,1:nump,:),3),1);

effcmp = zeros(1,40)
for i=1:40
    effcmp(i) = sum(e75(i,:)>e50(i,:))/length(find(e75(i,:)>0));
end
figure;plot(effcmp);
%% bin total angle by finding first index where total angle exceeds a bin edge
%for each simulation, and get the mean resources gathered by that point

% nump = 1304;

% pos_HMM_50 = pos_HMM_50(:,:,1:nump);
% pos_HMM_75 = pos_HMM_75(:,:,1:nump);
% efficiency = efficiency(:,1:nump,:);
% efficiency_50 = efficiency_50(:,1:nump,:);

            dY = diff(cat(1,zeros(1,size(pos_HMM_50,3)),squeeze(pos_HMM_50(:,2,:))));
            dX = diff(cat(1,zeros(1,size(pos_HMM_50,3)),squeeze(pos_HMM_50(:,1,:))));
            allTan = atan2(dY,dX);
            
            allTan = diff(cat(1,zeros(1,size(pos_HMM_50,3)),allTan));
            allTan = allTan.*(allTan < pi).*(allTan > -pi) + ... %fix pi flips
                (allTan >= pi) .* (allTan - 2*pi)  +  (allTan <= -pi) .* (allTan + 2*pi);

            allTan_50 = cumsum(abs(allTan),1);
            xx_50 = cumsum(sum(efficiency_50,3),1); 
            

binz = [0:pi:35];
angle_mat = [];
for i=1:length(binz)-1
    over_bin = allTan_50 >= binz(i) & allTan_50 < binz(i+1);
    [a b] = max(over_bin,[],1);
    xx_inds = xx_50(sub2ind(size(xx_50),b(a~=0),find(a)));
    angle_mat = [angle_mat mean(xx_inds(:))];
end

figure;plot(binz(1:end-1),angle_mat)

            dY = diff(cat(1,zeros(1,size(pos_HMM_75,3)),squeeze(pos_HMM_75(:,2,:))));
            dX = diff(cat(1,zeros(1,size(pos_HMM_75,3)),squeeze(pos_HMM_75(:,1,:))));
            allTan = atan2(dY,dX);
            
            allTan = diff(cat(1,zeros(1,size(pos_HMM_50,3)),allTan));
            allTan = allTan.*(allTan < pi).*(allTan > -pi) + ... %fix pi flips
                (allTan >= pi) .* (allTan - 2*pi)  +  (allTan <= -pi) .* (allTan + 2*pi);

            allTan_75 = cumsum(abs(allTan),1);
            xx_75 = cumsum(sum(efficiency,3),1); 
            

binz = [0:pi:35];
angle_mat = [];
for i=1:length(binz)-1
    over_bin = allTan_75 >= binz(i) & allTan_75 < binz(i+1);
    [a b] = max(over_bin,[],1);
    xx_inds = xx_75(sub2ind(size(xx_75),b(a~=0),find(a)));
    angle_mat = [angle_mat mean(xx_inds(:))];
end

hold on;plot(binz(1:end-1),angle_mat,'r')

%%  graph the stability of the mean for angle turned after 10 resources -- use this to calculate a separate N
angle_mat_50 = zeros(10,2e4);
fanothresh_50  = [];
for i = 1:10
cnt = 1;
for j = 1:10:2e5
    allTan_50_ = allTan_50(:,1:j);
     xx_50 = cumsum(sum(efficiency_50(:,1:j,:),3),1); 
     j
    xx_inds = xx_50 >= i;
    [a b] = max(xx_inds,[],1);
    xx_inds = allTan_50_(sub2ind(size(allTan_50_),b(a~=0),find(a)));
    angle_mat_50(i,cnt) = [mean(xx_inds)];

%angle_mat_50_sem = [angle_mat_50_sem std(xx_inds)/sqrt(nump)];
cnt = cnt+1;

%check if fano has dropped below threshold
cv_ind = find(diff(angle_mat_50(i,:))~=0 & ~isnan(angle_mat_50(i,1:end-1)));
this_amat = angle_mat_50(i,diff(angle_mat_50(i,:))~=0 & ~isnan(angle_mat_50(i,1:end-1)));
if length(this_amat)>40
    fano_fac = var(this_amat(end-40:end))/mean(this_amat(end-40:end));
    if fano_fac < 1e-5
        fanothresh_50(i) = cv_ind(end-20); 
        break;
    end
end
end

end

angle_mat_75 = zeros(10,2e4);
fanothresh_75  = [];
for i = 1:10
cnt = 1;
for j = 1:10:2e5
    allTan_75_ = allTan_75(:,1:j);
     xx_75 = cumsum(sum(efficiency(:,1:j,:),3),1); 
     j
    xx_inds = xx_75 >= i;
    [a b] = max(xx_inds,[],1);
    xx_inds = allTan_75_(sub2ind(size(allTan_75_),b(a~=0),find(a)));
    angle_mat_75(i,cnt) = [mean(xx_inds)];

%angle_mat_50_sem = [angle_mat_50_sem std(xx_inds)/sqrt(nump)];
cnt = cnt+1;

%check if fano has dropped below threshold
cv_ind = find(diff(angle_mat_75(i,:))~=0 & ~isnan(angle_mat_75(i,1:end-1)));
this_amat = angle_mat_75(i,diff(angle_mat_75(i,:))~=0 & ~isnan(angle_mat_75(i,1:end-1)));
if length(this_amat)>40
    fano_fac = var(this_amat(end-40:end))/mean(this_amat(end-40:end));
    if fano_fac < 1e-5
        fanothresh_75(i) = cv_ind(end-20); 
        break;
    end
end
end

end

%% Check p-values for i = 10
j = fanothresh_75(end)*10; %calculated this via skipping 10 for speed
 allTan_75_ = allTan_75(:,1:j);
     xx_75 = cumsum(sum(efficiency(:,1:j,:),3),1); 
     j
    xx_inds = xx_75 >= i;
    [a b] = max(xx_inds,[],1);
    xx_inds_ = allTan_75_(sub2ind(size(allTan_75_),b(a~=0),find(a)));

     allTan_50_ = allTan_50(:,1:j);
     xx_50 = cumsum(sum(efficiency_50(:,1:j,:),3),1); 
     j
    xx_inds = xx_50 >= i;
    [a b] = max(xx_inds,[],1);
    xx_inds = allTan_50_(sub2ind(size(allTan_50_),b(a~=0),find(a)));

%% Flagged for deletion -- but should double check that this is really unneeded
figure;%plot(angle_mat_50,'.-')
%errorbar(angle_mat_50,angle_mat_50_sem,'b');
plot(angle_mat_50);

angle_mat_75 = [];
angle_mat_75_sem = [];
for i=10
    xx_inds = xx_75 >= i;
    [a b] = max(xx_inds,[],1);
    xx_inds = allTan_75(sub2ind(size(allTan_50),b(a~=0),find(a)));
    angle_mat_75 = [angle_mat_75 mean(xx_inds)];
     angle_mat_75_sem = [angle_mat_75_sem std(xx_inds)/sqrt(length(xx_inds))];
%angle_mat_75_sem = [angle_mat_75_sem std(xx_inds)/sqrt(nump)];
end
hold on;%plot(angle_mat_75,'.-r')
errorbar(angle_mat_75,angle_mat_75_sem,'r');

%%  mean angle turned for 1, 2, 3, 4, resources gatehred, etc.
angle_mat_50 = [];
angle_mat_50_sem = [];
ind_mat_50 = [];
for i=1:7
    xx_inds = xx_50 >= i;
    [a b] = max(xx_inds,[],1);
    xx_inds = allTan_50(sub2ind(size(allTan_50),b(a~=0),find(a)));
    angle_mat_50 = [angle_mat_50 mean(xx_inds)];
     angle_mat_50_sem = [angle_mat_50_sem std(xx_inds)/sqrt(length(xx_inds))];
     ind_mat_50 = [ind_mat_50 length(xx_inds)/nump];
%angle_mat_50_sem = [angle_mat_50_sem std(xx_inds)/sqrt(nump)];
xx_inds_50 = xx_inds;
end
figure(99);%plot(angle_mat_50,'.-')
errorbar(angle_mat_50,angle_mat_50_sem,'b');
figure(100);plot(ind_mat_50,'.-');


angle_mat_75 = [];
angle_mat_75_sem = [];
ind_mat_75 = [];
for i=1:7
    xx_inds = xx_75 >= i;
    [a b] = max(xx_inds,[],1);
    xx_inds = allTan_75(sub2ind(size(allTan_50),b(a~=0),find(a)));
    angle_mat_75 = [angle_mat_75 mean(xx_inds)];
     angle_mat_75_sem = [angle_mat_75_sem std(xx_inds)/sqrt(length(xx_inds))];
     ind_mat_75 = [ind_mat_75 length(xx_inds)/nump];
%angle_mat_75_sem = [angle_mat_75_sem std(xx_inds)/sqrt(nump)];
end
figure(99);hold on;%plot(angle_mat_75,'.-r')
errorbar(angle_mat_75,angle_mat_75_sem,'r');
figure(100);hold on;plot(ind_mat_75,'.-r');

%% FOr elife revision ... make plots of p value of rep number
% use t-test2 because that is the distribution i use for calculating
% statistical power
for i =1:1500
[h p(i)] = ttest2(hold75(40,1:i),hold50(40,1:i));
end
figure;plot(p)
hold on;plot([0 1500],[0.05 0.05],'--k')
xlabel('# reps','fontsize',20)
ylabel('p-value','fontsize',20);
set(gca,'fontsize',10);
%% Plot number of sims vs. stat power based on the "population" mean and 
% std for both conditions. These will be then used to plot panels in S13
load 15Jan14_posHMM_50_and_real_basedon3_23_1_newangledistribution_expi_118
load 26Jan16_allCum_50_75_anglenorm_forstatpower;
rng(2); %randomize trajectories
pos_HMM_50 = single(pos_HMM_50); %reduce vector file size
pos_HMM_75 = single(pos_HMM_real);
randangs = rand(1,size(pos_HMM_50,3))*2*pi;
randangs = repmat(randangs,[size(pos_HMM_50,1) 1]);
pos_(:,1,:) = squeeze(pos_HMM_50(:,1,:)).*cos(randangs) - squeeze(pos_HMM_50(:,2,:)).*sin(randangs);
pos_(:,2,:) = squeeze(pos_HMM_50(:,1,:)).*sin(randangs) + squeeze(pos_HMM_50(:,2,:)).*cos(randangs);
pos_HMM_50 = single(pos_);

pos_(:,1,:) = squeeze(pos_HMM_75(:,1,:)).*cos(randangs) - squeeze(pos_HMM_75(:,2,:)).*sin(randangs);
pos_(:,2,:) = squeeze(pos_HMM_75(:,1,:)).*sin(randangs) + squeeze(pos_HMM_75(:,2,:)).*cos(randangs);
pos_HMM_75 = single(pos_);

            dY = diff(cat(1,zeros(1,size(pos_HMM_50,3)),squeeze(pos_HMM_50(:,2,:))));
            dX = diff(cat(1,zeros(1,size(pos_HMM_50,3)),squeeze(pos_HMM_50(:,1,:))));
            allTan = atan2(dY,dX);
            
            allTan = diff(cat(1,zeros(1,size(pos_HMM_50,3)),allTan));
            allTan = allTan.*(allTan < pi).*(allTan > -pi) + ... %fix pi flips
                (allTan >= pi) .* (allTan - 2*pi)  +  (allTan <= -pi) .* (allTan + 2*pi);

            allTan_50 = cumsum(abs(allTan),1);
            xx_50 = hold50; 
            
                        dY = diff(cat(1,zeros(1,size(pos_HMM_75,3)),squeeze(pos_HMM_75(:,2,:))));
            dX = diff(cat(1,zeros(1,size(pos_HMM_75,3)),squeeze(pos_HMM_75(:,1,:))));
            allTan = atan2(dY,dX);
            
            allTan = diff(cat(1,zeros(1,size(pos_HMM_50,3)),allTan));
            allTan = allTan.*(allTan < pi).*(allTan > -pi) + ... %fix pi flips
                (allTan >= pi) .* (allTan - 2*pi)  +  (allTan <= -pi) .* (allTan + 2*pi);

            allTan_75 = cumsum(abs(allTan),1);
            xx_75 = hold75; 
            
    xx_inds = xx_50 >= 1;
    [a b] = max(xx_inds,[],1);
    xx_inds = allTan_50(sub2ind(size(allTan_50),b(a~=0),find(a)));
    a50 = mean(xx_inds);
    sa50 = std(xx_inds);
    xx_inds = xx_75 >= 1;
    [a b] = max(xx_inds,[],1);
    xx_inds = allTan_75(sub2ind(size(allTan_75),b(a~=0),find(a)));
    a75 = mean(xx_inds);
    
m50 = mean(hold50(40,:));
std50 = std(hold50(40,:));
m75 = mean(hold75(40,:));

nn = 1:1500;
powout = sampsizepwr('t',[m50 std50],m75,[],nn);
figure;plot(nn,powout);
powout = sampsizepwr('t',[a50 sa50],a75,[],nn);
hold on;plot(nn,powout,'r');
%% for dist normalized

load 15Jan14_posHMM_50_and_real_basedon3_23_1_newangledistribution_expi_118
load 15Jan14_posHMM_50_swimdist_15081_good_angle_hist
load 26Jan16_allCum_50_75_distnorm_forstatpower;

pos_HMM_50 = single(pos_HMM_50); %reduce vector file size
pos_HMM_75 = single(pos_HMM_real);
randangs = rand(1,size(pos_HMM_50,3))*2*pi;
randangs = repmat(randangs,[size(pos_HMM_50,1) 1]);
pos_(:,1,:) = squeeze(pos_HMM_50(:,1,:)).*cos(randangs) - squeeze(pos_HMM_50(:,2,:)).*sin(randangs);
pos_(:,2,:) = squeeze(pos_HMM_50(:,1,:)).*sin(randangs) + squeeze(pos_HMM_50(:,2,:)).*cos(randangs);
pos_HMM_50 = single(pos_);

pos_(:,1,:) = squeeze(pos_HMM_75(:,1,:)).*cos(randangs) - squeeze(pos_HMM_75(:,2,:)).*sin(randangs);
pos_(:,2,:) = squeeze(pos_HMM_75(:,1,:)).*sin(randangs) + squeeze(pos_HMM_75(:,2,:)).*cos(randangs);
pos_HMM_75 = single(pos_);



            dY = diff(cat(1,zeros(1,size(pos_HMM_50,3)),squeeze(pos_HMM_50(:,2,:))));
            dX = diff(cat(1,zeros(1,size(pos_HMM_50,3)),squeeze(pos_HMM_50(:,1,:))));
            allTan = atan2(dY,dX);
            
            allTan = diff(cat(1,zeros(1,size(pos_HMM_50,3)),allTan));
            allTan = allTan.*(allTan < pi).*(allTan > -pi) + ... %fix pi flips
                (allTan >= pi) .* (allTan - 2*pi)  +  (allTan <= -pi) .* (allTan + 2*pi);

            allTan_50 = cumsum(abs(allTan),1);
            xx_50 = hold50; 
            
                        dY = diff(cat(1,zeros(1,size(pos_HMM_75,3)),squeeze(pos_HMM_75(:,2,:))));
            dX = diff(cat(1,zeros(1,size(pos_HMM_75,3)),squeeze(pos_HMM_75(:,1,:))));
            allTan = atan2(dY,dX);
            
            allTan = diff(cat(1,zeros(1,size(pos_HMM_50,3)),allTan));
            allTan = allTan.*(allTan < pi).*(allTan > -pi) + ... %fix pi flips
                (allTan >= pi) .* (allTan - 2*pi)  +  (allTan <= -pi) .* (allTan + 2*pi);

            allTan_75 = cumsum(abs(allTan),1);
            xx_75 = hold75; 
            
    xx_inds = xx_50 >= 1;
    [a b] = max(xx_inds,[],1);
    xx_inds = allTan_50(sub2ind(size(allTan_50),b(a~=0),find(a)));
    a50 = mean(xx_inds);
    sa50 = std(xx_inds);
    xx_inds = xx_75 >= 1;
    [a b] = max(xx_inds,[],1);
    xx_inds = allTan_75(sub2ind(size(allTan_75),b(a~=0),find(a)));
    a75 = mean(xx_inds);
    
m50 = mean(hold50(40,:));
std50 = std(hold50(40,:));
m75 = mean(hold75(40,:));

nn = 1:500;
powout = sampsizepwr('t',[m50 std50],m75,[],nn);
figure;plot(nn,powout);hold on;
powout = sampsizepwr('t',[a50 sa50],a75,[],nn);
plot(nn,powout,'r');
%% Calcuate # sims needed to reach 0.9 stat power given the difference 
%in mean/std over 1e6 sims that I can generate
load 26Jan16_allCum_50_75_anglenorm_forstatpower;
pows = zeros(1,1e4);
for i=1:length(pows)
    i
    m50 = mean(hold50(40,1+round(rand(10000,1)*(1e6-1))));
	s50 = std(hold50(40,1+round(rand(10000,1)*(1e6-1))));
    m75 = mean(hold75(40,1+round(rand(10000,1)*(1e6-1))));
	s75 = std(hold75(40,1+round(rand(10000,1)*(1e6-1))));
    pows(i) = sampsizepwr('t',[m50 mean([s50 s75])],m75,0.9,[],'tail','both');
end

figure;hist(pows,100);

%% Calculate stat power using monte carlo
load 26Jan16_allCum_50_75_anglenorm_forstatpower;

currp = zeros(1,1e4);
trackp = zeros(1,1e4);
for i=1000
    tic
    for j=1:10000
       
        currp(j) = ranksum(hold50(40,1+round(rand(i,1)*(1e6-1))), ...
            hold75(40,1+round(rand(i,1)*(1e6-1))));
    end
    toc
    trackp(i) = mean(currp<=0.05); %what ratio of random draws are significant?
end
figure;plot(trackp);
title('statpower, angle norm');
%%
load 26Jan16_allCum_50_75_distnorm_forstatpower;
currp = zeros(1,1e4);
trackp = zeros(1,1e4);
for i=100
    tic
    for j=1:10000
       
        currp(j) = ranksum(hold50(40,1+round(rand(i,1)*(1e6-1))), ...
            hold75(40,1+round(rand(i,1)*(1e6-1))));
    end
    toc
    trackp(i) = mean(currp<=0.05); %what ratio of random draws are significant?
end
figure;plot(trackp);
title('statpower, dist norm');
%% For elife revision...load in eff_struct and explore data

cd('D:\HBO-PAPER-EMOO-PC\HBO PAPER');
%load 20Dec15_effstruct_formarkov;
load 12Jan16_effstruct_formarkov_swimdist;


for i=1:length(eff_struct)
    figure(61);hold on;plot(eff_struct(i).eff_75-eff_struct(i).eff_50,'k');
    figure(63);hold on;plot(100*(eff_struct(i).eff_75./eff_struct(i).eff_50-1),'k');
end

std_50 = [];
std_75 = [];
for i=1:length(eff_struct)
    this_50 = [eff_struct(1:i).eff_50];
    std_50 = [std_50 std(this_50(end,:))];
    this_75 = [eff_struct(1:i).eff_75];
    std_75 = [std_75 std(this_75(end,:))];
end

figure;plot(std_50);
hold on;plot(std_75,'r');

stdplat = 20; %num of sims after which the total std stabilizes
this_50 = [eff_struct(1:20).eff_50];
figure;errorbar(mean(this_50,2),std(this_50,0,2)/sqrt(stdplat));
this_75 = [eff_struct(1:20).eff_75];
hold on;errorbar(mean(this_75,2),std(this_75,0,2)/sqrt(stdplat),'r');

figure;hold on;
for i=1:20
    plot([1 2],[this_75(end,i) this_50(end,i)],'.-')
end

errorbar([1 2],[mean(this_75(end,:)) mean(this_50(end,:))],[std(this_75(end,:)) std(this_50(end,:))]/sqrt(20),'k')

%For normalized curve with error bars:
stdplat = 20; %num of sims after which the total std stabilizes
this_50_ = [eff_struct(1:20).eff_50];
this_75_ = [eff_struct(1:20).eff_75];
this_50 = bsxfun(@rdivide,this_50_,max(this_75_,[],1));
this_75 = bsxfun(@rdivide,this_75_,max(this_75_,[],1));
figure;errorbar(mean(this_50,2),std(this_50,0,2)/sqrt(stdplat));

hold on;errorbar(mean(this_75,2),std(this_75,0,2)/sqrt(stdplat),'r');
%%
meanp = zeros(1,1e3);for i=1:1e3;meanp(i) = mean(sum(foundRes(:,round(1:i)+1),1));end
figure;plot(meanp)




% For Misha: plot p value over rep #

for i =1:2000
[p(i) h] = ranksum(hold75(40,1:i),hold50(40,1:i));
end

%% FOr more involved angle turned statistics, using the fanothresh variables calculated above
%load 11Feb16_intelli_critN_swimdist
load 11Feb16_intelli_critN_anglenorm %note, for this to work, must load in proper allTan_75 matrices using script above

mean_75 = [];
sem_75 = [];
mean_50 = [];
sem_50 = [];

succ_75 = [];
succ_50 = [];
for cc = 10%:10
    
j = min([fanothresh_75(cc)*10 fanothresh_50(cc)*10]); %calculated this via skipping 10 for speed
 allTan_75_ = allTan_75(:,1:j);
     xx_75 = cumsum(sum(efficiency(:,1:j,:),3),1); 
     j
    xx_inds = xx_75 >= cc;
    [a b] = max(xx_inds,[],1);
    xx_inds_ = allTan_75_(sub2ind(size(allTan_75_),b(a~=0),find(a)));
    mean_75(cc) = mean(xx_inds_);
    sem_75(cc) = std(xx_inds_)/sqrt(length(xx_inds_));
    succ_75(cc) = length(xx_inds_)/j;
    
     allTan_50_ = allTan_50(:,1:j);
     xx_50 = cumsum(sum(efficiency_50(:,1:j,:),3),1); 
     j
    xx_inds = xx_50 >= cc;
    [a b] = max(xx_inds,[],1);
    xx_inds = allTan_50_(sub2ind(size(allTan_50_),b(a~=0),find(a)));
    mean_50(cc) = mean(xx_inds);
    sem_50(cc) = std(xx_inds)/sqrt(length(xx_inds));
    succ_50(cc) = length(xx_inds)/j;

end
    figure(1);errorbar(mean_75,sem_75,'r');
    hold on;errorbar(mean_50,sem_50)
    
    figure(2);hold on;plot(succ_50,'.-');hold on;plot(succ_75,'.-r')
