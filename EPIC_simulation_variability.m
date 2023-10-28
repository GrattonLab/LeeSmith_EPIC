%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPIC
%
% To run this script:
% No dataset or script required
%
% What this script does:
% Runs a simulation of data with a pre-set between-subject mean and variation
% but varying levels of within-subject variability
% Sim 1. Pre-set between-subject variability of 30 ms
% Sim 2. Varying between-subject variability (5ms 10 20 30...)
% Sim 3. Different number of simulated participants (100 200 300 400 500)
%
% What this script outputs:
% Figure 6.
%
% Created on 02/09/2023 by Caterina and Hyejin
% Last modified on 07/13/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng('shuffle')  % reset the seed of the random number generator

% Parameter settings
n_subs = 100;  % number of simulated subjects
bs_dist_n = 100;  % number of data points of between-subject distribution
ws_dist_n = 100;  % number of data points of within-subject distribution
numTest = 100;  % 100 repetitions to get mean and its error bars as 95% CI
UpB = 97.5;  % CI upper bound
LwB = 2.5;  % lower bound

% BETWEEN SUBJECT 'TRUE' DISTRIBUTION (values from Robinson and Steyvers' [2023] data)
bs_mean = 92;  % approx congruency effect (in ms)
bs_std = 30;  % approx congruency effect var across participants with MANY trials (in ms)
bs_std_list = [5 10 20 30 40 50 60];  % if you want to loop through different true bs variances
% note, bs std is SET!!

% WITHIN SUBJECT VARIANCE ESTIMATES (values estimated from our study)
ws_std_list = [90 60 45 30 20 15 10 7]./2; % approx from plot, Method 2 (ms)
% correspond with ws 95% CI with 50, 100, 200, 400, 800, 1600, 3200, 6400
% trials; divided by 2 so only a single std rather than 2 away

% Varying number of simulated subjects
n_subs_list = [50 100 200 300 400 500 1000];

%% Simulation 1. Fixed between-subject variance
% Preassignment
ws_exper_mean = nan(length(ws_std_list),n_subs);  % mean CE of a subject (across varying ws SD)
bs_exper_mean = nan(numTest,length(ws_std_list));  % group mean CE
bs_exper_std = nan(numTest,length(ws_std_list));  % (observed) CE between-subject SD
for l = 1:numTest
    % generate a distribution of random values with the between-subject estimates above
    bs_dist = randn(1,bs_dist_n)*bs_std+bs_mean;

    % loop through different within-subject variances
    for ws = 1:length(ws_std_list)
        ws_std = ws_std_list(ws);

        for n = 1:n_subs
            disp(['Testing' num2str(l) ': Subject#' num2str(n) ' step' num2str(ws)])
            % simulate a session as a random draw from this distribution
            ws_mean = bs_dist(randi(bs_dist_n,1,1));
            ws_dist = randn(1,ws_dist_n)*ws_std+ws_mean;
            ws_exper_mean(ws,n) = ws_dist(randi(ws_dist_n,1,1));
        end
    end
    % now calculate the apparent between-subject variance
    bs_exper_mean(l,:) = mean(ws_exper_mean,2);
    bs_exper_std(l,:) = std(ws_exper_mean,0,2);
end
mean_bs_exper_std = mean(bs_exper_std,1);  % mean SD of 100 testings
CI_bs_exper_std = (prctile(bs_exper_std,UpB)-prctile(bs_exper_std,LwB))/2;  % 95% CI

% Plot the results
cmap = hsv(length(bs_std_list));  % colormap
figure
a = errorbar(ws_std_list,mean_bs_exper_std,CI_bs_exper_std,'LineWidth',1.5); hold on
a.Marker = 'o';
%a.MarkerSize = 4;
a.Color = cmap(4,:);  % 4th is 30 ms
a.CapSize = 15;
yline(bs_std,'r--','LineWidth',1.5);
set(gca,'FontSize',16)
xlabel('Within-subject standard deviation (ms)','FontSize',20)
ylabel('Apparent between-subject standard deviation (ms)','FontSize',20)
ylim([0 85])
legend('','True between-subject standard deviation','FontSize',18)

%% Simulation 2. Varying the size of between-subject variability
% Preassignment
bs_std_cell = nan(length(bs_std_list),length(ws_std_list));
conIntvlh_cell = nan(length(bs_std_list),length(ws_std_list));
% loop through different true bs variances
for b = 1:length(bs_std_list)
    bs_std = bs_std_list(b);

    % loop through different ws variances
    ws_exper_sim = nan(n_subs,numTest,length(ws_std_list));
    for l = 1:numTest
        bs_dist = randn(1,bs_dist_n)*bs_std+bs_mean;  % generate a distribution of random values with the between sub estimates above

        for ws = 1:length(ws_std_list)
            ws_std = ws_std_list(ws);
            for n = 1:n_subs  % sample a sub at random from the bs distribution with replacement
                disp(['BS step' num2str(b) ' - Testing' num2str(l) ': Subject#' num2str(n) ' WS step' num2str(ws)])
                % set that to the 'true' ws_mean for this sub
                ws_mean = bs_dist(randi(bs_dist_n,1,1));

                % generate a rand dist for this subject w a given amt of variance
                ws_dist = randn(1,ws_dist_n)*ws_std+ws_mean;

                ws_exper_sim(n,l,ws) = ws_dist(randi(ws_dist_n,1,1));
            end
        end
    end

    % now calculate the apparent bs variance from random draw
    sim_bs_std = squeeze(std(ws_exper_sim,'omitnan'));
    bs_exper_std_M = mean(sim_bs_std);

    % Confidence interval
    conIntvlh = (prctile(sim_bs_std,UpB)-prctile(sim_bs_std,LwB))/2;
    bs_std_cell(b,:) = bs_exper_std_M;
    conIntvlh_cell(b,:) = conIntvlh;
end

% Plot results
figure
for b = 1:length(bs_std_list)
    a = errorbar(ws_std_list,bs_std_cell(b,:),conIntvlh_cell(b,:),'LineWidth',1.5); hold on  % random draw 원하면 bs_std_cell 대신 bs_exper_std_RD
    a.Marker = 'o';
    %a.MarkerSize = 4;
    a.Color = cmap(b,:);
    a.CapSize = 15;
end
set(gca,'FontSize',16)
xlabel('Within-subject standard deviation (ms)','FontSize',20)
ylabel('Apparent between-subject standard deviation (ms)','FontSize',20);
ylim([0 85])
lgd = legend('5 ms','10 ms','20 ms','30 ms','40 ms','50 ms','60 ms','Location','southeast','FontSize',12);
title(lgd,'True between-subject standard deviation','FontSize',13)

%% Simulation 3. Varying the number of simulated participants
bs_std = 30;
% Preassignment
bs_std_cell2 = nan(length(n_subs_list),length(ws_std_list));
conIntvlh_cell2 = nan(length(n_subs_list),length(ws_std_list));
% loop through different true bs variances
for s = 1:length(n_subs_list)
    n_subs = n_subs_list(s);

    % loop through different ws variances
    ws_exper_sim = nan(n_subs,numTest,length(ws_std_list));

    for l = 1:numTest
        bs_dist = randn(1,bs_dist_n)*bs_std+bs_mean;  % generate a distribution of random values with the between sub estimates above
        for ws = 1:length(ws_std_list)
            ws_std = ws_std_list(ws);
            for n = 1:n_subs  % sample a sub at random from the bs distribution with replacement
                disp(['Sample size step' num2str(s) ' - Testing' num2str(l) ': Subject#' num2str(n) ' WS step' num2str(ws)])
                % set that to the 'true' ws_mean for this sub
                ws_mean = bs_dist(randi(bs_dist_n,1,1));

                % generate a rand dist for this subject w a given amt of variance
                ws_dist = randn(1,ws_dist_n)*ws_std+ws_mean;
                ws_exper_sim(n,l,ws) = ws_dist(randi(ws_dist_n,1,1));
            end
        end
    end

    % now calculate the apparent bs variance
    % BS calculated as the mean of 100 simulations
    sim_bs_std = squeeze(std(ws_exper_sim,'omitnan'));
    bs_exper_std_M = mean(sim_bs_std);
    bs_std_cell2(s,:) = bs_exper_std_M;

    % Confidence Interval
    conIntvlh2 = (prctile(sim_bs_std,UpB)-prctile(sim_bs_std,LwB))/2;
    conIntvlh_cell2(s,:) = conIntvlh2;
end

% Plot the results
figure
cmap = hsv(length(n_subs_list));  % colormap
for s = 1:length(n_subs_list)
    a = errorbar(ws_std_list,bs_std_cell2(s,:),conIntvlh_cell2(s,:),'LineWidth',1.5); hold on
    a.Marker = 'o';
    %a.MarkerSize = 4;
    a.Color = cmap(s,:);
    a.CapSize = 15;
end
set(gca,'FontSize',16)
xlabel('Within-subject standard deviation (ms)','FontSize',20)
ylabel('Apparent between-subject standard deviation (ms)','FontSize',20)
ylim([0 85])
lgd = legend('50','100','200','300','400','500','1,000','Location','southeast','FontSize',14);
title(lgd,'Sample size','FontSize',16)

save('sim_varResults.mat','mean_bs_exper_std','CI_bs_exper_std','bs_std_cell','conIntvlh_cell',...
    'bs_std_cell2','conIntvlh_cell2')