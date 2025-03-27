%% Script for averaging data generated using scriptA and scriptB
%% Author: Ilaria Carta
%% This code works with previously saved variables (perievents created with scriptB for each exp, then stacked and saved in a matfile)
% computes zscores and area under curve

clear all
close all
data_path = 'C:\Users\user_domain\datafolder';                             %specify data path
addpath = data_path;
behavior_name = 'tailsuspension';
savdir = data_path;                                                        %select folder where you want to save results

load('all_tailsusp.mat');                                                  %use sample data provided or replace with data in the same format
%% split male and female datasets
all_events_M = [grooming_perievents_G4M4_ts; grooming_perievents_G6M1_ts; grooming_perievents_G6M3_ts; grooming_perievents_G6M4_ts; grooming_perievents_G7M2_ts];  %list all the perievent matrices that you need to analyze
all_events_F = [grooming_perievents_G5F1_ts; grooming_perievents_G5F2_ts; grooming_perievents_G7F6_ts];  %list all the perievent matrices that you need to analyze

%% calculations around event

all_events = [all_events_M; all_events_F];

figure

perievent_duration = 201;
d=perievent_duration;
plot(1:d, all_events(:,1:d), 'k');
hold on
plot(1:d, mean(all_events(:,1:d)), 'r', 'Linewidth',2);
axis tight
box off
%ylim([0 0.05])
title('deltaF/F all events') 


%% calculate zscores using baseline zscore (all datapoints and 1s bins)

length_baseline = 100;
rate=20;
n_bins=floor(length(all_events)./rate);

%%calculate zscore for Females

[m, n] = size(all_events_F);
zscore_container_F = zeros(m,n);
baseline_5s = zeros(m,length_baseline);
binned_zscore_container_F = zeros(n,n_bins);

for i=1:m
baseline_5s(i,:) = all_events_F(i,1:length_baseline);
zscore_container_F(i,:) =(all_events_F(i,:)-mean(baseline_5s(i,:)))./std(baseline_5s(i,:));
end

zscore_container_F(any(isnan(zscore_container_F),2),:) = [];                   %removes rows with NaNs
clear m n;
[m, n] = size(zscore_container_F);

for j=1:m
    for k=1:20:(n-20)                                                      %this number 20 is to bin the data every 20 datapoints, so every second
    binned_zscore_container_F(j,k)=nanmean(zscore_container_F(j,(k:k+20)));
    end
end

bins_zscore_container_F = binned_zscore_container_F(1:m,1:20:end);
clear baseline_5s;

%%calculate zscore for Males

[p, q] = size(all_events_M);

zscore_container_M = zeros(p,q);
baseline_5s = zeros(p,length_baseline);
binned_zscore_container_M = zeros(q,n_bins);

for i=1:p
baseline_5s(i,:) = all_events(i,1:length_baseline);
zscore_container_M(i,:) =(all_events(i,:)-mean(baseline_5s(i,:)))./std(baseline_5s(i,:));
end

zscore_container_M(any(isnan(zscore_container_M),2),:) = [];               %removes rows with NaNs
clear p q;
[p, q] = size(zscore_container_M);

for j=1:p
    for k=1:20:(q-20)                                                      %this number 20 is to bin the data every 20 datapoints, so every second
    binned_zscore_container_M(j,k)=nanmean(zscore_container_M(j,(k:k+20)));
    end
end

bins_zscore_container_M = binned_zscore_container_M(1:p,1:20:end);


%%peak zscore calculations
pre_event_start = 1;
pre_event_end = 101;
post_event_end = 201;                                                      %201 for 5 s, 401 for 10s
z_post_F = zscore_container_F(1:end, pre_event_end:post_event_end);        %change to 2:end if need to separate later from 1st event
z_post_F = zscore_container_F(1:end, pre_event_end:post_event_end);        %change to 2:end if need to separate later from 1st event
max_z_F=max(z_post_F, [], 2);                                              %for the row maximums

z_post_M = zscore_container_M(1:end, pre_event_end:post_event_end);        %change to 2:end if need to separate later from 1st event
z_post_M = zscore_container_M(1:end, pre_event_end:post_event_end);        %change to 2:end if need to separate later from 1st event
max_z_M=max(z_post_M, [], 2);                                              %for the row maximums

%% plot zscore data

% female heatmap
fps=20;                                                                    %specify frame rate
figure
subplot (3,2,1)
im_data = imagesc(zscore_container_F);
h = im_data;
h.AlphaData = ones(size(h.CData));
cb=colorbar;
caxis([0 10]);                                                             %caxis decides the lower and upper limits of the colorbar. 
xticks(1:5*fps:40*fps);
xticklabels(-5:5:30);
xlim([1 301])
xlabel('time from event onset(s)');
set(get(cb,'Title'),'String','z-score'); 
title('z score - all F events') 
xlim([1 301])
ax = gca;
ax.FontSize = 8; 

%male heatmap
subplot (3,2,2)
im_data = imagesc(zscore_container_M);
h = im_data;
h.AlphaData = ones(size(h.CData));
cb=colorbar;
caxis([0 10]);                                                             %caxis decides the lower and upper limits of the colorbar. 
xticks(1:5*fps:40*fps);
xticklabels(-5:5:30);
xlabel('time from event onset(s)');
set(get(cb,'Title'),'String','z-score'); 
title('z score - all M events') 
xlim([1 301])                                                              %customize axis limits
ax = gca;
ax.FontSize = 8; 

%line plot for both sexes
subplot (3,2,3)
stdshade(zscore_container_F, 0.2, 'm');                                    %stdshade function was downloaded separately
hold on
stdshade(zscore_container_M, 0.2, 'b');
axis tight
xticks(1:5*fps:40*fps);
xticklabels(-5:5:30);
xlabel('time from approach(s)');
ylabel('z-score');
title('avg z score') 
xlim([1 301])
ylim([-5 12]);                                                             %customize axis limits
title('mean z score and s.e.m.') 

%plot binned zscore
subplot (3,2,4)

avg_z_F = mean(bins_zscore_container_F);
avg_z_M = mean(bins_zscore_container_M);

st_error_F= sem(bins_zscore_container_F);
st_error_M= sem(bins_zscore_container_M);

s=scatter(1:n_bins,avg_z_F, 'm');
s.Marker = 'none';
er_M = errorbar(1:n_bins,avg_z_F,st_error_F);    
er_M.Color = 'm'; 
er_M.LineStyleMode = 'auto'; 

hold on
s=scatter(1:n_bins,avg_z_M, 'b');
s.Marker = 'none';
er_M = errorbar(1:n_bins,avg_z_M,st_error_M);    
er_M.Color = 'b'; 
er_M.LineStyleMode = 'auto'; 

ylabel('z-score');
xlabel('bins');
axis tight
box off
title('binned z-score and s.e.m.') 


%%Display and calculate area under the curve

subplot (3,2,5)                                                             

avg_F = mean(zscore_container_F);
avg_M = mean(zscore_container_M);

plot (1:401, avg_F(1:401));hold on
area(avg_F(1:401), 'Facecolor', 'm' , 'Facealpha', 0.2)
plot (1:401, avg_M(1:401));hold on
area(avg_M(1:401), 'Facecolor', 'b' , 'Facealpha', 0.2)
xlim([1 301]);                                                             %customize axis limits
xticks(1:5*fps:40*fps);
xticklabels(-5:5:30);
xlabel('time from event onset(s)');
box off
title('AUC') 

%calculate area under the curve 5s before and 5s after onset

auc_raw_M = trapz(zscore_container_M');                                    %the trapz function computes the AUC, while the area function is only for plotting
auc_raw_M_pre = trapz(zscore_container_M(1:p,1:201)');                     %isolate pre-onset 
auc_raw_M_post = trapz(zscore_container_M(1:p,201:301)');                  %isolate post-onset(make sure the pre and post time windows have the same length for the area to be comparable)

auc_raw_F = trapz(zscore_container_F');
auc_raw_F_pre = trapz(zscore_container_F(1:m,1:201)');
auc_raw_F_post = trapz(zscore_container_F(1:m,201:301)');


%%Plot AUC pre/post with error bars

subplot (3,2,6)

max_size = p+m;
auc_rawdata = nan(2,max_size);
auc_rawdata(1,1:m)=auc_raw_F_post(1:m);
auc_rawdata(2,1:p)=auc_raw_M_post(1:p);

bar(1:2,nanmean(auc_rawdata'));                                            %plot AUC values in a bargraph
hold on
plot(1:2,auc_rawdata,'o');
xlab = {'females' 'males'};
set(gca,'xtick',1:2,'xticklabel',xlab, 'Fontsize' , 8);
ylabel('AUC');
title('average AUC and s.e.m.') 

%save this dataset
dataset_name = strcat(behavior_name, 'results_all_trials');                %combine strings to save results
savdir = data_path;
save(fullfile(savdir,dataset_name),'zscore_container_F', 'zscore_container_M', 'bins_zscore_container_F', 'bins_zscore_container_M', 'max_z_F', 'max_z_M', 'auc_raw_F_pre', 'auc_raw_F_post', 'auc_raw_M_pre', 'auc_raw_M_post');

