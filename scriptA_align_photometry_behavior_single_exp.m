%% Load photometry and behavior data, normalize to isosbestic, and align
%% Author: Ilaria Carta 
% process each experiment individually

clear all
close all
addpath 'C:\Users\user_domain\datafolder';                                 %specify data folder path

%% open 560 trace
filename = '041523_G7F2_03_410';                                           %use sample data provided or replace with data in the same format
delimiter = ' ';
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
data470 = dataArray{:, 1};
fclose(fileID);
clearvars filename delimiter formatSpec fileID dataArray ans;

%% open 470 trace
filename = '041523_G7F2_03_470';                                           %use sample data provided or replace with data in the same format
delimiter = ' ';
formatSpec = '%f%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
data410 = dataArray{:, 1};
fclose(fileID);
clearvars filename delimiter formatSpec fileID dataArray ans;

recording_duration_s = floor(length(data470)/20);                          %calculates the rec duration in s (if the framerate used is 40, every second we get 20 datapoints of 470 signal and 20 of 410 signal)
data470_cut = data470(1:recording_duration_s*20);                          %crop the first second, we will do this also for the behavior

%% open behavior analysis file and specify parameters

behav_scoring = load('G7F2_PUP4_last.txt');
frame_rate=30;                                                             %specify framerate for behavior video, if different than photometry sampling rate
video_duration_s = floor(length(behav_scoring)/frame_rate);
data_offset = video_duration_s - recording_duration_s;                     %in case the video is longer than the recording, crop the eccess
behav_cut=behav_scoring(1:end,:);
stimulus=behav_cut(:,11);
maximum = max(max(stimulus));
[x,y]=find(stimulus==maximum);

pup_removal= x(3,1);
pup_intro= x(2,1);

investig=behav_cut(:,1)';                                                  %from the behavior datasheet, specify the columns corresponding of the behavior of interest
pupgroom=behav_cut(:,2)';


%% Define time and time bins
s=40;                                                                      %number of total datapoints in 1 second
m=40*60;                                                                   %number of datapoints in one minute
d=recording_duration_s;                                                    %recording length abbreviation
fs=s/2;                                                                    %number of datapoints for each trace (20 for 470, 20 for 410)
time=[1:length(data470)];
round_duration=time(1:fs*d);                                               %round the duration to the final second so it's easier to crop from recording and behavior, if needed
data1=data410(1:fs*d,1);
data2=data470(1:fs*d,1);


%% plot photometry data only (1st figure)

subplot (3,1,1);
plot(round_duration,data1);                                                     %plot raw 470 signal
xlabel('time(s)')
ylabel('avg pixel int')
xticks(1*fs:30*fs:fs*d)
xticklabels(0:30:d)
%set(gca,'XTick',[]);
axis tight

subplot (3,1,2);                                                           %plot raw 410 signal
plot(round_duration, data2,'r');
xlabel('time(s)')
ylabel('avg pixel int')
xticks(1*fs:30*fs:fs*d)
xticklabels(0:30:d)
%set(gca,'XTick',[]);
axis tight


%%Code for normalizing (kindly provided by Neurophotometrics)
figure
subplot(2,1,1)
time_F_matrix=[round_duration' data1 data2];
temp_fit = fit(time_F_matrix(:,1),time_F_matrix(:,2),'exp2');
plot(round_duration,data1-temp_fit(data1));
title('Linearize by Fitting with Biexponential')
xlabel('Time of day in total ms')
ylabel('F')

subplot (2,1,2)
%fit isosbestic
temp_fit = fit(time_F_matrix(:,1),time_F_matrix(:,3),'exp2'); %note, in this case, I am ignoring the first 2000 points where there is this weird fast decay to get a better fit. experimentally, i normally set things up so this isn't an important time in the recording / animal is habituating.
%scale fit to calcium dependent data
fit2 = fitlm(temp_fit(time_F_matrix(:,1)),time_F_matrix(:,2));
%calculate a crude dF/F by subtracting and then dividing by the fit
plot(time_F_matrix(:,1),100*(time_F_matrix(:,2)-(fit2.Fitted))./(fit2.Fitted))
xlabel('Time(s)')
ylabel('crude dF/F (%)')
title('Linearizing + Normalizing Using Isosbestic')

%data correction
FP.fakebackground = 0.017;
FP.corrected_data(:,1) = time_F_matrix(:,1); %we'll keep the time vector in the first column to make things easier

figure
for i = 2:3
       temp_fit = fit(time_F_matrix(:,1),time_F_matrix(:,3),'exp2'); %note, in this case, I am ignoring the first 2000 points where there is this weird fast decay to get a better fit. experimentally, i normally set things up so this isn't an important time in the recording / animal is habituating.
%scale fit to calcium dependent data
    fit2 = fitlm(temp_fit(time_F_matrix(:,1)),time_F_matrix(:,i));
%calculate a crude dF/F by subtracting and then dividing by the fit
subplot(3,1,i);
%plot(fit2);
%calculate dF/F by subtracting the background
    FP.corrected_data(:,i) = (time_F_matrix(:,i)-((fit2.Fitted)-FP.fakebackground))./(fit2.Fitted-FP.fakebackground);
clear temp_fit fit2
end

f5 = figure;
for i = 1:2
    subplot(2,1,i)
   
    plot(FP.corrected_data(:,1),FP.corrected_data(:,i+1))
     title('Corrected Data')
     xlabel('time of day in total ms')
     ylabel('%dF/F')
end

%%end of code section from Neurophotometrics

GCamP_corrected_NPM_method = smooth(FP.corrected_data(:,2),0.005);                   %smoothing (optional)
norm_GCamP_trace = GCamP_corrected_NPM_method - min(GCamP_corrected_NPM_method);     %match lowest point with zero
start_trace = 1;
end_trace= d*fs;

%% plot data from photometry and behavior together (2nd figure)
figure
ax1=axes('units','inches','position',[2 2.8 12 4]);                        %in order to overlay two plots with different number of data points such as photometry and behavioral data, we create individual axes, this line specifies the first axis position.
plot(round_duration(start_trace:end_trace),norm_GCamP_trace(start_trace:end_trace),'k','LineWidth',1.5);                                         
hold on 
%the following 3 lines are optional in case you wish to plot the
%introduction of the stimulus
qx = [0 1 1 0]*5+pup_intro*0.6667;                                         %plot the time of pup introduction (optional). Make sure to use a converting factor if dealing with different framerates. In our case the pup introduction time is expressed based on behavior time (30 fps), but recording framerate for 470 is 20, so if multiplied by 0.6667 (20/30) it will match the corresponding time in the recording
qy = [0 0 1 1]*6.5+(-0.045);
patch(qx, qy, 'Red', 'FaceColor', [1,0,0], 'LineWidth', 0.05, 'LineStyle', '-', 'FaceAlpha', 0.8,'edgealpha',0.2);

ylim([0 0.2]);
xticks(1*fs:20*fs:fs*d);
xticklabels(0:20:d);
axis tight;
box off
ylabel('\DeltaF/F');
xlabel('time(s)');
ylim([-0.0009 0.007]);                                                     %manually customize axis based on min-max
%xlim([start_trace end_trace]);

%% plot behaviors underneath trace
%
start_ethogram = 1*frame_rate;                                             %crop video as done for recording 
end_ethogram = d*frame_rate;

ax2=axes('units','inches','position',[2 3 12 .2]);                         %this specifies the position of the second axis for behavior ethogram (1st behavior)
imagesc(investig(start_ethogram:end_ethogram),'AlphaData', 0.6);
colormap(ax2,[0.8 0.8 0.8; 0 1 0; 0 0 1]);  
set(gca,'YTick',[]);
set(gca,'XTick',[]);
box off
hold on


ax3=axes('units','inches','position',[2 2.8 12 .2]);                       %this specifies the position of the second axis for behavior ethogram (2nd behavior)
imagesc(pupgroom(start_ethogram:end_ethogram),'AlphaData', 0.6);
colormap(ax3,[0.8 0.8 0.8; 1 1 1; 0.5 0 0.8]);           
set(gca,'YTick',[]);
set(gca,'XTick',[]);
box off
hold on

%colors available for other behaviors, if needed:
%colormap(ax4,[0.8 0.8 0.8; 1 1 1;1 0.3 0]);  orange
%colormap(ax5,[0.8 0.8 0.8; 1 1 1; 1 0 0]);  red
%colormap(ax5,[0.8 0.8 0.8; 1 1 1; 0 0.5 0.5]); teal 
%colormap(ax5,[0.8 0.8 0.8; 1 1 1;  0.2 0.8 0.8]); cyan
%colormap(ax2,[0.8 0.8 0.8; 0 1 0; 0 0 1]); blue 
%colormap(ax2,[0.8 0.8 0.8; 0 1 0; 1 0 1]); magenta
%colormap(ax4,[0.8 0.8 0.8; 1 1 1;  0 1 0]); light green
%colormap(ax4,[0.8 0.8 0.8; 1 1 1;  0.2 0.5 0.2]);  dark green
%colormap(ax5,[0.8 0.8 0.8; 1 1 1; 0.5 0 0.8]); purple 
