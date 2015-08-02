% classification.m
%% import

clear all;
close all;
addpath /PaqTools
addpath /data

directory_name = '/data';
files = dir(fullfile(directory_name,'*.paq'));
% aside: 20100416_16_10_29.paq is very noisy; best to exclude it

%% store paq object data as voltage and current
%paqobj = paq('/home/lane/Desktop/Peter_Data_SpikeTests/20100419_17_45_39.paq');
% spike datasets:           % in /Peter_Data_SpikeTests/Practice_spiketests
% 080109-05_17_01PM.paq
% 080709-07_21_16PM.paq
% 20091106_14_42_48.paq
% 20091111_11_34_25.paq
% 20091113_14_51_42.paq

% 022409-12_47_12PM.paq     % in /data folder
voltage = cell(1,length(files));
current = cell(1,length(files));

for i = 1:length(files);
    paqobj(i) = paq(['/data/' files(i).name]);
    traces = paqobj(i).data('channels',1:4,[0,0]);
    % the spike test will have the maximum variance
    if max(var(traces(:,1:2))) > max(var(traces(:,3:4)));
        voltage{:,i} = traces(:,1);
        current{:,i} = traces(:,2);
    else
        voltage{:,i} = traces(:,3);
        current{:,i} = traces(:,4);
    end
    figure; plot(voltage{:,i})
    %figure; plot(current)
end


% Throw out Neuron 4 - cannot calculte R_input because all positive,
% spiking voltage depolarizations (not normal current steps)


%% determine idealized spike trains
close all
% samplingrate = 10000;  % Hz

% each cell of "spiketrains" will be a vector of spike times in seconds
spiketrains = cell(1,length(voltage));
% each cell of "amplitudes" will be a vector of actual voltage peaks of
% the idealized spike train
amplitudes = cell(1,length(voltage));
% each goes{1,neuron} will be 0's and 1's corresponding to the current step, 
% each goes{2,neuron} will be the current amplitude of that given step
% each goes{3,neuron} will be the first time of each current step
% each goes{4,neuron} will be the (normalized) current amplitude of each current step
goes = cell(4,length(voltage));
% for each neuron the time range will be different
time = cell(1,length(voltage));

for neuron = 1:length(voltage);
    samplingrate = paqobj(neuron).SampleRate;
    rawtrainV = voltage{:,neuron};
    rawtrainC = current{:,neuron};
    endtime = length(rawtrainV)/samplingrate;
    % create a vector of time in seconds that translates directly to the
    % index number of the voltage and current data
    times = 0:1/samplingrate:(endtime-1/samplingrate);
    % to find spikes, first take all time points when voltage was > 0
    abovezero = find(rawtrainV > 0);
    % initialize variables
    numspikes = 1; ends = cell(1,1); starts = cell(1,1); k=1;
    % take "mountain climbing" approach to find when spike peaks
    for instance = 2:length(abovezero);
        starts{1} = abovezero(1);
        if abovezero(instance)~=abovezero(instance-1)+1;
            if abovezero(instance)~=abovezero(instance-1)+2;
                numspikes = numspikes+1;
                ends{k} = abovezero(instance-1); 
                starts{k+1} = abovezero(instance);
                k=k+1;
            end
        end
        ends{k} = abovezero(length(abovezero));
    end
    amps = zeros(numspikes,1);
    spiketrain = zeros(numspikes,1);
    for spike = 1:numspikes;
        index = find(rawtrainV(starts{spike}:ends{spike})==max(rawtrainV(starts{spike}:ends{spike})));
        index = starts{spike}+index(1)-1;
        % the actual voltage amplitude of spike at time of idealized spike
        amps(spike) = rawtrainV(index);
        % the time stamp of the idealized spike
        spiketrain(spike) = times(index); % in secs
    end
    spiketrains{1,neuron} = spiketrain;
    amplitudes{1,neuron} = amps;
    
    % normalize the current trace steady state to zero
    normalized = rawtrainC - mode(rawtrainC);
    abovezero = find(abs(normalized) > 5);
    indices = zeros(length(times),1);
    indices(2:length(times),1) = find(times) - 1;
    gocues = ismember(indices,abovezero);
    % clean up single-point noise
    for go = 3:length(gocues)-2;
        if gocues(go+1) == 0;
            if gocues(go-1) == 0;
                gocues(go) = 0;
            elseif gocues(go-2) == 0;
                gocues(go) = 0;
            end
        elseif gocues(go+2) == 0;
            if gocues(go-1) == 0;
                gocues(go) = 0;
            end
        end
        if abs(gocues(go+1)) > 0;
            if abs(gocues(go-1)) > 0;
                gocues(go) = gocues(go-1);
            end
        end
    end
    cues = []; u = 1;
    camp = [];
    for go = 2:length(gocues);
        if gocues(go) > 0;
            if gocues(go-1) == 0;
                cues(u) = times(go);
                if length(gocues) > (go+samplingrate-1);
                    camp(u) = mean(normalized(go:go+samplingrate-1));
                end
                u=u+1;
            end
        end
    end
    goes{1,neuron} = gocues;
    goes{2,neuron} = goes{1,neuron}.*rawtrainC;
    goes{3,neuron} = cues;
    goes{4,neuron} = camp;
    time{1,neuron} = times;
end

% plot(time{1,1},goes{1,1},time{1,1},current{1,1})
% clear abovezero amps ans camp cues ends endtime go gocues index indices normalized rawtrainC rawtrainV traces starts spiketrain
clearvars -except spiketrains amplitudes goes time voltage current paqobj
%% Raster plot
figure
hold on    
k=0;
for j = 1:length(spiketrains);
    sptimes = spiketrains{j};
    for i = 1:length(sptimes)  	 	%Going through each spike time
        line([sptimes(i) sptimes(i)], [k k+1],'Color','k') 	%Create a tick mark at sp time
    end
    k=k+1;
end
%ylim([0 length(spike)])			
xlabel('Time (sec)')	
ylabel('Neurons')
title('Raster plot of all neurons')




%% determine ISI distributions
close all
isi_dist = cell(length(spiketrains),1);
last_isis = cell(length(spiketrains),1);
firingrates = zeros(length(spiketrains),17); % 17 chosen empirically
latency = cell(length(spiketrains),1);
curramp = cell(length(spiketrains),1);
for j = 1:length(spiketrains);
    sptimes = spiketrains{j};
    % initialize isi's for this neuron
    isis=[]; h=1;
    rate=[];
    % determine isi's
    for i = 1:length(sptimes)-1;
        isis(h,1) = sptimes(i+1)-sptimes(i);
        rate(h,1) = 1/isis(h,1);
        h=h+1;
    end
    % break down isi's for each current step
    h=1; start=1;
    breaks = find(isis>1);
    for stop = 1:length(breaks);
        isi_dist{j,h} = isis(start:breaks(stop)-1);
        if isfinite(mean(rate(start:breaks(stop)-1)));
            firingrates(j,h) = mean(rate(start:breaks(stop)-1));
        end
        h=h+1;
        start = breaks(stop)+1;
    end
    %figure; plot(isi_dist{j,h-1});
    last_isis{j,1} = isi_dist{j,h-1};
    
    % determine latency
    cues = goes{3,j}; % PROBLEM WITH GOES{3,28} -> 572 points  -> attempted noise reduction resolution, still VERY noisy; discard this neuron
    camp = goes{4,j};
    lat = zeros(length(cues),1);
    curr = zeros(length(cues),1);
    for go = 1:length(cues);
        distance = abs(cues(go)-sptimes);
        if min(distance) < 1;
            lat(go) = min(distance);
            if length(camp) >= go;
                curr(go) = camp(go);
            end
        end
    end
    lat = lat(lat~=0);
    curr = curr(curr~=0);
    latency{j,1} = lat;
    curramp{j,1} = curr;
    %figure; plot(latency{j}); title('Latency'); ylabel('Time (sec)')
    % most show an exponential decay in latency, although some (#34) show
    % logarithmic increase in latency, and some peak then decrease (#25)
    % also the incline of the curve varies somewhat
end
              

%% determine input resistance
R_input = zeros(length(spiketrains),23);
R_input1 = zeros(length(spiketrains),1);
Linearity = zeros(length(spiketrains),1);
for neuron = 1:length(spiketrains);
    V_ss = mode(voltage{neuron})*10^(-3);
    V_step = zeros(length(goes{4,neuron}),1);
    I = V_step;
    for step = 1:length(goes{4,neuron});
        nospike = find(goes{4,neuron}(goes{4,neuron}(step)==curramp{neuron}));
        if isempty(nospike);
            I(step) = goes{4,neuron}(step)*10^(-12);
            start = int32(goes{3,neuron}(step)*paqobj(neuron).SampleRate);
            ending = int32((goes{3,neuron}(step)+1)*paqobj(neuron).SampleRate);
            V_step(step) = mean(voltage{neuron}(start:ending));
            V_step(step) = V_step(step)*10^(-3);
            R_input(neuron,step) = (V_step(step) - V_ss)/I(step);
            if goes{4,neuron}(step) < 0 && goes{4,neuron}(step+1) > 0;
                R_input1(neuron) = mean(R_input(neuron,step-3:step-1));
            end
        end
    end
    V_step = V_step(V_step~=0);
    I = I(1:length(V_step));
    V_diff = V_step-V_ss;
    %{
    % Plot the V-I relationship
    figure; plot(I,V_diff);
    xlabel('Current')
    ylabel('Voltage difference')
    title(['V-I relationship, Neuron ' num2str(neuron)]);
    %}
    % Determine how linear the V-I relationship is:
    [P,S] = polyfit(I,V_diff,1);
    [Y,Delta] = polyval(P,I,S);
    Linearity(neuron) = mean(Delta);
    Linearity = Linearity(Linearity~=Inf);
    removeme = find(isnan(Linearity));
    for remove = 1:length(removeme);
        Linearity(removeme(remove)) = 0;
    end
end

figure; hist(Linearity);
xlabel('Degree of Linearity')
ylabel('# of Neurons')
title('Histogram of Linearity of V-I relationship for all neurons')

figure; hist(R_input1);
xlabel('Input Resistances (ohms)')
ylabel('# of Neurons')
title('Histogram of Average Input Resistance')

%% Determine voltage threshold for spiking
threshold = zeros(length(spiketrains),1);
for neuron = 1:length(spiketrains);
    threshold(neuron) = (curramp{neuron}(1)*10^(-12))*R_input1(neuron); % in Volts
    threshold(neuron) = threshold(neuron)*10^3; % convert to mV
end


figure; hist(threshold); xlabel('Spiking Threshold in mV'); 
ylabel('# of Neurons'); title('Histogram of Spiking Thresholds')

%% Latency metric
nneurons = length(voltage);
r_s = zeros(nneurons,1);
log_lat = r_s;
avg_lat = r_s;
for i = 1:nneurons;
    numlat = length(latency{i,:});
    t = 0:numlat-1;
    [P,S] = polyfit(t,log(latency{i,:})',1);
    [r,p] = corrcoef(P(1)*t-P(2),log(latency{i,:})');
    r_s(i) = r(2);
    [Y,Delta] = polyval(P,t,S);
    log_lat(i) = mean(Delta);
    avg_lat(i) = mean(latency{i});
end

%figure; hist(r_s,10)
figure; hist(log_lat)
xlabel('Wellness of Logarithmic Fit to Latencies')
ylabel('# of Neurons')
title('Histogram of how logarithmic a neuron"s latencies are')  

figure; hist(avg_lat,8)
title('Histogram of Average Latencies')
ylabel('# of Neurons')
xlabel('Latency of first spike in current step after current steps on')


%% Firing Rate metrics
fr = zeros(length(spiketrains),1);
fr_std = zeros(length(spiketrains),1);
for neuron = 1:length(spiketrains);
    frs = firingrates(neuron,:);
    frs = frs(frs~=0);
    fr(neuron) = mean(frs);
    fr_std(neuron) = std(frs);
end

figure; hist(fr,9);
xlabel('Mean Firing Rate'); ylabel('# of Neurons')
title('Histogram of Firing Rates')

figure; hist(fr_std);
xlabel('Standard Deviation of Firing Rates'); ylabel('# of Neurons')
title('Histogram of Standard Deviations of Firing Rates')

%% Spike Amplitude Metrics
% First, we need to add a "current step" dimension to our "amplitudes" cell
% array.
sp_amplitudes = cell(length(spiketrains),1);
spiketrains_perstep = cell(length(spiketrains),1);
avg_fits = zeros(length(spiketrains),1);
amp_stds = zeros(length(spiketrains),1);
for neuron = 1:length(spiketrains);
    fit = zeros(length(goes{4,neuron}));
    stds = zeros(length(goes{4,neuron}));
    i=1;
    for step = 1:length(goes{4,neuron});
        start = goes{3,neuron}(step);
        ending = (goes{3,neuron}(step)+1);
        indices = find(spiketrains{neuron} > start & spiketrains{neuron} < ending);
        if isempty(indices) == 0;
            sp_amplitudes{neuron,i} = amplitudes{neuron}(indices);
            spiketrains_perstep{neuron,i} = spiketrains{neuron}(indices);
            nspikes = length(indices);
            x = 1:nspikes-1; % let's model all spikes excluding the first of each current step
            logs = log(sp_amplitudes{neuron,i})';
            [P,S] = polyfit(x,logs(1:nspikes-1),1);
            [Y,Delta] = polyval(P,x,S);
            onlyfin = isfinite(Delta);
            fit(i) = mean(Delta(onlyfin));
            stds(i) = std(sp_amplitudes{neuron,i});
            i=i+1;
        end
    end
    fit = fit(fit~=0);
    onlyfin = isfinite(fit);
    if isnan(mean(fit(onlyfin)))
        avg_fits(neuron) = 0;
    else
        avg_fits(neuron) = mean(fit(onlyfin));
    end
    stds = stds(stds~=0);
    if isnan(mean(stds))
        amp_stds(neuron) = 0;
    else
        amp_stds(neuron) = mean(stds);
    end
end
        
figure; hist(avg_fits);
title('Histogram of how logarithmic subsequent spike amplitudes in a current step are')
ylabel('# of Neurons')
xlabel(' Wellness of Fitted Logarithmic Curve to Subsequent Spike Amplitudes')

figure; hist(amp_stds);
title('Histogram of Spike Amplitude Standard Deviations (avg. within current step std)')
xlabel('Average Within-Current Step Standard Deviation of Spike Amplitudes')
ylabel('# of Neurons')


%% Spike Width
% strategy is to take the width at half voltage between peak and baseline
% where baseline is minimum voltage between last two spikes
widths = zeros(length(isi_dist),1);
widths_std = widths;
for neuron = 1:length(isi_dist);
    spwidth = zeros(size(isi_dist,2),1);
    for step = 1:size(isi_dist,2);
        if isempty(spiketrains_perstep{neuron,step}) == 0;
            nisis = size(isi_dist{neuron,step},1); % nisis will be nspikes-1
            if nisis == 0;
                peak_index = int32(spiketrains_perstep{neuron,step}(1)*paqobj(neuron).SampleRate);
                baseline = min(voltage{neuron}(peak_index:(peak_index+(paqobj(neuron).SampleRate)/2)));
            else
                after = int32(((isi_dist{neuron,step}(nisis))/2)*paqobj(neuron).SampleRate);
                secondtolast = int32(spiketrains_perstep{neuron,step}(length(spiketrains_perstep{neuron,step})-1)*paqobj(neuron).SampleRate); % spike... should be equal to nisis
                baseline = min(voltage{neuron}(secondtolast:(secondtolast+after)));
            end
            spwidths = zeros(nisis+1,1);
            spheights = zeros(nisis+1,1);
            stds = zeros(nisis+1,1);
            for i = 1:length(spiketrains_perstep{neuron,step}); % should = nisis+1
                peak_index = int32(spiketrains_perstep{neuron,step}(i)*paqobj(neuron).SampleRate);
                peak = voltage{neuron}(peak_index);
                spheights(i) = (peak+baseline)/2; % note that these are actually voltages, not time widths
                [a firstside] = min(abs(spheights(i)-voltage{neuron}(peak_index-300:peak_index)));
                [a secondside] = min(abs(spheights(i)-voltage{neuron}(peak_index:peak_index+300)));
                spwidths(i) = ((300-firstside)+secondside)/paqobj(neuron).SampleRate; % in seconds
            end
            spwidth(step) = mean(spwidths);
            stds(step) = std(spwidths);
        end
    end
    spwidth = spwidth(spwidth~=0);
    widths(neuron) = mean(spwidth);
    widths_std(neuron) = mean(stds);
end
        
figure; hist(widths);
xlabel('Average Time durations of Action Potentials at Half-Width')
ylabel('# of Neurons')
title('Histogram of Spike Widths')

figure; hist(widths_std);
xlabel('Average Standard Deviations of Spike Widths')
ylabel('# of Neurons')
title('Histogram of Standard Deviations of Spike Widths')


%% Important metrics
%{
METRICS:
R_input1(neuron)                % Average Input Resistance across curr. steps
R_input(neuron,current step)    % Each calculated input resistance
Linearity(neuron)               % How linear each V-I relationship is
threshold(neuron)               % spiking threshold in voltage
latency{neuron,current step}    % Time (sec) until first spike after current steps on
avg_lat(neuron)                 % Average latencies
log_lat(neuron)                 % How logarithmic a neuron's latencies are across current steps
fr(neuron)                      % Average firing rate (averaging fr across current steps)
fr_std(neuron)                  % Average standard deviation of firing rates (just looking at 1/ISIs)
                                    % i.e., to capture adaptation
avg_fits(neuron)                % How (on avg.) subsequent spike amplitudes fit a log. curve
amp_stds(neuron)                % Average standard deviation of spike amplitudes within current step
widths(neuron)                  % Time (sec) duration of spike at half-width
widths_std(neuron)              % Average standard deviation of within-current step a.p widths


isi_dist{neuron,current step}  % each row of cells is a different neuron, 
                               % column 1 is the first current step that
                               % there was spiking for a given neuron
(last_isis{neuron,1})
latency{neuron,1}       % a vector of latencies (sec) across current steps
amplitudes{1,neuron}    % peak action potential amplitude
firingrates(neuron,:)   % same kind of structure as isi_dist but a matrix
                        % with zeros where there was 1 or no spikes.
                        % Calculated by taking mean(1/isi_i)
curramp{neuron,1}       % each cell is a vector of current amplitudes
                        % corresponding to the current steps associated
                        % with spiking; so the first number is the current
                        % in pA of the first step that caused the neuron to
                        % spike.


OTHER STUFF:
Whole trace:
current{1,neuron}
voltage{1,neuron}
time{1,neuron}  % for each neuron the time range (secs) will be different

Spike data:
spiketrains{1,neuron}   % vectors of spike times in seconds
amplitudes{1,neuron}    % vectors of voltage amplitudes for each spike

Current steps:
goes{1,neuron} % are 0's and 1's corresponding to the current step, 
goes{2,neuron} % are current amplitudes for each given step
goes{3,neuron} % is the first time of each current step
goes{4,neuron} % is the normalized (relative to base line) amplitude (pA) of each current step

%}

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACTUAL CLASSIFICATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Actual Classification

% Let's do PCA on the variables (metrics) we've collected thus far
% Let X be a matrix of size n_neurons by n_variables

% because of Matlab truncation, we'll need to shrink R_input1
R_input1 = R_input1/mean(R_input1);
X = [R_input1 Linearity threshold avg_lat log_lat fr fr_std avg_fits amp_stds widths widths_std];
[pc,score,eigvals,tsquared] = princomp(X);
per_variance = eigvals/sum(eigvals);

[IDX,Centers,SumD,D] = kmeans(X,3);
figure;
plot3(X(IDX==1,1),X(IDX==1,2),X(IDX==1,3),'r.', ...
    X(IDX==2,1),X(IDX==2,2),X(IDX==2,3),'b.', ...
    X(IDX==3,1),X(IDX==3,2),X(IDX==3,3),'g.', Centers(:,1),Centers(:,2),Centers(:,3),'kx');
title('k-Means Clustering of Neurons')

[IDX4,Centers4,SumD4,D4] = kmeans(X,4);


%% Cross Validation of k-means method
permutations = 100;
cross_validation = zeros(size(goes,2),permutations);
accuracies = zeros(3,permutations);
for perm = 1:permutations;
    random_samples = randsample(30,round(rand*25));
    while length(random_samples) < 3;
        random_samples = randsample(30,round(rand*25));
    end
    X_valid = zeros(length(random_samples),11);
    for samp = 1:length(random_samples);
        X_valid(samp,:) = [R_input1(samp) Linearity(samp) threshold(samp) avg_lat(samp) log_lat(samp) fr(samp) fr_std(samp) avg_fits(samp) amp_stds(samp) widths(samp) widths_std(samp)];
    end
    IDX_valid = kmeans(X_valid,3);
    for samp = 1:length(random_samples);
        cross_validation(random_samples(samp),perm) = IDX_valid(samp);
    end
    same_1 = find(cross_validation(:,perm)==1);
    same_2 = find(cross_validation(:,perm)==2);
    same_3 = find(cross_validation(:,perm)==3);
    
    % do they cluster the same?
    accuracy = zeros(3,1);
    numdifferent1 = length(find(IDX(same_1)~=mode(IDX(same_1))));
    accuracy(1) = (length(same_1)-numdifferent1)/length(same_1);
    numdifferent2 = length(find(IDX(same_2)~=mode(IDX(same_2))));
    accuracy(2) = (length(same_2)-numdifferent2)/length(same_2);
    numdifferent3 = length(find(IDX(same_3)~=mode(IDX(same_3))));
    accuracy(3) = (length(same_3)-numdifferent3)/length(same_3);
    accuracies(:,perm) = accuracy;
    for i = 1:length(accuracies);
        if isnan(mean(accuracies(:,i))) == 1;
            accuracies = [accuracies(:,1:(i-1)) accuracies(:,(i+1):length(accuracies))];
        end
    end
end

mean_acc3 = mean(mean(accuracies));
            

% let's see what it's like for various k
ks = 2:6;
mean_acc = zeros(1,length(ks));
for k = 1:length(ks);
permutations = 100;
clusters = ks(k);
cross_validation = zeros(size(goes,2),permutations);
accuracies = zeros(clusters,permutations);
for perm = 1:permutations;
    random_samples = randsample(30,round(rand*25));
    while length(random_samples) < clusters;
        random_samples = randsample(30,round(rand*25));
    end
    X_valid = zeros(length(random_samples),11);
    for samp = 1:length(random_samples);
        X_valid(samp,:) = [R_input1(samp) Linearity(samp) threshold(samp) avg_lat(samp) log_lat(samp) fr(samp) fr_std(samp) avg_fits(samp) amp_stds(samp) widths(samp) widths_std(samp)];
    end
    IDX_valid = kmeans(X_valid,clusters);
    for samp = 1:length(random_samples);
        cross_validation(random_samples(samp),perm) = IDX_valid(samp);
    end
    same = cell(1,clusters);
    for clusts = 1:clusters;
        same{clusts} = find(cross_validation(:,perm)==clusts);
    end
    
    % do they cluster the same?
    accuracy = zeros(clusters,1);
    for clus = 1:clusters
        % Note that when length(same{clus})==1, the accuracy will always be
        % 100%
        numdifferent = length(find(IDX(same{clus})~=mode(IDX(same{clus}))));
        accuracy(clus) = (length(same{clus})-numdifferent)/length(same{clus});
        accuracies(:,perm) = accuracy;
    end
    for i = 1:length(accuracies);
        if isnan(mean(accuracies(:,i))) == 1;
            accuracies = [accuracies(:,1:(i-1)) accuracies(:,(i+1):length(accuracies))];
        end
    end
end

mean_acc(k) = mean(mean(accuracies));
end            



figure; plot(ks,mean_acc);
xlabel('Number of Clusters k')
ylabel('Accuracy')
title('Cross-validation.  Similarity of Clustering')

close

% step-wise classification
% Step 1) Can RS and FS cells be pulled apart by a k-means of 2 along only
% planes of firing rate?

[IDX_fr Centers_fr] = kmeans(X(:,6),3);
plot3(X(IDX_fr==1,1),X(IDX_fr==1,2),X(IDX_fr==1,3),'r.', ...
    X(IDX_fr==2,1),X(IDX_fr==2,2),X(IDX_fr==2,3),'b.', ...
    X(IDX_fr==3,1),X(IDX_fr==3,2),X(IDX_fr==3,3),'g.', Centers_fr(1),Centers_fr(2),Centers_fr(3),'kx');
title('2-Means Clustering of Neurons by Firing Rate')


% close all

%% Neural Network Construction

% Train neural network
% targets = [0 0 0 1 1 1; 0 0 0 1 1 1];
targets = zeros(length(IDX),size(X,2));
targets(:,1) = IDX';
for var = 2:size(X,2);
    targets(:,var) = targets(:,1);
end
targets = targets';
input = X';
net=newff(input,targets,5);
net=train(net,input,targets);

%{
clearvars -except X IDX per_variance net paqobj input targets
save afewthings.mat
%}

%% Sanity Check
% The first 30 cells (out of 39) in the training set

%{
spiketest1
20100419_17_45_39.paq
RS adapting

spiketest2
20100419_17_47_17.paq
Rs adapting

spiktest1
20100416_16_10_29.paq
RS small adaptation

spiketest
220100414_12_31_04.paq
FS

spiketest2
20100405_15_15_45.paq
RS little or no adaptation
large AHP

spiketest2
20100405_15_17_20.paq
RS addapting almost FS
maybe a low thresh?
spikes eventually get exhausted

spiketest2
20100327_15_21_12.paq
-noise hump made cell spike
RS addapter
dend burst

spiketest2
20100223_10_32_51.paq
stutt

spiketest2 manual
20100223_10_34_38.paq
stutt / FS

spiketest2
20100219_10_31_39.paq
stutt FS

spiketest2
20100216_09_49_29.paq
only spiked twice for a long time
eventually burst

spiketest2
20100213_13_36_01.paq
huge inp res
bursting fs

spiketest2
20100210_12_23_00.paq
maybe LThresh spiker with addapting

spiketest2
20100209_11_53_12.paq
stutterer

spiketest1
20100208_18_24_28.paq
RS addapting

spiketest1
20100206_13_43_52.paq
RS addapting

spiketest1
20100204_15_34_45.paq
RS addapting

spiketest1:
\20100202_10_32_53.paq
regular/adapting.

spiketest2:
20100202_10_34_00.paq
regular/adapting. 

spiketest1
20100201_17_15_05.paq
RS addapting dend burst

spiketest1
20100128_12_26_57.paq
RS addapting
high thresh,...

spiketest2
20100127_16_02_16.paq
RS addapting (strong addap)

spiketest2
20100126_10_15_57.paq
Stutt FS

spiketest2
20100120_13_34_41.paq
FS maybe stutt

spiketest1
20100118_13_38_57.paq
crazy spiker

spiketest2
20100118_13_40_10.paq
crazy spiker

spiketest2
20100116_14_09_05.paq
very sparse spiking
dendritic burst

manual spiketest
20100116_14_10_58.paq
interesting spiking, not sure if I should call it RS

spiketest2
20100115_17_12_51.paq
RS addapting w/ burst

spiketest1
20100114_12_53_14.paq
RS addapting
bursts
%}
%{
%% Proper testing
close all
addpath /home/lane/Desktop/Peter_Data_SpikeTests/Practice_spiketests
directory_name = '/home/lane/Desktop/Peter_Data_SpikeTests/Practice_spiketests';
files_testing = dir(fullfile(directory_name,'*.paq'));

%% store paq object data as voltage and current
%paqobj = paq('/home/lane/Desktop/Peter_Data_SpikeTests/20100419_17_45_39.paq');
% spike datasets:           % in /Peter_Data_SpikeTests/Practice_spiketests
% 080109-05_17_01PM.paq
% 080709-07_21_16PM.paq
% 20091106_14_42_48.paq
% 20091111_11_34_25.paq
% 20091113_14_51_42.paq
close all
clear var

% 022409-12_47_12PM.paq     % in /data folder
voltage = cell(1,length(files_testing));
current = cell(1,length(files_testing));

for i = 1:length(files_testing);
    paqobj2(i) = paq(['/home/lane/Desktop/Peter_Data_SpikeTests/Practice_spiketests/' files_testing(i).name]);
    traces = paqobj2(i).data('channels',1:4,[0,0]);
    % the spike test will have the maximum variance
    if max(var(traces(:,2))) > max(var(traces(:,4)));
        voltage{:,i} = traces(:,1);
        current{:,i} = traces(:,2);
    else
        voltage{:,i} = traces(:,3);
        current{:,i} = traces(:,4);
    end
    figure; plot(voltage{i})
    title(['Test Neuron ' num2str(i)])
    %figure; plot(current)
end


% Throw out Neuron 4 - cannot calculte R_input because all positive,
% spiking voltage depolarizations (not normal current steps)


%% determine idealized spike trains
close all
% samplingrate = 10000;  % Hz

% each cell of "spiketrains" will be a vector of spike times in seconds
spiketrains = cell(1,length(voltage));
% each cell of "amplitudes" will be a vector of actual voltage peaks of
% the idealized spike train
amplitudes = cell(1,length(voltage));
% each goes{1,neuron} will be 0's and 1's corresponding to the current step, 
% each goes{2,neuron} will be the current amplitude of that given step
% each goes{3,neuron} will be the first time of each current step
% each goes{4,neuron} will be the (normalized) current amplitude of each current step
goes = cell(4,length(voltage));
% for each neuron the time range will be different
time = cell(1,length(voltage));

for neuron = 1:length(voltage);
    samplingrate = paqobj(neuron).SampleRate;
    rawtrainV = voltage{:,neuron};
    rawtrainC = current{:,neuron};
    endtime = length(rawtrainV)/samplingrate;
    % create a vector of time in seconds that translates directly to the
    % index number of the voltage and current data
    times = 0:1/samplingrate:(endtime-1/samplingrate);
    % to find spikes, first take all time points when voltage was > 0
    abovezero = find(rawtrainV > 0);
    % initialize variables
    numspikes = 1; ends = cell(1,1); starts = cell(1,1); k=1;
    % take "mountain climbing" approach to find when spike peaks
    for instance = 2:length(abovezero);
        starts{1} = abovezero(1);
        if abovezero(instance)~=abovezero(instance-1)+1;
            if abovezero(instance)~=abovezero(instance-1)+2;
                numspikes = numspikes+1;
                ends{k} = abovezero(instance-1); 
                starts{k+1} = abovezero(instance);
                k=k+1;
            end
        end
        ends{k} = abovezero(length(abovezero));
    end
    amps = zeros(numspikes,1);
    spiketrain = zeros(numspikes,1);
    for spike = 1:numspikes;
        index = find(rawtrainV(starts{spike}:ends{spike})==max(rawtrainV(starts{spike}:ends{spike})));
        index = starts{spike}+index(1)-1;
        % the actual voltage amplitude of spike at time of idealized spike
        amps(spike) = rawtrainV(index);
        % the time stamp of the idealized spike
        spiketrain(spike) = times(index); % in secs
    end
    spiketrains{1,neuron} = spiketrain;
    amplitudes{1,neuron} = amps;
    
    % normalize the current trace steady state to zero
    normalized = rawtrainC - mode(rawtrainC);
    abovezero = find(abs(normalized) > 5);
    indices = zeros(length(times),1);
    indices(2:length(times),1) = find(times) - 1;
    gocues = ismember(indices,abovezero);
    % clean up single-point noise
    for go = 3:length(gocues)-2;
        if gocues(go+1) == 0;
            if gocues(go-1) == 0;
                gocues(go) = 0;
            elseif gocues(go-2) == 0;
                gocues(go) = 0;
            end
        elseif gocues(go+2) == 0;
            if gocues(go-1) == 0;
                gocues(go) = 0;
            end
        end
        if abs(gocues(go+1)) > 0;
            if abs(gocues(go-1)) > 0;
                gocues(go) = gocues(go-1);
            end
        end
    end
    cues = []; u = 1;
    camp = [];
    for go = 2:length(gocues);
        if gocues(go) > 0;
            if gocues(go-1) == 0;
                cues(u) = times(go);
                if length(gocues) > (go+samplingrate-1);
                    camp(u) = mean(normalized(go:go+samplingrate-1));
                end
                u=u+1;
            end
        end
    end
    goes{1,neuron} = gocues;
    goes{2,neuron} = goes{1,neuron}.*rawtrainC;
    goes{3,neuron} = cues;
    goes{4,neuron} = camp;
    time{1,neuron} = times;
end

% plot(time{1,1},goes{1,1},time{1,1},current{1,1})
% clear abovezero amps ans camp cues ends endtime go gocues index indices normalized rawtrainC rawtrainV traces starts spiketrain
clearvars -except spiketrains amplitudes goes time voltage current paqobj
%% Raster plot
figure
hold on    
k=0;
for j = 1:length(spiketrains);
    sptimes = spiketrains{j};
    for i = 1:length(sptimes)  	 	%Going through each spike time
        line([sptimes(i) sptimes(i)], [k k+1],'Color','k') 	%Create a tick mark at sp time
    end
    k=k+1;
end
%ylim([0 length(spike)])			
xlabel('Time (sec)')	
ylabel('Neurons')
title('Raster plot of all neurons')




%% determine ISI distributions
close all
isi_dist = cell(length(spiketrains),1);
last_isis = cell(length(spiketrains),1);
firingrates = zeros(length(spiketrains),17); % 17 chosen empirically
latency = cell(length(spiketrains),1);
curramp = cell(length(spiketrains),1);
for j = 1:length(spiketrains);
    sptimes = spiketrains{j};
    % initialize isi's for this neuron
    isis=[]; h=1;
    rate=[];
    % determine isi's
    for i = 1:length(sptimes)-1;
        isis(h,1) = sptimes(i+1)-sptimes(i);
        rate(h,1) = 1/isis(h,1);
        h=h+1;
    end
    % break down isi's for each current step
    h=1; start=1;
    breaks = find(isis>1);
    for stop = 1:length(breaks);
        isi_dist{j,h} = isis(start:breaks(stop)-1);
        if isfinite(mean(rate(start:breaks(stop)-1)));
            firingrates(j,h) = mean(rate(start:breaks(stop)-1));
        end
        h=h+1;
        start = breaks(stop)+1;
    end
    %figure; plot(isi_dist{j,h-1});
    last_isis{j,1} = isi_dist{j,h-1};
    
    % determine latency
    cues = goes{3,j}; % PROBLEM WITH GOES{3,28} -> 572 points  -> attempted noise reduction resolution, still VERY noisy; discard this neuron
    camp = goes{4,j};
    lat = zeros(length(cues),1);
    curr = zeros(length(cues),1);
    for go = 1:length(cues);
        distance = abs(cues(go)-sptimes);
        if min(distance) < 1;
            lat(go) = min(distance);
            if length(camp) >= go;
                curr(go) = camp(go);
            end
        end
    end
    lat = lat(lat~=0);
    curr = curr(curr~=0);
    latency{j,1} = lat;
    curramp{j,1} = curr;
    %figure; plot(latency{j}); title('Latency'); ylabel('Time (sec)')
    % most show an exponential decay in latency, although some (#34) show
    % logarithmic increase in latency, and some peak then decrease (#25)
    % also the incline of the curve varies somewhat
end
              

%% determine input resistance
R_input = zeros(length(spiketrains),23);
R_input1 = zeros(length(spiketrains),1);
Linearity = zeros(length(spiketrains),1);
for neuron = 1:length(spiketrains);
    V_ss = mode(voltage{neuron})*10^(-3);
    V_step = zeros(length(goes{4,neuron}),1);
    I = V_step;
    for step = 1:length(goes{4,neuron});
        nospike = find(goes{4,neuron}(goes{4,neuron}(step)==curramp{neuron}));
        if isempty(nospike);
            I(step) = goes{4,neuron}(step)*10^(-12);
            start = int32(goes{3,neuron}(step)*paqobj(neuron).SampleRate);
            ending = int32((goes{3,neuron}(step)+1)*paqobj(neuron).SampleRate);
            V_step(step) = mean(voltage{neuron}(start:ending));
            V_step(step) = V_step(step)*10^(-3);
            R_input(neuron,step) = (V_step(step) - V_ss)/I(step);
            if goes{4,neuron}(step) < 0 && goes{4,neuron}(step+1) > 0;
                R_input1(neuron) = mean(R_input(neuron,step-3:step-1));
            end
        end
    end
    V_step = V_step(V_step~=0);
    I = I(1:length(V_step));
    V_diff = V_step-V_ss;
    %{
    % Plot the V-I relationship
    figure; plot(I,V_diff);
    xlabel('Current')
    ylabel('Voltage difference')
    title(['V-I relationship, Neuron ' num2str(neuron)]);
    %}
    % Determine how linear the V-I relationship is:
    [P,S] = polyfit(I,V_diff,1);
    [Y,Delta] = polyval(P,I,S);
    Linearity(neuron) = mean(Delta);
    Linearity = Linearity(Linearity~=Inf);
    removeme = find(isnan(Linearity));
    for remove = 1:length(removeme);
        Linearity(removeme(remove)) = 0;
    end
end

figure; hist(Linearity);
xlabel('Degree of Linearity')
ylabel('# of Neurons')
title('Histogram of Linearity of V-I relationship for all neurons')

figure; hist(R_input1);
xlabel('Input Resistances (ohms)')
ylabel('# of Neurons')
title('Histogram of Average Input Resistance')

%% Determine voltage threshold for spiking
threshold = zeros(length(spiketrains),1);
for neuron = 1:length(spiketrains);
    threshold(neuron) = (curramp{neuron}(1)*10^(-12))*R_input1(neuron); % in Volts
    threshold(neuron) = threshold(neuron)*10^3; % convert to mV
end


figure; hist(threshold); xlabel('Spiking Threshold in mV'); 
ylabel('# of Neurons'); title('Histogram of Spiking Thresholds')

%% Latency metric
nneurons = length(voltage);
r_s = zeros(nneurons,1);
log_lat = r_s;
avg_lat = r_s;
for i = 1:nneurons;
    numlat = length(latency{i,:});
    t = 0:numlat-1;
    [P,S] = polyfit(t,log(latency{i,:})',1);
    [r,p] = corrcoef(P(1)*t-P(2),log(latency{i,:})');
    r_s(i) = r(2);
    [Y,Delta] = polyval(P,t,S);
    log_lat(i) = mean(Delta);
    avg_lat(i) = mean(latency{i});
end

%figure; hist(r_s,10)
figure; hist(log_lat)
xlabel('Wellness of Logarithmic Fit to Latencies')
ylabel('# of Neurons')
title('Histogram of how logarithmic a neuron"s latencies are')  

figure; hist(avg_lat,8)
title('Histogram of Average Latencies')
ylabel('# of Neurons')
xlabel('Latency of first spike in current step after current steps on')


%% Firing Rate metrics
fr = zeros(length(spiketrains),1);
fr_std = zeros(length(spiketrains),1);
for neuron = 1:length(spiketrains);
    frs = firingrates(neuron,:);
    frs = frs(frs~=0);
    fr(neuron) = mean(frs);
    fr_std(neuron) = std(frs);
end

figure; hist(fr,9);
xlabel('Mean Firing Rate'); ylabel('# of Neurons')
title('Histogram of Firing Rates')

figure; hist(fr_std);
xlabel('Standard Deviation of Firing Rates'); ylabel('# of Neurons')
title('Histogram of Standard Deviations of Firing Rates')

%% Spike Amplitude Metrics
% First, we need to add a "current step" dimension to our "amplitudes" cell
% array.
sp_amplitudes = cell(length(spiketrains),1);
spiketrains_perstep = cell(length(spiketrains),1);
avg_fits = zeros(length(spiketrains),1);
amp_stds = zeros(length(spiketrains),1);
for neuron = 1:length(spiketrains);
    fit = zeros(length(goes{4,neuron}));
    stds = zeros(length(goes{4,neuron}));
    i=1;
    for step = 1:length(goes{4,neuron});
        start = goes{3,neuron}(step);
        ending = (goes{3,neuron}(step)+1);
        indices = find(spiketrains{neuron} > start & spiketrains{neuron} < ending);
        if isempty(indices) == 0;
            sp_amplitudes{neuron,i} = amplitudes{neuron}(indices);
            spiketrains_perstep{neuron,i} = spiketrains{neuron}(indices);
            nspikes = length(indices);
            x = 1:nspikes-1; % let's model all spikes excluding the first of each current step
            logs = log(sp_amplitudes{neuron,i})';
            [P,S] = polyfit(x,logs(1:nspikes-1),1);
            [Y,Delta] = polyval(P,x,S);
            onlyfin = isfinite(Delta);
            fit(i) = mean(Delta(onlyfin));
            stds(i) = std(sp_amplitudes{neuron,i});
            i=i+1;
        end
    end
    fit = fit(fit~=0);
    onlyfin = isfinite(fit);
    if isnan(mean(fit(onlyfin)))
        avg_fits(neuron) = 0;
    else
        avg_fits(neuron) = mean(fit(onlyfin));
    end
    stds = stds(stds~=0);
    if isnan(mean(stds))
        amp_stds(neuron) = 0;
    else
        amp_stds(neuron) = mean(stds);
    end
end
        
figure; hist(avg_fits);
title('Histogram of how logarithmic subsequent spike amplitudes in a current step are')
ylabel('# of Neurons')
xlabel(' Wellness of Fitted Logarithmic Curve to Subsequent Spike Amplitudes')

figure; hist(amp_stds);
title('Histogram of Spike Amplitude Standard Deviations (avg. within current step std)')
xlabel('Average Within-Current Step Standard Deviation of Spike Amplitudes')
ylabel('# of Neurons')


%% Spike Width
% strategy is to take the width at half voltage between peak and baseline
% where baseline is minimum voltage between last two spikes
widths = zeros(size(isi_dist,1),1);
widths_std = widths;
for neuron = 1:size(isi_dist,1);
    spwidth = zeros(size(isi_dist,2),1);
    for step = 1:size(isi_dist,2);
        if isempty(spiketrains_perstep{neuron,step}) == 0;
            nisis = size(isi_dist{neuron,step},1); % nisis will be nspikes-1
            if nisis == 0;
                peak_index = int32(spiketrains_perstep{neuron,step}(1)*paqobj(neuron).SampleRate);
                baseline = min(voltage{neuron}(peak_index:(peak_index+(paqobj(neuron).SampleRate)/2)));
            else
                after = int32(((isi_dist{neuron,step}(nisis))/2)*paqobj(neuron).SampleRate);
                secondtolast = int32(spiketrains_perstep{neuron,step}(length(spiketrains_perstep{neuron,step})-1)*paqobj(neuron).SampleRate); % spike... should be equal to nisis
                baseline = min(voltage{neuron}(secondtolast:(secondtolast+after)));
            end
            spwidths = zeros(nisis+1,1);
            spheights = zeros(nisis+1,1);
            stds = zeros(nisis+1,1);
            for i = 1:length(spiketrains_perstep{neuron,step}); % should = nisis+1
                peak_index = int32(spiketrains_perstep{neuron,step}(i)*paqobj(neuron).SampleRate);
                peak = voltage{neuron}(peak_index);
                spheights(i) = (peak+baseline)/2; % note that these are actually voltages, not time widths
                [a firstside] = min(abs(spheights(i)-voltage{neuron}(peak_index-300:peak_index)));
                [a secondside] = min(abs(spheights(i)-voltage{neuron}(peak_index:peak_index+300)));
                spwidths(i) = ((300-firstside)+secondside)/paqobj(neuron).SampleRate; % in seconds
            end
            spwidth(step) = mean(spwidths);
            stds(step) = std(spwidths);
        end
    end
    spwidth = spwidth(spwidth~=0);
    widths(neuron) = mean(spwidth);
    widths_std(neuron) = mean(stds);
end
        
figure; hist(widths);
xlabel('Average Time durations of Action Potentials at Half-Width')
ylabel('# of Neurons')
title('Histogram of Spike Widths')

figure; hist(widths_std);
xlabel('Average Standard Deviations of Spike Widths')
ylabel('# of Neurons')
title('Histogram of Standard Deviations of Spike Widths')


%% Neural Network Classification
load afewthings.mat;
% Let's do PCA on the variables (metrics) we've collected thus far
% Let X be a matrix of size n_neurons by n_variables

% because of Matlab truncation, we'll need to shrink R_input1
R_input1 = R_input1/mean(R_input1);
X_test = [R_input1 Linearity threshold avg_lat log_lat fr fr_std avg_fits amp_stds widths widths_std];
[pc,score,eigvals,tsquared] = princomp(X);
per_variance_test = eigvals/sum(eigvals);

[IDX_test,Centers,SumD,D] = kmeans(X_test,3);
figure;
plot3(X(IDX==1,1),X(IDX==1,2),X(IDX==1,3),'r.', ...
    X(IDX==2,1),X(IDX==2,2),X(IDX==2,3),'b.', ...
    X(IDX==3,1),X(IDX==3,2),X(IDX==3,3),'g.', Centers(:,1),Centers(:,2),Centers(:,3),'kx');
title('k-Means Clustering of Neurons')



% Neural Network Classification
input_testing = X_test';

trainOutputs = sim(net,X');
testOutputs = sim(net,input_testing);

%}
