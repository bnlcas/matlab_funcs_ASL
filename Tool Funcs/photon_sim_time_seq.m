function [ICR, OCR] = photon_sim(CountRate, pulse_time,clear_time)
%% Takes a Incident Count Rate (cps) and pulse shapping time (in us)
% and calculates the average OCR through repeated montecarlo sampling
rng('shuffle');

sample_rate = 10^6;
Iterations = 1000;

%% Loop through several simulated 1 second acuquisitions and take the ICR and OCR of each
ICRs = zeros(Iterations,1);
OCRs = zeros(Iterations,1);

parfor j = 1:Iterations

%% Generate Event Indecies
%event = sample_rate*sort(rand(CountRate,1));
event_ind = sample_rate*sort(rand(CountRate,1));

%% Clear out areas between cleartime and sample time
use_mc = true;
if use_mc
    clear_events = false(length(event_ind),1);
    for i = 2:length(event_ind)
        if (event_ind(i) < (event_ind(i-1)+pulse_time)) & (event_ind(i) > (event_ind(i-1)+clear_time))
            clear_events(i) = true;
        end
    end
    event_ind(clear_events) = [];
end

%% Plot PileUp
%ICR = length(event_ind); % number of events
%OCR = 0;
ICRs(j) = length(event_ind);
for i = 2:(ICRs(j)-1)
    
    if (event_ind(i)>(event_ind(i-1)+pulse_time) & ((event_ind(i)+pulse_time)<event_ind(i+1)))
        OCRs(j) = OCRs(j) +1;
    end
end
% ICRs(j) = ICR;
% OCRs(j) = OCR;



end
ICR = mean(ICRs);
OCR = mean(OCRs);

