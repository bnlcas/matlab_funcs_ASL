function Mean_OCR = photon_sim(ICR, pulse_time)
%% Takes a Incident Count Rate (cps) and pulse shapping time (in us)
% and calculates the average OCR through repeated montecarlo sampling
rng('shuffle');

samplepts = 10^6;

Iterations = 1000;
OCR = zeros(Iterations,1);
% Generate ICR number of random photon events in increasing time order

parfor j = 1:Iterations
    start_times = sort(samplepts*rand(ICR,1));
    for i =2:(ICR-1)
        if (start_times(i) > (start_times(i-1)+pulse_time)) & ((start_times(i)+pulse_time)<start_times(i+1))
            OCR(j) = OCR(j)+1;
        end
    end
end

Mean_OCR = mean(OCR);
