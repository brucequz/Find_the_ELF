function [estimation] = Viterbi_Decoder(trellis,Channel_observations)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
mapping_tab = [0 0; 0 1; 1 0; 1 1];
mapping_tab = [1 1; 1 -1; -1 1; -1 -1];
total_length = length(Channel_observations)/2;
metrics = 3000*ones(trellis.numStates, 1+total_length);
metrics(1,1)=0;
previous_stat = zeros(trellis.numStates, total_length);
previous_input = zeros(trellis.numStates, total_length);

for cur_stage = 1:total_length
    cur_channel_obs = Channel_observations(2*cur_stage-1:2*cur_stage).';
    for cur_state = 1:trellis.numStates
        % we first check 0
        output_0 = trellis.outputs(cur_state,1);
        next_state = trellis.nextStates(cur_state,1)+1;
        cur_metric = norm(mapping_tab(output_0+1,:)-cur_channel_obs)^2;
        total_metric = metrics(cur_state,cur_stage)+cur_metric;
        if (total_metric<=metrics(next_state,cur_stage+1))
            metrics(next_state,cur_stage+1)=total_metric;
            previous_input(next_state,cur_stage)=0;
            previous_stat(next_state,cur_stage)=cur_state;
        end
        % then we chekc 1
        output_1 = trellis.outputs(cur_state,2);
        next_state = trellis.nextStates(cur_state,2)+1;
        cur_metric = norm(mapping_tab(output_1+1,:)-cur_channel_obs)^2;
        total_metric = metrics(cur_state,cur_stage)+cur_metric;
        if (total_metric<=metrics(next_state,cur_stage+1))
            metrics(next_state,cur_stage+1)=total_metric;
            previous_input(next_state,cur_stage)=1;
            previous_stat(next_state,cur_stage)=cur_state;
        end
    end
end

%disp(metrics);

estimation = zeros(1,total_length);
for ii = total_length:-1:1
    if ii == total_length
        estimation(ii)=previous_input(1,ii);
        p_state = previous_stat(1,ii);
    else
        estimation(ii)=previous_input(p_state,ii);
        p_state = previous_stat(p_state,ii);
    end
end
end