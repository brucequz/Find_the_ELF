load("received_signal.mat", "received_signal");
load("messages.mat", "messages");
num_trials=size(received_signal,1);
trellis = poly2trellis(4,[13 17]);

d_min =6;

SNR = [1 2 2.5 3 3.5 4 4.5 5 6 7];
workers = 1;
%message = [1 1 0 1 ];
count = 1;
FER = zeros(1,7);
for cur_snr = SNR
    total_error = 0;
    total_count = 0;
    
    cur_error = 0;
    curtotal = 0;
    while cur_error < 200 
        for ii = 6:num_trials
        chan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)','SNR',cur_snr);
        convencoder = comm.ConvolutionalEncoder(trellis,"TerminationMethod","Terminated");
        modBPSK = comm.BPSKModulator;
               
            curtotal=curtotal+1;
            message = messages(ii,:);
            %message = (randi(2,[1,64])-1).';
            codeword = convencoder(message');
            %modSignal = real(modBPSK(codeword));
            % receivedSignal = chan(modSignal);
            receivedSignal = received_signal(ii,:)';
            estimation = Viterbi_Decoder(trellis,receivedSignal);
            if(~isequal(message,estimation(1:5)))
                cur_error=cur_error+1;
            end
        end
    end
    fprintf("%f \n", sum(cur_error)/sum(curtotal));
    FER(count)=sum(cur_error)/sum(curtotal);
    count=count+1;
end

