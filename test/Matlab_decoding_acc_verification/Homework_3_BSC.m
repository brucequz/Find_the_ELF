trellis = poly2trellis(4,[15 17]);


SNR = [0.002];
workers = 10;
%message = [1 1 0 1 ];
count = 1;
FER = zeros(1,7);
for cur_snr = SNR
    total_error = 0;
    total_count = 0;
    cur_error = zeros(1,workers);
    curtotal = zeros(1,workers);
    parfor ii = 1:workers
        chan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)','SNR',cur_snr);
        convencoder = comm.ConvolutionalEncoder(trellis,"TerminationMethod","Terminated");
        modBPSK = comm.BPSKModulator;
        while cur_error(ii)<10        
        curtotal(ii)=curtotal(ii)+1;
        message = (randi(2,[1,64])-1).';
        codeword = convencoder(message);
        modSignal = real(modBPSK(codeword));
        %%add BSC noise
        ran_err = binornd(1,cur_snr,1,134).';
        ran_err(ran_err==1)=-1;
        ran_err(ran_err==0)=1;
        receivedSignal=modSignal.*ran_err;
        estimation = Viterbi_Decoder(trellis,receivedSignal);
        if(~isequal(message.',estimation(1:64)))
            cur_error(ii)=cur_error(ii)+1;
        end
        end
    end
    fprintf("%f \n", sum(cur_error)/sum(curtotal));
    FER(count)=sum(cur_error)/sum(curtotal);
    count=count+1;
end


%message = zeros(1,64).';
message = [1 1 0 1 1 0 1].';


