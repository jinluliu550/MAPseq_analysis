%% Organize data from MAPseq/BARseq experiments
%make the projection matrices from the MAPseq barcodematrices
oldfolder=cd('XC9');
[BXC91,XC91spikes,XC91filtaligneddepth,XC91AllPN,XC91LocAllPN]=analyzebarcodematrix1('barcodematrixXC9.mat','individualdatawithdistancealigned.mat','spikes.mat');
%when prompted, input lower bound = 1100; upper bound = 1300;
cd('../XC14');

[BXC141,XC141spikes,~,~,~]=analyzebarcodematrix1('barcodematrixZL113XC.mat','none','spikes.mat');

cd('../XC28');

[BXC281,XC281spikes,XC281filtaligneddepth,XC281AllPN,XC281LocAllPN]=analyzebarcodematrix1('barcodematrix.mat','individualdatawithdistancealignedcombinedXC28.mat','spikes.mat');
%when prompted, input lower bound=950; upper bound=1150;

cd(oldfolder);

%Get individual spike-in normalized files
brain_C9=XC91AllPN./repmat(XC91spikes,length(XC91AllPN),1)*max(XC91spikes);
brain_C14=BXC141./repmat(XC141spikes,length(BXC141),1)*max(XC141spikes);
brain_C28=XC281AllPN./repmat(XC281spikes,length(XC281AllPN),1)*max(XC281spikes);

% Save individual dataframes
writematrix(brain_C9,'brain_c9.csv');
writematrix(brain_C14,'brain_c14.csv');
writematrix(brain_C28,'brain_c28.csv');