%% Test arrayDim Object

clear all; clear classes;
test = arrayDim('dimName',{'Time' 'Frequency' 'Channel'},...
                'dimSize',{1000 100 5},...
                'dimValues',{[1:1000] [] []},...
                'dimLabels',{[] [] {'A' 'B' 'C' 'D' 'E' }});
              
test2 = arrayDim('dimName',{'Time' 'Frequency' 'Channel'},...
                 'dimSize',{999 100 5},...
                 'dimValues',{[1001:1999] [] []});