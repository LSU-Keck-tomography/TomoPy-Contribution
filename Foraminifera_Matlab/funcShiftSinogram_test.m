%SHIFT GREATER THAN 0
test=[1:5
6:10
11:15
16:20
21:25
26:30
31:35
36:40]

shift=2;
numberColumns=5;
numberAngles=10;

transposetest=transpose(test)

newtest=transposetest(:,(shift:numberColumns))

[r,c]=size(newtest)

paddedtest=padarray(newtest,[numberAngles-r numberColumns-c],'post')

finaltest=transpose(paddedtest)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SHIFT LESS THAN 0
test=[1:5
6:10
11:15
16:20
21:25
26:30
31:35
36:40]

shift=-2;
numberColumns=5;
numberAngles=10;

transposetest=transpose(test)

newtest=transposetest(:,(1:numberColumns+shift))

[r,c]=size(newtest)

paddedtest=padarray(newtest,[numberAngles-r numberColumns-c],'post')

finaltest=transpose(paddedtest)