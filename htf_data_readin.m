%sam htf salt data input
clear variables
%xlsread
SolarSaltData = xlsread('C:\Users\sethb\Desktop\Energy Storage\Model Data\sam-htf-property-tables.xlsx',1);
save('SolarSaltData')
clear variables
HitecXLData = xlsread('C:\Users\sethb\Desktop\Energy Storage\Model Data\sam-htf-property-tables.xlsx',3);
save('HitecXLData')
clear variables
HitecData = xlsread('C:\Users\sethb\Desktop\Energy Storage\Model Data\sam-htf-property-tables.xlsx',5);
save('HitecData')
%define data
% T_SS = SolarSaltData(:,1);
% h_SS = SolarSaltData(:,7);
% T_HT = HitecData(:,1);
% h_HT = HitecData(:,7);
% T_HTXL = HitecXLData(:,1);
% h_HTXL = HitecXLData(:,7);