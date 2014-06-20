%clear all


filename='/home/wesley/ncfiles/smallcape_force_0001.nc';
ncid=netcdf.open(filename,'NC_NOWRITE');

x = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'x'));
y = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'y'));
lon = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lon'));
lat = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lat'));
lonc = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'lonc'));
latc = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'latc'));
ua = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'ua'),[0 0],[end 1])';
va = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'va'),[1 1],[end 1])';
time = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time'));
netcdf.close(ncid);

ua

time=time+datenum('17-Nov-1858');  %On this data, the MJD = 0

order={'M2','S2','N2','K2','K1','O1','P1','Q1'}

coef = ut_solv(time, double(ua(:,0)), double(va(:,0)), uvllnode(0), 'auto','NoTrend','Rmin','OLS','NoDiagn','LinCI');

save ../speedmatlabcoef.mat coef


coef = ut_solv(time, double(ua(:,0)), [], uvllnode(0), 'auto','NoTrend','Rmin','OLS','NoDiagn','LinCI');

save ../matlabcoef.mat coef
