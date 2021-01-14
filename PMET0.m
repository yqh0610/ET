%%%%ET0计算函数
%%根据彭曼(PM)公式计算日尺度参考作物蒸散发ET0 (FAO56)

% % % 示例数据 % % % 
%%%和田1999
% ET0=PMET0(1,-5.6,-0.5,-8.9,86,1.5,5.2,1375,37.133333);
% ET0=PMET0(2,-7.8,-2.7,-10.8,92,0.5,4,1375,37.133333);

%%%需要的参数依次是日序数，日平均温度(摄氏度)，日最高温度(摄氏度)，日最低温度(摄氏度)，平均相对湿度（%），风速(m/s)，日照时长(h)，站点海拔高度（m），站点纬度（度）
function ET0=PMET0(yearday,Tmean,Tmax,Tmin,RHmean,u,sunshineHour,DEM,latitude)	
P=101.3*((293-0.0065*DEM)/293)^5.26;%%大气压P （kPa) %%有实测还是应该直接用实测
gama=0.665/1000*P;%%湿度计常数γ(kPa℃-1)
T=(Tmin+Tmax)/2;
es=0.6108*exp(17.269*Tmax/(237.3+Tmax))/2+0.6108*exp(17.269*Tmin/(237.3+Tmin))/2;%%平均饱和水汽压es(kPa)
ea=RHmean/100*(es);%%实际水汽压ea(kPa)  实际水汽压有实测的用实测，有好几种估算算法，这是由平均相对湿度计算来的
VPD=es-ea;%饱和水汽压差（kPa）
deta = (4098*0.6108*exp((17.27*Tmean)/(Tmean+237.3)))/(Tmean+237.3)^2;%饱和水汽压曲线斜率Δ(kPa℃-1)
dr=1+0.033*cos(2*pi/365*yearday);%日地间相对距离的倒数dr
delta=0.409*sin(2*pi/365*yearday-1.39);%太阳磁偏角δ
fi=pi/180*latitude;%w纬度转换为弧度 弧度
omegas=acos(-tan(fi)*tan(delta));%%这个有两种计算方法 太阳时角ws
N=24/pi*omegas;%最大日照时数N
Ra=24*60/pi*0.082*dr*(omegas*sin(fi)*sin(delta)+cos(fi)*cos(delta)*sin(omegas));%天顶辐射Ra(MJm-2day-1)
Rs=(0.25+0.5*sunshineHour/N)*Ra;%太阳辐射Rs(MJm-2day-1)
Rso=(0.75+(2*10^(-5))*DEM)*Ra;%%晴空太阳辐射Rs0(MJm-2day-1)
Rns=(1-0.23)*Rs;%净短波辐射Rns(MJm-2day-1)
Rnl=4.903*10^(-9)*((Tmax+273.16)^4+(Tmin+273.16)^4)/2*(0.34-0.14*ea^0.5)*(1.35*Rs/Rso-0.35);%%净长波辐射Rnl(MJm-2day-1)
u2 = u*4.87/(log(67.8*10-5.42));%%风速（m/s）转换到两米 一般国家站风速测点10-12米，此处10米可根据实际情况更改
G = 0;%土壤热通量G（MJm-2day-1）日尺度上忽略，记为0
Rn=Rns-Rnl;%净辐射Rn(MJ m-2day-1)
ET0 = (0.408*deta*(Rn-G)+gama*(900/(Tmean+273.16))*u2*VPD)/(deta+gama+gama*0.34*u2);%%mm
end




