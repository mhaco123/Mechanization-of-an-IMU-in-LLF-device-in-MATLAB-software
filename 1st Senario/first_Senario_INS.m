% In the name of GOD 
% 990921_M.h.Aghanoori_982015004

clc,clear,close all
%% Earth components
a=6378137;
e=0.08181919;
we=7.292115*10^-5;
% e=0;we=0;
a1=9.7803267714;
a2=0.0052790414;
a3=0.0000232718;
a4=-0.0000030876910891;
a5=0.0000000043977311;
a6=0.0000000000007211;
%% normal test

%% 1st Senario جسم ساکن در ارتفاع
F = readtable('Accelerometer.csv','PreserveVariableNames',true);  % skips the first three rows of data
W = readtable('Orientation.csv','PreserveVariableNames',true);  % skips the first three rows of data
t(1)=0;d_t=0.01;time=100;
Ve(1)=0;Vn(1)=0;Vu(1)=0;
d_Ve(1)=0;d_Vn(1)=0;d_Vu(1)=0;
r(1)=0;p(1)=0;y(1)=0;
d_r(1)=0;d_p(1)=0;d_A(1)=0;
lat(1)=31;lon(1)=52;h(1)=1450;
d_lat(1)=0;d_lon(1)=0;d_h(1)=0;
N(1)=0;E(1)=0;U(1)=h(1);

% % Fb=zeros(3,time*100);
% Fb(:,1)=T(:,2);
Fx=(F.x)';Wx=(W.pitch)';
Fy=(F.y)';Wy=(W.roll)';
Fz=(F.z)';Wz=(W.yaw)';
Fb=[Fx;Fy;Fz];
beta_Fx=0;beta_Fy=0;beta_Fz=0;
var_Fx=0;var_Fy=0;var_Fz=0;
d_Fx(1)=0;d_Fy(1)=0;d_Fz(1)=0.00;
W=[Wx;Wy;Wz];
beta_Wx=0;beta_Wy=0;beta_Wz=0;
var_Wx=0;var_Wy=0;var_Wz=0;
d_Wx(1)=0;d_Wy(1)=0;d_Wz(1)=0;
White=0.1;

%% 3rd step of Mechnization calculate & update of Rotation matrix
Vl(:,1)=[Ve(1);Vn(1);Vu(1)];
R_b_l1=[cosd(y)*cosd(r)-sind(y)*sind(p)*sind(r),...
    -sind(y)*cosd(p),cosd(y)*sind(r)+sind(y)*sind(p)*cos(r);...
    sind(y)*cosd(r)+cosd(y)*sind(p)*sind(r),...
    cosd(y)*cosd(p),sind(y)*sind(r)-cosd(y)*sind(p)*cosd(r);...
    -cosd(p)*sind(r),sind(p),cosd(p)*cosd(r)];

q4=0.5*sqrt(1+R_b_l1(1,1)+R_b_l1(2,2)+R_b_l1(3,3));
q1=0.25*(R_b_l1(3,2)-R_b_l1(2,3))/q4;
q2=0.25*(R_b_l1(1,3)-R_b_l1(3,1))/q4;
q3=0.25*(R_b_l1(2,1)-R_b_l1(1,2))/q4;
q(:,1)=[q1;q2;q3;q4];
%% UPDATING
for i=1:length(F.x)


d_Fx(i+1)=d_Fx(i)+0.5*(-beta_Fx*d_Fx(i)+(sqrt(2*beta_Fx*var_Fx^2))*White)*d_t;
d_Fy(i+1)=d_Fy(i)+0.5*(-beta_Fy*d_Fy(i)+(sqrt(2*beta_Fy*var_Fy^2))*White)*d_t;
d_Fz(i+1)=d_Fz(i)+0.5*(-beta_Fz*d_Fz(i)+(sqrt(2*beta_Fz*var_Fz^2))*White)*d_t;

d_Wx(i+1)=d_Wx(i)+0.5*(-beta_Wx*d_Wx(i)+(sqrt(2*beta_Wx*var_Wx^2))*White)*d_t;
d_Wy(i+1)=d_Wy(i)+0.5*(-beta_Wy*d_Wy(i)+(sqrt(2*beta_Wy*var_Wy^2))*White)*d_t;
d_Wz(i+1)=d_Wz(i)+0.5*(-beta_Wz*d_Wz(i)+(sqrt(2*beta_Wz*var_Wz^2))*White)*d_t;

RN=a/((1-e^2*sind(lat(i))^2)^(1/2));
RM=(a*(1-e^2))/((1-e^2*sind(lat(i))^2)^(3/2));
Sb=[0   Wz(i) -Wy(i) Wx(i);
    -Wz(i) 0 Wx(i) Wy(i)
    Wy(i) -Wx(i) 0 Wz(i)
    -Wx(i) -Wy(i) -Wz(i) 0]*d_t;

teta=sqrt((Wx(i)*d_t)^2+(Wy(i)*d_t)^2+(Wz(i)*d_t)^2);

q(:,i+1)=q(:,i)+(0.5*(2*(cos(teta/2)-1)*eye(4)+(2/teta)*sin(teta/2)*Sb))*q(:,i);
q1=q(1,i+1);
q2=q(2,i+1);
q3=q(3,i+1);
q4=q(4,i+1);

R_b_l=[q1^2-q2^2-q3^2+q4^2 2*(q1*q2-q3*q4) 2*(q1*q3+q2*q4)
     2*(q1*q2+q3*q4) -q1^2+q2^2-q3^2+q4^2 2*(q2*q3-q1*q4)
     2*(q1*q3-q2*q4) 2*(q2*q3+q1*q4) -q1^2-q2^2+q3^2+q4^2];
%
d_p(i+1)=d_p(i)+0.5*(((d_Vn(i))/(RM+h(i)))+R_b_l(1,1)*(d_Wx(i)+d_Wx(i+1))...
    +R_b_l(1,2)*(d_Wy(i)+d_Wy(i+1)) +R_b_l(1,3)*(d_Wz(i)+d_Wz(i+1)))*d_t;
p(i+1)=round(asind(R_b_l(3,2)),1)-d_p(i+1);

d_r(i+1)=d_r(i)+0.5*((-(d_Ve(i))/(RN+h(i)))+R_b_l(2,1)*(d_Wx(i)+d_Wx(i+1))...
    +R_b_l(2,2)*(d_Wy(i)+d_Wy(i+1))+R_b_l(2,3)*(d_Wz(i)+d_Wz(i+1)))*d_t;
r(i+1)=round(-atand(R_b_l(3,1)/R_b_l(3,3)),1)-d_r(i+1);

d_A(i+1)=d_A(i)+0.5*(((-tand(lat(i))*(d_Ve(i)))/(RN+h(i)))+...
    R_b_l(3,1)*(d_Wx(i)+d_Wx(i+1))+R_b_l(3,2)*(d_Wy(i)+d_Wy(i+1))+...
    R_b_l(3,3)*(d_Wz(i)+d_Wz(i+1)))*d_t;
y(i+1)=round(-atand(R_b_l(1,2)/R_b_l(2,2)),1)-d_A(i+1);
%
omelel=[0 (-Ve(i)*tand(lat(i)))/(RN+h(i)) Ve(i)/(RN+h(i))
        (Ve(i)*tand(lat(i)))/(RN+h(i)) 0 Vn(i)/(RM+h(i))
        -Ve(i)/(RN+h(i)) -Vn(i)/(RM+h(i)) 0];

omeeil=[0 -we*sind(lat(i)) we*cosd(lat(i))
        we*sind(lat(i)) 0 0
        -we*cosd(lat(i)) 0 0];
    
g=a1*(1+a2*sind(lat(i))^2+a3*sind(lat(i))^4)+((a4+a5*sind(lat(i))^2)*h(i))+a6*h(i)^2;
gl=[0;0;-g];

dvbl=R_b_l*Fb(:,i)*d_t;

dvl(:,i+1)=dvbl-(2*omeeil+omelel)*Vl(:,i)*d_t+gl*d_t;
dvl(3,i+1)=0;
Vl(:,i+1)=Vl(:,i)+0.5*(dvl(:,i)+dvl(:,i+1));

Fe(i)=R_b_l(1,1)*Fx(i)+R_b_l(1,2)*Fy(i)+...
    R_b_l(1,3)*Fz(i);
Fn(i)=R_b_l(2,1)*Fx(i)+R_b_l(2,2)*Fy(i)+...
    R_b_l(2,3)*Fz(i);
Fu(i)=R_b_l(3,1)*Fx(i)+R_b_l(3,2)*Fy(i)+...
    R_b_l(3,3)*Fz(i);

d_Ve(i+1)=d_Ve(i)+0.5*(Fu(i)*(d_r(i)+d_r(i+1))-Fn(i)*(d_A(i)+d_A(i+1)))+...
    R_b_l(1,1)*d_Fx(i)+R_b_l(1,2)*d_Fy(i)+R_b_l(1,3)*d_Fz(i)*d_t;
Ve(i+1)=Vl(1,i+1)-d_Ve(i+1);

d_Vn(i+1)=d_Vn(i)+0.5*(-Fu(i)*(d_p(i)+d_p(i+1))+Fe(i)*(d_A(i)+d_A(i+1))+...
    R_b_l(2,1)*d_Fx(i)+R_b_l(2,2)*d_Fy(i)+R_b_l(2,3)*d_Fz(i))*d_t;
Vn(i+1)=Vl(2,i+1)-d_Vn(i+1);

d_Vu(i+1)=d_Vu(i)+0.5*(Fn(i)*(d_p(i)+d_p(i+1))-Fe(i)*(d_r(i)+d_r(i+1))+...
    R_b_l(3,1)*d_Fx(i)+R_b_l(3,2)*d_Fy(i)+R_b_l(3,3)*d_Fz(i))*d_t;
Vu(i+1)=Vl(3,i+1)-d_Vu(i+1);

d_lat(i+1)=d_lat(i)+0.5*((d_Vn(i)+Vn(i+1))/(RN+h(i)))*d_t;
lat(i+1)=(lat(i)+0.5*((Vn(i)+Vn(i+1))/(RN+h(i)))*d_t)-d_lat(i+1);

d_lon(i+1)=d_lon(i)+0.5*((d_Ve(i)+d_Ve(i+1))/((RN+h(i))*cos(d_lat(i))))*d_t;
lon(i+1)=(lon(i)+0.5*((Ve(i)+Ve(i+1))/((RN+h(i))*cos(lat(i))))*d_t)-d_lon(i+1);


d_h(i+1)=d_h(i)+0.5*(d_Vu(i)+d_Vu(i+1))*d_t;
h(i+1)=(h(i)+0.5*(Vu(i)+Vu(i+1))*d_t)-d_h(i+1);

E(i+1)=E(i)+0.5*(Ve(i)+Ve(i+1))*d_t;
N(i+1)=N(i)+0.5*(Vn(i)+Vn(i+1))*d_t;
U(i+1)=U(i)+0.5*(Vu(i)+Vu(i+1))*d_t;
t(i+1)=t(i)+d_t;
if h(i+1)<=0
    break
end
end
%%
% H = [eye(6) zeros(6,9)];
save('INS_LLF.mat','lat','lon','h','E','N','U','Ve','Vn','Vu','r','p','y','d_lat'...
,'d_lon','d_h','d_Fx','d_Fy','d_Fz','d_Wx','d_Wy','d_Wz','d_p','d_r','d_A','d_Ve','d_Vn','d_Vu','W'...
,'Fe','Fn','Fu','R_b_l');
subplot(3,4,1)
plot(t,lat,'g','linewidth',1)
xlabel('Time(s)');ylabel('Latitude');
grid
%%
subplot(3,4,5)
plot(t,lon,'g','linewidth',1)
xlabel('Time(s)');ylabel('longtitude');
grid
%%
subplot(3,4,9)
plot(t,h,'g','linewidth',1)
xlabel('Time(s)');ylabel('h');
grid
%%
subplot(3,4,2)
plot(t,E,'r','linewidth',1)
xlabel('Time(s)');ylabel('E');
grid
%%
subplot(3,4,6)
plot(t,N,'r','linewidth',1)
xlabel('Time(s)');ylabel('N');
grid
%%
subplot(3,4,10)
plot(t,U,'r','linewidth',1)
xlabel('Time(s)');ylabel('U');
grid
%%
subplot(3,4,3)
plot(t,Ve,'c','linewidth',1)
xlabel('Time(s)');ylabel('Ve');
grid
%%
subplot(3,4,7)
plot(t,Vn,'c','linewidth',1)
xlabel('Time(s)');ylabel('Vn');
grid
%%
subplot(3,4,11)
plot(t,Vu,'c','linewidth',1)
xlabel('Time(s)');ylabel('Vu');
grid
%%
subplot(3,4,4)
plot(t,r,'linewidth',1)
xlabel('Time(s)');ylabel('roll');
grid
%%
subplot(3,4,8)
plot(t,p,'linewidth',1)
xlabel('Time(s)');ylabel('pich');
grid
%%
subplot(3,4,12)
plot(t,y,'linewidth',1)
xlabel('Time(s)');ylabel('yaw');
grid