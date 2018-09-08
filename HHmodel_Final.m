% iap*0.01 is an input current to model
num_peaks=zeros();
for iap=1:60
[time, amp]=HHmodel(0.01*iap);  %function called

%% calculation using Findpeaks Function
[peaks, locs]=findpeaks(amp);   %finding peaks
for n=1:length(peaks)
if (peaks(n)>=8) % minimum amplitude to be called as AP
ppeaks(n)=peaks(n);
else
ppeaks(n)=0;
end
end
num_peaks(iap)=sum(ppeaks>0); %Finding number of peaks for one iap
X(iap)=iap*0.01;
%%     Calculation using Zero crossing    
%      lamp=length(amp);
%      newamp=amp-mean(amp);
%      zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);                    % Returns Zero-Crossing Indices Of Argument Vector
%      zx = zci(newamp);
%      nzx(iap,1)=iap*0.01;
%      p=zx(length(zx)/2,1);
%      q=zx(length(zx)/2-1,1);
%      nzx(iap,2)=(length(zx)-6)/2;
%      if (nzx(iap,1)>0.1)
%          if(sum(amp>0)<6000)
%           nzx(iap,2)=0;    
%      end
%      end

end

%% Plotting plots
Y=num_peaks;
plot(X,Y);
xlabel('I_ext');
ylabel('Spikes per second')
hold on;
for l=2:length(Y) %finding I1, I2, I3.
if(Y(l)>0 && Y(l-1)==0)
 I1=X(l-1);
end
if(Y(l)>Y(l-1)+5)
 I2=X(l-1);
end
if(Y(l)<Y(l-1)-5)
 I3=X(l-2);
end
end
I1=[I1,I1];
y=[0,120];
plot(I1,y,'r');
I2=[I2,I2];
plot(I2,y,'b');
I3=[I3,I3];
plot(I3,y,'g');
text(0.5,100,'I1 => .0200');
text(0.5,95,'I2 => .0600');
text(0.5,90,'I3 => .48');
%% Function made of HHmodel
function [t,vhist]= HHmodel(ImpCur)
gkmax=.36;
gnamax=1.20;
vna=50; 
vk=-77; 
gl=0.003;
vl=-54.387; 
cm=.01; 
disp(ImpCur);
dt=0.01;
niter=100000;
t=(1:niter)*dt;
iapp=ImpCur*ones(1,niter);
%for i=1:100
 %   iapp(1,i)=ImpCur;
 %end;
v=-64.9964;
m=0.0530;
h=0.5960;
n=0.3177;

gnahist=zeros(1,niter);
gkhist=zeros(1,niter);
vhist=zeros(1,niter);
mhist=zeros(1,niter);
hhist=zeros(1,niter);
nhist=zeros(1,niter);


for iter=1:niter
  gna=gnamax*m^3*h; 
  gk=gkmax*n^4; 
  gtot=gna+gk+gl;
  vinf = ((gna*vna+gk*vk+gl*vl)+ iapp(iter))/gtot;
  tauv = cm/gtot;
  v=vinf+(v-vinf)*exp(-dt/tauv);
  alpham = 0.1*(v+40)/(1-exp(-(v+40)/10));
  betam = 4*exp(-0.0556*(v+65));
  alphan = 0.01*(v+55)/(1-exp(-(v+55)/10));
  betan = 0.125*exp(-(v+65)/80);
  alphah = 0.07*exp(-0.05*(v+65));
  betah = 1/(1+exp(-0.1*(v+35)));
  taum = 1/(alpham+betam);
  tauh = 1/(alphah+betah);
  taun = 1/(alphan+betan);
  minf = alpham*taum;
  hinf = alphah*tauh;
  ninf = alphan*taun;
  m=minf+(m-minf)*exp(-dt/taum);
  h=hinf+(h-hinf)*exp(-dt/tauh);
  n=ninf+(n-ninf)*exp(-dt/taun);
  vhist(iter)=v; mhist(iter)=m; hhist(iter)=h; nhist(iter)=n;
end

% figure(1)
% %subplot(2,1,1)
% plot(t,vhist)
% title('voltage vs time')

% 
% figure(2)
% %subplot(2,1,2)
% plot(t,mhist,'y-', t,hhist,'g.',t,nhist,'b-')
% legend('m','h','n')
% 
% figure(3)
% gna=gnamax*(mhist.^3).*hhist; 
%   gk=gkmax*nhist.^4;
%   clf
%   plot(t,gna,'r');
%   hold on
%   plot(t,gk,'b');
%   legend('gna','gk')
%   hold off
end