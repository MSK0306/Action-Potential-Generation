%Muhammed Saadeddin Kocak 2232346
%EE416 Term Project
t=linspace(-99,30000,30100);%usec
Im(1:100)=0;%uA/sec^2
impulseduration=200;%usec
impulseamplitude=53;%uA/cm^2
numberofconsecutivestimuli=2;
timedelaybetweenstimuli=5000;%usec
Cm=1;%uF/cm^2
Vm(1:length(t))=-90;%mV
vm(1:length(t))=0;%Vmembrane-Vrest
gNarest=120;%mS/cm^2
gKrest=36;%mS/cm^2
gLrest=0.3;%mS/cm^2
gK(1:length(t))=gKrest;
gNa(1:length(t))=gNarest;
gL(1:length(t))=gLrest;

VNa=25;%mV
VK=-102;%mV
VL=-79.387;%mV
for i=100:length(t)
   if t(i)>-1 && t(i)<impulseduration
       Im(i)=impulseamplitude;%uA/cm^2
   else
       Im(i)=0;
   end
   %For the second stimuli
   if numberofconsecutivestimuli~=1
   if t(i)>timedelaybetweenstimuli && t(i)<timedelaybetweenstimuli+impulseduration
       Im(i)=impulseamplitude;%uA/cm^2
   end
   end
   %
   an(i)=((0.01*(10-vm(i)))/(exp((10-vm(i))/10)-1))/1000;%1/usec
   bn(i)=(0.125*exp((-vm(i))/80))/1000;%1/usec
   am(i)=((0.1*(25-vm(i)))/(exp(0.1*(25-vm(i)))-1))/1000;%1/usec
   bm(i)=(4*exp((-vm(i))/18))/1000;%1/usec
   ah(i)=(0.07*exp((-vm(i))/20))/1000;%1/usec
   bh(i)=(1/(exp((30-vm(i))/10)+1))/1000;%1/usec
   n(100)=an(100)/(an(100)+bn(100));
   m(100)=am(100)/(am(100)+bm(100));
   h(100)=ah(100)/(ah(100)+bh(100));
   gK(i)=gKrest*(n(i).^4);%mS/cm^2
   gNa(i)=gNarest*(m(i).^3)*h(i);%mS/cm^2
   IK(i)=gK(i)*(Vm(i)-VK);%uA/cm^2
   INa(i)=gNa(i)*(Vm(i)-VNa);%uA/cm^2
   IL(i)=gLrest*(Vm(i)-VL);%uA/cm^2
   Ic(i)=Im(i)-IK(i)-INa(i)-IL(i);
   %Estimating deltan deltam deltah
   deltan(i)=(an(i)*(1-n(i)))-(bn(i)*n(i));
   deltam(i)=(am(i)*(1-m(i)))-(bm(i)*m(i));
   deltah(i)=(ah(i)*(1-h(i)))-(bh(i)*h(i));
   deltaVm(i)=(Im(i)-IK(i)-INa(i)-IL(i))/(Cm*1000);%mV
   %Calculating next variables
   n(i+1)=n(i)+deltan(i);
   m(i+1)=m(i)+deltam(i);
   h(i+1)=h(i)+deltah(i);
   Vm(i+1)=Vm(i)+deltaVm(i);
   vm(i+1)=vm(i)+deltaVm(i);
end
figure
plot(t,Im(1:length(t)));%Both for applied and membrane current
title('Applied Input Current (Total Mebrane Current) vs Time');
xlabel('Time (usec)');
ylabel('Current(uA/cm^2)');
figure
hold on
plot(t,IK(1:length(t)));
plot(t,INa(1:length(t)));
plot(t,IL(1:length(t)));
plot(t,Ic(1:length(t)));
legend({'IK','INa','IL','Ic'});
title('Ionic Currrents and Capacitive Current vs Time');
xlabel('Time (usec)');
ylabel('Current(uA/cm^2)');
hold off
figure
hold on
plot(t,gK(1:length(t)));
plot(t,gNa(1:length(t)));
plot(t,gL(1:length(t)));
legend({'gK','gNa','gL'});
title('Channel Conductances vs Time');
xlabel('Time (usec)');
ylabel('Conductance(mS/cm^2)');
hold off
figure
plot(t,Vm(1:length(t)));
title('Membrane Voltage vs Time');
xlabel('Time (usec)');
ylabel('Voltage(mV)')
%Extra
figure
hold on
plot(t,n(1:length(t)));
plot(t,m(1:length(t)));
plot(t,h(1:length(t)));
legend({'n','m','h'});
title('Activation Parameters');
xlabel('Time (usec)');
ylabel('Gating')
hold off
