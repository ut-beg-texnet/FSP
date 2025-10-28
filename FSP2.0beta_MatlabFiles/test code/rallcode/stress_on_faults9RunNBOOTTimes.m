%%Coulomb failure criterion
function [resultMPa]=stress_on_faults9RunNBOOTTimes(S1,S2,S3,Pp,mu,a,b,c,faults,failureCriteria) 
%MohrFrac_arb(S1, S2, S3,Pp,coeff friction, a, b, c, ...
%               faults=2 column matrix with strike(0-360), dip(0-90) (use
%               right-hand rule)), e.g faults=[208,60;36,62;...]
%All angle information should be input as DEGREES, stresses as MPA (output
%will be the same)
% a=trend of S1...exception when S1 is vertical a=trend SHmax-90 degrees
% b=-plunge of S1 (plunge is angle down from horizontal)
% c=plunge of S2, 0 if S1 or S3 is vertical, 90 if S2 is vertical
% failureCriteria 1=PPf, 2=SCU (note this is not built!, 3=CFF

% format short

%converting fault and stress orientation info from degrees to radians 
str=(pi./180) .* faults(:,1);
dip=(pi./180) .* faults(:,2);
a=(pi./180) .* a;
b=(pi./180) .* b;
c=(pi./180) .* c;

%defines principal stress tensor
S=[S1 0 0;0 S2 0;0 0 S3];
%transformation from principal stress to geographic coordinates
R1=[cos(a)*cos(b) sin(a)*cos(b) -sin(b);
    cos(a)*sin(b)*sin(c)-sin(a)*cos(c) sin(a)*sin(b)*sin(c)+cos(a)*cos(c) cos(b)*sin(c);
    cos(a)*sin(b)*cos(c)+sin(a)*sin(c) sin(a)*sin(b)*cos(c)-cos(a)*sin(c) cos(b)*cos(c)];

Sg=R1'*S*R1;

%the following part of the code is dependent on fault orientation
%for loop makes calculations for each fault one at a time
fault_no=[1:1:length(faults(:,1))]';
Sn=zeros(length(fault_no),1);% preallocate
rake=Sn;% preallocate
tau=Sn;% preallocate
Sn_eff=Sn;% preallocate
SCU=Sn;
 CFF=Sn;% preallocate
PPF=Sn;% preallocate
 for k=1:length(fault_no) 

    % transformation from geographic to fault coordinate system
    R2=[cos(str(k)) sin(str(k)) 0;
         sin(str(k))*cos(dip(k)) -cos(str(k))*cos(dip(k)) -sin(dip(k));
         -sin(str(k))*sin(dip(k)) cos(str(k))*sin(dip(k)) -cos(dip(k))];

    Sf=R2*Sg*R2';
    %Solves for normal stress resolved on the fault surface
    Sn(k)=Sf(3,3);

    %Solve for rake of the slip vector
    if (Sf(3,2)>0)
        rake(k)=(atan(Sf(3,2)/(Sf(3,1))));
    elseif (Sf(3,2)<0)&&(Sf(3,1)>0)
        rake(k)=pi-atan(Sf(3,2)/(-Sf(3,1)));
    else
        rake(k)=atan((-Sf(3,2))/(-Sf(3,1)))-pi;
    end
    %Solve for shear stress resolved on the fault plane
    R3=[cos(rake(k)) sin(rake(k)) 0;
        -sin(rake(k)) cos(rake(k)) 0;
        0 0 1];

    Sr=R3*Sf*R3';

    tau(k,1)=Sr(3,1);
    %solve for Coulomb failure function (CFF) (proxy for likelihood of slip) of fault
    Sn_eff(k,1)=Sn(k)-Pp;
    CFF(k,1)=abs(tau(k)) - mu*Sn_eff(k);
% horizontalPorePressureToFailure=effective normal stress-tau/mu
    PPF(k,1)=(abs(tau(k))./mu)-Sn_eff(k)  ; 
%     rake(k)=rad2deg(rake(k));
 end
 
 switch failureCriteria
    case 1
        resultMPa=PPF;
     case 2
         disp('haven''t coded SCU yet!')
     case 3
            resultMPa=CFF;    
 end

%output data
%changing signs to make output easier to understand
% rake=(sign(tau).*rake');
% tau=abs(tau);

% FaultNo_Strike_Dip_ShearStress_EffNormalStress_Rake_CFF=[fault_no,faults,tau, Sn_eff, rake, CFF];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This part of the code is for visualizing the above results on 
%a 3D Mohr Diagram
if 0 % 1 is for debugging only, this plots mohr circle, 
%     don't run thousands of times! run no more than say 10 
S1_eff=S1-Pp;
S2_eff=S2-Pp;
S3_eff=S3-Pp;

%%%%Building the 3D diagram
format compact                    
format long e                    
eta=linspace(0,pi,50);
% close all
figure(9931)
% figure
set(gcf,'color','w')
subplot(2,1,2)
%Plotting S3-S2 cicle
x_32 = ((S2_eff-S3_eff)/2)*cos(eta)+((S2_eff-S3_eff)/2+S3_eff);                   % generate x-coordinate
y_32 = ((S2_eff-S3_eff)/2)*sin(eta);                   
plot(x_32,y_32);    

% if thisIteration==totalIterations  % only run these on last call of plotting, saves time
axis('equal');                    
% set(gca, 'Xlim', [-1, S1_eff+2], 'YLim', [0, (S1-S3)/2*1.1])
xlabel('Effective Normal Stress, MPa')
ylabel('Shear Stress, MPa')
title('3D Mohr Diagram')
% end
hold on

%Plotting S2-S1 cicle
x_21 = ((S1_eff-S2_eff)/2)*cos(eta)+((S1_eff-S2_eff)/2+S2_eff);                   % generate x-coordinate
y_21 = ((S1_eff-S2_eff)/2)*sin(eta);                   
plot(x_21,abs(y_21));   
% plot(x_21,y_21); 

%Plotting S3-S1 cicle
x_31 = ((S1_eff-S3_eff)/2)*cos(eta)+((S1_eff-S3_eff)/2+S3_eff);                   % generate x-coordinate
y_31 = ((S1_eff-S3_eff)/2)*sin(eta);                   
plot(x_31,y_31);                        

%plot failure line
y_failure=[0 (S1-S3)/2*1.1];
x_failure=y_failure/mu;
plot(x_failure, y_failure, 'r-')

%%%%Adding results from previous calculations


subplot(2,1,1)
scatter(faults(:,1), faults(:,2), 30, CFF, 'filled')
% title('Fault Orientations (labeled by input order & colorcoded by CFF)')

hold on
% for i=1:length(fault_no)
% str=num2str(i);
% % text(faults(i,1)+3, faults(i,2)+3,str,'fontsize',12)
% end
set(gca,'Xlim',[0 360], 'Ylim', [0 90],'fontsize',14)
% if thisIteration==totalIterations  % only run these on last call of plotting, saves time
title('Fault Orientations (colorcoded by CFF)')
xlabel('Strike')
ylabel('Dip')
colorbar
% end
%set(gca, 'Xlim', [0 180], 'Ylim', [-90 90])
subplot(2,1,2)
% 
scatter(Sn_eff, abs(tau), 30, CFF, 'filled')
% if thisIteration==totalIterations  % only run these on last call of plotting, saves time
title('3D Mohr Diagram (colorcoded by CFF)','fontsize',14)
    set(gca,'fontsize',14)
colorbar
% end
format short
end