

% [data_] = Spire_Pseudorange_Processing(2018,[9],[23:24],83,2,1,0);

%%% PLOT solution 
mon = data_.mon; dday = data_.dday; hh=data_.hh; 
mm = data_.mm; ss = data_.ss; 
dsp3_p = data_.dsp3_p;
dsp3_v = data_.dsp3_v;
b = data_.dsp3_b/c*1e9;

yyyy = data_.yyyy; 
dpost_pos = data_.dpost_pos; dpost_vel = data_.dpost_vel;
dpre_pos = data_.dpre_pos;

dates = datetime(yyyy,mon,dday,hh,mm,ss); 

[dsp3_p_1,ind1]=rmoutliers(dsp3_p(1,:));
dates_1 = dates(~ind1);
[dsp3_p_2,ind2]=rmoutliers(dsp3_p(2,:));
dates_2 = dates(~ind2);
[dsp3_p_3,ind3]=rmoutliers(dsp3_p(3,:));
dates_3 = dates(~ind3);

[b_3,ind4]=rmoutliers(b);
dates_43 = dates(~ind4);
%%
figure; 
subplot(3,1,1);
plot(dates_1,dsp3_p_1,'.'); grid on; ylabel('X (m)'); title('Difference with given POD position- Sat ID=083')
set(gca,'FontSize',14)
subplot(3,1,2);
plot(dates_2,dsp3_p_2,'.'); grid on; ylabel('Y (m)'); 
set(gca,'FontSize',14)
subplot(3,1,3);
plot(dates_3,dsp3_p_3,'.'); grid on; ylabel('Z (m)'); 
set(gca,'FontSize',14)


figure; 
plot(dates,dpost_pos,'.'); grid on; ylabel('(m)'); title('Postfit Residuals - position, Sat ID = 083')
set(gca,'FontSize',14)
ylim([-40 40])
figure; 
plot(dates,dpre_pos,'.'); grid on; ylabel('(m)'); title('Prefit Residuals - position, Sat ID = 083')
set(gca,'FontSize',14)
ylim([-100 100])

%%
figure;
plot(dates_41, b_1)
hold on
plot(dates_42, b_2)
plot(dates_43, b_3)
ylabel('Clock error (ns)')
title('Difference with given POD clock bias')
%% Plot the time difference

