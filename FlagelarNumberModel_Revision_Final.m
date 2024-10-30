%% DFFM dataset
%load('2024-02-20_DataRFT_ForSave.mat', 'DataDFFM');
DataDFFM = readtable('2024-10-01_DataDFFM_Final.txt');
%% Flagellar number model 

% Params Rod
lmbd  = 2.3; % um % Turner 2000
R = 0.2; %um % Turner 2000 %0.20 for indiv
r = 0.01; %d=0.020; %um  %https://www.sciencedirect.com/science/article/pii/S0022283685702823?via%3Dihub

dlmbd = 0;%0.1; % variation is [-delta/2 +delta/2]
dR = 0.0;%0.02;
dr = 0;

% eta = 1; % mPa.s 

% Params Rod
Lbody = 2.5;
dbody = 0.8;

dLb = 0.8;
ddb = 0.2;

nst = 5000;


% Setup
dataset = struct();
dataset.flagellaParams = {lmbd,R,r,dlmbd,dR,dr};
dataset.bodyParams = {Lbody,dbody,dLb,ddb};
dataset.model = cell(9,1);

dataset.Nf = (1:8);
dataset.Lf = 5.3983 + 3.6564 * log10(dataset.Nf);

dataset.model{1} = 'Lighthill';
dataset.model{2} = 'GrayHancock';
dataset.model{3} = 'RC';
for i=4:6
    dataset.model{i} = 'Lighthill-Constantgk';
    dataset.gk{i} = 0.5+0.1*(i-4);
end
for i=7:9
    dataset.model{i} = 'GrayHancock-Constantgk';
    dataset.gk{i} = 0.5+0.1*(i-7);
end

dataset.color = {[193 39 45], [46 49 146], [0 0 0], [237 28 36], [241 90 36], [247 147 30], [0 113 188], [41 171 226], [0 169 157]};
dataset.color =cellfun(@(x) x/255, dataset.color,'UniformOutput',false);
%% Compute
for i=1:3
    [Data,Fct] = generateModelFlagellarPropulsion(dataset.model{i},dataset.Lf,dataset.Nf,lmbd,R,r,Lbody,dbody,dlmbd,dR,dr,dLb,ddb,'Nstat',nst);
    dataset.data{i} = Data;
    dataset.fct{i} = Fct;
end
for i=4:9
    [Data,Fct] = generateModelFlagellarPropulsion(dataset.model{i},dataset.Lf,dataset.Nf,lmbd,R,r,Lbody,dbody,dlmbd,dR,dr,dLb,ddb,'Nstat',nst,'gkcst',dataset.gk{i});
    dataset.data{i} = Data;
    dataset.fct{i} = Fct;
end


%% errorbar plots for resp to ref

figure
subplot(2,3,1)
plot(DataDFFM.nflag,DataDFFM.v0,'ok','MarkerFaceColor','k','MarkerSize',10)
hold on
cellfun( @(x,y) errorbar(dataset.Nf, x.mean.U,x.std.U,'Color',y), dataset.data, dataset.color)
%legend({'Data',dataset.model{:}})
title('v0')
set(gca,'FontSize',14)

subplot(2,3,4)
plot(DataDFFM.nflag,DataDFFM.v0/mean(DataDFFM.v0),'ok','MarkerFaceColor','k','MarkerSize',10)
hold on
cellfun( @(x,y) errorbar(dataset.Nf, x.mean.U./mean(x.mean.U), x.std.U./mean(x.mean.U),'Color',y), dataset.data, dataset.color)
%legend([{'Data'},dataset.model{:}])
title('normalized v0')
set(gca,'FontSize',14)

%figure
subplot(2,3,2)
plot(DataDFFM.nflag,DataDFFM.v0wf,'ok','MarkerFaceColor','k','MarkerSize',10)
hold on
cellfun( @(x,y) errorbar(dataset.Nf, x.mean.U_wf,x.std.U_wf,'Color',y), dataset.data, dataset.color)
%legend({'Data',dataset.model{:}})
title('v0/wf')
set(gca,'FontSize',14)

subplot(2,3,5)
plot(DataDFFM.nflag,DataDFFM.v0wf/mean(DataDFFM.v0wf),'ok','MarkerFaceColor','k','MarkerSize',10)
hold on
cellfun( @(x,y) errorbar(dataset.Nf, x.mean.U_wf./mean(x.mean.U_wf), x.std.U_wf./mean(x.mean.U_wf),'Color',y), dataset.data, dataset.color)
%legend([{'Data'},dataset.model{:}])
title('normalized v0/wf')
set(gca,'FontSize',14)

%figure
subplot(2,3,3)
plot(DataDFFM.nflag,DataDFFM.wbwf,'ok','MarkerFaceColor','k','MarkerSize',10)
hold on
cellfun( @(x,y) errorbar(dataset.Nf, x.mean.wb_wf,x.std.wb_wf,'Color',y), dataset.data, dataset.color)
%legend({'Data',dataset.model})
title('wb/wf')
set(gca,'FontSize',14)

subplot(2,3,6)
plot(DataDFFM.nflag,DataDFFM.wbwf/mean(DataDFFM.wbwf),'ok','MarkerFaceColor','k','MarkerSize',10)
hold on
cellfun( @(x,y) errorbar(dataset.Nf, x.mean.wb_wf./mean(x.mean.wb_wf), x.std.wb_wf./mean(x.mean.wb_wf),'Color',y), dataset.data, dataset.color)
legend([{'Data'},dataset.model{:}])
title('normalized wb/wf')
set(gca,'FontSize',14)
%%
figure
hold on
cellfun( @(x,y) errorbar(dataset.Nf, 180/pi*x.mean.theta,180/pi*x.std.theta,'Color',y), dataset.data, dataset.color)
legend([dataset.model])
title('theta')
set(gca,'FontSize',14)

%% errorbar plots for Extended Data Fig 5

rng = [1,2,5,6,8,9];

figure
subplot(2,3,1)
plot(DataDFFM.nflag,DataDFFM.v0,'ok','MarkerFaceColor','k','MarkerSize',10)
hold on
cellfun( @(x,y) errorbar(dataset.Nf, x.mean.U,x.std.U,'Color',y), dataset.data(rng), dataset.color(rng))
%legend({'Data',dataset.model{:}})
title('v0')
set(gca,'FontSize',14,'YLim',[0 40],'XLim',[0 8])

subplot(2,3,4)
plot(DataDFFM.nflag,DataDFFM.v0/mean(DataDFFM.v0),'ok','MarkerFaceColor','k','MarkerSize',10)
hold on
cellfun( @(x,y) errorbar(dataset.Nf, x.mean.U./mean(x.mean.U), x.std.U./mean(x.mean.U),'Color',y), dataset.data(rng), dataset.color(rng))
%legend([{'Data'},dataset.model{:}])
title('normalized v0')
set(gca,'FontSize',14,'YLim',[0.75 1.2],'XLim',[0 8])

%figure
subplot(2,3,2)
plot(DataDFFM.nflag,DataDFFM.v0wf,'ok','MarkerFaceColor','k','MarkerSize',10)
hold on
cellfun( @(x,y) errorbar(dataset.Nf, x.mean.U_wf,x.std.U_wf,'Color',y), dataset.data(rng), dataset.color(rng))
%legend({'Data',dataset.model{:}})
title('v0/wf')
set(gca,'FontSize',14,'YLim',[0 0.25],'XLim',[0 8])

subplot(2,3,5)
plot(DataDFFM.nflag,DataDFFM.v0wf/mean(DataDFFM.v0wf),'ok','MarkerFaceColor','k','MarkerSize',10)
hold on
cellfun( @(x,y) errorbar(dataset.Nf, x.mean.U_wf./mean(x.mean.U_wf), x.std.U_wf./mean(x.mean.U_wf),'Color',y), dataset.data(rng), dataset.color(rng))
%legend([{'Data'},dataset.model{:}])
title('normalized v0/wf')
set(gca,'FontSize',14,'YLim',[0.75 1.2],'XLim',[0 8])

%figure
subplot(2,3,3)
plot(DataDFFM.nflag,DataDFFM.wbwf,'ok','MarkerFaceColor','k','MarkerSize',10)
hold on
cellfun( @(x,y) errorbar(dataset.Nf, x.mean.wb_wf,x.std.wb_wf,'Color',y), dataset.data(rng), dataset.color(rng))
%legend({'Data',dataset.model})
title('wb/wf')
set(gca,'FontSize',14,'YLim',[0 0.3],'XLim',[0 8])

subplot(2,3,6)
plot(DataDFFM.nflag,DataDFFM.wbwf/mean(DataDFFM.wbwf),'ok','MarkerFaceColor','k','MarkerSize',10)
hold on
cellfun( @(x,y) errorbar(dataset.Nf, x.mean.wb_wf./mean(x.mean.wb_wf), x.std.wb_wf./mean(x.mean.wb_wf),'Color',y), dataset.data(rng), dataset.color(rng))
legend([{'Data'},dataset.model{rng}])
title('normalized wb/wf')
set(gca,'FontSize',14,'YLim',[0.4 1.5],'XLim',[0 8])

%% theta for ED Fig 5
figure
hold on
cellfun( @(x,y) errorbar(dataset.Nf, 180/pi*x.mean.theta,180/pi*x.std.theta,'Color',y), dataset.data(rng), dataset.color(rng))
legend([dataset.model(rng)])
title('theta')
set(gca,'FontSize',14)
%% motor torque

figure
hold on
cellfun( @(x,y) errorbar(dataset.Nf, x.mean.GammaM, x.std.GammaM,'Color',y), dataset.data(rng), dataset.color(rng))
legend([dataset.model(rng)])
title('Motor torque')
set(gca,'FontSize',14)

%% Constant theta 

theta0 = [0,2.5,5,10,20,30];
for i=10:15
    dataset.model{i} = 'Lighthill-Constantgk';
    dataset.gk{i} = 0.7;
    
    [Data,Fct] = generateModelFlagellarPropulsion(dataset.model{i},dataset.Lf,dataset.Nf,lmbd,R,r,Lbody,dbody,dlmbd,dR,dr,dLb,ddb,'Nstat',nst,'gkcst',dataset.gk{i},'theta0',theta0(i-9)*pi/180);
    dataset.data{i} = Data;
    dataset.fct{i} = Fct;
end

figure
plot(DataDFFM.nflag,DataDFFM.v0,'ok','MarkerFaceColor','k','MarkerSize',10)
hold on
cellfun( @(x,y) errorbar(dataset.Nf, x.mean.U, x.std.U), dataset.data([6,10:15]))
legend(["Data";"Lh gk-0.7";num2str(theta0')])
title('v0 fixed theta')
set(gca,'FontSize',14,'YLim',[10 30],'XLim',[0 8])

%% Swimming speed for Fig 2d

figure
plot(DataDFFM.nflag,DataDFFM.v0,'ok','MarkerFaceColor','k','MarkerSize',10)
hold on
cellfun( @(x,y) errorbar(dataset.Nf, x.mean.U,x.std.U,'Color',y), dataset.data(6), dataset.color(6))
%legend({'Data',dataset.model{:}})
xlabel('N flagella')
ylabel('U (um/s)')
legend('Data','Model')
set(gca,'FontSize',14,'YLim',[0,30],'PlotBoxAspectRatio',[1 1 1])

%% Plot frequencies Fig 2c
figure
plot(DataDFFM.nflag,DataDFFM.wM,'ok','MarkerFaceColor','k','MarkerSize',10)
hold on
plot(DataDFFM.nflag,DataDFFM.wf,'ok','MarkerFaceColor','r','MarkerSize',10)
%plot(DataDFFM.nflag,DataDFFM.wb,'ok','MarkerFaceColor','b','MarkerSize',10)
legend('Of+Ob','Of')
xlabel('N flagella')
ylabel('Frequency (Hz)')
set(gca, 'FontSize',18,'PlotBoxAspectRatio',[1 1 1],'YLim',[160 230],'XLim',[0 8])

%% Exporting main Figure data
model_mainfigure = [dataset.Nf', dataset.data{6}.mean.U', dataset.data{6}.std.U'];
save('2024-10-01_ModelSpeed.txt','model_mainfigure','-ascii','-tabs')
%% Exporting ED Fig 5
dataToExport = dataset.data(rng);
gk=dataset.gk(rng);

TablesToExport = cell(size(dataToExport));
CellsToExport = cell(2,length(dataToExport));
L = table((1:8)','VariableNames',{'l'});

for i=1:length(dataToExport)
    Mean = struct2table(structfun(@(x) x',dataToExport{i}.mean,'UniformOutput',false));
    Std = struct2table(structfun(@(x) x',dataToExport{i}.std,'UniformOutput',false));
    Mean.Properties.VariableNames = cellfun(@(x)['mean_',x],Mean.Properties.VariableNames,'UniformOutput',false);
    Std.Properties.VariableNames = cellfun(@(x)['std_',x],Std.Properties.VariableNames,'UniformOutput',false);
    
    TablesToExport{i} = [L,Mean,Std];
    writetable(TablesToExport{i},['model_',dataToExport{i}.model,num2str(gk{i}),'.xlsx']);
    CellsToExport{2,i} = TablesToExport{i};
    CellsToExport{1,i} = dataToExport{i}.model;
end
