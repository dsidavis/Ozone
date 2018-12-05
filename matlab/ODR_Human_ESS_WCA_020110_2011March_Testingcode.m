clc, clear all, close all

global taxi
%
prompt = {'Use current parameters as Seed for FEV1? (yes,no)','Current parameters for FVC? (yes,no)',...
    'FEV1 only or FEV1 & FVC?(FEV1=1,FEV1&FVC=2)','Include Filtered Air?(yes=1,no=2)',...
    'Subtract Filtered Air?(yes=1,no=2)','Subtract Filtered Air Mean or Protocol(s)?(Mn=1,Prot=2)','Group Protocols?(yes=1,no=2)'};
dlgTitle = 'Demographics';
def = {'yes','yes','1','1','2','2','2'}';
lineNo = 1;
file = inputdlg(prompt,dlgTitle,lineNo,def)
drawnow;

state        = 'fat';
curseed{1}   = char(file{1})
curseed{2}   = char(file{2})
pftest       = str2num(file{3})
taxi.FiltAir = str2num(file{4})
subtractFA   = str2num(file{5})
subtractMnOrProt = str2num(file{6})
group_prot  = str2num(file{7})
% num_grps    = str2num(file{8})

samplerate   = 200;       %200 Hz 

% xlsfiles = {'WCA1997';'WCA2000b';'WCA2002';'WCA2003ab';'WCA2006ab';'ESS2009';'ESS2009_Prot2to5'}
% FAprots  = [2,5,2,1,1,1,1]
xlsfiles = {'WCA1997';'WCA2000b';'WCA2002';'WCA2003ab';'WCA2006ab';'ESS2009'}
FAprots  = [2,5,2,1,1,1]
num_grps_all = [4,5,3,5,6]

prot_grpnm = {'A','B','C','D','E','F'};

prot_grps.WCA1997.A = [1,3,4];
prot_grps.WCA1997.B = [5,6,7];
prot_grps.WCA1997.C = [1];
prot_grps.WCA1997.D = [5];

prot_num.WCA1997 = {'prot134','prot567','prot1','prot5'};

prot_grps.WCA2000b.A = [1:4];
prot_grps.WCA2000b.B = [3,4];
prot_grps.WCA2000b.C = [1];
prot_grps.WCA2000b.D = [2];
prot_grps.WCA2000b.E = [3];

prot_num.WCA2000b = {'prot1to4','prot34','prot1','prot2','prot3'};

prot_grps.WCA2002.A = [1,3];
prot_grps.WCA2002.B = [1,3,5,6];
prot_grps.WCA2002.C = [4];

prot_num.WCA2002 = {'prot13','prot1356','prot4'};

prot_grps.WCA2003ab.A = [3:6];
prot_grps.WCA2003ab.B = [3,4];
prot_grps.WCA2003ab.C = [5,6];
prot_grps.WCA2003ab.D = [7];
prot_grps.WCA2003ab.E = [8,9];

prot_num.WCA2003ab = {'prot3to6','prot34','prot56','prot7','prot89'};

prot_grps.WCA2006ab.A = [2:6];
prot_grps.WCA2006ab.B = [2,4];
prot_grps.WCA2006ab.C = [3,5,6];
prot_grps.WCA2006ab.D = [8];
prot_grps.WCA2006ab.E = [8,9];
prot_grps.WCA2006ab.F = [9];

prot_num.WCA2006ab = {'prot2to6','prot24','prot356','prot8','prot89','prot9'};

for study = 1:6
    clear rsqr_dpft_a adjrsqr_dpft_a param_dpft mn_1st2hr_dpft std_1st2hr_dpft stdBelmn_2hr_dpft underBaseline_2hr 
    clear under_std_1st2hr_dpft above_std_1st2hr_dpft mn_FA_dpft std_FA_dpft stdBelmn_FA_dpft underBaseline_FA
    clear under_std_FA_dpft above_std_FA_dpft baseline_dpft dpft_modxl dpft_msrxl rsqr0cept_dpft_a incl_samples
   
    FAprot  = FAprots(study)
    
    if study <6
        num_grps = num_grps_all(study)
    end

    xlsfile   = xlsfiles{study}
    paramfile = strcat('ODR_',xlsfile,'_optparam.mat');      
    
    if group_prot == 1 && study<6
        num_grps = num_grps_all(study)
        subgroup = '_subgrps';
    elseif group_prot == 2
        num_grps = 1
        subgroup = '_all';
    end
    
    if taxi.FiltAir == 1
        if subtractFA == 2 
            savename = strcat('ODR_',xlsfile,'_optparam.mat')
            savename2 = strcat('ODR_',xlsfile,subgroup,'_optparam.mat')
            xlsfile_b = xlsfile
        elseif subtractFA == 1 
            if subtractMnOrProt == 1
                savename = strcat('ODR_',xlsfile,'_PFTminusMnFA_optparam.mat')
                savename2 = strcat('ODR_',xlsfile,subgroup,'_PFTminusMnFA_optparam.mat')
                xlsfile_b = strcat(xlsfiles{study},'_PFTminusMnFA')
            elseif subtractMnOrProt == 2
                savename = strcat('ODR_',xlsfile,'_PFTminusFAprot_optparam.mat')
                savename2= strcat('ODR_',xlsfile,subgroup,'_PFTminusFA_optparam.mat')
                xlsfile_b = strcat(xlsfiles{study},'_PFTminusFAprot')
            end
        end
    elseif taxi.FiltAir == 2
        if subtractFA == 2    
            savename = strcat('ODR_',xlsfile,'_NoFAinMinErr_optparam.mat')
            savename2 = strcat('ODR_',xlsfile,subgroup,'_NoFAinMinErr_optparam.mat')
            xlsfile_b = strcat(xlsfiles{study},'_NoFAInMinErr')
        elseif subtractFA == 1 
            if subtractMnOrProt == 1
                savename = strcat('ODR_',xlsfile,'_NoFAinMinErr_PFTminusMnFA_optparam.mat')
                savename2 = strcat('ODR_',xlsfile,subgroup,'_NoFAinMinErr_PFTminusMnFA_optparam.mat')
                xlsfile_b = strcat(xlsfiles{study},'_NoFAInMinErr_PFTminusMnFA')
            elseif subtractMnOrProt == 2
                savename = strcat('ODR_',xlsfile,'_NoFAinMinErr_PFTminusFAprot_optparam.mat')
                savename2 = strcat('ODR_',xlsfile,subgroup,'_NoFAinMinErr_PFTminusFAprot_optparam.mat')
                xlsfile_b = strcat(xlsfiles{study},'_NoFAInMinErr_PFTminusFAprot')               
            end                
        end
    end
    
        % taxi.incl_prot

%     cd C:\LungMatLab\sigmoid
% %     load(paramfile)

cd C:\Users\Susan\Documents\Sue\ODR
%     cd C:\Sue\ODR

    load(savename)

%     [jnk,tx]       = xlsread(strcat('ozonedata\ODR_',xlsfile,'_091006\ODR_',xlsfile,'_FEV1_091006.xls'), 'A7:A100');
% 
%     subject        = xlsread(strcat('ozonedata\ODR_',xlsfile,'_091006\ODR_',xlsfile,'_RVE_091006.xls'), 'B1:AZ1');
%     VErest_msrd    = xlsread(strcat('ozonedata\ODR_',xlsfile,'_091006\ODR_',xlsfile,'_RVE_091006.xls'), 'B7:AZ100');
%     VEex_msrd      = xlsread(strcat('ozonedata\ODR_',xlsfile,'_091006\ODR_',xlsfile,'_EXVE_091006.xls'), 'B7:AZ100');
%     O3_msrd        = xlsread(strcat('ozonedata\ODR_',xlsfile,'_091006\ODR_',xlsfile,'_O3_091006.xls'), 'B7:AZ100');
%     pft_msrd.FEV1  = xlsread(strcat('ozonedata\ODR_',xlsfile,'_091006\ODR_',xlsfile,'_FEV1_091006.xls'), 'B7:AZ100');
%     pft_msrd.FVC   = xlsread(strcat('ozonedata\ODR_',xlsfile,'_091006\ODR_',xlsfile,'_FVC_091006.xls'), 'B7:AZ100');

    date = '_100121.xls';

    
    [jnk,tx_orig]       = xlsread(strcat('ODR_workbooks_100121\ODR_',xlsfile,date),'FEV1','A7:A100');

    subject        = xlsread(strcat('ODR_workbooks_100121\ODR_',xlsfile,date),'RVE', 'B1:AZ1');
    VErest_msrd_orig    = xlsread(strcat('ODR_workbooks_100121\ODR_',xlsfile,date),'RVE', 'B7:AZ100');
    VEex_msrd_orig      = xlsread(strcat('ODR_workbooks_100121\ODR_',xlsfile,date),'EXVE', 'B7:AZ100');
    O3_msrd_orig        = xlsread(strcat('ODR_workbooks_100121\ODR_',xlsfile,date),'O3', 'B7:AZ100');
    pft_msrd.FEV1_orig  = xlsread(strcat('ODR_workbooks_100121\ODR_',xlsfile,date),'FEV1', 'B7:AZ100');
    pft_msrd.FVC_orig   = xlsread(strcat('ODR_workbooks_100121\ODR_',xlsfile,date),'FVC', 'B7:AZ100');

%     sampleperiod_VErest = VErest_msrd_orig(:,1);
%     sampleperiod_VEex   = VEex_msrd_orig(:,1);
%     sampleperiod_O3     = O3_msrd_orig(:,1);
%     sampleperiod_FEV1   = pft_msrd.FEV1_orig(:,1);
%     sampleperiod_FVC    = pft_msrd.FVC_orig(:,1);

    protocol.VErest_orig = VErest_msrd_orig(:,2);
    protocol.VEex_orig   = VEex_msrd_orig(:,2);
    protocol.O3_orig     = O3_msrd_orig(:,2);
    protocol.FEV1_orig   = pft_msrd.FEV1_orig(:,2);
    protocol.FVC_orig    = pft_msrd.FVC_orig(:,2);

    if protocol.VErest_orig(1) == 2
        protocol.VErest_orig = VErest_msrd_orig(:,2)-1;
        protocol.VEex_orig   = VEex_msrd_orig(:,2)-1;
        protocol.O3_orig     = O3_msrd_orig(:,2)-1;
        protocol.FEV1_orig   = pft_msrd.FEV1_orig(:,2)-1;
        protocol.FVC_orig    = pft_msrd.FVC_orig(:,2)-1;
    end

    time.VErest_orig = VErest_msrd_orig(:,3:4);  %time is in minutes
    time.VEex_orig   = VEex_msrd_orig(:,3:4);
    time.O3_orig     = O3_msrd_orig(:,3:4);
    time.FEV1_orig   = pft_msrd.FEV1_orig(:,3:4);
    time.FVC_orig    = pft_msrd.FVC_orig(:,3:4);

    nonNaN_sI = find(subject>0);            %MatLab R2007a replaces Excel sheet blanks with NaN. More recent versions of MatLab do not do this.
    subject = subject(nonNaN_sI);

    % nonNaN_OI = find(O3_msrd_orig(1,:)~=NaN);
    msr.VErest_orig = VErest_msrd_orig(:,nonNaN_sI(1:end));
    msr.VEex_orig = VEex_msrd_orig(:,nonNaN_sI(1:end));
    msr.O3_orig   = O3_msrd_orig(:,nonNaN_sI(1:end));
    msr.FEV1_orig = pft_msrd.FEV1_orig(:,nonNaN_sI(1:end));
    msr.FVC_orig  = pft_msrd.FVC_orig(:,nonNaN_sI(1:end));

    rand('state',sum(100*clock));

    OptionsDiagnostics          = 'off';
    OptionsParameterDisplay     = 'final';
    OptionsJacobian             = 'off';
    OptionsLargeScale           = 'on';
    OptionsMaxFunEval           = 10000;
    OptionsMaxIter              = 10000;
    OptionsTolFun               = 1e-004;
    OptionsTolX                 = 1e-004;              
    OptionsLSHessian            = 'off';
    OptionsMSLevenbergMarquardt = 'on';
    OptionsDiffMinChange        = 2.0e-10;

    opts=optimset('DiffMinChange',OptionsDiffMinChange,'TolX',OptionsTolX,'Diagnostics',OptionsDiagnostics,'MaxFunEvals',OptionsMaxFunEval,'MaxIter',...
        OptionsMaxIter,'Display',OptionsParameterDisplay,'LevenbergMarquardt',OptionsMSLevenbergMarquardt,'Jacobian',...
        OptionsJacobian,'LargeScale',OptionsLargeScale);
    refine  = 1;
    reltol  = 8e-004

    % x_L = [-0.15       1     10    10 ]
    % x_U = [ 0       2500    800    15 ]


    % x_L = [-0.25    10      3   3.99 ];
    % x_U = [ 0     2500   1000    4.01 ];


    x_L = [-0.15       5       1    18 ]   %[A,DOS,tau,slope]
    x_U = [ 0       2500    1000    20 ]


    minmax   = 0;
    mvden    = 6; 
    ps       = 60;
    modl     =  2;
    errgoal  =  0; % min
    ac       = [2,2];% acceleration constants, only used for modl=0
    Iwt      = [0.9,0.4];  % intertia weights, only used for modl=0
    shw      =  100;   % how often to update display
    epoch    = 2000; % max iterations
    wt_end   = 1000; % iterations it takes to go from Iwt(1) to Iwt(2), only for modl=0
    errgrad  = 1e-004;   % lowest error gradient tolerance
    errgraditer = 7; % max # of epochs without error change >= errgrad
    PSOseed  = 1;    % if=1 then can input particle starting positions, if= 0 then all random
    %   % starting particle positions (first 20 at zero, just for an example)

    % good_s = [1,5,9,10,11,12,13,16,17,20,22,24,25,26,27,28,30,31];

    linrfit  = [1 0];
    whichstats = {'yhat','rsquare','adjrsquare','beta'};
    PFTS  = {'FEV1','FVC'};
    dPFTS = {'dFEV1','dFVC'};
    measure = {'VErest','VEex','O3','FEV1','FVC'};

    for grp = 1:num_grps 
%     for grp = 5:num_grps 
        clear incl_samples VE O3
        clear rsqr_dpft_a adjrsqr_dpft_a param_dpft mn_1st2hr_dpft std_1st2hr_dpft stdBelmn_2hr_dpft underBaseline_2hr 
        clear under_std_1st2hr_dpft above_std_1st2hr_dpft mn_FA_dpft std_FA_dpft stdBelmn_FA_dpft underBaseline_FA
        clear under_std_FA_dpft above_std_FA_dpft baseline_dpft dpft_modxl dpft_msrxl rsqr0cept_dpft_a
        
        if study <6
            Prot_tested = prot_num.(xlsfile)(grp)
        end
        
    for pft = 1:pftest
        PFT = char(PFTS(pft))
        PFT2 = strcat(PFT,'_orig')
        
        clear dpft_msrd TimePt dpftmod_all dpftmod_all_a dpftmod_all_xls doserate cumDose Mndoserate cnt_sub DR_2D MnDR_2D 
        clear cumDose_subset deltapft cnt avg_deltapft avg_sig1 R_sq_0cept_all R_sq_0cept_avg R_sq_avg R_sq_all
        clear RSQR_dpft_s1i0 yhat_dpft_s1i0 RSQR_dpft yhat_dpft R_sq_0cept R_sq RSQR_dpft_int0 yhat_dpft_int0 num_samples
        cnt = 1;
        
        for w = 1:length(measure)
            msrmt= char(measure(w));
            msrmt2 = strcat(msrmt,'_orig');
                        
            if taxi.FiltAir == 1                                %This redefines incl_prot and incl_samples if FA is included
                if group_prot == 2
                    incl_prot = [1:max(protocol.(msrmt2))];
                    incl_samples.(msrmt) = 1:length(protocol.(msrmt2));
                else 
                    incl_samples.(msrmt) = [];
                    prot_letr = char(prot_grpnm);
                    prot_letr2 = prot_letr(grp);
                    incl_prot = prot_grps.(xlsfile).(prot_letr2);
                    for n = 1:length(incl_prot)                        
                        incl_samp  = find(protocol.(msrmt2)==incl_prot(n));
                        incl_samples.(msrmt) = [incl_samples.(msrmt);incl_samp];                    
                    end
                end 
            elseif taxi.FiltAir == 2
                if FAprot == 1 
                    if study == 4
                        incl_prot     = [3:max(protocol.(msrmt2))];
                        incl_samples.(msrmt)  = find(protocol.(msrmt2)>2);
                    elseif study == 5
                        incl_prot     = [2:6,8:max(protocol.(msrmt2))];
                        incl_samples.(msrmt)  = [find(protocol.(msrmt2)>1 & protocol.(msrmt2)<7);find(protocol.(msrmt2)>7)];
                    elseif study == 6
                        incl_prot     = [2:max(protocol.(msrmt2))];
                        incl_samples.(msrmt)  = [find(protocol.(msrmt2)>1)];
                    end
                elseif FAprot == 2 && (study==1 || study==3) 
                    incl_prot         = [1,3:max(protocol.(msrmt2))];
                    incl_samples.(msrmt)      = [find(protocol.(msrmt2)<2);find(protocol.(msrmt2)>2)];
                elseif FAprot == 5 && study==2
                    incl_prot         = [1:4,6:max(protocol.(msrmt2))];
                    incl_samples.(msrmt)      = [find(protocol.(msrmt2)<5);find(protocol.(msrmt2)>5)];
                end 
            end
                protocol.(msrmt)   = protocol.(msrmt2)(incl_samples.(msrmt),:);
                time.(msrmt)       = time.(msrmt2)(incl_samples.(msrmt),:);  %time is in minutes
                msr.(msrmt)        = msr.(msrmt2)(incl_samples.(msrmt),:);
                num_samples.(msrmt)= length(incl_samples.(msrmt)); 
        end
        
        tx  = tx_orig(incl_samples.FEV1,:);
        
        protocol.VErest = protocol.VErest_orig(incl_samples.(msrmt),:);
        protocol.VEex   = protocol.VEex_orig(incl_samples.(msrmt),:);
        protocol.O3     = protocol.O3_orig(incl_samples.(msrmt),:);
        protocol.FEV1   = protocol.FEV1_orig(incl_samples.(msrmt),:);
        protocol.FVC    = protocol.FVC_orig(incl_samples.(msrmt),:);

        time.VErest = time.VErest_orig(incl_samples.(msrmt),:);  %time is in minutes
        time.VEex   = time.VEex_orig(incl_samples.(msrmt),:);
        time.O3     = time.O3_orig(incl_samples.(msrmt),:);
        time.FEV1   = time.FEV1_orig(incl_samples.(msrmt),:);
        time.FVC    = time.FVC_orig(incl_samples.(msrmt),:);

        msr.VErest = msr.VErest_orig(incl_samples.(msrmt),:);
        msr.VEex_msrd_orig = msr.VEex_orig(incl_samples.(msrmt),:);
        msr.O3_msrd_orig   = msr.O3_orig(incl_samples.(msrmt),:);
        msr.FEV1 = msr.FEV1_orig(incl_samples.(msrmt),:);
        msr.FVC  = msr.FVC_orig(incl_samples.(msrmt),:);
        
        [prot,Iprot,junk] = unique(protocol.VErest(:,1));
        if prot(1)>1
            prot = 1:length(prot);
        end
        
%         for s = 1:size(msr.VEex,2)
        for s = 1:2
            clear VE O3 taxi.msrd_dpft taxi.dpft taxi.MnDR taxi.DR taxi.pred_dpft stats msrmod DR_2D_a MnDR_2D_a cumD
            clear TimePt FA_samples taxi.incl_prot Coefficients_0cept S_err_0cept XTXI_0cept R_sq_0cept_a
            clear Coefficients S_err XTXI R_sq_a 
            
            if FAprot == 1 
                if study == 4
                    FA_samples     = find(protocol.FEV1_orig<3);
                    FA_samples_a   = find(protocol.FEV1_orig==1);
                    FA_samples_b   = find(protocol.FEV1_orig==2);

                    mn_FA_dpft(s,pft)  = mean(msr.(PFT2)(FA_samples,s),1);
                    std_FA_dpft(s,pft) = std(msr.(PFT2)(FA_samples,s),0,1);
                elseif study == 5
                    FA_samples     = [find(protocol.FEV1_orig==1); find(protocol.FEV1_orig==7)];
                    FA_samples_a   = find(protocol.FEV1_orig==1);
                    FA_samples_b   = find(protocol.FEV1_orig==7);

                    mn_FA_dpft(s,pft)  = mean(msr.(PFT2)(FA_samples,s),1);
                    std_FA_dpft(s,pft) = std(msr.(PFT2)(FA_samples,s),0,1);
                elseif study == 6
                    FA_samples     = find(protocol.FEV1_orig==1);

                    mn_FA_dpft(s,pft)  = mean(msr.(PFT2)(FA_samples,s),1)
                    std_FA_dpft(s,pft) = std(msr.(PFT2)(FA_samples,s),0,1);
                end
            elseif FAprot == 2 && (study==1 || study==3) 
                FA_samples     = find(protocol.FEV1_orig==2);

                mn_FA_dpft(s,pft)  = mean(msr.(PFT2)(FA_samples,s),1);
                std_FA_dpft(s,pft) = std(msr.(PFT2)(FA_samples,s),0,1);
            elseif FAprot == 5 && study==2
                FA_samples     = find(protocol.FEV1_orig==5);

                mn_FA_dpft(s,pft)  = mean(msr.(PFT2)(FA_samples,s),1);
                std_FA_dpft(s,pft) = std(msr.(PFT2)(FA_samples,s),0,1);
            end                         
            
            sub = strcat('sub',num2str(subject(s))) 
            for n = 1:size(msr.VEex,1)
                VE(time.VEex(n,1)+1:time.VEex(n,2)+1,protocol.VEex(n,1))        = msr.VEex(n,s);
            end
            for n = 1:size(msr.VErest,1)
                VE(time.VErest(n,1)+1:time.VErest(n,2)+1,protocol.VErest (n,1)) = msr.VErest(n,s);
            end

            for n = 1:size(msr.O3,1)
                O3(time.O3(n,1)+1:time.O3(n,2)+1,protocol.O3(n,1))  = msr.O3(n,s);
                if msr.O3(n,s) == 0
                    O3(time.O3(n,1)+1:time.O3(n,2)+1,protocol.O3(n,1)) = 0.00001;
                end
            end

%             VE(1,incl_prot) = msr.VErest([1;Iprot(1:end-1)+1],s); 
%             O3(1,incl_prot) = 0.00001; 
        %         O3(:,prot) = O3(:,prot) - O3(1,prot);
        
%             VE = VE(:,incl_prot);
%             O3 = O3(:,incl_prot);  

            [prot,Iprot,junk] = unique(protocol.VErest(:,1));
            
            VE = VE(:,incl_prot);
            O3 = O3(:,incl_prot);
            
            VE(1,:) = msr.VErest([1;Iprot(1:end-1)+1],s); 
            O3(1,:) = 0.00001; 

            [max_t,Imx_t] = max(time.O3(:,2));
            t = 0:max_t;
            taxi.t = repmat(t',1,length(incl_prot));

            doserate(:,:,s)  = VE.*O3*1.96;          %1.96 = ppm to ug/liter ozone conversion
            cumDose(:,:,s)   = cumsum(doserate(:,:,s).*[ones(1,size(doserate,2));diff(taxi.t)]);
            Mndoserate(2:size(taxi.t,1),:,s)= cumDose(2:end,:,s)./taxi.t(2:end,:);
        %     Mndoserate(1,:,s)  = 0;
            Mndoserate(1,:,s)= doserate(1,:,s);

            x_U(2) = max(max(cumDose(:,:,s)))

            time_pft_I      = find(time.(PFT)(:,2) == 0);
            difftime        = diff(time_pft_I);       
            [hrs,Ihrs,junk] = unique(difftime);
        %     [minIhrs,Indx] = min(Ihrs);
        
            if subtractFA == 1
                if subtractMnOrProt == 1
                    msr.(PFT)(:,s)          = msr.(PFT)(:,s) - mn_FA_dpft(s,pft);
                    msr.(PFT)(time_pft_I,s) = 0
                elseif subtractMnOrProt== 2
                    if study == 4
                        for n = 1:length(difftime)
                            if difftime(n)<=length(FA_samples_a)
                                msr.(PFT)(time_pft_I(n):time_pft_I(n+1)-1,s)= msr.(PFT)(time_pft_I(n):time_pft_I(n+1)-1,s) - mean([msr.(PFT2)(FA_samples_a(1:difftime(n)),s),msr.(PFT2)(FA_samples_b(1:difftime(n)),s)],2);
                            elseif difftime(n)>length(FA_samples_a)
                                diffbet = difftime(n)-length(FA_samples_a)
                                FAmn = repmat(mn_FA_dpft(s,pft),diffbet,1);
                                msr.(PFT)(time_pft_I(n):time_pft_I(n+1)-1,s)= msr.(PFT)(time_pft_I(n):time_pft_I(n+1)-1,s) - mean([[msr.(PFT2)(FA_samples_a,s);FAmn],[msr.(PFT2)(FA_samples_b,s);FAmn]],2); 
                            end
                        end
                        msr.(PFT)(time_pft_I(end):end,s)= msr.(PFT)(time_pft_I(end):end,s) - mean([msr.(PFT2)(FA_samples_a(1:length(msr.(PFT)(time_pft_I(end):end,s))),s),msr.(PFT2)(FA_samples_b(1:length(msr.(PFT)(time_pft_I(end):end,s))),s)],2);
                    else
                        for n = 1:length(difftime)
                            if study==5
%                                 if time_pft_I(n) <7
                                if difftime(n) <7
                                    FA_samples = FA_samples_a;
%                                 elseif time_pft_I(n) >6
                                elseif difftime(n) >6
                                    FA_samples = FA_samples_b;
                                end
                            end
                            if difftime(n)<=length(FA_samples)
                                msr.(PFT)(time_pft_I(n):time_pft_I(n+1)-1,s)= msr.(PFT)(time_pft_I(n):time_pft_I(n+1)-1,s) - msr.(PFT2)(FA_samples(1:difftime(n)),s);
                            elseif difftime(n)>length(FA_samples)
                                diffbet = difftime(n)-length(FA_samples)
                                FAmn = repmat(mn_FA_dpft(s,pft),diffbet,1);
                                msr.(PFT)(time_pft_I(n):time_pft_I(n+1)-1,s)= msr.(PFT)(time_pft_I(n):time_pft_I(n+1)-1,s) - [msr.(PFT2)(FA_samples,s);FAmn]; 
                            end
                        end
                        msr.(PFT)(time_pft_I(end):end,s)= msr.(PFT)(time_pft_I(end):end,s) - msr.(PFT2)(FA_samples(1:length(msr.(PFT)(time_pft_I(end):end,s))),s);
                    end
                end
            end       
        
            if length(incl_prot) == 1
                TimePt           = zeros(1,length(incl_prot));
                dpft_msrd(:,:,s) = zeros(1,length(incl_prot));               
            else
                TimePt           = zeros(max(hrs),length(incl_prot));
                dpft_msrd(:,:,s) = zeros(max(hrs),length(incl_prot));
            end
            for n = 1:length(time_pft_I)
                if n<length(time_pft_I)
                    TimePt(1:time_pft_I(n+1)-time_pft_I(n),n)      = time.(PFT)(time_pft_I(n):time_pft_I(n+1)-1,2);
                    dpft_msrd(1:time_pft_I(n+1)-time_pft_I(n),n,s) = msr.(PFT)(time_pft_I(n):time_pft_I(n+1)-1,s);
                    cumD(:,n)                                      = cumDose(TimePt(:,n)+1,n,s);
                    DR_2D_a(:,n)                                   = doserate(TimePt(:,n)+1,n,s);
                    MnDR_2D_a(:,n)                                 = Mndoserate(TimePt(:,n)+1,n,s);    
                else
                    TimePt(1:length(time.(PFT)(time_pft_I(n):end,2)),n)     = time.(PFT)(time_pft_I(n):end,2);
                    dpft_msrd(1:length(time.(PFT)(time_pft_I(n):end,2)),n,s)= msr.(PFT)(time_pft_I(n):end,s);
                    cumD(:,n)                                       = cumDose(TimePt(:,n)+1,n,s);
                    DR_2D_a(:,n)                                      = doserate(TimePt(:,n)+1,n,s);
                    MnDR_2D_a(:,n)                                    = Mndoserate(TimePt(:,n)+1,n,s);    
                end           
            end
            
            taxi.incl_prot =incl_prot;
            taxi.DRindx = 1+TimePt;

            taxi.DR   = doserate(:,:,s);
            taxi.MnDR = Mndoserate(:,:,s);

            taxi.msrd_dpft = dpft_msrd(:,:,s);

            taxi.err = 6000;
            taxi.fit = 0;

            switch curseed{pft}       
                case 'yes'               
    % %                x0_pft = [ODR.(dPFT).A(s),ODR.(dPFT).DOS(s),ODR.(dPFT).tau(s),ODR.(dPFT).slope(s)]

                    x0_pft = [ODR.(PFT).A(s),ODR.(PFT).DOS(s),ODR.(PFT).tau(s),ODR.(PFT).slope(s)]
                otherwise
                    if pft == 1
                        x0_pft = [-0.05,1200,100,15]  
                    else
                        x0_pft = [ODR.FEV1.A(s),ODR.FEV1.DOS(s),ODR.FEV1.tau(s),ODR.FEV1.slope(s)]
                    end                    
            end 

            y_0  = taxi.msrd_dpft(1,1); 

            taxi.pft = 1;    
            taxi.pso = 1;    

            dims     = length(x0_pft);
            varrange = [x_L', x_U'];
            mv       =[(varrange(:,2)-varrange(:,1))./mvden];
            PSOseedValue = repmat(x0_pft,round(0.5*ps),1);
        %     PSOseedValue = repmat(x0_pft,ps,1);            
            psoparams    =[shw epoch ps ac(1) ac(2) Iwt(1) Iwt(2) wt_end errgrad errgraditer errgoal modl PSOseed];

        %  % run pso
            [pso_out,te,tr]=pso_Trelea_vectorized('ModODR_020110', dims,mv, varrange, minmax, psoparams,[],PSOseedValue);
        %   
        % %-------------------------------------------------------------------------- 
        %     disp('Best fit parameters: cost = 0.5*sum((taxi.expForce-modForce).^2)');

            disp([' cost = ',num2str(pso_out(end))]);
            disp(['max(te) = ',num2str(max(te))]);

            x0_pft = pso_out(1:end-1)'
            cost  = pso_out(end);

            mindiff = repmat(5e-4,1,length(x0_pft)); 
            taxi.err = 6000;
            count = 1;
            cycle = 0;
            taxi.pso = 0;
            while count == 1 && cycle<5      
                [estPara,resnorm,residual,exitflag,output]= lsqnonlin(@ModODR_020110,x0_pft,x_L,x_U,opts);        
                if abs(estPara-x0_pft) < mindiff
                    count = 0
                else
                    count = 1 
                    x0_pft   = estPara
                    cycle = cycle + 1;
                end        
            end     

            x0_pft = estPara
            cost    = 0.5*resnorm

            ODR.(PFT).A(s)    = x0_pft(1);
            ODR.(PFT).DOS(s)  = x0_pft(2);
            ODR.(PFT).tau(s)  = x0_pft(3);
            ODR.(PFT).slope(s)= x0_pft(4);
        %     ODR.(PFT).BS(s)   = x0_pft(5);

%             taxi.DOS = x0_pft(2);
            taxi.tau = x0_pft(3);
            taxi.slope = x0_pft(4);
            taxi.fit = 1;
            taxi.err = 6000;
            count = 1;
            cycle = 0;

            taxi.pso = 1;    
            dims     = length(x0_pft(1:2));
            varrange = [x_L(1:2)', x_U(1:2)'];
            mv       =[(varrange(:,2)-varrange(:,1))./mvden];
%             PSOseedValue = repmat(x0_pft(1),round(0.5*ps),1);
        %     PSOseedValue = repmat(x0_pft,ps,1);
            PSOseedValue = repmat(1.4*x0_pft(1:2),round(0.3*ps),1);
            psoparams    =[shw epoch ps ac(1) ac(2) Iwt(1) Iwt(2) wt_end errgrad errgraditer errgoal modl PSOseed];

        %  % run pso
            [pso_out,te,tr]=pso_Trelea_vectorized('ModODR_020110', dims,mv, varrange, minmax, psoparams,[],PSOseedValue);
        %   
        % %-------------------------------------------------------------------------- 
        %     disp('Best fit parameters: cost = 0.5*sum((taxi.expForce-modForce).^2)');

            disp([' cost = ',num2str(pso_out(end))]);
            disp(['max(te) = ',num2str(max(te))]);

            x0_pft = pso_out(1:end-1)'
            cost  = pso_out(end);
            
            taxi.pso = 0;
            
            while count == 1 && cycle<5      
                [estPara,resnorm,residual,exitflag,output]= lsqnonlin(@ModODR_020110,x0_pft(1:2),x_L(1:2),x_U(1:2),opts);        
                if abs(estPara-x0_pft(1:2)) < mindiff(1:2)
                    count = 0
                else
                    count = 1 
                    x0_pft   = estPara
                    cycle = cycle + 1;
                end        
            end     

            x0_pft = estPara
            cost    = 0.5*resnorm

            ODR.(PFT).A(s)    = x0_pft(1);
            ODR.(PFT).DOS(s)  = x0_pft(2);
            
            pred_PFT     = taxi.pred_dpft(:,:);
            msrd_PFT     = taxi.msrd_dpft(:,:);

            linrfit  = polyfit(msrd_PFT(:),pred_PFT(:),1);
            [RSQR_dpft(s),yhat_dpft(:,s)]=rsqr_101907(linrfit,msrd_PFT(:),pred_PFT(:));
            linrfit_s1i0      = [1 0];
            [RSQR_dpft_s1i0(s),yhat_dpft_s1i0(:,s)]=rsqr_101907(linrfit_s1i0,msrd_PFT(:),pred_PFT(:));
            linrfit_int0 = [linrfit(1),0];
            [RSQR_dpft_int0(s),yhat_dpft_int0(:,s)]=rsqr_101907(linrfit_int0,msrd_PFT(:),pred_PFT(:));
%             [RSQR_dpft(s),yhat_dpft(:,s)]=rsqr_101907(linrfit,pred_PFT(:),msrd_PFT(:));
        
            stats            = regstats(pred_PFT(:),msrd_PFT(:),'linear',whichstats);
            yhat_dpft_a      = stats.yhat;
            beta_dpft_a      = stats.beta;
            RSQR_dpft_a(s)   = stats.rsquare;
            adjRSQR_dpft_a(s)= stats.adjrsquare;
            
            [Coefficients_0cept, S_err_0cept, XTXI_0cept, R_sq_0cept_a] = mregress(pred_PFT(:), msrd_PFT(:), 0);
            [Coefficients, S_err, XTXI, R_sq_a]       = mregress(pred_PFT(:), msrd_PFT(:), 1);
%             [Coefficients, S_err, XTXI, R_sq, F_val, Coef_stats, Y_hat, residuals, covariance] = mregress(msrd_PFT(:),pred_PFT(:),1);

            R_sq_0cept(s) = R_sq_0cept_a;
            R_sq(s)       = R_sq_a;
            
%             msrmod = [msrd_PFT(:),pred_PFT(:)];
%             RSQR_dpft_b(s) = apR2p(msrmod);

        %     save (paramfile, 'ODR')
            save (savename2, 'ODR')

            % set(0,'DefaultAxesColorOrder',[1 0 0;0 1 0;0 0 1;.5 .5 .5]);

            t60 = TimePt;
            dpftmod_all_a(:,s) = zeros(size(msr.(PFT),1),1);
            DR_2D(:,s) = zeros(size(msr.(PFT),1),1);
            MnDR_2D(:,s) = zeros(size(msr.(PFT),1),1);
            cumDose_subset(:,s) = zeros(size(msr.(PFT),1),1);

            for n = 1:length(time_pft_I)
                if n<length(time_pft_I)
                    dpftmod_all_a(time_pft_I(n):time_pft_I(n+1)-1,s) = taxi.pred_dpft(1:time_pft_I(n+1)-time_pft_I(n),n);
                    DR_2D(time_pft_I(n):time_pft_I(n+1)-1,s)        = DR_2D_a(1:time_pft_I(n+1)-time_pft_I(n),n);
                    MnDR_2D(time_pft_I(n):time_pft_I(n+1)-1,s)      = MnDR_2D_a(1:time_pft_I(n+1)-time_pft_I(n),n);
                    cumDose_subset(time_pft_I(n):time_pft_I(n+1)-1,s)= cumD(1:time_pft_I(n+1)-time_pft_I(n),n);
                else
                    dpftmod_all_a(time_pft_I(n):end,s)  = taxi.pred_dpft(1:length(time.(PFT)(time_pft_I(n):end,2)),n);
                    DR_2D(time_pft_I(n):end,s)         = DR_2D_a(1:length(time.(PFT)(time_pft_I(n):end,2)),n);
                    MnDR_2D(time_pft_I(n):end,s)       = MnDR_2D_a(1:length(time.(PFT)(time_pft_I(n):end,2)),n);
                    cumDose_subset(time_pft_I(n):end,s)=cumD(1:length(time.(PFT)(time_pft_I(n):end,2)),n);
                end           
            end

        %     cumDose_subset(:,s) = cumD(:);
        %     DR_2D(:,s)          = DR_2D_a(:);
        %     MnDR_2D(:,s)        = MnDR_2D_a(:);

            deltapft(:,:,s) = taxi.pred_deltapft;
        %    sig1(:,:,s)    = taxi.sig1_dpft;

            mn_1st2hr_dpft(s,pft)  = mean(mean(dpft_msrd(1:3,:,s),1));
            std_1st2hr_dpft(s,pft) = mean(std(dpft_msrd(1:3,:,s),0,1));

            stdBelmn_2hr_dpft(s,pft) = mn_1st2hr_dpft(s,pft)-std_1st2hr_dpft(s,pft);

            underBaseline_2hr(s,pft)     = 100*length(find(dpft_msrd(:,:,s)<stdBelmn_2hr_dpft(s)))./num_samples.FEV1;
            under_std_1st2hr_dpft(s,pft) = 100*length(find(dpft_msrd(:,:,s)<-std_1st2hr_dpft(s)))./num_samples.FEV1;

            above_std_1st2hr_dpft(s,pft) = 100*length(find(dpft_msrd(:,:,s)>std_1st2hr_dpft(s,pft)))./num_samples.FEV1;
            
%             if FAprot == 1 
%                 if study == 4
%                     mn_FA_dpft(s,pft)  = mean(mean(dpft_msrd(:,1:2,s),1));
%                     std_FA_dpft(s,pft) = mean(std(dpft_msrd(:,1:2,s),0,1));
%                 elseif study == 5
%                     mn_FA_dpft(s,pft)  = mean(mean(dpft_msrd(:,[1,7],s),1));
%                     std_FA_dpft(s,pft) = mean(std(dpft_msrd(:,[1,7],s),0,1));
%                 elseif study == 6
%                     mn_FA_dpft(s,pft)  = mean(dpft_msrd(:,1,s),1);
%                     std_FA_dpft(s,pft) = std(dpft_msrd(:,1,s),0,1);
%                 end
%             elseif FAprot == 2 && (study==1 || study==3) 
%                 mn_FA_dpft(s,pft)  = mean(dpft_msrd(:,2,s),1);
%                 std_FA_dpft(s,pft) = std(dpft_msrd(:,2,s),0,1);
%             elseif FAprot == 5 && study==2
%                 mn_FA_dpft(s,pft)  = mean(dpft_msrd(:,5,s),1);
%                 std_FA_dpft(s,pft) = std(dpft_msrd(:,5,s),0,1);
%             end          
            
%             if taxi.FiltAir == 1
                stdBelmn_FA_dpft(s,pft) = mn_FA_dpft(s,pft)-std_FA_dpft(s,pft);
                underBaseline_FA(s,pft) = 100*length(find(dpft_msrd(:,:,s)<stdBelmn_FA_dpft(s,pft)))./num_samples.FEV1;
                under_std_FA_dpft(s,pft)= 100*length(find(dpft_msrd(:,:,s)<-std_FA_dpft(s,pft)))./num_samples.FEV1;
                above_std_FA_dpft(s,pft)= 100*length(find(dpft_msrd(:,:,s)>std_FA_dpft(s,pft)))./num_samples.FEV1;
%             else
%                 stdBelmn_FA_dpft = NaN;
%                 underBaseline_FA = NaN;
%                 under_std_FA_dpft= NaN;
%                 above_std_FA_dpft= NaN;                
%             end

            if pft == 1 && s == 12 && grp == 2
                figure
                 subplot(2,1,1)
                 plot(cumD(:),taxi.msrd_dpft(:),'b*',cumD(:),taxi.pred_dpft(:), 'ro',...
                     cumD(:),repmat(mn_1st2hr_dpft(s,pft),length(cumD(:))),'g--',cumD(:),repmat(stdBelmn_2hr_dpft(s,pft),length(cumD(:))),'c'),...
                     grid on, title(sub),xlabel('Cumulative Dose (ppm?)'),ylabel([PFT,' (% of initial)']),legend('Meas','Mod','Mean 2h','Std 2h' )
                 subplot(2,1,2)
                 plot(t60(:),taxi.msrd_dpft(:),'b*',t60(:),taxi.pred_dpft(:), 'ro'),...
                     grid on, title(sub),xlabel('Time (min)'),ylabel([PFT,' (% of initial)']),legend('Meas','Mod')
            end

            cnt_sub(cnt) = s;
            cnt = cnt +1;
        end

        clear VE O3 taxi.msrd_dpft taxi.dFVC taxi.MnDR taxi.DR taxi.pred_dpft taxi.pred_dFVC stats msrmod

        ODR.(PFT).A_all_avg     = mean(ODR.(PFT).A(1:cnt-1));
        ODR.(PFT).DOS_all_avg   = mean(ODR.(PFT).DOS(1:cnt-1));
        ODR.(PFT).tau_all_avg   = mean(ODR.(PFT).tau(1:cnt-1));
        ODR.(PFT).slope_all_avg = mean(ODR.(PFT).slope(1:cnt-1));
        
        dpftmsr_all_xls = msr.(PFT)(:,1:cnt-1);
        dpftmod_all_xls = dpftmod_all_a(:,1:cnt-1);

        dpftmsr_all = msr.(PFT)(:,1:cnt-1);
        dpftmod_all = dpftmod_all_a(:,1:cnt-1);
        
        linrfit_s1i0      = [1 0];
        linrfit  = polyfit(dpftmsr_all(:),dpftmod_all(:),1);
        [RSQR_dpft_all,yhat_dpft_all]=rsqr_101907(linrfit_s1i0,dpftmsr_all(:),dpftmod_all(:));
        [RSQR_dpft_all,yhat_dpft_all]=rsqr_101907(linrfit,dpftmsr_all(:),dpftmod_all(:));

        stats             = regstats(dpftmod_all(:),dpftmsr_all(:),'linear',whichstats);
        yhat_dpft_all_a   = stats.yhat;
        beta_dpft_all_a   = stats.beta;
        RSQR_dpft_all_a   = stats.rsquare;
        adjRSQR_dpft_all_a= stats.adjrsquare;
        
        [Coefficients_all, S_err_all, XTXI_all, R_sq_0cept_all, F_val_all, Coef_stats_all, Y_hat_all, residuals_all, covariance_all] = mregress(dpftmod_all(:),dpftmsr_all(:), 0);
        [Coefficients_all, S_err_all, XTXI_all, R_sq_all, F_val_all, Coef_stats_all, Y_hat_all, residuals_all, covariance_all] = mregress(dpftmod_all(:),dpftmsr_all(:), 1);

%         msrmod_all = [dpftmsr_all(:),dpftmod_all(:)];
        % RSQR_dpft_all = apR2p(msrmod_all);

        x_U(2) = mean(max(cumDose_subset))

        taxi.err = 1e+006;
        taxi.fit = 0;

        switch curseed{pft}       
            case 'yes'    
    % %            x0_pft = [ODR.(dPFT).A_avg,ODR.(dPFT).DOS_avg,ODR.(dPFT).tau_avg,ODR.(dPFT).slope_avg]

               x0_pft = [ODR.(PFT).A_avg,ODR.(PFT).DOS_avg,ODR.(PFT).tau_avg,ODR.(PFT).slope_avg];
            otherwise
                if pft == 1
                   x0_pft = [-0.05,1200,100,15]; 
                else
                   x0_pft = [ODR.FEV1.A_avg,ODR.FEV1.DOS_avg,ODR.FEV1.tau_avg,ODR.FEV1.slope_avg]
                end                
        end 

        taxi.DR    = mean(doserate(:,:,:),3);
        taxi.MnDR  = mean(Mndoserate(:,:,:),3);

        % taxi.msrd_dpft = mean(msr.(PFT)(:,:),2);
        taxi.msrd_dpft = mean(dpft_msrd,3);

        y_0  = taxi.msrd_dpft(1,1); 

        taxi.pft = 1;
        taxi.pso = 1;

        dims     = length(x0_pft);
        varrange = [x_L', x_U'];
        mv       =[(varrange(:,2)-varrange(:,1))./mvden];
        %    PSOseedValue = repmat([0],ps-10,1);
        PSOseedValue = repmat(x0_pft,round(0.5*ps),1);
        % PSOseedValue = repmat(x0_pft,ps,1);            
        psoparams    =[shw epoch ps ac(1) ac(2) Iwt(1) Iwt(2) wt_end errgrad errgraditer errgoal modl PSOseed];

        %  % run pso
        [pso_out,te,tr]=pso_Trelea_vectorized('ModODR_020110', dims,mv, varrange, minmax, psoparams,[],PSOseedValue);
        %   
        % %-------------------------------------------------------------------------- 
        % disp('Best fit parameters: cost = 0.5*sum((taxi.expForce-modForce).^2)');

        disp([' cost = ',num2str(pso_out(end))]);
        disp(['max(te) = ',num2str(max(te))]);

        x0_pft = pso_out(1:end-1)'
        cost  = pso_out(end)

        mindiff = repmat(5e-4,1,length(x0_pft));
        taxi.err = 1e+006;
        taxi.pso = 0;
        count = 1;
        cycle = 0;
        while count == 1 && cycle<5      
            [estPara,resnorm,residual,exitflag,output]= lsqnonlin(@ModODR_020110,x0_pft,x_L,x_U,opts);        
            if abs(estPara-x0_pft) < mindiff
                count = 0
            else
                count = 1 
                x0_pft   = estPara
                cycle = cycle + 1;
            end        
        end     

        x0_pft = estPara
        cost    = 0.5*resnorm

        ODR.(PFT).A_avg    = x0_pft(1);
        ODR.(PFT).DOS_avg  = x0_pft(2);
        ODR.(PFT).tau_avg  = x0_pft(3);
        ODR.(PFT).slope_avg= x0_pft(4);

%         taxi.DOS = x0_pft(2);
        taxi.tau = x0_pft(3);
        taxi.slope = x0_pft(4);
        taxi.fit = 1;
        taxi.err = 1e+006;
        taxi.pso = 0;
        count = 1;
        cycle = 0;
        while count == 1 && cycle<5      
            [estPara,resnorm,residual,exitflag,output]= lsqnonlin(@ModODR_020110,x0_pft(1:2),x_L(1:2),x_U(1:2),opts);        
            if abs(estPara-x0_pft(1:2)) < mindiff(1:2)
                count = 0
            else
                count = 1 
                x0_pft   = estPara
                cycle = cycle + 1;
            end        
        end     

        x0_pft = estPara
        cost    = 0.5*resnorm

        ODR.(PFT).A_avg    = x0_pft(1);
        ODR.(PFT).DOS_avg  = x0_pft(2);

        avgPred_dpft = zeros(size(msr.(PFT),1),1);

        for n = 1:length(time_pft_I)
            if n<length(time_pft_I)
                avgPred_dpft(time_pft_I(n):time_pft_I(n+1)-1) = taxi.pred_dpft(1:time_pft_I(n+1)-time_pft_I(n),n);
            else
                avgPred_dpft(time_pft_I(n):end)  = taxi.pred_dpft(1:length(time.(PFT)(time_pft_I(n):end,2)),n);
            end           
        end
        avgMsr_dpft_xls = mean(msr.(PFT),2)
        avgMsr_dpft = mean(dpftmsr_all,2)
%         avgMsr_dpft = mean(msr.(PFT),2);
        msrmod      = [avgMsr_dpft,avgPred_dpft(:,:)];

        avgCumDose = mean(cumDose_subset,2);

        stats              = regstats(avgPred_dpft(:,:),avgMsr_dpft,'linear',whichstats);
        yhat_dpft_avg_a    = stats.yhat;
        beta_dpft_avg_a    = stats.beta;
        RSQR_dpft_avg_a    = stats.rsquare;
        adjRSQR_dpft_avg_a = stats.adjrsquare;
        
        [Coefficients_avg, S_err_avg, XTXI_avg, R_sq_0cept_avg, F_val_avg, Coef_stats_avg, Y_hat_avg, residuals_avg, covariance_avg] = mregress(avgPred_dpft(:,:),avgMsr_dpft, 0);
        [Coefficients_avg, S_err_avg, XTXI_avg, R_sq_avg, F_val_avg, Coef_stats_avg, Y_hat_avg, residuals_avg, covariance_avg] = mregress(avgPred_dpft(:,:),avgMsr_dpft, 1);

        avg_deltapft(:,:,s)= taxi.pred_deltapft;
        avg_sig1(:,:,s)    = taxi.sig1_dpft;

        % save (paramfile, 'ODR')
        save (savename2, 'ODR')

        % t60 = TimePt;
        t60 = time.(PFT)(:,2);
        % figure
        %  subplot(2,1,1)
        %  plot(t60,avgMsr_dpft,'b*',t60,avgPred_dpft, 'ro'),...
        %      grid on, title('Mean'),xlabel('Time (min)'),ylabel('dpft'),legend('Meas','Mod')
        %  subplot(2,1,2)
        %  plot(avgMsr_dpft,avgPred_dpft, 'k+'),...
        %      grid on, title('Mean'),xlabel('Measured'),ylabel('Modeled')%,legend('pft','FVC')

        % rsqr_dpft = [RSQR_dpft,NaN,RSQR_dpft_avg,NaN,RSQR_dpft_all];
        rsqr_dpft_a(:,:,pft)    = [RSQR_dpft_a,NaN,RSQR_dpft_avg_a,NaN,RSQR_dpft_all_a];
        adjrsqr_dpft_a(:,:,pft) = [adjRSQR_dpft_a,NaN,adjRSQR_dpft_avg_a,NaN,adjRSQR_dpft_all_a];
        rsqr0cept_dpft_a(:,:,pft) = [R_sq_0cept,NaN,R_sq_0cept_avg,NaN,R_sq_0cept_all];                

        param_dpft(:,:,pft) = [ODR.(PFT).A,NaN,ODR.(PFT).A_avg,NaN,ODR.(PFT).A_all_avg;ODR.(PFT).DOS,NaN,ODR.(PFT).DOS_avg,NaN,ODR.(PFT).DOS_all_avg;...
            ODR.(PFT).tau,NaN,ODR.(PFT).tau_avg,NaN,ODR.(PFT).tau_all_avg;ODR.(PFT).slope,NaN,ODR.(PFT).slope_avg,NaN,ODR.(PFT).slope_all_avg]; 

        avg_DR     = mean(DR_2D,2);
        avg_MnDR   = mean(MnDR_2D,2);

        avg_mn_1st2hr_dpft  = mean(mn_1st2hr_dpft(:,pft));
        avg_std_1st2hr_dpft = mean(std_1st2hr_dpft(:,pft));

        avg_stdBelmn_2hr_dpft = avg_mn_1st2hr_dpft-avg_std_1st2hr_dpft;

        avg_underBaseline_2hr     = 100*length(find(avgMsr_dpft<avg_stdBelmn_2hr_dpft))./num_samples.FEV1;
        avg_under_std_1st2hr_dpft = 100*length(find(avgMsr_dpft<-avg_std_1st2hr_dpft))./num_samples.FEV1;

        avg_above_std_1st2hr_dpft = 100*length(find(avgMsr_dpft>avg_std_1st2hr_dpft))./num_samples.FEV1;

%         if taxi.FiltAir == 1
            avg_mn_FA_dpft  = mean(mn_FA_dpft(:,pft));
            avg_std_FA_dpft = std(std_FA_dpft(:,pft));

            avg_stdBelmn_FA_dpft = avg_mn_FA_dpft-avg_std_FA_dpft;

            avg_underBaseline_FA  = 100*length(find(avgMsr_dpft<avg_stdBelmn_FA_dpft))./num_samples.FEV1;
            avg_under_std_FA_dpft = 100*length(find(avgMsr_dpft<-avg_std_FA_dpft))./num_samples.FEV1;

            avg_above_std_FA_dpft = 100*length(find(avgMsr_dpft>avg_std_FA_dpft))./num_samples.FEV1;
%         else
%             avg_mn_FA_dpft  = NaN;
%             avg_std_FA_dpft = NaN;
%             avg_stdBelmn_FA_dpft = NaN;
%             avg_underBaseline_FA  = NaN;
%             avg_under_std_FA_dpft = NaN;
%             avg_above_std_FA_dpft = NaN;
%         end

        baseline_dpft(:,:,pft) = [mn_1st2hr_dpft(:,pft)',NaN,avg_mn_1st2hr_dpft;std_1st2hr_dpft(:,pft)',NaN,avg_std_1st2hr_dpft;...
            mn_FA_dpft(:,pft)',NaN,avg_mn_FA_dpft;std_FA_dpft(:,pft)',NaN,avg_std_FA_dpft;...
            underBaseline_2hr(:,pft)',NaN,avg_underBaseline_2hr;underBaseline_FA(:,pft)',NaN,avg_underBaseline_FA;...
            repmat(NaN,1,length(under_std_1st2hr_dpft(:,pft)')+2);under_std_1st2hr_dpft(:,pft)',NaN,avg_under_std_1st2hr_dpft;under_std_FA_dpft(:,pft)',NaN,avg_under_std_FA_dpft;...
            repmat(NaN,1,length(above_std_1st2hr_dpft(:,pft)')+2);above_std_1st2hr_dpft(:,pft)',NaN,avg_above_std_1st2hr_dpft;above_std_FA_dpft(:,pft)',NaN,avg_above_std_FA_dpft];

        figure
         subplot(2,1,1)
         plot(avgCumDose(:,:),avgMsr_dpft,'b*',avgCumDose(:,:),avgPred_dpft(:,:), 'ro',...
             avgCumDose,repmat(avg_mn_1st2hr_dpft,length(avgCumDose)),'g--',avgCumDose,repmat(avg_stdBelmn_2hr_dpft,length(avgCumDose)),'c'),...
             grid on, title('Mean'),xlabel('Cumulative Dose (ppm?)'),ylabel([PFT,' (% of initial)']),legend('Meas','Mod')
         subplot(2,1,2)
         plot(t60,avgMsr_dpft,'b*',t60,avgPred_dpft(:,:), 'ro'),...
             grid on, title('Mean'),xlabel('Time (min)'),ylabel([PFT,' (% of initial)']),legend('Meas','Mod')

        
         
        h_sub        = subject';
        sub_cell     = num2cell([h_sub,NaN]);
        twomeans     = {'MeanData','','MeanParam'};
        onemean      = {'MeanData'};
        h_twomeans   = [sub_cell,twomeans];
        h_onemean    = [sub_cell,onemean];         
         
%         h_sub        = 1:cnt-1;
%         h_all        = {'MeanParam'};
%         h_avg        = {'MeanData'};
        % h_param  = {'A'; 'DOS';'tau'; 'slope'; 'BS'};
        h_param      = {'A';'DOS';'tau';'slope'};
        h_RSQR_yint0  = {'RSQR_yint0'};
        h_RSQR_a     = {'RSQR_trueCoeff'};
        h_RSQR_a_adj = {'RSQRadj_trueCoeff'};
%         h_pftRSQR{pft}= strcat(PFT,'_RSQR_1');
        h_pftRSQR_a{pft}= strcat(PFT,'_RSQR');
        h_pft{pft}   = PFT;
        h_MnStd      = {'Mean 1st 2h';'Std 1st 2h';'Mean FA';'Std FA';'%<(Mn-Std) 2h';'%<(Mn-Std) FA';...
            ' ';'%<-Std 1st 2h';'%<-Std FA';' ';'%>Std 1st 2h';'%>Std FA'}; 

        h_protocol = tx;
        h_model   = {'Modeled'};
        h_measure = {'Measured'};
        h_DR   = {'DR'};
        h_mnDR = {'mnDR'};
        h_cumDose = {'Cumulative Dose'};

        h_time = {'Time(min)'};

        NaN_fill = repmat(NaN,size(avg_MnDR,1),1);
        NaN_fill_2 = repmat(NaN,size(dpftmod_all_xls,1),1);

        DRxl = [DR_2D,NaN_fill,avg_DR];
        MnDRxl = [MnDR_2D,NaN_fill,avg_MnDR];
        cDosex1 = [cumDose_subset,NaN_fill,avgCumDose];

        dpft_modxl(:,:,pft) = [dpftmod_all_xls,NaN_fill_2,avgPred_dpft];
        dpft_msrxl(:,:,pft) = [msr.(PFT)(:,1:cnt-1),NaN_fill_2,avgMsr_dpft_xls];

        xlmnDR_hrow = num2str(size(DRxl,1)+4);
        xlmnDRrow = num2str(size(DRxl,1)+5);

        xlmsFEV_hrow = num2str(size(DRxl,1)+26);
        xlmsFEVrow = num2str(size(DRxl,1)+27);
        
        xlmnPFT_hrow = num2str(length(subject)+4);
        xlmnPFTrow = num2str(length(subject)+6);
        
        if length(incl_prot)>1
            oneprot = find(protocol.(PFT)==incl_prot(2));
        end
    end

    datime = datestr(now,'yymmddTHHMM');
    colon = find(datime == ':');
    space = find(datime == ' ');
    datime(colon) = '_';
    datime(space) = '_';
    
if group_prot == 2
    DRnm = 'DR';
    cumDosenm = 'cumDose';
    dpftnm  = 'dpft';
    dFEV1nm = 'dFEV1';
    dFVCnm  = 'dFVC';
elseif group_prot == 1
    DRnm = char(strcat('DR_',prot_num.(xlsfile)(grp)));
    cumDosenm = char(strcat('cumDose_',prot_num.(xlsfile)(grp)));
    dpftnm  = char(strcat('dpft_',prot_num.(xlsfile)(grp)));
    dFEV1nm = char(strcat('dFEV1_',prot_num.(xlsfile)(grp)));
    dFVCnm  = char(strcat('dFVC_',prot_num.(xlsfile)(grp)));
end

    % prompt={'Do you want the newly created Excel file to go here?                   If not, please enter the correct path'};
    % dlgTitle='Output Path';
    % numlines=1;
    default={strcat('C:\Sue\ODR\ODR_ModelParamRSQR_',xlsfile_b,'_',datime,'.xls')};

    % outputfile = inputdlg(prompt,dlgTitle,numlines,default);
    % File=char(outputfile);

    File=char(default);

    Excel = actxserver ('Excel.Application');
    if ~exist(File,'file')
        ExcelWorkbook = Excel.workbooks.Add;
        ExcelWorkbook.SaveAs(File)
        ExcelWorkbook.Close(false);
    end
    ExcelWorkbook = Excel.workbooks.Open(File);

%     xlswrite1(File,h_sub,DRnm,'B1');
    xlswrite1(File,h_onemean,DRnm,'B1');
%     xlswrite1(File,h_avg,DRnm,'AG1');
    xlswrite1(File,h_DR,DRnm,'A2');
    xlswrite1(File,h_protocol,DRnm,'A3');
    xlswrite1(File,DRxl,DRnm,'B3');
    xlswrite1(File,h_mnDR,DRnm,strcat('A',xlmnDR_hrow));
    xlswrite1(File,h_protocol,DRnm,strcat('A',xlmnDRrow));
    xlswrite1(File,MnDRxl,DRnm,strcat('B',xlmnDRrow));

%     xlswrite1(File,h_sub,cumDosenm,'B1');
    xlswrite1(File,h_onemean,cumDosenm,'B1');
%     xlswrite1(File,h_avg,cumDosenm,'AG1');
    xlswrite1(File,h_cumDose,cumDosenm,'A2');
    xlswrite1(File,h_protocol,cumDosenm,'A3');
    xlswrite1(File,cDosex1,cumDosenm,'B3');

    if length(incl_prot)>1
        xlswrite1(File,h_measure,dpftnm,'A1');
        xlswrite1(File,h_protocol(oneprot,:)',dpftnm,'D2');
        xlswrite1(File,h_sub',dpftnm,'A3');
        xlswrite1(File,{h_pft{1}},dpftnm,'D1');
        xlswrite1(File,dpft_msrxl(oneprot,:,1)',dpftnm,'D3');
        xlswrite1(File,h_avg,dpftnm,strcat('A',xlmnPFT_hrow));
        xlswrite1(File,h_all,dpftnm,strcat('A',xlmnPFTrow));
        xlswrite1(File,{h_pftRSQR_a{1}},dpftnm,'B2');
        xlswrite1(File,rsqr_dpft_a(:,:,1)',dpftnm,'B3');
        if pft ==2
            xlswrite1(File,h_protocol(oneprot,:)',dpftnm,'N2');
            xlswrite1(File,{h_pft{2}},dpftnm,'N1');
            xlswrite1(File,dpft_msrxl(oneprot,:,2)',dpftnm,'N3');
            xlswrite1(File,{h_pftRSQR_a{2}},dpftnm,'M2');
            xlswrite1(File,rsqr_dpft_a(:,:,2)',dpftnm,'M3');
        end
    end

%     xlswrite1(File,h_sub,dFEV1nm,'C1');
    xlswrite1(File,h_twomeans,dFEV1nm,'C1');
%     xlswrite1(File,h_avg,dFEV1nm,'AI1'); 
%     xlswrite1(File,h_all,dFEV1nm,'AK1');
    xlswrite1(File,h_RSQR_a,dFEV1nm,'A2');
    xlswrite1(File,rsqr_dpft_a(:,:,1),dFEV1nm,'C2');
    xlswrite1(File,h_RSQR_a_adj,dFEV1nm,'A3');
    xlswrite1(File,adjrsqr_dpft_a(:,:,1),dFEV1nm,'C3');
%     xlswrite1(File,h_RSQR_yint0,dFEV1nm,'A4');
%     xlswrite1(File,rsqr0cept_dpft_a(:,:,1),dFEV1nm,'C4');
 
    xlswrite1(File,h_param,dFEV1nm,'A6');
    xlswrite1(File,param_dpft(:,:,1),dFEV1nm,'C6');

    xlswrite1(File,h_MnStd,dFEV1nm,'A11');
    xlswrite1(File,baseline_dpft(:,:,1),dFEV1nm,'C11');

    xlswrite1(File,h_time,dFEV1nm,'B24');
    xlswrite1(File,time.FEV1(:,2),dFEV1nm,'B25');
    xlswrite1(File,h_model,dFEV1nm,'A24');
    xlswrite1(File,h_protocol,dFEV1nm,'A25');
    xlswrite1(File,dpft_modxl(:,:,1),dFEV1nm,'C25');
    xlswrite1(File,h_time,dFEV1nm,strcat('B',xlmsFEV_hrow));
    xlswrite1(File,time.FEV1(:,2),dFEV1nm,strcat('B',xlmsFEVrow));
    xlswrite1(File,h_measure,dFEV1nm,strcat('A',xlmsFEV_hrow));
    xlswrite1(File,h_protocol,dFEV1nm,strcat('A',xlmsFEVrow));
    xlswrite1(File,dpft_msrxl(:,:,1),dFEV1nm,strcat('C',xlmsFEVrow));

    if pft ==2
%         xlswrite1(File,h_sub,dFVCnm,'C1');
        xlswrite1(File,h_twomeans,dFVCnm,'C1');
%         xlswrite1(File,h_avg,dFVCnm,'AI1'); 
%         xlswrite1(File,h_all,dFVCnm,'AK1');
        xlswrite1(File,h_RSQR_a,dFVCnm,'A2');
        xlswrite1(File,rsqr_dpft_a(:,:,2),dFVCnm,'C2');
        xlswrite1(File,h_RSQR_a_adj,dFVCnm,'A3');
        xlswrite1(File,adjrsqr_dpft_a(:,:,2),dFVCnm,'C3');
%         xlswrite1(File,h_RSQR_yint0,dFVCnm,'A4');
%         xlswrite1(File,rsqr0cept_dpft_a(:,:,2),dFVCnm,'C4');

        xlswrite1(File,h_param,dFVCnm,'A6');
        xlswrite1(File,param_dpft(:,:,2),dFVCnm,'C6');

        xlswrite1(File,h_MnStd,dFVCnm,'A11');
        xlswrite1(File,baseline_dpft(:,:,2),dFVCnm,'C11');

        xlswrite1(File,h_time,dFVCnm,'B24');
        xlswrite1(File,time.FVC(:,2),dFVCnm,'B25');
        xlswrite1(File,h_model,dFVCnm,'A24');
        xlswrite1(File,h_protocol,dFVCnm,'A25');
        xlswrite1(File,dpft_modxl(:,:,2),dFVCnm,'C25');
        xlswrite1(File,h_time,dFVCnm,strcat('B',xlmsFEV_hrow));
        xlswrite1(File,time.FVC(:,2),dFVCnm,strcat('B',xlmsFEVrow));
        xlswrite1(File,h_measure,dFVCnm,strcat('A',xlmsFEV_hrow));
        xlswrite1(File,h_protocol,dFVCnm,strcat('A',xlmsFEVrow));
        xlswrite1(File,dpft_msrxl(:,:,2),dFVCnm,strcat('C',xlmsFEVrow));
    end
    % % 
    % 
    ExcelWorkbook.Save
    ExcelWorkbook.Close(false)  % Close Excel workbook.
    Excel.Quit;
    delete(Excel); 
end
end
