clc, clear all, close all

global taxi
%
prompt = {'Include Filtered Air?(yes=1,no=2)','Subtract Filtered Air?(yes=1,no=2)','Subtract Filtered Air Mean or Protocol(s)?(Mn=1,Prot=2)'};
dlgTitle = 'Demographics';
def = {'1','2','2','2'}';
lineNo = 1;
file = inputdlg(prompt,dlgTitle,lineNo,def)
drawnow;

pftest      = '1';
taxi.FiltAir= str2num(file{1})
subtractFA  = str2num(file{2})
subtractMnOrProt = str2num(file{3})

xls_study   = {'AROPREDv1';'IMMPREDv1';'KOZPREDv1';'OSSPREDv1';'HR7LONGv1';'HR7PREDv1'}
xls_substudy = {'ARO','','','','','';'IMM','OZI','POPCOM','RESTOZ','','';'KOZ','','','','','';'OSS-2','OTH','','','','';...
    'PROLNGOZ','','','','','';'OZI2','PINOKIOZ1','PINOKIOZ2','POKOZ','REDUNDOZ2','RESTOZ2'}

FAprots  = [1,1,1,1,1,1]

pft = 1;

datime = datestr(now,'yymmddTHHMM');
colon = find(datime == ':');
space = find(datime == ' ');
datime(colon) = '_';
datime(space) = '_';

for study = 1:1    
    FAprot  = FAprots(study)

    xlsfile   = strcat('z4901_',xls_study{study})
    paramfile = strcat('UCD-EPAModDataV2.xls');      

    subgroup = '_all';

    cd C:\Sue\ODR
 
    BSA       = xlsread(strcat('C:\Sue\ODR\EPA\Mod_EPA_Data_fxdESS\',paramfile),'Est Parameters', 'I2:I1002');
    Time_orig = xlsread(strcat('EPA\Mod_EPA_Data_fxdESS\',xlsfile),'TIME', 'D2:AI290');
    
    O3_orig   = 
    C:\Sue\ODR\EPA\Mod_EPA_Data_fxdESS
    
    if taxi.FiltAir == 1
        if subtractFA == 2 
            xlsfile_b = xlsfile
            savename3 = strcat('ODR_',xls_study{study},subgroup,'_optparam.mat')
        elseif subtractFA == 1 
            if subtractMnOrProt == 1
                xlsfile_b = strcat(xls_study{study},'_PFTminusMnFA')
                savename3 = strcat('ODR_',xls_study{study},subgroup,'_PFTminusMnFA_optparam.mat')
            elseif subtractMnOrProt == 2
                xlsfile_b = strcat(xls_study{study},'_PFTminusFAprot')
                savename3= strcat('ODR_',xls_study{study},subgroup,'_PFTminusFA_optparam.mat')                    
            end
        end
    elseif taxi.FiltAir == 2
        if subtractFA == 2    
            xlsfile_b = strcat(xls_study{study},'_NoFAInMinErr')
            savename3 = strcat('ODR_',xlsfile,subgroup,'_NoFAinMinErr_optparam.mat')
        elseif subtractFA == 1 
            if subtractMnOrProt == 1
                xlsfile_b = strcat(xls_study{study},'_NoFAInMinErr_PFTminusMnFA')
                savename3 = strcat('ODR_',xlsfile,subgroup,'_NoFAinMinErr_PFTminusMnFA_optparam.mat')                    
            elseif subtractMnOrProt == 2
                xlsfile_b = strcat(xls_study{study},'_NoFAInMinErr_PFTminusFAprot') 
                savename3 = strcat('ODR_',xlsfile,subgroup,'_NoFAinMinErr_PFTminusFAprot_optparam.mat')                    
            end                
        end
    end
      
    default={strcat('C:\Sue\ODR\ODR_ModelParamRSQR_',xlsfile_b,'_',datime,'.xls')};

    File=char(default);

    Excel = actxserver ('Excel.Application');
    if ~exist(File,'file')
        ExcelWorkbook = Excel.workbooks.Add;
        ExcelWorkbook.SaveAs(File)
        ExcelWorkbook.Close(false);
    end
    ExcelWorkbook = Excel.workbooks.Open(File);     
    
    savename3 = char(savename3);
    load(savename3)  

    
%     if study < 6
%         substudy_I = strmatch(xls_substudy_b(study), Studyname_orig, 'exact');
%     end
            
    for substudy = 1:1        
        clear rsqr_dpft_a adjrsqr_dpft_a param_dpft mn_1st2hr_dpft std_1st2hr_dpft stdBelmn_2hr_dpft underBaseline_2hr VE O3
        clear under_std_1st2hr_dpft above_std_1st2hr_dpft mn_FA_dpft std_FA_dpft stdBelmn_FA_dpft underBaseline_FA pftmod_all_xls_subcol_c
        clear under_std_FA_dpft above_std_FA_dpft baseline_dpft dpft_modxl dpft_msrxl rsqr0cept_dpft_a incl_samples cumDose_subset 
        
%         if study == 6
            substudy_I = strmatch(xls_substudy(study,substudy), Studyname_orig, 'exact');
%         end 
        
        Baseline  = Baseline_orig(substudy_I,:);
        Studyname = Studyname_orig(substudy_I,:);
        VEex      = VEex_orig(substudy_I,:);
        O3        = O3_orig(substudy_I,:);
        FEV1      = FEV1_orig(substudy_I,:);
        Time      = Time_orig(substudy_I,:);

        protocol_numb = Baseline(:,2);
        incl_prot = unique(protocol_numb);
        subject = unique(Baseline(:,1));

        TIME   = Time(:,4:end);
        
        if taxi.FiltAir == 1
            if subtractFA == 2 
                savename  = strcat('ODR_',xls_study{study},'_',Studyname(1),'_optparam.mat')
                savename2 = strcat('ODR_',xls_study{study},subgroup,'_',Studyname(1),'_optparam.mat')
            elseif subtractFA == 1 
                if subtractMnOrProt == 1
                    savename  = strcat('ODR_',xls_study{study},'_',Studyname(1),'_PFTminusMnFA_optparam.mat')
                    savename2 = strcat('ODR_',xls_study{study},subgroup,'_',Studyname(1),'_PFTminusMnFA_optparam.mat')
                elseif subtractMnOrProt == 2
                    savename = strcat('ODR_',xls_study{study},'_',Studyname(1),'_PFTminusFAprot_optparam.mat')
                    savename2= strcat('ODR_',xls_study{study},subgroup,'_',Studyname(1),'_PFTminusFA_optparam.mat')
                end
            end
        elseif taxi.FiltAir == 2
            if subtractFA == 2    
                savename  = strcat('ODR_',xlsfile,'_',Studyname(1),'_NoFAinMinErr_optparam.mat')
                savename2 = strcat('ODR_',xlsfile,subgroup,'_',Studyname(1),'_NoFAinMinErr_optparam.mat')
            elseif subtractFA == 1 
                if subtractMnOrProt == 1
                    savename  = strcat('ODR_',xlsfile,'_',Studyname(1),'_NoFAinMinErr_PFTminusMnFA_optparam.mat')
                    savename2 = strcat('ODR_',xlsfile,subgroup,'_',Studyname(1),'_NoFAinMinErr_PFTminusMnFA_optparam.mat')
                elseif subtractMnOrProt == 2
                    savename  = strcat('ODR_',xlsfile,'_',Studyname(1),'_NoFAinMinErr_PFTminusFAprot_optparam.mat')
                    savename2 = strcat('ODR_',xlsfile,subgroup,'_',Studyname(1),'_NoFAinMinErr_PFTminusFAprot_optparam.mat')
                end                
            end
        end            
                
        rand('state',sum(100*clock));

        whichstats = {'yhat','rsquare','adjrsquare','beta'};
        PFTS  = {'FEV1'};
        dPFTS = {'dFEV1'};
        measure = {'VEex','O3','FEV1'}; 

        clear rsqr_dpft_a adjrsqr_dpft_a param_dpft mn_1st2hr_dpft std_1st2hr_dpft stdBelmn_2hr_dpft underBaseline_2hr incl_samples
        clear under_std_1st2hr_dpft above_std_1st2hr_dpft mn_FA_dpft std_FA_dpft stdBelmn_FA_dpft underBaseline_FA
        clear under_std_FA_dpft above_std_FA_dpft baseline_dpft dpft_modxl dpft_msrxl rsqr0cept_dpft_a

        clear dpft_msrd TimePt dpftmod_all dpftmod_all_a dpftmod_all_xls doserate cumDose Mndoserate cnt_sub DR_2D MnDR_2D 
        clear deltapft cnt avg_deltapft avg_sig1 R_sq_0cept_all R_sq_0cept_avg R_sq_avg R_sq_all VE_A 
        clear RSQR_dpft_s1i0 yhat_dpft_s1i0 RSQR_dpft yhat_dpft R_sq_0cept R_sq RSQR_dpft_int0 yhat_dpft_int0 num_samples
        
        clear taxi.msrd_dpft taxi.dpft taxi.MnDR taxi.DR taxi.pred_dpft stats msrmod DR_2D_a MnDR_2D_a cumD
        clear TimePt FA_sample taxi.incl_prot Coefficients_0cept S_err_0cept XTXI_0cept R_sq_0cept_a
        clear Coefficients S_err XTXI R_sq_a O3_A cumDose_subset SubProt RSQR_dpft_a adjRSQR_dpft_a ODR_b
        
        PFT = char(PFTS)
        PFT2 = strcat(PFT)

        msr.FEV1 = []
        msr.O3   = [];
        msr.VEex = [];

        cnt = 1;

        for w = 1:length(measure)
            msrmt= char(measure(w));
            msr1.(msrmt) = eval(msrmt);
            if w==3
                for x = 4:size(msr1.(msrmt),2)
                    msr.(msrmt)(:,x-3) = 100*(msr1.(msrmt)(:,x)-msr1.(msrmt)(:,4))./msr1.(msrmt)(:,4);   %Percent change FEV1 from baseline
                end
            else
                msr.(msrmt)= msr1.(msrmt)(:,4:end);
            end                
            num_samples.(msrmt)= length(protocol_numb); 
        end

        msr.VErest = 7.6119*Baseline(:,6);  %Schelegle constant 7.6119 = slope of BSA vs VErest for Schelegle test subjects

        msr_subcol_c.(PFT)=[]
        
        O3_A(:,1:TIME(:,1)+1) = repmat(msr.O3(:,1),1,length(1:TIME(:,1)+1));

        startI = find(TIME(1,1:20)==0); 
        
        if TIME(:,5)>0
            if startI == 5
                for n = 1:size(msr.VEex,2)
                    if n<size(msr.VEex,2)
                        VE_A(:,TIME(:,n+4)+1:TIME(:,n+12)) = repmat(msr.VEex(:,n),1,length(TIME(:,n+4)+1:TIME(:,n+12)));
                        VE_A(:,TIME(:,n+12)+1:TIME(:,n+5)) = repmat(msr.VErest,1,length(TIME(:,n+12)+1:TIME(:,n+5))); 
    %                     O3_A(:,TIME(:,n+4)+1:TIME(:,n+12)) = repmat(msr.O3(:,n),1,length(TIME(:,n+4)+1:TIME(:,n+12)));
                    elseif n==size(msr.VEex,2)
                        VE_A(:,TIME(:,n+4)+1:TIME(:,n+12)) = repmat(msr.VEex(:,n),1,length(TIME(:,n+4)+1:TIME(:,n+12)));
                        VE_A(:,TIME(:,n+12)+1:TIME(:,2)+1) = repmat(msr.VErest,1,length(TIME(:,n+12)+1:TIME(:,2)+1));
                        if TIME(:,2) > TIME(:,1)
                            O3_A(:,TIME(:,1)+2:TIME(:,2)+1)  = repmat(msr.O3(:,n),1,length(TIME(:,1)+2:TIME(:,2)+1));
                        end
                    end
                end
            elseif startI == 13
                for n = 1:size(msr.VEex,2)
                    VE_A(:,TIME(:,n+4)+1:TIME(:,n+13)) = repmat(msr.VEex(:,n),1,length(TIME(:,n+4)+1:TIME(:,n+13)));
                    VE_A(:,TIME(:,n+12)+1:TIME(:,n+4)) = repmat(msr.VErest,1,length(TIME(:,n+12)+1:TIME(:,n+4))); 
    %                 O3_A(:,TIME(:,n+12)+1:TIME(:,n+4)) = repmat(msr.O3(:,n),1,length(TIME(:,n+12)+1:TIME(:,n+4)));
                end
                VE_A(:,TIME(:,n+13)+1:TIME(:,2)+1) = repmat(msr.VErest,1,length(TIME(:,n+13)+1:TIME(:,2)+1));
                if TIME(:,2) > TIME(:,1)
                    O3_A(:,TIME(:,1)+2:TIME(:,2)+1)  = repmat(msr.O3(:,end),1,length(TIME(:,1)+2:TIME(:,2)+1));
                end
            end
        else
            VE_A(:,1:TIME(:,1)+1) = repmat(msr.VEex(:,1),1,length(1:TIME(:,1)+1)); 
%             if TIME(:,2) > TIME(:,1)
%                 O3_A(:,TIME(:,1)+2:TIME(:,2)+1)  = repmat(msr.O3(:,n),1,length(TIME(:,1)+2:TIME(:,2)+1));
%             end
        end

%         TimePt = [0;TIME(1,:)'];
        FEV1indx = TIME(1,21:end)'+1;
        FEV1indx = FEV1indx(~isnan(FEV1indx));

        msr_subcol_b.(PFT)= []; 
           
        for s = 1:length(subject)

            clear VE O3 taxi.msrd_dpft taxi.dpft taxi.MnDR taxi.DR taxi.pred_dpft stats msrmod DR_2D_a MnDR_2D_a cumD
            clear FA_sample taxi.incl_prot Coefficients_0cept S_err_0cept XTXI_0cept R_sq_0cept_a
            clear Coefficients S_err XTXI R_sq_a 

            msr_subcol_a.(PFT) = [];
            msr_subcol_b.(PFT) = [];
            dpftmod_all_xls_subcol_a.(PFT) = [];
            dpftmod_all_xls_subcol_b.(PFT) = [];    

            if cnt ==1
                s1 = s;
            end

            sub = strcat('sub',num2str(subject(s))) 
            sub_I = find(Baseline(:,1) == subject(s));

            if FAprot == 1 
                FA_sample    = find(protocol_numb==1 & Baseline(:,1) == subject(s));
                if isempty(FA_sample)
%                     FA_surrogate = 0;
                    mn_FA_dpft(s,pft)  = NaN;
                    std_FA_dpft(s,pft) = NaN;
                else
                    mn_FA_dpft(s,pft)  = nanmean(msr.(PFT2)(FA_sample,:),2);
                    std_FA_dpft(s,pft) = nanstd(msr.(PFT2)(FA_sample,:),0,2);
                end
            end 

            VE = VE_A(sub_I,:)';
            O3 = O3_A(sub_I,:)';
            
            dpft_msrd(:,:,s) = repmat(NaN,length(FEV1indx),length(incl_prot));
            doserate(:,:,s)  = repmat(NaN,size(O3_A,2),length(incl_prot));
            cumDose(:,:,s)   = repmat(NaN,size(O3_A,2),length(incl_prot));
            Mndoserate(:,:,s)= repmat(NaN,size(O3_A,2),length(incl_prot));          
                                   
            dpft_msrd(1:length(FEV1indx),protocol_numb(sub_I),s) = msr.(PFT)(sub_I,:)';

            t = 0:TIME(:,2);
            
            doserate(1:size(VE,1),protocol_numb(sub_I),s)             = VE.*O3*1.96;          %1.96 = ppm to ug/liter ozone conversion
            taxi.t = repmat(t',1,size(doserate,2));
            cumDose(1:size(doserate(:,:,s),1),:,s) = cumsum(doserate(:,:,s).*[ones(1,size(doserate(:,:,s),2));diff(taxi.t)]);
            Mndoserate(2:size(taxi.t,1),:,s)       = cumDose(2:end,:,s)./taxi.t(2:end,:);

            Mndoserate(1,:,s)= doserate(1,:,s);

            x_U(2) = max(max(cumDose(:,:,s)))

            if x_U(2) < x_L(2)
                x_U(2)= x_L(2) + 0.001;
            end

            if subtractFA == 1
                if subtractMnOrProt == 1
                    msr.(PFT)(sub_I,:)           = msr.(PFT)(sub_I,:) - mn_FA_dpft(s,pft);
%                     msr.(PFT)(sub_I,1) = 0
                elseif subtractMnOrProt== 2 
                    if isempty(FA_sample)==0
                       for rowI = 1:length(sub_I)
                           msr_b.(PFT)(sub_I(rowI),:)= msr.(PFT)(sub_I(rowI),:) - msr.(PFT)(sub_I(1),:);
                       end
                       msr.(PFT)(sub_I,:) = msr_b.(PFT)(sub_I,:);
                    end
                end
            end 
            
            for n = 1:size(cumDose,2)
%                 if n<size(cumDose,2)
                    cumD(:,n)      = cumDose(FEV1indx,n,s);
                    DR_2D_a(:,n)   = doserate(FEV1indx,n,s);
                    MnDR_2D_a(:,n) = Mndoserate(FEV1indx,n,s);    
%                 else
%                     cumD(:,n)      = cumDose(FEV1indx,n,s);
%                     DR_2D_a(:,n)   = doserate(FEV1indx,n,s);
%                     MnDR_2D_a(:,n) = Mndoserate(FEV1indx,n,s);    
%                 end           
            end
            taxi.FEV1indx = repmat(FEV1indx,1,size(doserate,2));

            taxi.DR   = doserate(:,:,s);
            taxi.MnDR = Mndoserate(:,:,s);

            taxi.msrd_dpft = dpft_msrd(:,:,s);

            taxi.err = 6000;
            taxi.fit = 0;
            
            substudyname = char(Studyname(1));
            dash         = find(substudyname == '-');
            space        = find(substudyname == ' ');
            substudyname(dash) = '_';
            substudyname(space) = '_';

                        
            y_0  = taxi.msrd_dpft(1,1);
            
            %--------------------------------------

            FEV1indx = taxi.FEV1indx;

            A     = taxi.A;
            DOS   = taxi.DOS;
            tau   = taxi.tau;
            slope = taxi.slope;
            DR    = taxi.DR;
            mnDR  = taxi.MnDR;
            t     = taxi.t;

            y = zeros(size(DR));
            ycum = y;
            for n =2:length(t)
                sig1(n,:)     = DR(n,:)./(1+exp(-slope*(t(n,:)-(DOS./mnDR(n,:)))));
                ceoss         = squeeze(sig1(n,:)).*tau/log(2);
                y(n,:)        = (squeeze(ceoss) - squeeze(ycum(n-1,:))).*(1-exp(-log(2)./tau));
                ycum(n,:)     = ycum(n-1,:) + y(n,:);
                deltapft(n,:) = A*ycum(n,:);
            end
            for m = 1:size(FEV1indx,2)
                dpft_min(:,m) = deltapft(FEV1indx(:,m),m);
            end

            if taxi.FiltAir == 1
                ymod   = squeeze(dpft_min(1:end,1:end,:));
            elseif taxi.FiltAir == 2
                ymod   = squeeze(dpft_min(1:end,:,:));
            end    
            
            %--------------------------------------------------------------

            pred_PFT     = taxi.pred_dpft(:,:);
            msrd_PFT     = taxi.msrd_dpft(:,:);


        %     save (paramfile, 'ODR')
%             save (savename2, 'ODR')
            save (savename3, 'ODR')

            t60 = taxi.FEV1indx;

            dpftmod_all_a(sub_I,:) = taxi.pred_dpft(:,protocol_numb(sub_I))';

            allsub = find(Baseline(:,1)>= subject(s1) & Baseline(:,1));
            numbprot = unique(Baseline(allsub,2));
            subprotI = 1;

            for p = 1:numbprot(end)              
                SubProt_a       = find(Baseline(sub_I,2) == p); 
                if isempty(SubProt_a)==1
                    dpftmod_all_xls_subcol_a.(PFT) = repmat(NaN,size(msr.(PFT),2),1);
                    msr_subcol_a.(PFT)             = repmat(NaN,size(dpftmod_all_a,2),1);
                else
                    dpftmod_all_xls_subcol_a.(PFT) = [dpftmod_all_a(sub_I(subprotI),:)]';
                    msr_subcol_a.(PFT)             = [msr.(PFT)(sub_I(subprotI),:)]';  
                    subprotI = subprotI + 1;
                end
                dpftmod_all_xls_subcol_b.(PFT) = [dpftmod_all_xls_subcol_b.(PFT);dpftmod_all_xls_subcol_a.(PFT)];
                msr_subcol_b.(PFT)             = [msr_subcol_b.(PFT);msr_subcol_a.(PFT)];                    
            end

            pftmod_all_xls_subcol_c.(PFT)(:,s) = dpftmod_all_xls_subcol_b.(PFT);
            msr_subcol_c.(PFT)(:,s) = msr_subcol_b.(PFT);                

            DR_2D(:,cnt)           = DR_2D_a(:);
            MnDR_2D(:,cnt)         = MnDR_2D_a(:);
            cumDose_subset(:,cnt)  = cumD(:);

            deltapft(:,:,s) = taxi.pred_deltapft;

            mn_1st2hr_dpft(cnt,pft)  = nanmean(nanmean(dpft_msrd(1:3,:,s),1));
            std_1st2hr_dpft(cnt,pft) = nanmean(nanstd(dpft_msrd(1:3,:,s),0,1));

            stdBelmn_2hr_dpft(s,pft) = mn_1st2hr_dpft(cnt,pft)-std_1st2hr_dpft(cnt,pft);

            underBaseline_2hr(cnt,pft)     = 100*length(find(dpft_msrd(:,:,s)<stdBelmn_2hr_dpft(s)))./num_samples.FEV1;
            under_std_1st2hr_dpft(cnt,pft) = 100*length(find(dpft_msrd(:,:,s)<-std_1st2hr_dpft(cnt)))./num_samples.FEV1;

            above_std_1st2hr_dpft(cnt,pft) = 100*length(find(dpft_msrd(:,:,s)>std_1st2hr_dpft(cnt,pft)))./num_samples.FEV1;

            stdBelmn_FA_dpft(s,pft) = mn_FA_dpft(s,pft)-std_FA_dpft(s,pft);
            underBaseline_FA(cnt,pft) = 100*length(find(dpft_msrd(:,:,s)<stdBelmn_FA_dpft(cnt,pft)))./num_samples.FEV1;
            under_std_FA_dpft(cnt,pft)= 100*length(find(dpft_msrd(:,:,s)<-std_FA_dpft(s,pft)))./num_samples.FEV1;
            above_std_FA_dpft(cnt,pft)= 100*length(find(dpft_msrd(:,:,s)>std_FA_dpft(s,pft)))./num_samples.FEV1;

            if (s < 6 || s>90) || (study~=6 && study~=2 && study~=4)
                figure
                 subplot(2,1,1)
                 plot(cumD(:),taxi.msrd_dpft(:),'b*',cumD(:),taxi.pred_dpft(:), 'ro',...
                     cumD(:),repmat(mn_1st2hr_dpft(cnt,pft),length(cumD(:))),'g--',cumD(:),repmat(stdBelmn_2hr_dpft(s,pft),length(cumD(:))),'c'),...
                     grid on, title(sub),xlabel('Cumulative Dose (ppm?)'),ylabel([PFT,' (% of initial)']),legend('Meas','Mod','Mean 2h','Std 2h' )
                 subplot(2,1,2)
                 plot(t60(:),taxi.msrd_dpft(:),'b*',t60(:),taxi.pred_dpft(:), 'ro'),...
                     grid on, title(sub),xlabel('Time (min)'),ylabel([PFT,' (% of initial)']),legend('Meas','Mod')
            end
            all_s(cnt) = s;
            cnt = cnt +1;
        end

        clear VE O3 taxi.msrd_dpft taxi.dFVC taxi.MnDR taxi.DR taxi.pred_dpft taxi.pred_dFVC stats msrmod
        clear dpftmsr_all dpftmod_all

        sub_I_all = find(Baseline(:,1)>= subject(s1) & Baseline(:,1) <= subject(s));

        ODR.(PFT).(substudyname).A_all_avg     = nanmean(ODR.(PFT).(substudyname).A(1:cnt-1));
        ODR.(PFT).(substudyname).DOS_all_avg   = nanmean(ODR.(PFT).(substudyname).DOS(1:cnt-1));
        ODR.(PFT).(substudyname).tau_all_avg   = nanmean(ODR.(PFT).(substudyname).tau(1:cnt-1));
        ODR.(PFT).(substudyname).slope_all_avg = nanmean(ODR.(PFT).(substudyname).slope(1:cnt-1));

        dpftmsr_all(sub_I_all,:) = msr.(PFT)(sub_I_all,:);
        dpftmod_all(sub_I_all,:) = dpftmod_all_a(sub_I_all,:);




        ODR.(PFT).(substudyname).A_avg    = x0_pft(1);
        ODR.(PFT).(substudyname).DOS_avg  = x0_pft(2);
        ODR.(PFT).(substudyname).tau_avg  = x0_pft(3);
        ODR.(PFT).(substudyname).slope_avg= x0_pft(4);

        
        save (savename3, 'ODR')

        
        avg_DR     = nanmean(DR_2D,2);
        avg_MnDR   = nanmean(MnDR_2D,2);

        avg_mn_1st2hr_dpft  = nanmean(mn_1st2hr_dpft(:,pft));
        avg_std_1st2hr_dpft = nanmean(std_1st2hr_dpft(:,pft));

        avg_stdBelmn_2hr_dpft = avg_mn_1st2hr_dpft-avg_std_1st2hr_dpft;

        avg_underBaseline_2hr     = 100*length(find(avgMsr_dpft<avg_stdBelmn_2hr_dpft))./num_samples.FEV1;
        avg_under_std_1st2hr_dpft = 100*length(find(avgMsr_dpft<-avg_std_1st2hr_dpft))./num_samples.FEV1;

        avg_above_std_1st2hr_dpft = 100*length(find(avgMsr_dpft>avg_std_1st2hr_dpft))./num_samples.FEV1;

        avg_mn_FA_dpft  = nanmean(mn_FA_dpft(:,pft));
        avg_std_FA_dpft = nanstd(std_FA_dpft(:,pft));

        avg_stdBelmn_FA_dpft = avg_mn_FA_dpft-avg_std_FA_dpft;

        avg_underBaseline_FA  = 100*length(find(avgMsr_dpft<avg_stdBelmn_FA_dpft))./num_samples.FEV1;
        avg_under_std_FA_dpft = 100*length(find(avgMsr_dpft<-avg_std_FA_dpft))./num_samples.FEV1;

        avg_above_std_FA_dpft = 100*length(find(avgMsr_dpft>avg_std_FA_dpft))./num_samples.FEV1;

        baseline_dpft(:,:,pft) = [mn_1st2hr_dpft(:,pft)',NaN,avg_mn_1st2hr_dpft;std_1st2hr_dpft(:,pft)',NaN,avg_std_1st2hr_dpft;...
            mn_FA_dpft(:,pft)',NaN,avg_mn_FA_dpft;std_FA_dpft(:,pft)',NaN,avg_std_FA_dpft;...
            underBaseline_2hr(:,pft)',NaN,avg_underBaseline_2hr;underBaseline_FA(:,pft)',NaN,avg_underBaseline_FA;...
            repmat(NaN,1,length(under_std_1st2hr_dpft(:,pft)')+2);under_std_1st2hr_dpft(:,pft)',NaN,avg_under_std_1st2hr_dpft;under_std_FA_dpft(:,pft)',NaN,avg_under_std_FA_dpft;...
            repmat(NaN,1,length(above_std_1st2hr_dpft(:,pft)')+2);above_std_1st2hr_dpft(:,pft)',NaN,avg_above_std_1st2hr_dpft;above_std_FA_dpft(:,pft)',NaN,avg_above_std_FA_dpft];

        figure
         subplot(2,1,1)
         plot(avgCumDose(:),avgMsr_dpft(:),'b*',avgCumDose(:),avgPred_dpft(:), 'ro'),...
%              avgCumDose(:),repmat(avg_mn_1st2hr_dpft,length(avgCumDose(:)),1),'g--',avgCumDose(:),repmat(avg_stdBelmn_2hr_dpft,length(avgCumDose(:)),1),'c'),...
             grid on, title('Mean'),xlabel('Cumulative Dose (ppm?)'),ylabel([PFT,' (% of initial)']),legend('Meas','Mod')
         subplot(2,1,2)
         plot(t60(:),avgMsr_dpft','b*',t60(:),avgPred_dpft','ro'),...
             grid on, title('Mean'),xlabel('Time (min)'),ylabel([PFT,' (% of initial)']),legend('Meas','Mod')

        if length(subject) <21
            col_havg = strcat(sprintf('%c','A'+length(subject)+2));
            col_havg_b = strcat(sprintf('%c','A'+length(subject)+3));
            col_hall  = strcat(sprintf('%c','A'+length(subject)+5));
        elseif length(subject)>20 && length(subject)<23
            col_havg = strcat(sprintf('%c','A'+length(subject)+2));
            col_havg_b = strcat(sprintf('%c','A'+length(subject)+3));
            col_hall  = strcat(sprintf('%c','A'),sprintf('%c','A'+length(subject)-21)); 
        elseif length(subject)==23
            col_havg = strcat(sprintf('%c','A'+length(subject)+2));
            col_havg_b = strcat(sprintf('%c','A'),sprintf('%c','A'+length(subject)-23));
            col_hall  = strcat(sprintf('%c','A'),sprintf('%c','A'+length(subject)-21)); 
        elseif length(subject)<46
            col_havg = strcat(sprintf('%c','A'),sprintf('%c','A'+length(subject)-24));
            col_havg_b = strcat(sprintf('%c','A'),sprintf('%c','A'+length(subject)-23));
            col_hall  = strcat(sprintf('%c','A'),sprintf('%c','A'+length(subject)-21)); 
        elseif length(subject)>45 && length(subject)<48
            col_havg = strcat(sprintf('%c','A'),sprintf('%c','A'+length(subject)-24));
            col_havg_b = strcat(sprintf('%c','A'),sprintf('%c','A'+length(subject)-23));
            col_hall  = strcat(sprintf('%c','B'),sprintf('%c','A'+length(subject)-47)); 
        elseif length(subject)==48
            col_havg = strcat(sprintf('%c','A'),sprintf('%c','A'+length(subject)-24));
            col_havg_b = strcat(sprintf('%c','B'),sprintf('%c','A'+length(subject)-48));
            col_hall  = strcat(sprintf('%c','B'),sprintf('%c','A'+length(subject)-46)); 
        elseif length(subject)<71
            col_havg = strcat(sprintf('%c','B'),sprintf('%c','A'+length(subject)-49));
            col_havg_b = strcat(sprintf('%c','B'),sprintf('%c','A'+length(subject)-48));
            col_hall  = strcat(sprintf('%c','B'),sprintf('%c','A'+length(subject)-46)); 
        elseif length(subject)>70 && length(subject)<73
            col_havg = strcat(sprintf('%c','B'),sprintf('%c','A'+length(subject)-49));
            col_havg_b = strcat(sprintf('%c','B'),sprintf('%c','A'+length(subject)-48));
            col_hall  = strcat(sprintf('%c','C'),sprintf('%c','A'+length(subject)-71)); 
        elseif length(subject)==73
            col_havg = strcat(sprintf('%c','B'),sprintf('%c','A'+length(subject)-49));
            col_havg_b = strcat(sprintf('%c','C'),sprintf('%c','A'+length(subject)-73));
            col_hall  = strcat(sprintf('%c','C'),sprintf('%c','A'+length(subject)-71)); 
        else 
            col_havg = strcat(sprintf('%c','C'),sprintf('%c','A'+length(subject)-74));
            col_havg_b = strcat(sprintf('%c','C'),sprintf('%c','A'+length(subject)-73));
            col_hall  = strcat(sprintf('%c','C'),sprintf('%c','A'+length(subject)-71));           
        end

        h_sub        = subject';
        sub_cell     = num2cell([h_sub,NaN]);
        twomeans     = {'MeanData','','MeanParam'};
        onemean      = {'MeanData'};
        h_twomeans   = [sub_cell,twomeans];
        h_onemean    = [sub_cell,onemean];
        h_param      = {'A';'DOS';'tau';'slope'};
        h_RSQR_yint0  = {'RSQR_yint0'};
        h_RSQR_a     = {'RSQR_trueCoeff'};
        h_RSQR_a_adj = {'RSQRadj_trueCoeff'};
        h_pftRSQR_a{pft}= strcat(PFT,'_RSQR');
        h_pft{pft}   = PFT;
        h_MnStd      = {'Mean 1st 2h';'Std 1st 2h';'Mean FA';'Std FA';'%<(Mn-Std) 2h';'%<(Mn-Std) FA';...
            ' ';'%<-Std 1st 2h';'%<-Std FA';' ';'%>Std 1st 2h';'%>Std FA'};

        times       = repmat([0:length(FEV1indx)-1]',length(substudy_I),1);
        studyname_a = Studyname_orig(substudy_I);
        studyname_b = repmat(studyname_a(1),length(times),1);
        time_FEV1   = repmat((FEV1indx-1),size(TIME,1),1);
        
        h_protocol = strcat(studyname_b,'_',num2str(times));
        h_model   = {'Modeled'};
        h_measure = {'Measured'};
        h_DR   = {'DR'};
        h_mnDR = {'mnDR'};
        h_cumDose = {'Cumulative Dose'};

        h_time = {'Time(min)'};
        
        h_msr = {'%change FEV1'};

        NaN_fill = repmat(NaN,size(avg_MnDR,1),1);
        NaN_fill_2 = repmat(NaN,length(avgMsr_dpft_xls),1);

        DRxl = [DR_2D,NaN_fill,avg_DR];
        MnDRxl = [MnDR_2D,NaN_fill,avg_MnDR];
        cDosex1 = [cumDose_subset,NaN_fill,avgCumDose(:)];

        dpft_modxl(:,:,pft) = [pftmod_all_xls_subcol_c.(PFT),NaN_fill_2,avgPred_dpft];
        dpft_msrxl(:,:,pft) = [msr_subcol_c.(PFT),NaN_fill_2,avgMsr_dpft_xls];

        xlmnDR_hrow = num2str(size(DRxl,1)+4);
        xlmnDRrow = num2str(size(DRxl,1)+5);

        xlmsFEV_hrow = num2str(size(DRxl,1)+26);
        xlmsFEVrow = num2str(size(DRxl,1)+27);

        xlmnPFT_hrow = num2str(length(subject)+4);
        xlmnPFTrow = num2str(length(subject)+6);

        if length(incl_prot)>1
             oneprot = find(dpft_msrxl(:,1)==0);
        end

        studyname_c = char(studyname_b(1));

        DRnm = strcat('DR_',studyname_c);
        cumDosenm = strcat('cumDose_',studyname_c);
        dpftnm  = strcat('dpft_',studyname_c);
        dFEV1nm = strcat('dFEV1_',studyname_c);
        dFVCnm  = strcat('dFVC_',studyname_c);

        xlswrite1(File,h_onemean,DRnm,'B1');        
        xlswrite1(File,h_DR,DRnm,'A2');
        xlswrite1(File,h_protocol(1:size(DRxl,1)),DRnm,'A3');
        xlswrite1(File,DRxl,DRnm,'B3');
        xlswrite1(File,h_mnDR,DRnm,strcat('A',xlmnDR_hrow));
        xlswrite1(File,h_protocol(1:size(DRxl,1)),DRnm,strcat('A',xlmnDRrow));
        xlswrite1(File,MnDRxl,DRnm,strcat('B',xlmnDRrow));

        xlswrite1(File,h_onemean,cumDosenm,'B1');        
        xlswrite1(File,h_cumDose,cumDosenm,'A2');
        xlswrite1(File,h_protocol(1:size(cDosex1,1)),cumDosenm,'A3');
        xlswrite1(File,cDosex1,cumDosenm,'B3');

        xlswrite1(File,h_msr,dFEV1nm,'A1');
        xlswrite1(File,h_twomeans,dFEV1nm,'C1');

        xlswrite1(File,h_param,dFEV1nm,'A6');
        xlswrite1(File,param_dpft(:,:,1),dFEV1nm,'C6');

        xlswrite1(File,h_MnStd,dFEV1nm,'A11');
        xlswrite1(File,baseline_dpft(:,:,1),dFEV1nm,'C11');

        xlswrite1(File,h_time,dFEV1nm,'B24');
        xlswrite1(File,time_FEV1(1:size(dpft_modxl,1)),dFEV1nm,'B25');
        xlswrite1(File,h_model,dFEV1nm,'A24');
        xlswrite1(File,h_protocol(1:size(dpft_modxl,1)),dFEV1nm,'A25');
        xlswrite1(File,dpft_modxl(:,:,1),dFEV1nm,'C25');
        xlswrite1(File,h_time,dFEV1nm,strcat('B',xlmsFEV_hrow));
        xlswrite1(File,time_FEV1(1:size(dpft_modxl,1)),dFEV1nm,strcat('B',xlmsFEVrow));
        xlswrite1(File,h_measure,dFEV1nm,strcat('A',xlmsFEV_hrow));
        xlswrite1(File,h_protocol(1:size(dpft_modxl,1)),dFEV1nm,strcat('A',xlmsFEVrow));
        xlswrite1(File,dpft_msrxl(:,:,1),dFEV1nm,strcat('C',xlmsFEVrow));

        ExcelWorkbook.Save
    end
    ExcelWorkbook.Close(false)  % Close Excel workbook.
    Excel.Quit;
    delete(Excel);
end
