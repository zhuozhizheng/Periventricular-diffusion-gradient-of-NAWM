clear all;clc;
tic;
datapath='G:\Life_span_gradient\T1_NODDI';
savepath='G:\Life_span_gradient\results\Diffusion_value';
cd(datapath);f1=dir;f=f1(3:end);

for i=1:length(f)
    i
    mask_path=[datapath filesep f(i).name];
    cd(mask_path);
    
    fmask=dir('Ring_distance_map*.nii');
    if length(fmask)>0
        tem=spm_vol([mask_path filesep fmask.name]);
        mask=spm_read_vols(tem);
        
        lesion=dir('wFLAIR_2_T1_ples*.nii');
        
        if length(lesion)>0
        tem=spm_vol([mask_path filesep lesion.name]);
        mask_exclude=spm_read_vols(tem);
        mask_exclude(mask_exclude>0.01)=1;
        SE=[0 1 0;1 1 1;0 1 0];
        mask_exclude=imdilate(mask_exclude,SE);
        mask_exclude=imdilate(mask_exclude,SE);
        mask(mask_exclude>=0)=0;
        end
        
        [m1 m2 m3]=size(mask);
        tem_mask=mask;tem_mask=reshape(mask,m1*m2*m3,1);
        
        NODDI_path=mask_path;
        cd(NODDI_path);nii=dir('wB0_2_T1_FIT*nii');
        if length(nii)>0
        for j1=1:length(nii)
            tem=spm_vol([NODDI_path filesep nii(j1).name]);
            data=spm_read_vols(tem);
            tem_data=data;tem_data=reshape(data,m1*m2*m3,1);
            tem_mask(isnan(tem_data))=0;
            
            for j2=1:max(max(tem_mask)) 
                value=mean(tem_data(tem_mask==j2));
                list{i,1,j1}=f(i).name;
                list{i,j2+1,j1}=value;
            end
        end
        end
    else
        continue;
    end
end

cd(savepath);
for i=1:length(nii);
    tem_list=list(:,:,i);
    xlswrite('New_NAWM_Distance_Results.xlsx',tem_list,i);
end

           