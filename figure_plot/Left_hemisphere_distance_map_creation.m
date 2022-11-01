clear all;clc;

path='D:\ZZZ\Manuscripts\Gradient_study\lifespan_gradient\codes\code_github\example\SX0004_CUIXIAOHUA_20200823_R_CUIXIAOHUA'
cd(path)
tem=spm_vol('Ring_distance_map.nii');

data=spm_read_vols(tem);

for i=1:size(data,1)/2
    for j=1:size(data,2)
        for k=1:size(data,3)
            
            data(i,j,k)=0;
        end
    end
end

tem.fname=['Left_',tem.fname];
spm_write_vol(tem,data);