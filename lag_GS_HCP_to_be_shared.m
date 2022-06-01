function lag_GS_HCP_to_be_shared(subno,prepost) %subno=subject no; i=day1 or day2
% upsampling method is adapted from Yunjie Tong's code:
% https://github.com/TonglabPurdue/Systematic-low-frequency-oscillation-in-fMRI/blob/master/Delay_Map_fMRI/Delay_map.m
 
    dname1=('directory_name');
    % cd (dname1);
for phase=1:2
    mname0= 'path_to_whole_brain_mask';
    mname= 'path_to_reference_signal_mask';
    mask0 = spm_read_vols(spm_vol(mname0));
    mask = spm_read_vols(spm_vol(mname));
    [dx, dy, dz] = size(mask0);
    if phase==1
        LR='LR';
    else
        LR='RL';
    end
    fname='path_to_data_file';
    data = spm_read_vols(spm_vol(fname));
    n = size(data,4); % time points
    hdr1 = spm_vol(fname);
    hdr1 = hdr1(1);
    masked =data.* (mask0>0);
    rdata = reshape(masked, [dx*dy*dz,n]);
    maskreshaped = reshape(mask0, [dx*dy*dz,1]);
    s = sum(maskreshaped>0);
    ts = rdata(maskreshaped>0,:);
    dts = (detrend(ts.','linear')).';
    masked1 =data.* (mask>0);
    rdata1 = reshape(masked1, [dx*dy*dz,n]);
    maskreshaped1 = reshape(mask, [dx*dy*dz,1]);
    ts2 = rdata1(maskreshaped1>0,:);
    dts2 = (detrend(ts2.','linear')).';
    meants = detrend(mean(dts2),'linear');
%% Compute TD and R
    TR=0.72;
    upsampleratio=0.15;  % adapted from Yunjie Tong's Delay_map.m
    ratioNtoOld=floor(TR/upsampleratio); %TR/ratioNtoOld = new temporal resolution
    xcorr_range=floor(6/(TR/ratioNtoOld)); % +/- 6 sec 
    ref_ts=detrend(oversample_ts(meants.',ratioNtoOld),'linear');
parfor j=1:s
     ts0 = (dts(j,:)).';
     if (ts0(1,1)~=0)
     ts1=detrend(oversample_ts(ts0,ratioNtoOld),'linear');
     cor_gs(j,1)=corr(ts1,ref_ts); % zerolag Pearson correlation coeffecients R
     [r,p]=xcorr(ts1,ref_ts,xcorr_range,'coeff');
     T=p(find(r==max(r))); %T = corresponding time delay 
     lag_gs(j,1)= T*(TR/ratioNtoOld); % Time Delay TD         
     else
     lag_gs(j,1)= 0;
     cor_gs(j,1)= 0;
     end
end
%% Output Images      
    data= zeros([dx, dy, dz]);
    data(mask0>0) = lag_gs;
    fname1 = sprintf('%s/%03d_%01d_%s_TD.nii', dname1,subno,prepost,LR);
    hdr1.fname = fname1;
    spm_write_vol(hdr1, data);
    
    data(mask0>0) = cor_gs;
    fname1 = sprintf('%s/%03d_%01d_%s_R_.nii', dname1,subno,prepost,LR);
    hdr1.fname = fname1;
    spm_write_vol(hdr1, data);
end
end
