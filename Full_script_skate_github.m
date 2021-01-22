% This script will find bottom temperatures of climate models and add
% twenty year averages in the  2050 and 2100 time periods to the 
% high-resolution ROMS model to show climate model outputs in 
% higher resolution using the delta method. Assuming model data
% is separated into historical and RCP scenarios, need to input file
% path and beginning and end target years. 
% Will then send climate variables to R, run script, and returen thorny
% skate distributions and abundances

% Specify RCP Scenario and target variable.
RCP = RCP85; %RCP45 RCP85 TCS 
target_var = 1; % 1=thetao, 2=so

% Specify historical(target1&2) and RCP(target3&4) years.
target1 = 1955;
target2 = 2005;
target3 = 2006;
target4 = 2100;

% Climate models used
models = {'CCCma' ,'CMCC','LASG',     'MIROC', 'MPI',       'NOAA',      'ensemble_average';
          'CanESM','CMS', 'FGOALS-g2','MIROC5','MPI-ESM-MR','GFDL_CM2.6','ensemble'};
      
mn = 1:5; % Model number. GFDL_CM2.6 is a different resolution. Should run separately with TCS
%
q = 146; o = 331; p = 401; % Dimensions of years, WOA data
%
longabun_all = nan([length(mn) 4 q o p]);
shortabun_all = nan([length(mn) 4 q o p]);
llabun_all = nan([length(mn) 4 q o p]);
shortpa_all = nan([length(mn) 4 q o p]);
longpa_all = nan([length(mn) 4 q o p]);
%}

for modelnumber = mn
    inst = models{1,modelnumber};
    model = models{2,modelnumber}; 
    histo = 'historical';
    
% This will need to be run to extract climate model variables. Can be
% skipped if those are saved later (much faster to save and reload later)
    for variable = target_var 
        if variable == 1
            vari = 'thetao';
        else
            vari = 'so';
        end
       
        load('C:/Users/brian.grieve/Documents/Model_Data/WOA_CLIM/NWA/clim85-13'); % WOA climatology
        load('C:/Users/brian.grieve/Documents/Model_Data/ETOPO2v2g/etopo10.mat'); % ETOPO upscaled to 10 degree resolution
        filepath =  sprintf('cd E:/pcmdi/%s/%s/%s/%s',inst,model,histo,vari);
        eval(filepath);
        files = dir;

        % Identify files that contain target years
        for h = 3:length(files)
            startyeardex = strfind(files(h).name, 'r1i1p1_')+7;
            endyeardex = strfind(files(h).name, 'r1i1p1_')+14;
            firstyears1(h-2) = str2num(files(h).name(startyeardex:(startyeardex+3)));
            lastyears1(h-2) =  str2num(files(h).name(endyeardex:(endyeardex+3)));
        end;
        fileyears = 1+max(lastyears1-firstyears1);
        targetindex1 = find(target1>=firstyears1 & target1<=lastyears1);
        targetindex2 = find(target2>=firstyears1 & target2<=lastyears1);
        cutoff1 = target1-firstyears1(targetindex1); % Years in file before target date
        if strcmp(inst,'MOHC')==1
            cutoff1 = cutoff1-1; 
        end
        cutoff2 = fileyears-(lastyears1(targetindex2)-firstyears1(targetindex2)+1); % Number of years shorter last file is
        cutoff3 = fileyears-(lastyears1(targetindex1)-firstyears1(targetindex1)+1); % Number of years shorter first file is
        cutoff4 = abs(target2-lastyears1(targetindex2)); % Years in file after end date

        % Preallocation
        y = targetindex1+2:targetindex2+2;
        yy = max(y);
        x = (length(files)-length(y)-(length(files)-max(y,2)));
        xx = max(x);
        time = ncread(files(3).name,'time');
        m = size(time,1);
        lat = ncread(files(3).name,'lat');
        lon = ncread(files(3).name,'lon');
        if size(lat,2)>1
            [o p] = size(lat);
        else
            o = length(lat);
            p = length(lon);
        end
        allbots = nan([(yy-xx) m o p]);
        allsurf = allbots;

        for f = y
            filename = files(f).name;
            display(filename)
            data = ncread(filename,vari);
            time = ncread(filename,'time');
            m = size(time,1);
            lev = ncread(filename,'lev');
            data(data==0)=NaN;
            f1 = find(size(data)==length(time));
            f2 = find(size(data)==length(lev));
            f3 = find(size(data)==o);
            f4 = find(size(data)==p);
            data = permute(data,[f1 f2 f3 f4]);
            if f==y(1) & cutoff1>0
               data(1:(12*cutoff1),:,:,:) = [];
            elseif f==y(yy-xx) & cutoff4>0
                data((m-12*cutoff4+1):end,:,:,:) = [];
            end
            [m n o p] = size(data);
            mat = zeros([n 1]);
            bottomvar = nan([m o p]);

            if f==y(1)
                for j = 1:o
                    for k = 1:p
                        for i = 1:n
                            A = data(1,i,j,k);
                            B = isfinite(A);
                            mat(i,:) = B;
                            maxlayer = find(mat==1,1,'last');
                            if isempty(maxlayer)
                                maxlayer=NaN;
                            end
                            bottomlayer(j,k) = maxlayer;
                        end
                    end
                end
            end

            for t = 1:m
                for j = 1:o
                    for k = 1:p
                        A = bottomlayer(j,k);
                        if isnan(A)
                            continue
                        end
                        bottomvar(t,j,k) = data(t,A,j,k);
                        allbots(f-xx,t,j,k) = bottomvar(t,j,k);
                        allsurf(f-xx,t,j,k) = squeeze(data(t,1,j,k));
                    end
                end
            end
        end

        clear data

        % Reshape into [months lon lat]
        allbots = permute(allbots,[2 1 3 4]);
        bvarhist = reshape(allbots, [(size(allbots,1)*size(allbots,2)) size(allbots,3) size(allbots,4)]);
        allsurf = permute(allsurf,[2 1 3 4]);
        svarhist = reshape(allsurf, [(size(allsurf,1)*size(allsurf,2)) size(allsurf,3) size(allsurf,4)]);
        q = size(bvarhist,1);

        % Remove NaN layers (e.g. a result of fitting 5 years of data into ten year windows)
        % Ensures there is data for all years in climate model
        nanlayer = nan([q 1]);
        for i = 1:q
            A = squeeze(bvarhist(i,:,:));
            nanlayer(i) = sum(isnan(A(:)))==numel(A);
        end
        nanlayer = logical(nanlayer);
        nanlayer_hist = nanlayer;
        if nansum(nanlayer(1:fileyears*12))~= 12*cutoff1 & nansum(q-nanlayer(fileyears*12:end))~= 12*cutoff2
            figure; bar(nanlayer); title('NaN layers, historical');
            error('Historical NaN layers are not equal to removed files. Check nanlayer.') % Check if there is an equal number of nan years as expected
        elseif sum(nanlayer(fileyears*12+1:q-fileyears*12))>0
            warning('No data in some historical years. Check nanlayer')
            nanlayer(fileyears*12+1:q-fileyears*12)=false;
            bvarhist(nanlayer,:,:)=[];
            svarhist(nanlayer,:,:)=[];
            figure; bar(nanlayer_hist); title('NaN layers, historical');
        else
            bvarhist(nanlayer,:,:)=[];
            svarhist(nanlayer,:,:)=[];
        end
        %
        if strcmp(inst,'MOHC')
            bvarhist(1,:,:) = [];
            svarhist(1,:,:) = [];
        end
        %
        %
        clear A B f filename files fileyears i j k m mat maxlayer bottomvar t data allbot* z firstyears change change2 lastyears startyeardex endyeardex x y

        % Repeat for projected data
        filepath =  sprintf('cd E:/pcmdi/%s/%s/%s/%s',inst,model,RCP,vari);
        eval(filepath);
        files = dir;

        for h = 3:length(files)
            startyeardex = strfind(files(h).name, 'r1i1p1_')+7;
            endyeardex = strfind(files(h).name, 'r1i1p1_')+14;
            firstyears2(h-2) = str2num(files(h).name(startyeardex:(startyeardex+3)));
            lastyears2(h-2) =  str2num(files(h).name(endyeardex:(endyeardex+3)));
        end;
        fileyears = 1+max(lastyears2-firstyears2);
        targetindex3 = find(target3>=firstyears2 & target3<=lastyears2);
        targetindex4 = find(target4>=firstyears2 & target4<=lastyears2);
        cutoff1 = target3-firstyears2(targetindex3);
        cutoff2 = fileyears-(lastyears2(targetindex4)-firstyears2(targetindex4)+1);
        cutoff3 = fileyears-(lastyears2(targetindex3)-firstyears2(targetindex3)+1);
        cutoff4 = abs(target4-lastyears2(targetindex4));
        if strcmp(inst,'MOHC')==1
            cutoff1 = cutoff1-1; 
            cutoff2 = fileyears-(lastyears2(targetindex4(1))-firstyears2(targetindex4(1))+1);
            cutoff3 = 0;
            cutoff4 = 2100-target4; 
            fileyears=10; 
            targetindex4 = max(targetindex4); 
        end

        y = targetindex3+2:targetindex4+2;
        yy = max(y);
        x = (length(files)-length(y)-(length(files)-max(y))); 
        xx = max(x);
        time = ncread(files(4).name,'time');
        m = size(time,1);
        allbots = nan([(yy-xx) m o p]);
        allsurf = allbots;

        for f = y
            filename = files(f).name;
            display(filename)
            data = ncread(filename,vari);
            time = ncread(filename,'time');
            lev = ncread(filename,'lev');
            data(data==0)=NaN;
            f1 = find(size(data)==length(time));
            f2 = find(size(data)==length(lev));
            f3 = find(size(data)==o);
            f4 = find(size(data)==p);
            data = permute(data,[f1 f2 f3 f4]);
            m = size(data,1);
            if f==y(1) & cutoff1>0
                data(1:(12*cutoff1),:,:,:) = [];
            elseif f==y(yy-xx) & cutoff4>0
                data((m-12*cutoff4+1):end,:,:,:) = [];
            end
            [m n o p] = size(data);
            bottomvar = nan([m o p]);
            mat = nan([n 1]);


            if f==y(1)
                for j = 1:o
                    for k = 1:p
                        for i = 1:n
                            A = data(1,i,j,k);
                            B = isfinite(A);
                            mat(i,:) = B;
                            maxlayer = find(mat==1,1,'last');
                            if isempty(maxlayer)
                                maxlayer=NaN;
                            end
                            bottomlayer(j,k) = maxlayer;
                        end
                    end
                end
            end

            for t = 1:m
                for j = 1:o
                    for k = 1:p
                        A = bottomlayer(j,k);
                        if isnan(A)
                            continue
                        end
                        bottomvar(t,j,k) = data(t,A,j,k);
                        allbots(f-xx,t,j,k) = bottomvar(t,j,k);
                        allsurf(f-xx,t,j,k) = squeeze(data(t,1,j,k));
                    end
                end
            end
        end

        clear bottom* lev mat maxlayer temp time data

        allbots = permute(allbots,[2 1 3 4]);
        bvarproj = reshape(allbots, [(size(allbots,1)*size(allbots,2)) size(allbots,3) size(allbots,4)]);
        allsurf = permute(allsurf,[2 1 3 4]);
        svarproj = reshape(allsurf, [(size(allsurf,1)*size(allsurf,2)) size(allsurf,3) size(allsurf,4)]);

        q = size(bvarproj,1);

        clear allbots allsurf

        % Remove NaN Layers
        %
        disp('NaN removal');
        nanlayer = nan([q 1]);
        for i = 1:q
            A = squeeze(bvarproj(i,:,:));
            nanlayer(i) = sum(isnan(A(:)))==numel(A);
        end
        nanlayer_proj = nanlayer; 
        nanlayer = logical(nanlayer);
        if strcmp(inst,'MOHC')==0
            if nansum(nanlayer(1:fileyears*12))~= 12*cutoff3 & nansum(nanlayer(q-fileyears*12:end))~= 12*cutoff2
                figure; bar(nanlayer); title('NaN layers, projected');
                error('Projected NaN layers are not equal to removed files. Check nanlayer.') % Check if there is an equal number of nan years as expected
            elseif sum(nanlayer(fileyears*12+1:q-fileyears*12))>0
                warning('No data in some projected years. Check nanlayer')
                nanlayer(fileyears*12+1:q-fileyears*12)=false;
                bvarproj(nanlayer,:,:)=[];
                svarproj(nanlayer,:,:)=[];
                figure; bar(nanlayer_proj); title('NaN layers, projected');
            else
                bvarproj(nanlayer,:,:)=[];
                svarproj(nanlayer,:,:)=[];
            end
        elseif strcmp(inst,'MOHC')==1
             if nansum(nanlayer(1:fileyears*12))~= 12*cutoff3 & nansum(nanlayer(q-2*fileyears*12:end))~= 12*cutoff2+(12*fileyears-1)
                figure; bar(nanlayer); title('NaN layers, projected');
                error('Projected NaN layers are not equal to removed files. Check nanlayer.') % Check if there is an equal number of nan years as expected
            elseif sum(nanlayer(fileyears*12+1:q-2*fileyears*12))>0
                warning('No data in some projected years. Check nanlayer')
                nanlayer(fileyears*12+1:q-fileyears*12)=false;
                bvarproj(nanlayer,:,:)=[];
                svarproj(nanlayer,:,:)=[];
                figure; bar(nanlayer_proj); title('NaN layers, projected');
            else
                bvarproj(nanlayer,:,:)=[];
                svarproj(nanlayer,:,:)=[];
             end
        end
        %


        % Concat all years into one file with correct year
        bvarfull = cat(1,bvarhist,bvarproj);
        lengthhist = size(bvarhist,1); lengthproj = size(bvarproj,1); 
        clear bvarhist bvarproj
        svarfull = cat(1,svarhist,svarproj);
        clear svarhist svarproj
        if size(bvarfull,1)/12 ~=target4-target1+1
            error('Timeseries incorrect length')
        end
        %
        disp('Separate by year and month');
        qq = size(bvarfull,1);
        q = length([1:12:qq]);
        bvarfull2 = nan([q 12 o p]);
        for i = 1:12
            disp(i)
            bvarfull2(:,i,:,:) = bvarfull([i:12:qq],:,:); 
        end
        clear bvarfull
        svarfull2 = nan([q 12 o p]);
        for i = 1:12
            disp(i)
            svarfull2(:,i,:,:) = svarfull([i:12:qq],:,:);
        end
        %
        clear bvarfull svarfull A B cutoff* endyeardex f* i j k last* startyeardex t x* y*

        % Pick out delta of initial conditions. WOA dataset uses 1955-2012 mean conditions
        disp('Delta Clim'); 
        bvarclim = nan([12 o p]);
        bvardelta = nan([q 12 o p]);
        for i = 1:12
            disp(i)
            bvarclim(i,:,:) = nanmean(bvarfull2(1:2013-target1,i,:,:),1);
            for t = 1:q
                A = squeeze(bvarfull2(t,i,:,:));
                bvardelta(t,i,:,:) = A-squeeze(bvarclim(i,:,:));
            end
        end;
        clear bvarclim bvarfull2
        svarclim = nan([12 o p]);
        svardelta = nan([q 12 o p]);
        for i = 1:12
            disp(i)
            svarclim(i,:,:) = nanmean(svarfull2(1:2013-target1,i,:,:),1);
            for t = 1:q
                A = squeeze(svarfull2(t,i,:,:));
                svardelta(t,i,:,:) = A-squeeze(svarclim(i,:,:));
            end
        end
        clear svarclim svarfull2 

        if size(lat,2)==1
            lon2 = meshgrid(double(lon),double(lat));
            lat2 = meshgrid(double(lat),double(lon))';
        else
            lon2 = double(lon);
            lat2 = double(lat);
        end
        lon10 = lon10+360;
        if min(lon(:))<0
            lon2 = lon2+360;
        end

        [x(1) y(1)] = ind2sub([o p],dsearchn([lon2(:),lat2(:)],[280 30]));
        [x(2) y(2)] = ind2sub([o p],dsearchn([lon2(:),lat2(:)],[310 55]));
        svarreg = svardelta(:,:,x(1):x(2),y(1):y(2));
        bvarreg = bvardelta(:,:,x(1):x(2),y(1):y(2));

        % Interpolate temperature where coastline and climate model clash
        cd C:/Users/brian.grieve/Documents/Scripts
        bvarregfill = nan([q 12 x(2)-x(1)+1 y(2)-y(1)+1]);
        svarregfill = bvarregfill;
        for i = 1:q
            if ismember(i,[1:10:q]);
                disp(sprintf('fillnans_%d',i));
            end
            for t = 1:12
                bvarregfill(i,t,:,:) = fillnans(squeeze(bvarreg(i,t,:,:)),'radius',2);
                bvarregfill(i,t,:,:) = fillnans(squeeze(bvarregfill(i,t,:,:)),'radius',2);
                svarregfill(i,t,:,:) = fillnans(squeeze(svarreg(i,t,:,:)),'radius',2);
                svarregfill(i,t,:,:) = fillnans(squeeze(svarregfill(i,t,:,:)),'radius',2);
            end
        end
        bvardeltafill = bvardelta;
        bvardeltafill(:,:,x(1):x(2),y(1):y(2)) = bvarregfill;
        clear bvarregfill bvarreg bvardelta
        svardeltafill = svardelta;
        svardeltafill(:,:,x(1):x(2),y(1):y(2)) = svarregfill;

        clear lat252 lat lon x y t svarregfill svardelta svarreg

        %




        % Interpolate model resolution to WOA resolution
        o = size(lon10,1);
        p = size(lon10,2);
        vqsdelta = nan([q 12 o p]);
        vqbdelta = nan([q 12 o p]);
        warning('off');
        for i = 1:q
            disp(sprintf('griddata_%d',i));
            for mon = 1:12
                A = squeeze(svardeltafill(i,mon,:,:));
                vqsdelta(i,mon,:,:) = griddata(lon2,lat2,A,lon10,lat10,'nearest');
                A = squeeze(bvardeltafill(i,mon,:,:));
                vqbdelta(i,mon,:,:) = griddata(lon2,lat2,A,lon10,lat10,'nearest');
            end
        end
        warning('on'); 
        clear A svardeltafill bvardeltafill

        % Save these 3D (Time x 2Darea) ocean temperatures for later
        filepath =  sprintf('cd E:/pcmdi/%s/%s/%s',inst,model,RCP);
        eval(filepath);

        if variable==1
            disp('thetao')
            sst_ts = nan([q 12 o p]);
            bwt_ts = nan([q 12 o p]);
            for i = 1:q
                sst_ts(i,:,:,:) = squeeze(vqsdelta(i,:,:,:))+sst_mon;
                bwt_ts(i,:,:,:) = squeeze(vqbdelta(i,:,:,:))+bwt_mon;
            end
            a = isfinite(sst_ts); 
            sst_delta = nan(size(sst_ts));
            bwt_delta = sst_delta; 
            sst_delta(a) = vqsdelta(a);
            bwt_delta(a) = vqbdelta(a);
            save('temperatures_10.mat','sst_ts','bwt_ts','sst_delta','bwt_delta');
            clear sst_ts bwt_ts sst_delta bwt_delta vqsdelta vqbdelta
        elseif variable==2
            disp('so');
            sss_ts = nan([q 12 o p]);
            bws_ts = nan([q 12 o p]);
            for i = 1:q
                sss_ts(i,:,:,:) = squeeze(vqsdelta(i,:,:,:))+sss_mon;
                bws_ts(i,:,:,:) = squeeze(vqbdelta(i,:,:,:))+bws_mon;
            end
            a = isfinite(sss_ts); 
            sss_delta = nan(size(sss_ts));
            bws_delta = sss_delta; 
            sss_delta(a) = vqsdelta(a);
            bws_delta(a) = vqbdelta(a);
            save('salinities_10.mat','sss_ts','bws_ts','sss_delta','bws_delta');
            clear sss_ts bws_ts sss_delta bws_delta
        end

        clear A bvarreg* bvardelta bvardeltafill data filepath i m mon n nanlayer* qq svardelta* svarreg* 
    %
    end
    %}

%}
%
    %
    % Export to R. Run statistical models and predict skate abundance for
    % each climate model 
    disp('Export')
    filepath =  sprintf('cd D:/pcmdi/%s/%s/%s',inst,model,RCP);
    eval(filepath);
    load('salinities_10.mat','bws_ts','sss_ts'); load('temperatures_10.mat','bwt_ts','sst_ts');
    load('D:/Backup/2018_6_01/Model_Data/ETOPO2v2g/ETOPO_NWA.mat','etopo','latE','lonE'); 
    load('D:/Backup/2018_6_01/Model_Data/ETOPO2v2g/etopo10.mat')
    load('D:/Backup/2018_6_01/Model_Data/ETOPO2v2g/etopo_us.mat')
    load('D:/Backup/2018_6_01/Thorny_Skate/latlon10.mat'); 
    load('D:/Backup/2018_6_01/Calanus/Data/latlon.mat'); 

    if strcmp(model,'GFDL_CM2.6')==1
        lon = lon_gfdl; 
        lat = lat_gfdl;
        region = region_gfdl; 
        etopo = etopo_us; 
        load('C:/Users/brian.grieve/Documents/Thorny_Skate/Data/blank_varexp_gfdl.mat'); 
    else 
        lon = lon10; 
        lat = lat10;
        region = region10; 
        etopo = etopo10; 
        load('D:/Backup/2018_6_01/Thorny_Skate/Data/blank_varexp.mat'); 
    end
    
    a = nan(size(lon));
    a(etopo<=1000)=1;
    ind_fin = find(isfinite(a(:)) & isfinite(region(:)));
    %
    
    if size(sst_ts,1)==145
        sst_ts(146,:,:,:) = nan;
        sss_ts(146,:,:,:) = nan; 
        bwt_ts(146,:,:,:) = nan; 
        bws_ts(146,:,:,:) = nan; 
    end
    sst = nan([4 size(sst_ts,1) size(lon)]);
    sss = sst; 
    bwt = sst; 
    bws = sst; 
    for i = 1:4
        sst(i,:,:,:) = squeeze(nanmean(sst_ts(:,3*i-2:3*i,:,:),2));
        sss(i,:,:,:) = squeeze(nanmean(sss_ts(:,3*i-2:3*i,:,:),2));
        bwt(i,:,:,:) = squeeze(nanmean(bwt_ts(:,3*i-2:3*i,:,:),2));
        bws(i,:,:,:) = squeeze(nanmean(bws_ts(:,3*i-2:3*i,:,:),2));
    end

    cd D:/Backup/2018_6_01/Thorny_skate/Conditions_2D/DataIntoR
    var_exp = blank_varexp(ind_fin,:); 
    [q m o p] = size(sst_ts); 

    for i = 1:q
        for m = 1:4
            A = squeeze(sst(m,i,:,:));
            var_exp(:,16) = A(ind_fin);
            A = squeeze(sss(m,i,:,:));
            var_exp(:,18) = A(ind_fin);
            A = squeeze(bwt(m,i,:,:));
            var_exp(:,17) = A(ind_fin);
            A = squeeze(bws(m,i,:,:));
            var_exp(:,19) = A(ind_fin);
            var_exp(:,15) = m;
            if m==1
                var_exp_all = var_exp;
            else
                var_exp_all = [var_exp_all; var_exp];
            end
        end
        csvwrite(sprintf('skate_var_exp_all_%d.csv',i),var_exp_all);
    end
    seas_vec = var_exp_all(:,15);
    
    clear filepath sst* sss* bwt* bws* i var_exp* blank_varexp
%
%
    % Send data to R, run species distribution model Prediction_Script.R in R
    disp('R'); 
    cd 'C:/Program Files/R/R-4.0.1/bin'
    system('Rscript "C:\Users\brian.grieve\Documents\thornyskate\skate_pred3_short.R"');
    
    %
    % Import Data from R
    disp('Import from R'); 
    shortabun = nan([4 q o p]);
    longabun = nan([4 q o p]); 
    shortpa = nan([4 q o p]); 
    longpa = nan([4 q o p]); 
    llabun = nan([4 q o p]); 
    %
    filename = 'D:/Backup/2018_6_01/Thorny_skate/Conditions_2D/ExportFromR2/gamshort_pred.csv';
    gamshort = importdata(filename);
    gamshort = gamshort.data;
    gamshort(gamshort==-999)=nan;
    %
    filename = sprintf('D:/Backup/2018_6_01/Thorny_skate/Conditions_2D/ExportFromR2/gamlong_pred.csv');
    gamlong = importdata(filename);
    gamlong = gamlong.data;
    gamlong(gamlong==-999)=nan;
    %
    filename = sprintf('D:/Backup/2018_6_01/Thorny_skate/Conditions_2D/ExportFromR2/gamshort_pa.csv');
    gamshort_pa = importdata(filename);
    gamshort_pa = gamshort_pa.data;
    gamshort_pa(gamshort_pa==-999)=nan;
    %
    filename = sprintf('D:/Backup/2018_6_01/Thorny_skate/Conditions_2D/ExportFromR2/gamlong_pa.csv');
    gamlong_pa = importdata(filename);
    gamlong_pa = gamlong_pa.data;
    gamlong_pa(gamlong_pa==-999)=nan;
    %}
    filename = sprintf('D:/Backup/2018_6_01/Thorny_skate/Conditions_2D/ExportFromR2/gamll_pred.csv');
    gamll_abun = importdata(filename);
    gamll_abun = gamll_abun.data;
    gamll_abun(gamll_abun==-999)=nan;
    %}
     for i = 1:q
         for m = 1:4
            %
            map_g = nan(o,p); 
            map_g(ind_fin) = gamshort(seas_vec==m,i); 
            shortabun(m,i,:,:) = map_g;
             %
            map_g = nan(o,p); 
            map_g(ind_fin) = gamlong(seas_vec==m,i); 
            longabun(m,i,:,:) = map_g;
            map_g = nan(o,p); 
            map_g(ind_fin) = gamshort_pa(seas_vec==m,i); 
            shortpa(m,i,:,:) = map_g;
            %
            map_g = nan(o,p); 
            map_g(ind_fin) = gamlong_pa(seas_vec==m,i); 
            longpa(m,i,:,:) = map_g;
            %
            map_g = nan(o,p); 
            map_g(ind_fin) = gamll_abun(seas_vec==m,i); 
            llabun(m,i,:,:) = map_g;
            %}
         end
     end

    if modelnumber==2 % No 2100 here
        longabun(:,146,:,:) = nan; 
        shortabun(:,146,:,:) = nan; 
        shortpa(:,146,:,:) = nan; 
        longpa(:,146,:,:) = nan;
        llabun(:,146,:,:) = nan; 
    end
    
    clear A clim delta shortprop longprop llprop gammap short_abun_pred long_abun_pred 
    
     longabun_all(modelnumber-min(mn)+1,:,:,:,:) = longabun; 
     shortabun_all(modelnumber-min(mn)+1,:,:,:,:) = shortabun;
     llabun_all(modelnumber-min(mn)+1,:,:,:,:) = llabun; 
     shortpa_all(modelnumber-min(mn)+1,:,:,:,:) = shortpa; 
     longpa_all(modelnumber-min(mn)+1,:,:,:,:) = longpa;
end
%clear longabun shortabun llabun shortpa longpa
%}

    
