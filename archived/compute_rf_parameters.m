
a0 = [2 2 2 2 2 2];
options = optimset('Display','iter');
% Due to the image coordinates, x and y indices are swapped
f = @(a) a(1).*y.^2 + a(2).*x.*y + a(3).*x.^2+ a(4).*y +a(5).*x + a(6);
af = lsqnonlin(f, a0, [], [], options);

load('/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Analysis/20220607C_20noise.mat');
stixelSize = 20;

for jj = 1 : size(significance_maps,1)
    M = squeeze(significance_maps(jj,:,:));
    stmp = squeeze(spatial_maps(jj,:,:,1));
    
    if ~isempty(M)
        % Get the signficant points
        [row,col] = find(M > 0);
        P = [row,col];
        % Remove outliers.
        [B,TF] = rmoutliers(P);
        
        %------------------------------------------------------------------
        % Fit a 2D Gaussian to the points.
        try
            [z, a, b, alpha] = fitellipse(B');
            gauss_params(jj,:) = [1, alpha/pi*180,b,a,z(1),z(2),0];
%             gauss_params(jj,:) = fmgaussfit(1:size(stmp,1), 1:size(stmp,2), stmp.*M);
        catch
            if isempty(B)
                gauss_params(jj,[5,6]) = mean(P); % Just use the centroid points.
                gauss_params(jj,[3,4]) = [std(P(:,1)),std(P(:,2))];
            else
                gauss_params(jj,[5,6]) = mean(B); % Just use the centroid points.
                gauss_params(jj,[3,4]) = [std(B(:,1)),std(B(:,2))];
            end
        end
        % Check the params.
        gauss_params(jj,3:4) = abs(gauss_params(jj,3:4));
        if gauss_params(jj,4) > 3*gauss_params(jj,3)
            gauss_params(jj,4) = gauss_params(jj,3);
        end

        %------------------------------------------------------------------
%         % Fit an ellipse to the points.
%         x = B(:,1); y = B(:,2);
%         f = @(a) a(1).*y.^2 + a(2).*x.*y + a(3).*x.^2+ a(4).*y +a(5).*x + a(6);
%         af = lsqnonlin(f, a0, [], [], options);
        %------------------------------------------------------------------
    end
end

