clear all; clc; close all;

global CineFile RunProps ImageRange CropWindow ROI
CineFile = ''
ImageRange = [1 40000] + 0 ; %Two element vector [start end]. If end=inf; goes until it runs out of images.
r=128; %Reduced order model size; also size of the batch

dt = 800e-9;

PerformCrop = 0; %Flag to enable cropping
PerformROI = 0; %Flag to enable roi polygon
DoublePulse = true;

%% 2. Loads images
LoadPhantomLibraries();
RegisterPhantom(true);

RunProps=GetCineProperties(CineFile);

I = double(ReadCineFileImage(CineFile,RunProps.firstImage+10,0));
I2 = double(ReadCineFileImage(CineFile,RunProps.firstImage+11,0));

ImageRange(2) = min(ImageRange(2), RunProps.numberOfImages);

%Defines crop window
if PerformCrop
    figure; imagesc(I); daspect([1 1 1]); colormap gray
    title('Select cropping window');
    R=getrect();
    CropWindow=round(R);

    
    close all;

    I=imcrop(I,CropWindow);
else
    CropWindow=0;
end

%Defines ROI mask
if PerformROI
    figure; imagesc(I); daspect([1 1 1]); colormap gray
    title('Select ROI for masking');
    ROI=roipoly();
    
    close all;
else
    ROI=0;
end

%% 2. Build per-pulse means

n          = size(I,1) * size(I,2);  % flattened image dimension
total_imgs = ImageRange(2) - ImageRange(1) + 1;

mu1 = zeros(n, 1);
mu2 = zeros(n, 1);
count1 = 0;
count2 = 0;

disp('Calculating mean...')
for i = 1:total_imgs
    I = getImages(ImageRange(1) + i - 1);
    x = I(:);
    
    if mod(i, 2) == 1   % pulse 1 (odd)
        mu1    = mu1 + x;
        count1 = count1 + 1;
    else                 % pulse 2 (even)
        mu2    = mu2 + x;
        count2 = count2 + 1;
    end
    
    if mod(i, 1000) == 0
        disp(['Mean pass: ' num2str(i) '/' num2str(total_imgs)]);
    end
end

mu1 = mu1 / count1;
mu2 = mu2 / count2;

% Second loop for variance (needs means first)
var1 = zeros(n, 1);
var2 = zeros(n, 1);

for i = 1:total_imgs
    I = getImages(ImageRange(1) + i - 1);
    x = I(:);
    
    if mod(i, 2) == 1
        var1 = var1 + (x - mu1).^2;
    else
        var2 = var2 + (x - mu2).^2;
    end
    
    if mod(i, 1000) == 0
        disp(['Variance pass: ' num2str(i) '/' num2str(total_imgs)]);
    end
end

var1 = var1 / count1;
var2 = var2 / count2;

gain_ratio = sqrt(mean(var1) / mean(var2));
disp(['Pulse gain ratio: ' num2str(gain_ratio)]);

%% 3. Burn-in: accumulate SVD subspace U via batched incremental SVD
batch_size = r; %Batch equal to reduced dimension

% --- Initialise with first r snapshots ---
AllImages = zeros(n,r);
for i = 1:r
    I = getImages(ImageRange(1) + i - 1);
    if mod(i,1)==0
        AllImages(:,i) = I(:) - mu1;
    else
        AllImages(:,i) = I(:) - mu2;
    end
    disp(['First batch: ' num2str(i)]);
end
[U, S_mat, ~] = svd(AllImages, 'econ');
sigma = diag(S_mat);

% --- Preallocate convergence log ---
n_batches    = ceil((total_imgs - r) / batch_size);
conv_angles  = zeros(n_batches, 1);
batch_idx    = 0;
converged    = false;
tol          = 1e-3;  % set your tolerance here

% --- Stream in batches ---
i = r + 1;
while i <= total_imgs && ~converged

    % Load next batch
    i_end      = min(i + batch_size - 1, total_imgs);
    b          = i_end - i + 1;          % actual batch size (last batch may be smaller)
    X_batch    = zeros(n, b);

    for j = 1:b
        I = getImages(ImageRange(1) + i + j - 2);
        X_batch(:,j) = I(:);
    end

    % 1. Project batch onto current U
    P_coeff = U' * X_batch;              % r x b
    P       = X_batch - U * P_coeff;     % n x b  residual

    % 2. Thin QR of residual to get new orthogonal directions
    [Q, R] = qr(P, 0);                   % Q: n x b,  R: b x b

    % 3. Form (r+b) x (r+b) core matrix
    K = [ diag(sigma),  P_coeff ;
          zeros(b, r),  R       ];

    % 4. SVD of core, truncate to r
    [Uk, Sk, ~] = svd(K, 'econ');
    Uk    = Uk(:, 1:r);
    sigma = diag(Sk(1:r, 1:r));          % keep as vector

    % 5. Update U
    U_prev = U;
    U      = [U, Q] * Uk;               % n x r

    % 6. Convergence metric
    overlap      = svd(U_prev' * U);
    angle_metric = norm(ones(r,1) - overlap);

    batch_idx               = batch_idx + 1;
    conv_angles(batch_idx)  = angle_metric;

    disp(['Batch ending at image ' num2str(i_end) ...
          '  (' num2str(b) ' imgs)  angle metric: ' num2str(angle_metric)]);

    if angle_metric < tol
        disp(['Converged at image ' num2str(i_end)]);
        converged = true;
    end

    i = i_end + 1;
end

conv_angles = conv_angles(1:batch_idx);  % trim unused

% --- Plot ---
figure;
semilogy(conv_angles, 'b-o', 'MarkerSize', 3, 'LineWidth', 1.5);
yline(tol, 'r--', sprintf('tol = %.0e', tol), 'LineWidth', 1.2);
xlabel('Batch index');
ylabel('Subspace angle metric');
title(sprintf('Batched incremental SVD convergence  (r=%d, batch=%d)', r, batch_size));
grid on;

%% Eliminate the modes related to the mean and do the DMD now
badModes = 1:2;
U(:,badModes) = [];
r=size(U,2);


%% Perform Online DMD for entire dataset
P = zeros(r, r);
Q = zeros(r, r);
check_idx = 0;
lambda_prev   = nan(r, 1);
conv_metric=[];
eig_history=[];

for i = 1:2:total_imgs-1
    % Pulse 1 snapshot (odd)
    I1 = getImages(ImageRange(1) + i - 1);
    z1 = U' * (I1(:) - mu1);

    % Pulse 2 snapshot (even)
    I2 = getImages(ImageRange(1) + i);
    z2 = U' * (I2(:) - mu2) * gain_ratio;

    % Accumulate: x = z1, y = z2
    P = P + z1 * z1';
    Q = Q + z2 * z1';

    if mod(round(i/2), 200) == 0
        A_tilde        = Q / P;
        lambda         = sort(eig(A_tilde), 'ComparisonMethod', 'real');
        check_idx      = check_idx + 1;
        eig_history(:, check_idx) = lambda;

        if ~any(isnan(lambda_prev))
            % Match eigenvalues by nearest neighbour in complex plane
            D = abs(lambda - lambda_prev.');   % r x r complex distance matrix
            
            % Hungarian assignment — minimize total matching distance
            [assignment, ~] = munkres(D); 
            [~, assignment_idx] = max(assignment, [], 2);
            lambda_matched  = lambda(assignment_idx);
            
            conv_metric(check_idx) = mean(abs(lambda_matched - lambda_prev));

            cla;
            plot(lambda,'k*'); daspect([1 1 1]);
            tt = linspace(0,2*pi,200);hold on;
            
            plot(cos(tt),sin(tt),'b:');drawnow; pause(0.01)

            disp(['DMD check at image ' num2str(i) ...
                  '  eig metric: ' num2str(conv_metric(check_idx), '%.4e')]);
        end

        lambda_prev = lambda;

    end
end

% Final solve for reduced operator
A_tilde = Q / P;

% Eigendecomposition
[W, Lambda] = eig(A_tilde);
lambda      = diag(Lambda);

% Continuous time eigenvalues
omega       = log(lambda) / dt;
freq        = imag(omega) / (2*pi);  % Hz

% Full spatial modes
Phi = U * W;                         % n x r

disp('DMD accumulation complete.');


%% Final pass: find energy spectrum (now that eigenvalues are almost unitary)
W_inv      = pinv(W);          % r x r
E          = zeros(r, 1);      % accumulate |a_i|^2
count      = 0;

for i = 1:2:total_imgs-1       % same cross-pulse pairs

    % Pulse 1
    I1 = getImages(ImageRange(1) + i - 1);
    z1 = U' * (I1(:) - mu1);

    % Pulse 2
    I2 = getImages(ImageRange(1) + i);
    z2 = U' * (I2(:) - mu2) * gain_ratio;

    % Project both onto DMD basis and accumulate
    a1 = W_inv * z1;
    a2 = W_inv * z2;

    E     = E + abs(a1).^2 + abs(a2).^2;
    count = count + 2;

    if mod(round(i/2), 500) == 0
        disp(['Mode energy pass: ' num2str(i) '/' num2str(total_imgs)]);
    end
end

E = E / count;   % mean power per mode

%Reorders the frequencies and energies
[~,freqOrder] = sort(freq,'ascend');
frequencies=freq(freqOrder);
E=E(freqOrder);

negativeFreqs = frequencies<0;
frequencies(negativeFreqs)=[];
E(negativeFreqs)=[];

%%
%==========Mode Amplitudes===========
figure;
stem(frequencies, E, 'k^');
xlabel('Frequency (Hz)');
ylabel('Mean mode energy (Count^2)');
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
title('Shear Layer DMD power spectrum');





stop;
save('ODMD.mat','-v7.3')


%%
%==========Plot Mode Shape===========
freqPlot = 32e3; %Hz; will find the closest one
saveVid = 0;
freqIdx = find(min(abs(freqPlot-freq)) == abs(freqPlot-freq));
% freqIdx = 6;
freqMatch = freq(freqIdx);
ModeShape = Phi(:,freqIdx);
ModeShape = reshape(ModeShape,size(I));
cBound = [-1 1] * 0.3 * max(abs(ModeShape(:))); %Symmetric color axis bounds given the mode shape
theta=linspace(0,6*pi,100);
if saveVid
    vid=VideoWriter(['ODMD f=' num2str(freqMatch,'%0.1f') ' Hz'],'MPEG-4');
    vid.FrameRate = 30;
    vid.Quality = 90;
    open(vid)
end
figure('color', 'w');
 for i=1:length(theta)
     II=squeeze(ModeShape*exp(1i*theta(i)));
%     II=squeeze(ModeShape);
     imagesc(-X_DMD_mm, Y_DMD_mm,real(II)); daspect([1 1 1]);
    
    set(gca, 'YDir', 'normal');
%     caxis(cBound);
colormap redblue;
    xlabel('X [mm]'); ylabel('Y [mm]');
    title(['St = ' num2str(freqMatch*L/U,'%0.3f') ]);
    drawnow;
     pause(0.01);
if saveVid
        writeVideo(vid,getframe(gcf));
end
end


