[sample, freq] = audioread('xgrigo02.wav')
% 2:
ft = fft(sample, freq)
norm_ft = ft / freq
half_freq = freq/2
spectrum = abs(norm_ft(1:half_freq))
freq_axis = [1:half_freq]

plot(freq_axis, spectrum)
xlabel('f[Hz]')
ylabel('|ck|')
grid()
pause

% 3:
[val_max, arg_max] = max(abs(norm_ft))
disp(['Max spectrum module [Hz]: ' num2str(arg_max)])
pause

% 4:
b = [0.2324, -0.4112, 0.2324] 
a = [1, 0.2289, 0.4662]
zplane(b,a)
p = roots(a) 
if (isempty(p) | abs(p) < 1) 
  disp('stabilní')
else
  disp('nestabilní')
end
pause

% 5:
H = freqz(b,a,256); f=(0:255) / 256 * freq / 2
plot (f,abs(H)); grid; xlabel('f'); ylabel('|H(f)|')
disp('Horní propust')

pause

% 6:
res = filter(b, a, sample)

res_ft = fft(res, freq) / freq

plot(freq_axis, abs(res_ft(1:half_freq)))
pause

% 7:
[val_max, arg_max] = max(abs(res_ft))
disp(['Max filtered spectrum module [Hz]: ' num2str(arg_max)])
pause
 
% 9:

d = 10
samples_del10 = [zeros(d, 1); sample(1:length(sample)-d)]
[acorrel, d] = xcorr(samples_del10, sample, 'unbiased')
fifty_fifty = acorrel(length(acorrel)/2-49: length(acorrel)/2+50)
plot([-49: 50], fifty_fifty)
pause

% 10:
disp(['Hodnota koeficentu R[10]: ' num2str(fifty_fifty(length(fifty_fifty)/2+10))])
pause


% 11:
gmin = min(min(sample)); gmax = max(max(sample));
% budeme chtit 50 bodu:
kolik = 50;
g = linspace(gmin,gmax,kolik);

% 12:
% probability calculation from hist2opt.m
n = 10;
delayed_n = zeros(size(sample))
samp_len = length(sample)   

% delaying signal by 10 samples
delayed_n(n:samp_len) = sample(1:samp_len-n+1)

xmin = min(sample); xmax = max(sample);

kolik = 50;
x = linspace(xmin,xmax,kolik);

L = length(x);
N = length(sample);

h = zeros(L,L)

xcol = x(:); bigx = repmat(xcol,1,N); 

yr = sample(:)'; bigy = repmat(yr,L,1);
[dummy,ind1] = min(abs(bigy - bigx)); 
yr = delayed_n(:)'; bigy = repmat(yr,L,1);
[dummy,ind2] = min(abs(bigy - bigx)); 

for ii=1:N,
  d1 = ind1(ii);   d2 = ind2(ii); 
  h(d1,d2) = h(d1,d2) + 1; 
end

surf = (x(2) - x(1))^2;
p = h / N / surf;  

% image for ex. 11
imagesc (x,x,p); axis xy; colorbar; xlabel('x2'); ylabel('x1');

% Integral result check for ex. 12
disp(['Integral = ' num2str(sum(sum(p)) * surf)])

% 13:
x_col = x(:); X1 = repmat(x_col,1,L);
x_row = x_col'; X2 = repmat(x_row,L,1); 


auto_c = sum(sum (X1 .* X2 .* p)) * surf;
disp(['R[10] = ' num2str(auto_c)]);


