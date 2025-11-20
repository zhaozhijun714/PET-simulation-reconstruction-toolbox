function [x, H] = fbp2d_simulation(G, imgsiz, voxel_size, proj_size, yi, ni, ri, fbptype, filter, frequency_scaling, interp_method, use_psf, FWHM)
if isempty(fbptype)
    fbptype = 'siddon';
end
if isempty(filter)
    filter = 'linear';
end
if isempty(interp_method)
    interp_method = 'ram-lak';
end
if isempty(frequency_scaling)
    frequency_scaling = 1;
end
yi = reshape(yi, proj_size);
ni = reshape(ni, proj_size);
ri = reshape(ri, proj_size);
proj = (yi - ri) ./ ni;
proj(proj < 0) = 0;
theta = linspace(0, 179, proj_size(2));
if strcmp(fbptype, 'siddon')
    [proj_filter, H] = filterProjections_siddon_simulation(proj, filter, frequency_scaling);
    x = proj_back_simulation(G, proj_filter(:), imgsiz, voxel_size, use_psf, FWHM);
else strcmp(fbptype, 'standard')
    x = FilterBackProjection_simulation(proj, theta, interp_method, filter, frequency_scaling, double(imgsiz(1)));
    x = x(:);
end
end
function [p,H] = filterProjections_siddon_simulation(p_in, filter, d)
p = p_in;
len = size(p,1);
H = designFilter_siddon(filter, len, d);
if strcmpi(filter, 'none')
    return;
end
p(length(H),1)=0;
p = fft(p);
p = bsxfun(@times, p, H);
p = ifft(p,'symmetric');
p(len+1:end,:) = [];
end
function filt = designFilter_siddon(filter, len, d)
order = max(64,2^nextpow2(2*len));
if strcmpi(filter, 'none')
    filt = ones(1, order);
    return;
end
n = 0:(order/2);
filtImpResp = zeros(1,(order/2)+1);
filtImpResp(1) = 1/4;
filtImpResp(2:2:end) = -1./((pi*n(2:2:end)).^2);
filtImpResp = [filtImpResp filtImpResp(end-1:-1:2)];
filt = 2*real(fft(filtImpResp));
filt = filt(1:(order/2)+1);
w = 2*pi*(0:size(filt,2)-1)/order;
switch filter
    case 'ram-lak'
    case 'shepp-logan'
        filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));
    case 'cosine'
        filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d));
    case 'hamming'
        filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/d));
    case 'hann'
        filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./d)) / 2;
    otherwise
        error(message('images:iradon:invalidFilter'))
end
filt(w>pi*d) = 0;
filt = [filt' ; filt(end-1:-1:2)'];
end
function [img, H] = FilterBackProjection_simulation(varargin)
[p,theta,filter,d,interp,N] = parse_inputs(varargin{:});
len=size(p,1);
H = designFilter(filter, len, d);
p(length(H),1)=0;
p = fft(p);
for i = 1:size(p,2)
    p(:,i) = p(:,i).*H;
end
p = real(ifft(p));
p(len+1:end,:) = [];
img = zeros(N);
xax = (1:N)-ceil(N/2);
x = repmat(xax, N, 1);
y = rot90(x);
costheta = cos(theta);
sintheta = sin(theta);
ctrIdx = ceil(len/2);
imgDiag = 2*ceil(N/sqrt(2))+1;
if size(p,1) < imgDiag
    rz = imgDiag - size(p,1);
    p = [zeros(ceil(rz/2),size(p,2)); p; zeros(floor(rz/2),size(p,2))];
    ctrIdx = ctrIdx+ceil(rz/2);
end
if strcmp(interp, 'nearest neighbor')
    for i=1:length(theta)
        proj = p(:,i);
        t = round(x*costheta(i) + y*sintheta(i));
        img = img + proj(t+ctrIdx);
    end
elseif strcmp(interp, 'linear')
    for i=1:length(theta)
        proj = p(:,i);
        t = x.*costheta(i) + y.*sintheta(i);
        a = floor(t);
        img = img + (t-a).*proj(a+1+ctrIdx) + (a+1-t).*proj(a+ctrIdx);
    end
elseif strcmp(interp, 'spline')
    for i=1:length(theta)
        proj = p(:,i);
        taxis = (1:size(p,1)) - ctrIdx;
        t = x.*costheta(i) + y.*sintheta(i);
        projContrib = interp1(taxis,proj,t(:),'*spline');
        img = img + reshape(projContrib,N,N);
    end
end
img = img*pi/(2*length(theta));
end

function filt = designFilter(filter, len, d)
order = max(64,2^nextpow2(2*len));
filt = 2*( 0:(order/2) )./order;
w = 2*pi*(0:size(filt,2)-1)/order;
switch filter
    case 'ram-lak'
    case 'shepp-logan'
        filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));
    case 'cosine'
        filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d));
    case 'hamming'
        filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/d));
    case 'hann'
        filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./d)) / 2;
    otherwise
        error('Invalid filter selected.');
end
filt(w>pi*d) = 0;
filt = [filt' ; filt(end-1:-1:2)'];
end

function [p,theta,filter,d,interp,N] = parse_inputs(varargin)
if nargin<2
    error('Invalid input arguments.');
end
p = varargin{1};
theta = pi*varargin{2}/180;
N = 0;
d = 1;
filter = 'ram-lak';
interp = 'linear';
string_args = {'nearest neighbor', 'linear', 'spline', ...
    'ram-lak','shepp-logan','cosine','hamming', 'hann'};
for i=3:nargin
    arg = varargin{i};
    if ischar(arg)
        idx = strmatch(lower(arg),string_args);
        if isempty(idx)
            error(['Unknown input string: ' arg '.']);
        elseif prod(size(idx)) > 1
            error(['Ambiguous input string: ' arg '.']);
        elseif prod(size(idx)) == 1
            if idx <= 3   % It is the interpolatio
                interp = string_args{idx};
            elseif (idx > 3) & (idx <= 8)
                filter = string_args{idx};
            end
        end
    elseif prod(size(arg))==1
        if arg <=1
            d = arg;
        else
            N = arg;
        end
    else
        error('Invalid input parameters');
    end
end
if N==0
    N = 2*floor( size(p,1)/(2*sqrt(2)) );
end
if isempty(theta)
    theta = pi / size(p,2);
end
if prod(size(theta))==1
    theta = (0:(size(p,2)-1))* theta;
end
if length(theta) ~= size(p,2)
    error('THETA does not match the number of projections.');
end
end

