r = 2



r = max(sqrt(eps(class(r))), abs(r));
w = (abs(r)<pi) .* sin(r) ./ r

w = (abs(r)<1) .* (1 - r.^2).^2

w = 1 ./ (1 + r.^2)


w = 1 ./ (1 + abs(r))


w = 1 ./ max(1, abs(r))


r = max(sqrt(eps(class(r))), abs(r));
w = tanh(r) ./ r

w = ones(size(r))

w = 1 * (abs(r)<1)

w = exp(-(r.^2))