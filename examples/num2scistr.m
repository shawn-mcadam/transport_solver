function s = num2scistr(x)
%NUM2SCISTR Convert a number x to a LaTeX string of its base 10 scientific
% notation: $sign mantissa \cdot 10^{exponent}$
    if x < 0, sign = "-"; else, sign = ""; end
    x = abs(x);
    flr = floor(log10(x));
    frac = log10(x) - flr;
    s = sign + num2str(10^frac) + "\cdot10^{" + num2str(flr) + "}";
end