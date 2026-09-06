function label = classify_amplitude(amplitudeDeg)
%CLASSIFY_AMPLITUDE Apply the scope labels used by the LATO study.
%   LABEL = lato.scan.classify_amplitude(AMPLITUDEDEG) returns a string
%   array. Values greater than 45 deg are labeled large_amplitude and
%   values greater than 90 deg are labeled out_of_scope. Missing values
%   are labeled unavailable.

label = repmat("moderate_amplitude", size(amplitudeDeg));
label(isnan(amplitudeDeg)) = "unavailable";
label(amplitudeDeg > 45 & amplitudeDeg <= 90) = "large_amplitude";
label(amplitudeDeg > 90) = "out_of_scope";
end
