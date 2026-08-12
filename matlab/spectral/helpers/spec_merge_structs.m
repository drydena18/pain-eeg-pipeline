function base = spec_merge_structs(base, varargin)
% SPEC_MERGE_STRUCTS Copy all fields from one or more structs into base
% V 1.0.0
%
%   base = spec_merge_structs(base, s1, s2, ...)
%
% Fields in later inputs overwrite same-named fields from earlier ones
% (including base). Empty imputs ([], struct()) are skipped harmlessly.
%
% Used to assemble featChan/featGA in spectral_core.m from the separate
% whole_/pre_/post_/delta_ feature structs without repeating the
% fieldnames()-loop pattern at every call site.
%
% Inputs:
%   base        : struct to merge into (may be struct() to start fresh)
%   varargin    : any number of additional structs
%
% Output:
%   base : struct with all fields from base and varargin{:} combined

for i = 1:numel(varargin)
    s = varargin{i};
    if isempty(s) || ~isstruct(s)
        continue;
    end
    fn = fieldnames(s);
    for k = 1:numel(fn)
        base.(fn{k}) = s.(fn{k});
    end
end
end