function tak_sim_plot_alg_result(output)
% tak_function
%=========================================================================%
% - Comments
%=========================================================================%
% (06/29/2014)
%%
figure,imexpb
if isfield(output, 'fval')
    subplot(131),tplot(log10(output.fval(2:end))), title('log10(function value)')
end
if isfield(output, 'fval')
    subplot(132),tplot(log10(output.wdist)), title('log10(||wtrue-west||)')
end
subplot(133),tplot(log10(output.rel_changevec)), title('log10(wnew-wold)')
drawnow