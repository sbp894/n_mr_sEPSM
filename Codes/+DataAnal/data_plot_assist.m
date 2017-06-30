function plot_assist=data_plot_assist(spike_data,StimData,condition_var)

plot_assist.Fs=spike_data(condition_var).Fs;
plot_assist.paramsIN=spike_data(condition_var).paramsIN;
plot_assist.CurStims=StimData{condition_var};
