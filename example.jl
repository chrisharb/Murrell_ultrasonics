using TDMSReader, PyPlot
include("base.jl")
path = "C:\\Users\\cwaha\\Dropbox\\My PC (DESKTOP-8JF2H49)\\Documents\\Mistras_testing\\Run0130_part3.tdms"
TDMS = readtdms(path)
traces, D1, D2, D3 = read_US(TDMS,1)
plot(D1(:t_us),traces[1])
