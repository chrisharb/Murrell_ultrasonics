"""
Read ultrasonic scan data
"""
function read_US(path::String, Ind::Int)
        chan = "BR_Rx"
        BR_Rx = Dict(
        :t_stamp        => t_conv(fid.groups[chan][Ind].props["Time (s)"]),
        :t_us           => (0:size(fid.groups[chan][Ind].data,1)-1)./fid.groups[chan][Ind].props["Fs (Hz)"],
        :width_us       => fid.groups[chan][Ind].props["Width (us)"],
        :delay_us       => fid.groups[chan][Ind].props["Delay (us)"],
        :gain_dB        => fid.groups[chan][Ind].props["Gain (dB)"],
        :n_stack        => fid.groups[chan][Ind].props["Stacking"],
        :F_pulse_MHz    => fid.groups[chan][Ind].props["Frequency (MHz)"],
        :V_pulse        => fid.groups[chan][Ind].props["Voltage"],
        :PRF_us         => fid.groups[chan][Ind].props["PRF (us)"]
        )
        trace1 = fid.groups[chan][Ind].data
        chan = "TR_Rx"
        TR_Rx = Dict(
        :t_stamp        => t_conv(float(fid.groups[chan][Ind].props["Time (s)"])),
        :t_us           => (0:size(fid.groups[chan][Ind].data,1)-1)./fid.groups[chan][Ind].props["Fs (Hz)"],
        :width_us       => fid.groups[chan][Ind].props["Width (us)"],
        :delay_us       => fid.groups[chan][Ind].props["Delay (us)"],
        :gain_dB        => fid.groups[chan][Ind].props["Gain (dB)"],
        :n_stack        => fid.groups[chan][Ind].props["Stacking"],
        :F_pulse_MHz    => fid.groups[chan][Ind].props["Frequency (MHz)"],
        :V_pulse        => fid.groups[chan][Ind].props["Voltage"],
        :PRF_us         => fid.groups[chan][Ind].props["PRF (us)"])
        trace2 = fid.groups[chan][Ind].data
        chan = "TR_Rx"
        TR_Tx = Dict(:t_stamp => t_conv(float(fid.groups[chan][Ind].props["Time (s)"])),
        :t_us          => (0:size(fid.groups[chan][Ind].data,1)-1)./fid.groups[chan][Ind].props["Fs (Hz)"],
        :width_us       => fid.groups[chan][Ind].props["Width (us)"],
        :delay_us       => fid.groups[chan][Ind].props["Delay (us)"],
        :gain_dB        => fid.groups[chan][Ind].props["Gain (dB)"],
        :n_stack        => fid.groups[chan][Ind].props["Stacking"],
        :F_pulse_MHz    => fid.groups[chan][Ind].props["Frequency (MHz)"],
        :V_pulse        => fid.groups[chan][Ind].props["Voltage"],
        :PRF_us         => fid.groups[chan][Ind].props["PRF (us)"])
        trace3 = fid.groups[chan][Ind].data
        return [trace1, trace2, trace3], BR_Rx, TR_Rx, TR_Tx
end
"""
Convert labview time to Julia time
"""
function t_conv(t_s::Float64)
        t_ref = datetime2unix(DateTime(1904, 1, 1, 0, 0, 0))
        t_out = unix2datetime(t_s+t_ref)
end
