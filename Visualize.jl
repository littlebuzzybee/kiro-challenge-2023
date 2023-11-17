using Plots;

function visualize_overall(instance_name, path)
    plotd = plot(substation_coords[:,1], substation_coords[:,2], seriestype=:scatter, label=string(Vs)*" substations", dpi=300)
    plot!(wind_turbine_coords[:,1], wind_turbine_coords[:,2], seriestype=:scatter, label=string(Vt)*" turbines")
    plot!(main_land_station_coords[1,:], main_land_station_coords[2,:], seriestype=:scatter, label="land station")

    nums = [(substation_coords[s,1], substation_coords[s,2], (string(s), 8, 2.0, :bottom, :red)) for s in 1:Vs];
    annotate!(nums);


    savefig(plotd, path*instance_name*"_config.png")
end


function visualize_turbines(instance_name, path)
    plotd = plot(wind_turbine_coords[:,1], wind_turbine_coords[:,2], seriestype=:scatter, label=string(Vt)*" turbines", dpi=300)

    nums = [(wind_turbine_coords[t,1], wind_turbine_coords[t,2], (string(t), 8, 2.0, :bottom, :red)) for t in 1:Vt];
    annotate!(nums)

    savefig(plotd, path*instance_name*"_turbines.png")
end


function visualize_result(x, y_ss, z, y_ls, path, optimizer_name, loss)
    x = value.(x);
    y_ss = convert(y_ss);
    y_ls = value.(y_ls);
    z = value.(z);

    # visualize_overall

    plotd = plot(substation_coords[:,1], substation_coords[:,2], seriestype=:scatter, label=string(Vs)*" substations", dpi=300, title=optimizer_name*" loss="*string(loss))
    plot!(wind_turbine_coords[:,1], wind_turbine_coords[:,2], seriestype=:scatter, label=string(Vt)*" turbines")
    plot!(main_land_station_coords[1,:], main_land_station_coords[2,:], seriestype=:scatter, label="land station")

    nums = [(substation_coords[s,1], substation_coords[s,2], (string(s), 8, 2.0, :bottom, :red)) for s in 1:Vs];
    annotate!(nums)

    x_c = sum(x, dims=2); # sum over substation types
    for v=1:Vs
        # replot substation in grey id it was not built
        if sum(x_c[v]) < 0.5
            plot!(substation_coords[v,1:1], substation_coords[v,2:2], seriestype=:scatter, label="", color=:grey)
        end
    end
    



    y_ss_c = sum(y_ss, dims=3); # sum over cable types
    for v1=1:Vs
        for v2=1:Vs
            if y_ss_c[v1,v2] > 0
                plot!([substation_coords[v1,1], substation_coords[v2,1]], [substation_coords[v1,2], substation_coords[v2,2]], color=:red, arrow=false, linewidth=2, label="")
            end
        end
    end

    for v=1:Vs
        for t=1:Vt
            if z[v,t] > 0.5
                plot!([substation_coords[v,1], wind_turbine_coords[t,1]], [substation_coords[v,2], wind_turbine_coords[t,2]], color=:blue, arrow=false, linewidth=1, label="")
            end
        end
    end

    for v=1:Vs
        for q=1:Q0
            if y_ls[v,q] > 0.5
                plot!([main_land_station_coords[1], substation_coords[v,1]], [main_land_station_coords[2], substation_coords[v,2]], color=:green, arrow=false, linewidth=1, label="")
            end
        end
    end

    savefig(plotd, path*instance_name*"_result_"*optimizer_name*".png")
end

