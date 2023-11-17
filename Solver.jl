using JSON;
using LinearAlgebra;
using JuMP, SCIP, Gurobi # HiGHS, GLPK

include("Import.jl");



path = "/home/matthias/Documents/GitHub/kiro-challenge-2023/";
instance_name = "huge";


file = path * "instances/"*instance_name*".json";


instance, Q0, Qs, Vs, Ss, Nws, Vt,
c0, cp, ctf, ctl, main_land_station_coords, Cmax, Pmax,
cqf_ls, cql_ls, rq_ls, pq,
cqf_ss, cql_ss, rq_ss,
substation_coords,
cs, rs, ps,
dsw, dss, dls,
πw, pw, wind_turbine_coords  = import_init(file);


if false
    visualize_overall(instance_name, path*"visualisations/")
    visualize_turbines(instance_name, path*"visualisations/")
end










    return model, total_cost, x, y_ss, z, y_ls, onesSs, onesVs, onesQ0, onesVt, used_substations
end

# include("Models.jl");


@time model, total_cost, x, y_ss, z, y_ls, onesSs, onesVs, onesQ0, onesVt, used_substations = get_model2(SCIP.Optimizer);

set_time_limit_sec(model, 50.0)

@objective(model, Min, total_cost); # fonction objectif

optimize!(model);

include("Metrics.jl");

@time cost_loss = loss(x, y_ss, z, y_ls)



include("Visualize.jl")
visualize_result(x, y_ss, z, y_ls, path*"visualisations/", "scip", round(cost_loss[4], digits=2))