using JSON;
using LinearAlgebra;
using DataStructures;
# using DataFrames;





function import_init(file::String)

    instance = JSON.parsefile(file)



    ## GENERAL PARAMETERS
    curtailing_cost    = instance["general_parameters"]["curtailing_cost"];
    curtailing_penalty = instance["general_parameters"]["curtailing_penalty"];
    fixed_cost_cable   = instance["general_parameters"]["fixed_cost_cable"];


    main_land_station_coords = zeros(2);
    main_land_station_coords[1]  =  instance["general_parameters"]["main_land_station"]["x"];
    main_land_station_coords[2]  =  instance["general_parameters"]["main_land_station"]["y"];

    maximum_curtailing  = instance["general_parameters"]["maximum_curtailing"];
    maximum_power       = instance["general_parameters"]["maximum_power"];
    variable_cost_cable = instance["general_parameters"]["variable_cost_cable"];


    Q0       = length(instance["land_substation_cable_types"]);
    Qs = length(instance["substation_substation_cable_types"]);
    Vs            = length(instance["substation_locations"]);
    Ss                  = length(instance["substation_types"]);
    N_wind_scenarios                    = length(instance["wind_scenarios"]);
    Nwt                     = length(instance["wind_turbines"]);



    ## LAND-SUBSTATION CABLE DATA
    land_sub_cable_costs         = zeros(Float64, Q0, 2);
    land_sub_cable_ratings       = zeros(Float64, Q0);
    land_sub_cable_proba_failure = zeros(Float64, Q0);

    for i=1:Q0
        land_sub_cable_costs[i,1]              = Float64(instance["land_substation_cable_types"][i]["fixed_cost"]);
        land_sub_cable_costs[i,2]              = Float64(instance["land_substation_cable_types"][i]["variable_cost"]);
        land_sub_cable_proba_failure[i]        = Float64(instance["land_substation_cable_types"][i]["probability_of_failure"]);
        land_sub_cable_ratings[i]              = Float64(instance["land_substation_cable_types"][i]["rating"]);
    end 



    ## SUBSTATION-SUBSTATION CABLE DATA
    sub_sub_cable_costs  = zeros(Float64, Qs, 2);
    sub_sub_cable_rating            = zeros(Float64, Qs);

    for i=1:Qs
        sub_sub_cable_costs[i,1]   = Float64(instance["substation_substation_cable_types"][i]["fixed_cost"]);
        sub_sub_cable_costs[i,2]   = Float64(instance["substation_substation_cable_types"][i]["variable_cost"]);
        sub_sub_cable_rating[i]               = Float64(instance["substation_substation_cable_types"][i]["rating"]);
    end 


    ## SUBSTATION LOCATIONS COORDS

    substation_location_coords = zeros(Float64, Vs, 2);

    for i=1:Vs
        substation_location_coords[i,1]  = Float64(instance["substation_locations"][i]["x"]);
        substation_location_coords[i,2]  = Float64(instance["substation_locations"][i]["y"]);
    end 



    ## SUBSTATION TYPES DATA
    substation_cost   = zeros(Float64, Ss);
    substation_ratings = zeros(Float64, Ss);
    substation_proba_failure = zeros(Float64, Ss);

    for i=1:Ss
        substation_cost[i]            = Float64(instance["substation_types"][i]["cost"]);
        substation_proba_failure[i]   = Float64(instance["substation_types"][i]["probability_of_failure"]);
        substation_ratings[i]         = Float64(instance["substation_types"][i]["rating"]);
    end 

    ## WIND SCENARIOS DATA

    wind_scenario_pow_gen = zeros(Float64, N_wind_scenarios);
    wind_scenario_proba   = zeros(Float64, N_wind_scenarios);

    for i=1:N_wind_scenarios
        wind_scenario_pow_gen[i] = Float64(instance["wind_scenarios"][i]["power_generation"]);
        wind_scenario_proba[i]   = Float64(instance["wind_scenarios"][i]["probability"]);
    end



    ## WIND TURBINES DATA

    wind_turbine_coords = zeros(Float64, Nwt, 2);

    for i=1:Nwt
        wind_turbine_coords[i,1] = Float64(instance["wind_turbines"][i]["x"]);
        wind_turbine_coords[i,2] = Float64(instance["wind_turbines"][i]["y"]);
    end

return Q0, Qs, Vs, Ss, N_wind_scenarios, Nwt,
 
curtailing_cost, curtailing_penalty, fixed_cost_cable, main_land_station_coords, maximum_curtailing, maximum_power, variable_cost_cable,

land_sub_cable_costs, land_sub_cable_ratings, land_sub_cable_proba_failure,

sub_sub_cable_costs, sub_sub_cable_rating,

substation_location_coords,

substation_cost, substation_ratings, substation_proba_failure,

wind_scenario_pow_gen, wind_scenario_proba, wind_turbine_coords

end




path = "/home/matthias/Documents/GitHub/kiro-challenge-2023/";
file = path * "instances/toy.json";




Q0, Qs, Vs, Ss, N_wind_scenarios, Nwt,
 
curtailing_cost, curtailing_penalty, fixed_cost_cable, main_land_station_coords, maximum_curtailing, maximum_power, variable_cost_cable,

land_sub_cable_costs, land_sub_cable_ratings, land_sub_cable_proba_failure,

sub_sub_cable_costs, sub_sub_cable_rating,

substation_location_coords,

substation_cost, substation_ratings, substation_proba_failure,

wind_scenario_pow_gen, wind_scenario_proba, wind_turbine_coords  = import_init(file);


# create adjacency matrix between turbines substations and turbines
d_sub_turb = zeros(Vs, Nwt);
for i=1:Vs
    for j=1:Nwt
        d_sub_turb[i,j] = norm(substation_location_coords[i,:] - wind_turbine_coords[j,:]);
    end
end

# create adjacency matrix between substations and substations
d_sub_sub = zeros(Vs, Vs);
for i=1:Vs
    for j=1:Vs
        d_sub_sub[i,j] = norm(substation_location_coords[i,:] - substation_location_coords[j,:]);
    end
end


# create adjacency matrix between substations and main land station
d_sub_main = zeros(Vs);
for i=1:Vs
    d_sub_main[i] = norm(substation_location_coords[i,:] - main_land_station_coords);
end

using JuMP, SCIP, HiGHS, GLPK

function get_model()
    model = Model(SCIP.Optimizer);
    # substations
    @variable(model, x[1:Vs, 1:Ss], Bin); # substation présente (emplacement, type)



    # câbles
    @variable(model, z[1:Vs, 1:Nwt], Bin) # câble présent (substation <--> éolienne)
    @variable(model, y[1:Vs, 1:Vs, 1:Qs],  Bin)  # câble présent (substation <--> substation, type)
    @variable(model, u[1:Vs, 1:Q0], Bin)         # câble présent ([v0] <--> substation, type)


    # contraintes
    onesSs = ones(Int, Ss);
    onesVs = ones(Int, Vs);
    @constraint(model, x * onesSs .<= onesVs) #(1) au plus une sous-station par emplacement
    for i in 1:Vs # at most one cable type per edge between substations
        for j in 1:Vs
            @constraint(model, sum(y[i, j, k] for k in 1:Qs) <= 1)
        end
    end

    for v=1:Vs
        @constraint(model, sum(u[]))
    end

    return model
end


model = get_model();