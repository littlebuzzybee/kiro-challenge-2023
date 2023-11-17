function import_init(file::String)

    instance = JSON.parsefile(file)

    ## GENERAL PARAMETERS
    c0    = instance["general_parameters"]["curtailing_cost"];
    cp    = instance["general_parameters"]["curtailing_penalty"];
    ctf   = instance["general_parameters"]["fixed_cost_cable"]; # for a substation-turbine cable
    ctl   = instance["general_parameters"]["variable_cost_cable"]; # for a substation-turbine cable

    ## MAIN LAND STATION
    main_land_station_coords = zeros(2);
    main_land_station_coords[1]  =  instance["general_parameters"]["main_land_station"]["x"];
    main_land_station_coords[2]  =  instance["general_parameters"]["main_land_station"]["y"];

    Cmax  = instance["general_parameters"]["maximum_curtailing"];
    Pmax  = instance["general_parameters"]["maximum_power"];



    Q0   = length(instance["land_substation_cable_types"]); # number of land-substation cable types
    Qs   = length(instance["substation_substation_cable_types"]); # number of substation-substation cable types
    Vs   = length(instance["substation_locations"]); # number of substations
    Ss   = length(instance["substation_types"]); # number of substation types
    Nws  = length(instance["wind_scenarios"]);   # number of wind scenarios
    Vt   = length(instance["wind_turbines"]);    # number of wind turbines



    ## LAND-SUBSTATION CABLES DATA
    cqf_ls  = zeros(Float64, Q0);
    cql_ls  = zeros(Float64, Q0);
    pq      = zeros(Float64, Q0);
    rq_ls   = zeros(Float64, Q0);
    
    for i=1:Q0
        cqf_ls[i]   = Float64(instance["land_substation_cable_types"][i]["fixed_cost"]);
        cql_ls[i]   = Float64(instance["land_substation_cable_types"][i]["variable_cost"]);
        pq[i]       = Float64(instance["land_substation_cable_types"][i]["probability_of_failure"]);
        rq_ls[i]    = Float64(instance["land_substation_cable_types"][i]["rating"]);
    end 



    ## SUBSTATION-SUBSTATION CABLES DATA
    cqf_ss  = zeros(Float64, Qs);
    cql_ss  = zeros(Float64, Qs);
    rq_ss   = zeros(Float64, Qs);

    for i=1:Qs
        cqf_ss[i]   = Float64(instance["substation_substation_cable_types"][i]["fixed_cost"]);
        cql_ss[i]   = Float64(instance["substation_substation_cable_types"][i]["variable_cost"]);
        rq_ss[i]    = Float64(instance["substation_substation_cable_types"][i]["rating"]);
    end 


    ## SUBSTATION LOCATIONS COORDS
    substation_coords = zeros(Float64, Vs, 2);

    for i=1:Vs
        substation_coords[i,1]  = Float64(instance["substation_locations"][i]["x"]);
        substation_coords[i,2]  = Float64(instance["substation_locations"][i]["y"]);
    end 

    

    ## SUBSTATION TYPES DATA
    cs = zeros(Float64, Ss);
    rs = zeros(Float64, Ss);
    ps = zeros(Float64, Ss);

    for i=1:Ss
        cs[i]   = Float64(instance["substation_types"][i]["cost"]);
        ps[i]   = Float64(instance["substation_types"][i]["probability_of_failure"]);
        rs[i]   = Float64(instance["substation_types"][i]["rating"]);
    end 

    ## WIND SCENARIOS DATA
    πw = zeros(Float64, Nws);
    pw = zeros(Float64, Nws);

    for i=1:Nws
        πw[i] = Float64(instance["wind_scenarios"][i]["power_generation"]);
        pw[i] = Float64(instance["wind_scenarios"][i]["probability"]);
    end



    ## WIND TURBINES DATA
    wind_turbine_coords = zeros(Float64, Vt, 2);

    for i=1:Vt
        wind_turbine_coords[i,1] = Float64(instance["wind_turbines"][i]["x"]);
        wind_turbine_coords[i,2] = Float64(instance["wind_turbines"][i]["y"]);
    end

    ## DISTANCE MATRICES
    # create distance matrix between turbines substations and turbines
    dsw = zeros(Vs, Vt);
    for i=1:Vs
        for j=1:Vt
            dsw[i,j] = norm(substation_coords[i,:] - wind_turbine_coords[j,:]);
        end
    end

    # create distance matrix between substations and substations
    dss = zeros(Vs, Vs);
    for i=1:Vs
        for j=1:Vs
            dss[i,j] = norm(substation_coords[i,:] - substation_coords[j,:]);
        end
    end


    # create distance matrix between substations and main land station
    dls = zeros(Vs);
    for i=1:Vs
        dls[i] = norm(substation_coords[i,:] - main_land_station_coords);
    end

return instance, Q0, Qs, Vs, Ss, Nws, Vt,
c0, cp, ctf, ctl, main_land_station_coords, Cmax, Pmax, 
cqf_ls, cql_ls, rq_ls, pq,    # land <--> substations
cqf_ss, cql_ss, rq_ss,        # substations <--> substations
substation_coords,
cs, rs, ps,                   # substations
dsw, dss, dls,                # distances entre les différents nœuds
πw, pw, wind_turbine_coords
end