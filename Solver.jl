using JSON;
using LinearAlgebra;
using JuMP, SCIP, Gurobi # HiGHS, GLPK

include("Import.jl");



path = "/home/matthias/Documents/GitHub/kiro-challenge-2023/";
instance_name = "large";


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


function get_model2(optimizer)
    
    ## MODEL & VARIABLES
    model = Model(optimizer);

    ### substations
    @variable(model, x[1:Vs, 1:Ss], Bin); # substation présente (emplacement, type)
    ### câbles
    @variable(model, z[1:Vs, 1:Vt], Bin)  # câble présent (substation <--> éolienne)
    @variable(model, y_ss[1:Vs, 1:Vs, 1:Qs],  Bin) # câble présent (substation <--> substation, type)
    @variable(model, y_ls[1:Vs, 1:Q0], Bin) # câble présent ([land] <--> substation, type)
    
    for q=1:Qs
        @constraint(model, y_ss[:,:,q] .== y_ss[:,:,q]') 
        # symétrie pour les câbles entre sous-stations
    end
    for v1=1:Vs
        @constraint(model, y_ss[v1,v1,:] .== 0)
        # pas d'arête boucle sur le graphe des sous-stations
    end
    

       
    
    # get dimensions of the turbine array, and define the minimum number of substations and links to use
    if instance_name == "toy"
        dim_x = 1
        dim_y = 1
        Nmin_used_substations = 1;
        Nmin_links = 2*0;
        used_substations = sort(Vector(1:Vs));
    elseif instance_name == "small"
        dim_x = 5
        dim_y = 15
        Nmin_used_substations = 5;
        Nmin_links = 2*2;
        used_substations = sort(Vector(1:Vs));
    elseif instance_name == "medium"
        dim_x = 9
        dim_y = 15
        Nmin_used_substations = 3;
        Nmin_links = 2*1;
        used_substations = sort([3,8,13,18,23,28]);
    elseif instance_name == "large"
        dim_x = 8
        dim_y = 16
        Nmin_used_substations = 6;
        Nmin_links = 2*3;
        used_substations = sort([29,22,15,8,1, 35,28,21,14,7, 32,25,18,11,4, 34,27,20,13,6]); #31,24,17,10,3]);
    elseif instance_name == "huge"
        dim_x = 5
        dim_y = 25
        Nmin_used_substations = 8;
        Nmin_links = 1*3;
        used_substations = sort([4,13,22,31,40,49,58,67,76]);
    end

    used_substations_except_Vs = deleteat!(copy(used_substations), findall(x->x==Vs, used_substations));

    # limit the number of used substations to a subset to cut through the search space
    unused_substations = deleteat!(collect(1:Vs), used_substations);
    @constraint(model, x[unused_substations,:] .== 0);
    @constraint(model, z[unused_substations,:] .== 0);
    @constraint(model, y_ss[unused_substations,:,:] .== 0);
    @constraint(model, y_ls[unused_substations,:] .== 0);


    # force all wind turbines on a same line to be connected to the same substation
    if instance_name in ["large", "huge", "medium"]
        for k=0:dim_y-1 
            for t=1:dim_x-1
                @constraint(model, z[:,k*dim_x+t] .== z[:,k*dim_x+t+1]);
            end
        end
    end


    ### CONSTANTES AUXILIAIRES POUR LA RÉDUCTION DES SOMMES
    onesSs = ones(Int, Ss);
    onesVs = ones(Int, Vs);
    onesQ0 = ones(Int, Q0);
    onesVt = ones(Int, Vt);


    ## CONSTRAINTS
    @constraint(model, x * onesSs .<= onesVs) #(1) au plus une sous-station par emplacement


    ### Cables between land and substations
    for v=1:Vs
        @constraint(model, sum(y_ls[v, q] for q in 1:Q0) == sum(x[v, s] for s in 1:Ss)); #  y_ls * onesQ0 .== x * onesSs
        #(2) une sous-station est présente <=> un câble doit la relier à la station principale
    end

    ### Cables between substations and turbines
    for t=1:Vt
        @constraint(model, sum(z[v, t] for v in 1:Vs) == 1);
        #(3) chaque éolienne doit être reliée à exactement une sous-station
    end
    
    # Cables between substations and substations
    for v=1:Vs
        @constraint(model, sum(y_ss[v, vprime, q] for vprime in (v+1):Vs, q in 1:Qs) <= sum(x[v, s] for s in 1:Ss));
        #(4) chaque sous-station est reliée à au plus une autre sous-station
    end

    ### Additional constraints
    
    # (condition non écrite dans le sujet mais essentielle): si un câble relie une sous-station à une éolienne, alors la sous-station doit être présente
    # contrainte supprimée car contrainte pénible pour le solveur, qui est déjà assurée par la contrainte sur la charge (chaque station doit pouvoir supporter les éoliennes connectées donc le produit de la capacité par la présence est minoré donc  présence > 0)
#=     for v=1:Vs
        for t=1:Vt
            @constraint(model, z[v,t] => {sum(x[v,:]) >= 1});
        end 
    end =#

    Ω = argmax(πw);
    pow_max = πw[Ω];

    @constraint(model, z * onesVt .>= x * onesSs); # si une sous-station est présente, elle doit être reliée à des turbines

    
    @constraint(model, x * rs .>= pow_max * z * onesVt); # empêcher tout curtailing sans failure de sous-station
    # remplace le coût éq. (5) -> contraindre à 0:  i.e. chaque station doit pouvoir supporter la charge de toutes les éoliennes qui lui sont connectées
    
    @constraint(model, y_ls * rq_ls .>= pow_max * z * onesVt) # chaque câble entre la station principale et une sous-station doit pouvoir supporter la charge de toutes les éoliennes qui lui sont connectées
    
    

    # curt_v_under_ω_fail_v = @expression(model, pow_max * (z[v,:]'*onesVt) - sum(y_ss[v,vbar,:]' * rq_ss for vbar in (v+1):Vs)); # curtailing of v under ω and failure of v [summed over vbar]
    # pow_gen_turb_link_v = @expression(model, πw[Ω] *  sum(z[:,:])); # power generated by turbines linked to vbar [summed over vbar]
    # pow_sent_v_to_vbar = @expression(model, (0.5*sum(y_ss[v,:,:]*rq_ss) +  πw[Ω]*sum(z[v,:]))); # power sent from v to vbar [summed over vbar]
    # capa_cabl_subs_link_v = @expression(model, sum(x[:,:]*rs) + sum(y_ls[:,:]*rq_ls)); # capacity of the cable/substation linking vbar [summed over vbar]

    activate_no_curtail = false; # medium: non, huge: oui, large: non
    if activate_no_curtail
        bound = Cmax;
        curt = QuadExpr();
        # remplace le coût éq. (6): empêcher tout curtailing supérieur à Cmax de v en tout scénario avec failure de v
        for v in used_substations
            if v < Vs
                add_to_expression!(curt, pow_max * (z[v,:]'*onesVt) - sum(y_ss[v,vbar,:]' * rq_ss for vbar in (v+1):Vs)); # curtailing of v under ω and failure of v
            end
            for vbar in used_substations
                if vbar > v
                    add_to_expression!(curt, πw[Ω] *  sum(z[vbar,:])); # power generated by turbines linked to vbar 
                    add_to_expression!(curt,   0.5*(sum(y_ss[v,vbar,:]'*rq_ss) +  πw[Ω]*sum(z[v,vbar]))); # power sent from v to vbar
                    add_to_expression!(curt, - 0.5*(sum(x[vbar,:]'*rs) + sum(y_ls[vbar,:]'*rq_ls))); # capacity of the cable/substation linking vbar
                    end
            end
        end
        drop_zeros!(curt);
        @constraint(model, curt <=  bound)
    end

    


    @constraint(model, sum(y_ss) >= Nmin_links); # forcer des liaisons entre sous stations /!\ attention, bien se rappeler que y_ss est symétrique (imposer deux fois le nombre voulu de liaisons)
    @constraint(model, sum(x) >= Nmin_used_substations); # forcer l'utilisation d'au moins x substations

    
    for v=1:Vs
        @constraint(model, sum(y_ss[v,:,:]) .<= sum(x[v,:])); # interdire une liaison entre une sous-station et une autre si la première n'est pas présente i.e. si une liaison entre deux sous-stations est présente, alors les deux sous-stations doivent être présentes
    end
    
    
    ## CONSTRUCTION COSTS
    construction_cost = AffExpr();
    add_to_expression!(construction_cost, onesVs' * y_ls * cqf_ls + dls' * y_ls * cql_ls); # coûts des câbles entre land et sous-stations

    for q=1:Qs
        add_to_expression!(construction_cost, .5 *  sum(cql_ss[q] * y_ss[:,:,q] .* dss)); # coûts linéaires des câbles entre sous-stations
        add_to_expression!(construction_cost, .5 *  sum(cqf_ss[q] * y_ss[:,:,q]));        # coûts fixes des câbles entre sous-stations
        # facteur 1/2 car la matrice est symétrique donc on compte deux fois chaque câble
    end

    add_to_expression!(construction_cost, onesVs' * x * cs); # coûts de construction des sous-stations
    add_to_expression!(construction_cost, sum(ctf * z + ctl * dsw .* z)); # coûts des câbles entre land et éoliennes




    # OPERATIONAL COSTS

    activate_cost_curtailing = true; # medium: non, huge: non, large: non
    if activate_cost_curtailing # enlever le operational_cost pour des résultats meilleurs avec large, medium et small

        operational_cost = QuadExpr();
        cst_mlt = c0; # constant multiplier
        
        for v in used_substations_except_Vs # to avoid MethodError of empty ranges
            # boucle auparavant faite sur un cropped_range = deleteat!(copy(used_substations), findall(x->x==Vs,used_substations)) 
            # i.e.  used substations dans 1:Vs-1
            
            # premier terme de (6) sans le relu
            add_to_expression!(operational_cost,   cst_mlt * πw[Ω] * (z[v,:]'*onesVt));

            add_to_expression!(operational_cost, - cst_mlt * sum(y_ss[v,vbar,:]' * rq_ss for vbar in (v+1):Vs));


        for vbar in used_substations
                if vbar != v && vbar > v
                    # second terme de (6) sans le relu ni les min
                    # pas de facteur 1/2 ici puisque la somme se fait déjà exactement sur les couples uniques dans (v+1):Vs
                    add_to_expression!(operational_cost,  cst_mlt * πw[Ω] *  sum(z[vbar,:])); # power generated by turbines linked to v0
                    add_to_expression!(operational_cost,  cst_mlt * (y_ss[v,vbar,:]'*rq_ss +  πw[Ω]*sum(z[v,:])));  # power sent from v to vbar
                    add_to_expression!(operational_cost, -cst_mlt * (x[vbar,:]'*rs +  y_ls[vbar,:]'*rq_ls)); # capacity of the cable/subsattion linking vbar
                end
            end
        end

        total_cost = construction_cost   + operational_cost;
    else
        total_cost = construction_cost; 
    end


    return model, total_cost, x, y_ss, z, y_ls, onesSs, onesVs, onesQ0, onesVt, used_substations
end



@time model, total_cost, x, y_ss, z, y_ls, onesSs, onesVs, onesQ0, onesVt, used_substations = get_model2(Gurobi.Optimizer);

set_time_limit_sec(model, 50.0)

@objective(model, Min, total_cost); # fonction objectif

optimize!(model);

include("Metrics.jl");

@time cost_loss = loss(x, y_ss, z, y_ls)



include("Visualize.jl")
visualize_result(x, y_ss, z, y_ls, path*"visualisations/", "gurobi_temp", round(cost_loss[4], digits=2))