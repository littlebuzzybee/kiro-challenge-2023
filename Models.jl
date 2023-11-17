
function get_model1(optimizer)
    ## MODEL & VARIABLES
    model = Model(optimizer);

    ### substations
    @variable(model, x[1:Vs, 1:Ss], Bin); # substation présente (emplacement, type)

    ### câbles
    @variable(model, z[1:Vs, 1:Vt], Bin)         # câble présent (substation <--> éolienne)
    #@variable(model, y_ss[v1=1:Vs, v2=(v1+1):Vs, 1:Qs],  Bin)  
    @variable(model, y_ss[1:Vs, 1:Vs, 1:Qs],  Bin) # câble présent (substation <--> substation, type)
    for q=1:Qs
        @constraint(model, y_ss[:,:,q] .== y_ss[:,:,q]') # symétrie
    end
    for v1=1:Vs
        @constraint(model, y_ss[v1,v1,:] .== 0) # pas d'arête boucle sur le graphe des sous-stations
    end
    @variable(model, y_ls[1:Vs, 1:Q0], Bin) # câble présent ([land] <--> substation, type)

    @show size(z)

    ### Building substations
    onesSs = ones(Int, Ss);
    onesVs = ones(Int, Vs);
    onesQ0 = ones(Int, Q0);
    onesVt = ones(Int, Vt);


    ## CONSTRAINTS
    @constraint(model, x * onesSs .<= onesVs) #(1) au plus une sous-station par emplacement

    #=     for v1 in 1:Vs # redondant / inutile
        for v2 in 1:Vs
            @constraint(model, sum(y_ss[v1, v2, q] for q in 1:Qs) <= 1)
            # au plus un seul (type de) câble entre deux sous-stations
        end
    end =#

    ###Cables between land and substations
    for v=1:Vs
        @constraint(model, sum(y_ls[v, q] for q in 1:Q0) == sum(x[v, s] for s in 1:Ss));
        #(2) une sous-station est présente <=> un câble doit la relier à la station principale
    end

    ### Cables between substations and turbines
    for t=1:Vt
        @constraint(model, sum(z[v, t] for v in 1:Vs) == 1);
        #(3) chaque éolienne doit être reliée à exactement une sous-station
    end

    
    # Cables between substations and substations
    for v=1:Vs
        @constraint(model, sum(y_ss[v1, v2, q] for v1 in 1:Vs, v2 in (v1+1):Vs, q in 1:Qs) <= sum(x[v, s] for s in 1:Ss));
        #(4) chaque sous-station est reliée à au plus une autre sous-station
    end

    ### Additional constraints
    
    for v=1:Vs
        for t=1:Vt
            @constraint(model, z[v,t] => {sum(x[v,:]) >= 1});
            # (condition non écrite dans le sujet mais essentielle): si un câble relie une sous-station à une éolienne, alors la sous-station doit être présente
            # est prise en compte par les suivantes (désactivées, explorées dans une version ultérieure de cette fonction)
        end 
    end
    
    
    @constraint(model, x * rs .>= Pmax * z * onesVt);
    # chaque station doit pouvoir supporter la charge de toutes les éoliennes qui lui sont connectées
    @constraint(model, y_ls * rq_ls .>= Pmax * z * onesVt)
    # chaque câble entre la station principale et une sous-station doit pouvoir supporter la charge de toutes les éoliennes qui lui sont connectées
    
    # @constraint(model, sum(x) >= 2); # forcer l'utilisation d'au moins 5 substations
    # @constraint(model, sum(y_ss) >= 1); # forcer au moins 4 câbles les reliant entre elles
    # @constraint(model, x * onesSs  .<= z * onesVt) # présence station => présence éiloenne connectée à cette station

    ## OBJECTIVES & COSTS

    cost_ls = onesVs' * y_ls * cqf_ls + dls' * y_ls * cql_ls; # coûts des câbles entre land et sous-stations


    cost_ss = AffExpr();
    for q=1:Qs
        cost_ss += .5 *  sum(cql_ss[q] * y_ss[:,:,q] .* dss); # coûts linéaires des câbles entre sous-stations
        cost_ss += .5 *  sum(cqf_ss[q] * y_ss[:,:,q]);        # coûts fixes des câbles entre sous-stations
        # facteur 1/2 car la matrice est symétrique donc on compte deux fois chaque câble
    end



    cost_build_stat = onesVs' * x * cs; # coûts de construction des sous-stations

    cost_sw = sum(ctf * z + ctl * dsw .* z); # coûts des câbles entre land et éoliennes

    construction_cost = cost_ls + cost_ss +  cost_build_stat + cost_sw;


     
    return model, construction_cost, x, y_ss, z, y_ls, onesSs, onesVs, onesQ0, onesVt
end
