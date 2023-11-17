
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
    for v=1:Vs
        @constraint(model, y_ss[v,v,:] .== 0) # pas d'arête boucle sur le graphe des sous-stations
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
        Nmin_used_substations = 4;
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
        Nmin_used_substations = 3;
        Nmin_links = 2*2;
        used_substations = sort([29,22,15,8,1, 35,28,21,14,7, 32,25,18,11,4, 34,27,20,13,6]); #31,24,17,10,3]);
    elseif instance_name == "huge"
        dim_x = 5
        dim_y = 25
        Nmin_used_substations = 8;
        Nmin_links = 2*4;
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

    
    @constraint(model, x * rs .>= pow_max * z * onesVt); # empêcher tout curtailing SANS failure de sous-station
    # remplace le coût éq. (5) -> contraindre à 0:  i.e. chaque station doit pouvoir supporter la charge de toutes les éoliennes qui lui sont connectées
    
    @constraint(model, y_ls * rq_ls .>= pow_max * z * onesVt) # chaque câble entre la station principale et une sous-station doit pouvoir supporter la charge de toutes les éoliennes qui lui sont connectées

    # for v in used_substations
    #     @constraint(model, pow_max * (z[v,:]'*onesVt) - sum(y_ss[v,vbar,:]' * rq_ss for vbar in (v+1):Vs) <= 5000); # empêcher tout curtailing AVEC failure de sous-station
    #     # remplace le coût éq. (6) -> contraindre à 0:  i.e. chaque station doit pouvoir être prise en charge par les autres qui lui sont connectées si elle casse
    # end
    

    # curt_v_under_ω_fail_v = @expression(model, pow_max * (z[v,:]'*onesVt) - sum(y_ss[v,vbar,:]' * rq_ss for vbar in (v+1):Vs)); # curtailing of v under ω and failure of v [summed over vbar]
    # pow_gen_turb_link_v = @expression(model, πw[Ω] *  sum(z[:,:])); # power generated by turbines linked to vbar [summed over vbar]
    # pow_sent_v_to_vbar = @expression(model, (0.5*sum(y_ss[v,:,:]*rq_ss) +  πw[Ω]*sum(z[v,:]))); # power sent from v to vbar [summed over vbar]
    # capa_cabl_subs_link_v = @expression(model, sum(x[:,:]*rs) + sum(y_ls[:,:]*rq_ls)); # capacity of the cable/substation linking vbar [summed over vbar]

    activate_no_fail_all_avgsum_curtail = true; # true pour tous SAUF medium, avec la pondération égale (1/2 1/2)
    if activate_no_fail_all_avgsum_curtail
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

    activate_op_cost = true; # huge: oui, le reste non
    if activate_op_cost

        operational_cost = QuadExpr();
        cst_mlt = 5*c0; # constant multiplier
        
        for v in used_substations_except_Vs # to avoid MethodError of empty ranges
            # boucle auparavant faite sur un cropped_range = deleteat!(copy(used_substations), findall(x->x==Vs,used_substations)) 
            # i.e.  used substations dans 1:Vs-1
            
            # premier terme de (6) sans le relu
            add_to_expression!(operational_cost,   cst_mlt * πw[Ω] * (z[v,:]'*onesVt));

            add_to_expression!(operational_cost, - cst_mlt * sum(y_ss[v,vbar,:]' * rq_ss for vbar in (v+1):Vs));


        # for vbar in used_substations
        #     if vbar != v && vbar > v
        #         # second terme de (6) sans le relu ni les min
        #         # pas de facteur 1/2 ici puisque la somme se fait déjà exactement sur les couples uniques dans (v+1):Vs
        #         add_to_expression!(operational_cost,  cst_mlt * πw[Ω] *  sum(z[vbar,:])); # power generated by turbines linked to v0
        #         add_to_expression!(operational_cost,  cst_mlt * (y_ss[v,vbar,:]'*rq_ss +  πw[Ω]*sum(z[v,:])));  # power sent from v to vbar
        #         add_to_expression!(operational_cost, -cst_mlt * (x[vbar,:]'*rs +  y_ls[vbar,:]'*rq_ls)); # capacity of the cable/subsattion linking vbar
        #     end
        # end
    end

        total_cost = construction_cost   + operational_cost;
    else
        total_cost = construction_cost; 
    end











