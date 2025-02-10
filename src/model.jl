function GTAPModel(data; resName = :USA, capName = :cap, savfClosure = :capFix)

    @extract(data,
        acts, endw, comm, reg,
        A_qtm, alphac, au, aug, aus,
        conshr0,
        depr, dpgov, dppriv, dpsave,
        ESUBC, ESUBD, ESUBG, ESUBI, ESUBM, ESUBQ, ESUBS, ESUBT, ESUBVA, ETAFF, ETRAE, ETRAQ, evfb, evfp, evos,
        fincome0,
        gdpfc0,
        ifRes, INCPAR, indtax0, invwgt,
        kb0, ke0,
        makb, maks, mobFact,
        pca0, pcif0, pds0, pe0, peb0, pes0, pfa0, pfact0, pfactw0, pfd0, pfe0, pfm0, pfob0, pga0, pgd0, pgm0,
        pgov0, pia0, pid0, pim0, pint0, pinv0, pinvLag, pmds0, pms0, po0, pop0, ppa0, ppd0, ppm0, ps0, psaveLag,
        pt0, ptax, ptax0, ptrans0, pva0,
        qc0, qca0, qds0, qe0, qfa0, qfd0, qfe0, qfm0, qga0, qgd0, qgm0, qia0, qid0, qim0, qint0, qinv0, qms0,
        qo0, qpa0, qpd0, qpm0, qst0, qtm0, qtmfsd0, qva0, qxs0,
        rental0, risk, rorc0, rore0, RORFLEX, rorg0,
        save, save0, savf0, savwgt, shr_gov, shr_inv, shr_qcad, shr_qcas, shr_qfa, shr_qfd, shr_qfe, shr_qfes, 
        shr_qfm, shr_qgd, shr_qgm, shr_qid, shr_qim, shr_qint, shr_qpd, shr_qpm, shr_qst, shr_qva, shr_qxs, 
        shr_vtwr, SUBPAR,
        taxrexp0, taxrfu0, taxrgc0, taxric0, taxrimp0, taxriu0, taxrout0, taxrpc0, tfd0, tfe0, tfm0, tgd0,
        tgm0, tid0, tim0, tinc0, tmarg0, tms0, tpd0, tpm0, txs0, 
        u0, uelas0, uepriv0, ug0, up0, us0, 
        vcif, vdep, vdfb, vdfp, vdgb, vdgp, vdib, vdip, vdpb, vdpp, vfob, vkb, vmfb, vmfp, vmgb, vmgp, vmib,
        vmip, vmpb, vmpp, vmsb, vst, vtwr, vxsb, 
        xg0, 
        y0, yg0, yi0, yp0,
        zshr0,
    )

    GTAPMod = JuMP.Model(PATHSolver.Optimizer)

    @variables GTAPMod begin
        qint[a = acts, r = reg], (start = 1)                                    # Aggregate intermediate demand
        qva[a = acts, r =  reg], (start = 1)                                    # Demand for aggregate value added bundle
        po[a = acts, r = reg], (start = 1)                                      # Unit cost of production
        qfa[i = comm, a = acts, r = reg], (start = 1)                           # (Armington) demand for inputs
        pint[a = acts, r = reg], (start = 1)                                    # Aggregate cost of inputs
        qfe[f = endw, a = acts, r = reg], (start = 1)                           # Factor demand
        pva[a = acts, r = reg], (start = 1)                                     # Aggregate price of value added
        qca[i = comm, a = acts, r = reg], (start = 1)                           # Commodity i supplied by a
        qo[a = acts, r = reg], (start = 1)                                      # Output by activity a
        pca[i = comm, a = acts, r = reg], (start = 1)                           # Tax-inclusive price of commodity i supplied by a
        pds[i = comm, r = reg], (start = 1)                                     # Price of supply of commodity i
        ps[i = comm, a = acts, r = reg], (start = 1)                            # Tax-exclusive price of commodity i supplied by a
        taxrout[r = reg], (start = 1)                                           # Value of output tax revenues
        taxrfu[r = reg], (start = 1)                                            # Value of factor-use tax revenues
        taxriu[r = reg], (start = 1)                                            # Value of tax revenues on intermediate demand
        taxrpc[r = reg], (start = 1)                                            # Value of tax revenues on private demand
        taxrgc[r = reg], (start = 1)                                            # Value of tax revenues on government demand
        taxric[r = reg], (start = 1)                                            # Value of tax revenues on investment demand
        taxrexp[r = reg], (start = 1)                                           # Value of tax revenues on exports
        taxrimp[r = reg], (start = 1)                                           # Value of tax revenues on imports
        indtax[r = reg], (start = 1)                                            # Total value of indirect tax revenues
        fincome[r = reg], (start = 1)                                           # Factor income at basic prices net of depreciation
        y[r = reg], (start = 1)                                                 # Regional income
        uepriv[r = reg], (start = uepriv0[r])                                   # Elasticity of cost wrt utility from private consumption
        uelas[r = reg], (start = uelas0[r])                                     # Elasticity of cost of utility wrt utility
        yp[r = reg], (start = 1)                                                # Private consumption expenditure
        yg[r = reg], (start = 1)                                                # Government consumption expenditure
        save[r = reg], (start = 1)                                              # Regional supply of nominal savings
        u[r = reg], (start = u0[r])                                             # Total utility (per capita)
        us[r = reg], (start = us0[r])                                           # Utility derived from savings
        zshr[i = comm, r = reg], (start = 1)                                    # Auxiliary share variable for private consumption
        conshr[i = comm, r = reg], (start = 1)                                  # Private budget shares
        qpa[i = comm, r = reg], (start = 1)                                     # Private consumption (at the Armington level)
        up[r = reg], (start = 1)                                                # Private utility
        qga[i = comm, r = reg], (start = 1)                                     # Government consumption (at the Armington level)
        pgov[r = reg], (start = 1)                                              # Government expenditure price deflator
        xg[r = reg], (start = 1)                                                # Aggregate volume of government expenditures
        ug[r = reg], (start = 1)                                                # Per capita utility derived from public expenditures
        qia[i = comm, r = reg], (start = 1)                                     # Investment expenditure (at the Armington level)
        pinv[r = reg], (start = 1)                                              # Investment expenditure price deflator
        qinv[r = reg], (start = 1)                                              # Aggregate volume of investment expenditures
        pfd[i = comm, a = acts, r = reg], (start = 1)                           # End user price for domestic intermediate goods
        pfm[i = comm, a = acts, r = reg], (start = 1)                           # End user price for imported intermediate goods
        ppd[i = comm, r = reg], (start = 1)                                     # End user price for domestic private goods
        ppm[i = comm, r = reg], (start = 1)                                     # End user price for imported private goods
        pgd[i = comm, r = reg], (start = 1)                                     # End user price for domestic government goods
        pgm[i = comm, r = reg], (start = 1)                                     # End user price for imported government goods
        pid[i = comm, r = reg], (start = 1)                                     # End user price for domestic investment goods
        pim[i = comm, r = reg], (start = 1)                                     # End user price for imported investment goods
        pfa[i = comm, a = acts, r = reg], (start = 1)                           # Armington price for intermediate goods
        ppa[i = comm, r = reg], (start = 1)                                     # Armington price for private goods
        pga[i = comm, r = reg], (start = 1)                                     # Armington price for government goods
        pia[i = comm, r = reg], (start = 1)                                     # Armington price for investment goods
        qfd[i = comm, a = acts, r = reg], (start = 1)                           # Demand for domestic intermediate goods
        qfm[i = comm, a = acts, r = reg], (start = 1)                           # Demand for imported intermediate goods
        qpd[i = comm, r = reg], (start = 1)                                     # Demand for domestic private goods
        qpm[i = comm, r = reg], (start = 1)                                     # Demand for imported private goods
        qgd[i = comm, r = reg], (start = 1)                                     # Demand for domestic government goods
        qgm[i = comm, r = reg], (start = 1)                                     # Demand for imported government goods
        qid[i = comm, r = reg], (start = 1)                                     # Demand for domestic investment goods
        qim[i = comm, r = reg], (start = 1)                                     # Demand for imported investment goods
        qms[i = comm, r = reg], (start = 1)                                     # Aggregate import demand
        qxs[i = comm, s = reg, d = reg], (start = 1)                            # Export/import from region s towards region d
        pms[i = comm, d = reg], (start = 1)                                     # Aggregate import price
        pfob[i = comm, s = reg, d = reg], (start = 1)                           # FOB export price
        pcif[i = comm, s = reg, d = reg], (start = 1)                           # CIF import price
        pmds[i = comm, s = reg, d = reg], (start = 1)                           # Tariff inclusive bilateral import price
        qds[i = comm, r = reg], (start = 1)                                     # Aggregate domestic demand for domestic goods (/x qst)
        qc[i = comm, r = reg], (start = 1)                                      # Market equilibrium for domestically produced goods
        pes[f = endw, a = acts, r = reg], (start = 1)                           # Factor market allocation and equilibrium
        pe[f = endw, r = reg], (start = 1)                                      # Aggregate factor price
        pfe[f = endw, a = acts, r = reg], (start = 1)                           # Factor returns at purchasers prices
        peb[f = endw, a = acts, r = reg], (start = 1)                           # Factor returns at basic prices
        qe[f = endw, r = reg], (start = 1)                                      # Total factor stock
        qtmfsd[m = comm, i = comm, s = reg, d = reg], (start = 1)               # Demand for service m to transport i from s to d
        ptrans[i = comm, s = reg, d = reg], (start = 1)                         # Average price to transport i from s to d
        qtm[m = comm], (start = 1)                                              # Total demand for service m
        qst[m = comm, r = reg], (start = 1)                                     # Regional supply of service m
        pt[m = comm], (start = 1)                                               # Aggregate world price of service m
        gdpfc[r = reg], (start = 1)                                             # Nominal GDP at factor cost
        kb[r = reg], (start = 1)                                                # Non-normalized aggregate capital stock
        ke[r = reg], (start = 1)                                                # End-of-period capital stock
        rental[r = reg], (start = 1)                                            # After-tax gross rate of return to capital
        rorc[r = reg], (start = 1)                                              # Net rate of return to capital
        rore[r = reg], (start = 1)                                              # Expected rate of return to capital
        savf[r = reg], (start = savf0[r])                                       # Default capital account closure--determination of capital flows
        rorg, (start = 1)                                                       # Default capital account closure--determination of the global rate of return
        yi[r = reg], (start = 1)                                                # Nominal gross investment
        chisave, (start = 1)                                                    # Savings adjustment factor
        psave[r = reg], (start = 1)                                             # Price index for regional savings
        pfact[r = reg], (start = 1)                                             # Average regional factor price
        pfactw, (start = 1)                                                     # Average global factor price
        Walras, (start = 0)                                                     # Check on Walras' Law

        #   (Normally) exogenous variables
        ptax[i = comm, a = acts, r = reg], (start = ptax0[i,a,r])               # Output tax
        tfe[f = endw, a = acts, r = reg], (start = tfe0[f,a,r])                 # Tax rate on factor use
        tinc[f = endw, a = acts, r = reg], (start = tinc0[f,a,r])               # Tax on factor returns
        tfd[i = comm, a = acts, r = reg], (start = tfd0[i,a,r])                 # Tax on intermediate use of domestic goods
        tfm[i = comm, a = acts, r = reg], (start = tfm0[i,a,r])                 # Tax on intermediate use of imported goods
        tpd[i = comm, r = reg], (start = tpd0[i,r])                             # Tax on private use of domestic goods
        tpm[i = comm, r = reg], (start = tpm0[i,r])                             # Tax on private use of imported goods
        tgd[i = comm, r = reg], (start = tgd0[i,r])                             # Tax on public use of domestic goods
        tgm[i = comm, r = reg], (start = tgm0[i,r])                             # Tax on public use of imported goods
        tid[i = comm, r = reg], (start = tid0[i,r])                             # Tax on investment use of domestic goods
        tim[i = comm, r = reg], (start = tim0[i,r])                             # Tax on investment use of imported goods
        txs[i = comm, s = reg, d = reg], (start = txs0[i,s,d])                  # Tax on exports by source and destination
        tmarg[i = comm, s = reg, d = reg], (start = tmarg0[i,s,d])              # Total transport margin by source and destination
        tms[i = comm, s = reg, d = reg], (start = tms0[i,s,d])                  # Tax on imports by source and destination
        pop[r = reg], (start = pop0[r])                                         # Population
        chi_qe[f = endw, r = reg], (start = 1)                                  # Shifter on aggregate factor stock
        pnum, (start = 1)                                                       # Model numeraire--fixed outside of model
    end
    #   Production module --------------------------------------------------------------------------

    #   Intermediate demand bundle
    @constraint(GTAPMod, eq_qint[a in acts, r in reg; qint0[a,r] != 0],
        qint[a,r] - qo[a,r]*(po[a,r]/pint[a,r])^ESUBT[a,r] ⟂ qint[a,r])
    for a ∈ acts, r ∈ reg
        if qint0[a,r] == 0
            fix(qint[a,r], 0, force=true)
        end
    end

    #   Value added bundle
    @constraint(GTAPMod, eq_qva[a in acts, r in reg; qva0[a,r] != 0],
        qva[a,r] - qo[a,r]*(po[a,r]/pva[a,r])^ESUBT[a,r] ⟂ qva[a,r])
    for a ∈ acts, r ∈ reg
        if qva0[a,r] == 0
            fix(qva[a,r], 0, force=true)
        end
    end
    
    #   Unit cost of production
    @constraint(GTAPMod, eq_po[a in acts, r in reg; qo0[a,r] != 0],
        po[a,r]^(1 - ESUBT[a,r]) - (shr_qint[a,r]*pint[a,r]^(1 - ESUBT[a,r]) + shr_qva[a,r]*pva[a,r]^(1 - ESUBT[a,r])) ⟂ po[a,r])
    for a ∈ acts, r ∈ reg
        if qo0[a,r] == 0
            fix(po[a,r], 1, force=true)
        end
    end

    #   Armington intermediate demand
    @constraint(GTAPMod, eq_qfa[i in comm, a in acts, r in reg; qfa0[i,a,r] != 0],
        qfa[i,a,r] - qint[a,r]*(pint[a,r]/pfa[i,a,r])^ESUBC[a,r] ⟂ qfa[i,a,r])
    for i ∈ comm, a ∈ acts, r ∈ reg
        if qfa0[i,a,r] == 0
            fix(qfa[i,a,r], 0, force=true)
        end
    end

     #   Aggregate price of inputs
     @constraint(GTAPMod, eq_pint[a in acts, r in reg; qint0[a,r] != 0],
        pint[a,r]^(1 - ESUBC[a,r]) - (sum(shr_qfa[i,a,r]*pfa[i,a,r]^(1 - ESUBC[a,r]) for i in comm)) ⟂ pint[a,r])
    for a ∈ acts, r ∈ reg
        if qint0[a,r] == 0
            fix(pint[a,r], 1, force=true)
        end
    end

    #   Factor demand
    @constraint(GTAPMod, eq_qfe[f in endw, a in acts, r in reg; qfe0[f,a,r] != 0],
        qfe[f,a,r] - qva[a,r]*(pva[a,r]/pfe[f,a,r])^ESUBVA[a,r] ⟂ qfe[f,a,r])
    for f ∈ endw, a ∈ acts, r ∈ reg
        if qfe0[f,a,r] == 0
            fix(qfe[f,a,r], 0, force=true)
        end
    end

     #   Aggregate price of value added
    @constraint(GTAPMod, eq_pva[a in acts, r in reg; qva0[a,r] != 0],
        pva[a,r]^(1 - ESUBVA[a,r]) - (sum(shr_qfe[f,a,r]*pfe[f,a,r]^(1 - ESUBVA[a,r]) for f in endw)) ⟂ pva[a,r])
    for a ∈ acts, r ∈ reg
        if qva0[a,r] == 0
            fix(pva[a,r], 1, force=true)
        end
    end

     #   Make module --------------------------------------------------------------------------------

    #   Supply of commodity i by activity a
    @constraint(GTAPMod, eq_qca[i in comm, a in acts, r in reg; qca0[i,a,r] != 0],
        if ETRAQ[a,r] != Inf
            qca[i,a,r] - qo[a,r]*(ps[i,a,r]/po[a,r])^ETRAQ[a,r]
        else
            ps[i,a,r] - po[a,r]
        end
        ⟂ qca[i,a,r])
    for i ∈ comm, a ∈ acts, r ∈ reg
        if qca0[i,a,r] == 0
            fix(qca[i,a,r], 0, force=true)
        end
    end

    #   Defines output by activity a
    @constraint(GTAPMod, eq_qo[a in acts, r in reg; qo0[a,r] != 0],
        if ETRAQ[a,r] != Inf
            po[a,r]^(1 + ETRAQ[a,r]) - (sum(shr_qcas[i, a, r]*ps[i,a,r]^(1 + ETRAQ[a,r]) for i in comm)) 
        else
            qo[a,r] - sum(shr_qcas[i, a, r]*qca[i,a,r] for i in comm)
        end
        ⟂ qo[a,r])
    for a ∈ acts, r ∈ reg
        if qo0[a,r] == 0
            fix(qo[a,r], 0, force=true)
        end
    end

    #  Demand for commodity i produced by activity a (inverse CES)
    @constraint(GTAPMod, eq_pca[i in comm, a in acts, r in reg; qca0[i,a,r] != 0],
        if ESUBQ[i,r] != Inf
            qca[i,a,r] - qc[i,r]*(pds[i,r]/pca[i,a,r])^ESUBQ[i,r] 
        else
            pca[i,a,r] - pds[i,r]
        end    
        ⟂ pca[i,a,r])
    for i ∈ comm, a ∈ acts, r ∈ reg
        if qca0[i,a,r] == 0
            fix(pca[i,a,r], 1, force=true)
        end
    end

    #   Defines supply price of commodity i
    @constraint(GTAPMod, eq_pds[i in comm, r in reg; qc0[i,r] != 0],
        if ESUBQ[i,r] != Inf
            pds[i,r]^(1 - ESUBQ[i,r]) - (sum(shr_qcad[i,a,r]*pca[i,a,r]^(1 - ESUBQ[i,r]) for a in acts))
        else
            qc[i,r] - sum(shr_qcad[i,a,r]*qca[i,a,r] for a in acts)    
        end
        ⟂ pds[i,r])
    for i ∈ comm, r ∈ reg
        if qc0[i,r] == 0
            fix(pds[i,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_ps[i in comm, a in acts, r in reg; qca0[i,a,r] != 0],
        pca[i,a,r] - ps[i,a,r]*(1 + ptax[i,a,r])*(ps0[i,a,r]/pca0[i,a,r]) ⟂ ps[i,a,r])
    for i ∈ comm, a ∈ acts, r ∈ reg
        if qca0[i,a,r] == 0
            fix(ps[i,a,r], 1, force=true)
        end
    end

    #   Income and its distribution ----------------------------------------------------------------

    #   Tax revenues from intermediate demand
    @constraint(GTAPMod, eq_taxriu[r in reg],
        taxriu[r]*taxriu0[r] - (sum(tfd[i,a,r]*pds[i,r]*qfd[i,a,r]*pds0[i,r]*qfd0[i,a,r]
                             +      tfm[i,a,r]*pms[i,r]*qfm[i,a,r]*pms0[i,r]*qfm0[i,a,r] for i in comm, a in acts)) ⟂ taxriu[r])

    #   Tax revenues from private demand
    @constraint(GTAPMod, eq_taxrpc[r in reg],
        taxrpc[r]*taxrpc0[r] - (sum(tpd[i,r]*pds[i,r]*qpd[i,r]*pds0[i,r]*qpd0[i,r]
                             +      tpm[i,r]*pms[i,r]*qpm[i,r]*pms0[i,r]*qpm0[i,r] for i in comm)) ⟂ taxrpc[r])

    #   Tax revenues from government demand
    @constraint(GTAPMod, eq_taxrgc[r in reg],
        taxrgc[r]*taxrgc0[r] - (sum(tgd[i,r]*pds[i,r]*qgd[i,r]*pds0[i,r]*qgd0[i,r]
                             +      tgm[i,r]*pms[i,r]*qgm[i,r]*pms0[i,r]*qgm0[i,r] for i in comm)) ⟂ taxrgc[r])

    #   Tax revenues from investment demand
    @constraint(GTAPMod, eq_taxric[r in reg],
        taxric[r]*taxric0[r] - (sum(tid[i,r]*pds[i,r]*qid[i,r]*pds0[i,r]*qid0[i,r]
                             +      tim[i,r]*pms[i,r]*qim[i,r]*pms0[i,r]*qim0[i,r] for i in comm)) ⟂ taxric[r])

    #   Output tax revenues
    @constraint(GTAPMod, eq_taxrout[r in reg],
       taxrout[r]*taxrout0[r] - (sum((ps0[i,a,r]*qca0[i,a,r])*ptax[i,a,r]*ps[i,a,r]*qca[i,a,r] for i in comm, a in acts)) ⟂ taxrout[r])

    #   Factor use tax revenues
    @constraint(GTAPMod, eq_taxrfu[r in reg],
       taxrfu[r]*taxrfu0[r] - (sum(tfe[f,a,r]*peb[f,a,r]*qfe[f,a,r]*peb0[f,a,r]*qfe0[f,a,r] for f in endw, a in acts)) ⟂ taxrfu[r])

    #   Export tax revenues
    @constraint(GTAPMod, eq_taxrexp[r in reg],
        taxrexp[r]*taxrexp0[r] - (sum(txs[i,r,d]*pds[i,r]*qxs[i,r,d]*pds0[i,r]*qxs0[i,r,d] for i in comm, d in reg)) ⟂ taxrexp[r])

    #   Import tax revenues
    @constraint(GTAPMod, eq_taxrimp[r in reg],
        taxrimp[r]*taxrimp0[r] - (sum(tms[i,s,r]*pcif[i,s,r]*qxs[i,s,r]*pcif0[i,s,r]*qxs0[i,s,r] for i in comm, s in reg)) ⟂ taxrimp[r])

    #   Total indirect taxes
    @constraint(GTAPMod, eq_indtax[r in reg],
        indtax[r] - (taxrpc[r]*(taxrpc0[r]/indtax0[r]) + taxrgc[r]*(taxrgc0[r]/indtax0[r]) + taxric[r]*(taxric0[r]/indtax0[r]) + taxriu[r]*(taxriu0[r]/indtax0[r])
            + taxrfu[r]*(taxrfu0[r]/indtax0[r]) + taxrout[r]*(taxrout0[r]/indtax0[r]) + taxrexp[r]*(taxrexp0[r]/indtax0[r]) + taxrimp[r]*(taxrimp0[r]/indtax0[r])) ⟂ indtax[r])

    @constraint(GTAPMod, eq_fincome[r in reg],
       fincome[r] - (sum(peb[f,a,r]*qfe[f,a,r]*(peb0[f,a,r]*qfe0[f,a,r]/fincome0[r]) for f in endw, a in acts) - depr[r]*pinv[r]*kb[r]*(pinv0[r]*kb0[r]/fincome0[r])) ⟂ fincome[r])

    @constraint(GTAPMod, eq_y[r in reg],
       y[r] - (fincome[r]*(fincome0[r]/y0[r]) + indtax[r]*(indtax0[r]/y0[r])) ⟂ y[r])

    #   Disposition of income and final demand -----------------------------------------------------

    #   Elasticity of expenditure wrt utility from private consumption
    @constraint(GTAPMod, eq_uepriv[r in reg],
        uepriv[r] - sum(conshr[i,r]*conshr0[i,r]*INCPAR[i,r] for i in comm) ⟂ uepriv[r])

    #   Elasticity of total expenditure wrt to utility
    @constraint(GTAPMod, eq_uelas[r in reg],
        uelas[r]*(dppriv[r]/uepriv[r] + dpgov[r] + dpsave[r]) - 1 ⟂ uelas[r])

    #   Aggregate nominal private consumption
    @constraint(GTAPMod, eq_yp[r in reg],
        yp[r] - dppriv[r]*(uelas[r]/uepriv[r])*y[r]*(y0[r]/yp0[r]) ⟂ yp[r])

    #   Aggregate nominal government consumption
    @constraint(GTAPMod, eq_yg[r in reg],
        yg[r] - dpgov[r]*uelas[r]*y[r]*(y0[r]/yg0[r]) ⟂ yg[r])

    #   Domestic supply of savings
    @constraint(GTAPMod, eq_qsave[r in reg],
        save[r]*(save0[r]/y0[r]) - dpsave[r]*uelas[r]*y[r] ⟂ save[r])

    #   Utility function
    @constraint(GTAPMod, eq_u[r in reg],
        log(u[r]) - (log(au[r]) +  dppriv[r]*log(up[r]) + dpgov[r]*log(ug[r]) + dpsave[r]*log(us[r])) ⟂ u[r])

    #   Utility from savings
    @constraint(GTAPMod, eq_us[r in reg],
        us[r] - (save0[r]*aus[r]/pop[r])*(save[r]/psave[r]) ⟂ us[r])

    #   Auxiliary consumption variable
    @constraint(GTAPMod, eq_zshr[i in comm, r in reg; qpa0[i,r] != 0],
        zshr[i,r] - (alphac[i,r]*SUBPAR[i,r]*(((ppa0[i,r]/yp0[r])^SUBPAR[i,r])/zshr0[i,r]))*(up[r]^(INCPAR[i,r]*SUBPAR[i,r]))*((pop[r]*ppa[i,r]/yp[r])^SUBPAR[i,r]) ⟂ zshr[i,r])
    for i ∈ comm, r ∈ reg
        if qpa0[i,r] == 0
            fix(zshr[i,r], 0, force=true)
        end
    end

    #   Private consumption budget shares
    @constraint(GTAPMod, eq_conshr[i in comm, r in reg; qpa0[i,r] != 0],
        conshr[i,r] - (zshr[i,r]*(zshr0[i,r]/conshr0[i,r]))/(sum(zshr[j,r]*zshr0[j,r] for j in comm)) ⟂ conshr[i,r])
    for i ∈ comm, r ∈ reg
        if qpa0[i,r] == 0
            fix(conshr[i,r], 0, force=true)
        end
    end

    #   Household demand for goods and services
    @constraint(GTAPMod, eq_qpa[i in comm, r in reg; qpa0[i,r] != 0],
        ppa[i,r]*qpa[i,r] - conshr[i,r]*yp[r]*((conshr0[i,r]*yp0[r])/(ppa0[i,r]*qpa0[i,r])) ⟂ qpa[i,r])
    for i ∈ comm, r ∈ reg
        if qpa0[i,r] == 0
            fix(qpa[i,r], 0, force=true)
        end
    end

    #   Private utility
    @constraint(GTAPMod, eq_up[r in reg],
        1 - sum(ifelse(zshr0[i,r] != 0.0, (zshr0[i,r]/SUBPAR[i,r])*zshr[i,r], 0.0) for i in comm) ⟂ up[r])

    #   Government demand for goods and services
    @constraint(GTAPMod, eq_qga[i in comm, r in reg; qga0[i,r] != 0],
        qga[i,r] - (xg[r] * (pgov[r] / pga[i,r])^ESUBG[r]) ⟂ qga[i,r])
    for i ∈ comm, r ∈ reg
        if qga0[i,r] == 0
            fix(qga[i,r], 0, force=true)
        end
    end

    #   Government expenditure price deflator
    @constraint(GTAPMod, eq_pgov[r in reg],
        pgov[r]*xg[r] - (sum(shr_gov[i,r]*pga[i,r]*qga[i,r] for i in comm)) ⟂ pgov[r])

    #   Aggregate volume of government expenditures
    @constraint(GTAPMod, eq_xg[r in reg],
        pgov[r]*xg[r] - yg[r]*(yg0[r]/(pgov0[r]*xg0[r])) ⟂ xg[r])

    @constraint(GTAPMod, eq_ug[r in reg],
        ug[r] - (aug[r]*xg0[r]/pop[r])*xg[r] ⟂ ug[r])

    #   Investment demand for goods and services
    @constraint(GTAPMod, eq_qia[i in comm, r in reg; qia0[i,r] != 0],
        qia[i,r] - (qinv[r] * (pinv[r] / pia[i,r])^ESUBI[r]) ⟂ qia[i,r])
    for i ∈ comm, r ∈ reg
        if qia0[i,r] == 0
            fix(qia[i,r], 0, force=true)
        end
    end

    #   Investment expenditure price deflator
    @constraint(GTAPMod, eq_pinv[r in reg],
        pinv[r]*qinv[r] - (sum(shr_inv[i,r]*pia[i,r]*qia[i,r] for i in comm)) ⟂ pinv[r])

    #   Aggregate volume of investment expenditures
    @constraint(GTAPMod, eq_qinv[r in reg],
        pinv[r]*qinv[r] - yi[r]*(yi0[r]/(pinv0[r]*qinv0[r])) ⟂ qinv[r])

    #   Armington module ---------------------------------------------------------------------------

    #   End use prices for domestic and imported goods
    #   Assumes all users face same basic prices: pds for domestic goods and pms for imported goods

    @constraint(GTAPMod, eq_pfd[i in comm, a in acts, r in reg; qfd0[i,a,r] != 0],
        pfd[i,a,r] - (1 + tfd[i,a,r]) * pds[i,r]*(pds0[i,r]/pfd0[i,a,r]) ⟂ pfd[i,a,r])
    for i ∈ comm, a ∈ acts, r ∈ reg
        if qfd0[i,a,r] == 0
            fix(pfd[i,a,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_pfm[i in comm, a in acts, r in reg; qfm0[i,a,r] != 0],
        pfm[i,a,r] - (1 + tfm[i,a,r]) * pms[i,r]*(pms0[i,r]/pfm0[i,a,r]) ⟂ pfm[i,a,r])
    for i ∈ comm, a ∈ acts, r ∈ reg
        if qfm0[i,a,r] == 0
            fix(pfm[i,a,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_ppd[i in comm, r in reg; qpd0[i,r] != 0],
        ppd[i,r] - (1 + tpd[i,r]) * pds[i,r]*(pds0[i,r]/ppd0[i,r]) ⟂ ppd[i,r])
    for i ∈ comm, r ∈ reg
        if qpd0[i,r] == 0
            fix(ppd[i,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_ppm[i in comm, r in reg; qpm0[i,r] != 0],
        ppm[i,r] - (1 + tpm[i,r]) * pms[i,r]*(pms0[i,r]/ppm0[i,r]) ⟂ ppm[i,r])
    for i ∈ comm, r ∈ reg
        if qpm0[i,r] == 0
            fix(ppm[i,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_pgd[i in comm, r in reg; qgd0[i,r] != 0],
        pgd[i,r] - (1 + tgd[i,r]) * pds[i,r]*(pds0[i,r]/pgd0[i,r]) ⟂ pgd[i,r])
    for i ∈ comm, r ∈ reg
        if qgd0[i,r] == 0
            fix(pgd[i,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_pgm[i in comm, r in reg; qgm0[i,r] != 0],
        pgm[i,r] - (1 + tgm[i,r]) * pms[i,r]*(pms0[i,r]/pgm0[i,r]) ⟂ pgm[i,r])
    for i ∈ comm, r ∈ reg
        if qgm0[i,r] == 0
            fix(pgm[i,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_pid[i in comm, r in reg; qid0[i,r] != 0],
        pid[i,r] - (1 + tid[i,r]) * pds[i,r]*(pds0[i,r]/pid0[i,r]) ⟂ pid[i,r])
    for i ∈ comm, r ∈ reg
        if qid0[i,r] == 0
            fix(pid[i,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_pim[i in comm, r in reg; qim0[i,r] != 0],
        pim[i,r] - (1 + tim[i,r]) * pms[i,r]*(pms0[i,r]/pim0[i,r]) ⟂ pim[i,r])
    for i ∈ comm, r ∈ reg
        if qim0[i,r] == 0
            fix(pim[i,r], 1, force=true)
        end
    end

    #   Armington price

    @constraint(GTAPMod, eq_pfa[i in comm, a in acts, r in reg; qfa0[i,a,r] != 0],
        pfa[i,a,r]^(1 - ESUBD[i,r]) - (shr_qfd[i,a,r]*pfd[i,a,r]^(1 - ESUBD[i,r]) + shr_qfm[i,a,r]*pfm[i,a,r]^(1 - ESUBD[i,r])) ⟂ pfa[i,a,r])
    for i ∈ comm, a ∈ acts, r ∈ reg
        if qfa0[i,a,r] == 0
            fix(pfa[i,a,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_ppa[i in comm, r in reg; qpa0[i,r] != 0],
        ppa[i,r]^(1 - ESUBD[i,r]) - (shr_qpd[i,r]*ppd[i,r]^(1 - ESUBD[i,r]) + shr_qpm[i,r]*ppm[i,r]^(1 - ESUBD[i,r])) ⟂ ppa[i,r])
    for i ∈ comm, r ∈ reg
        if qpa0[i,r] == 0
            fix(ppa[i,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_pga[i in comm, r in reg; qga0[i,r] != 0],
        pga[i,r]^(1 - ESUBD[i,r]) - (shr_qgd[i,r]*pgd[i,r]^(1 - ESUBD[i,r]) + shr_qgm[i,r]*pgm[i,r]^(1 - ESUBD[i,r])) ⟂ pga[i,r])
    for i ∈ comm, r ∈ reg
        if qga0[i,r] == 0
            fix(pga[i,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_pia[i in comm, r in reg; qia0[i,r] != 0],
        pia[i,r]^(1 - ESUBD[i,r]) - (shr_qid[i,r]*pid[i,r]^(1 - ESUBD[i,r]) + shr_qim[i,r]*pim[i,r]^(1 - ESUBD[i,r])) ⟂ pia[i,r])
    for i ∈ comm, r ∈ reg
        if qia0[i,r] == 0
            fix(pia[i,r], 1, force=true)
        end
    end

    #   Top level sourcing
    @constraint(GTAPMod, eq_qfd[i in comm, a in acts, r in reg; qfd0[i,a,r] != 0],
        qfd[i,a,r] - (qfa[i,a,r]  * (pfa[i,a,r] / pfd[i,a,r])^ESUBD[i,r]) ⟂ qfd[i,a,r])
    for i ∈ comm, a ∈ acts, r ∈ reg
        if qfd0[i,a,r] == 0
            fix(qfd[i,a,r], 0, force=true)
        end
    end

    @constraint(GTAPMod, eq_qfm[i in comm, a in acts, r in reg; qfm0[i,a,r] != 0],
        qfm[i,a,r] - (qfa[i,a,r]  * (pfa[i,a,r] / pfm[i,a,r])^ESUBD[i,r]) ⟂ qfm[i,a,r])
    for i ∈ comm, a ∈ acts, r ∈ reg
        if qfm0[i,a,r] == 0
            fix(qfm[i,a,r], 0, force=true)
        end
    end

    @constraint(GTAPMod, eq_qpd[i in comm, r in reg; qpd0[i,r] != 0],
        qpd[i,r] - (qpa[i,r]  * (ppa[i,r] / ppd[i,r])^ESUBD[i,r]) ⟂ qpd[i,r])
    for i ∈ comm, r ∈ reg
        if qpd0[i,r] == 0
            fix(qpd[i,r], 0, force=true)
        end
    end

    @constraint(GTAPMod, eq_qpm[i in comm, r in reg; qpm0[i,r] != 0],
        qpm[i,r] - (qpa[i,r]  * (ppa[i,r] / ppm[i,r])^ESUBD[i,r]) ⟂ qpm[i,r])
    for i ∈ comm, r ∈ reg
        if qpm0[i,r] == 0
            fix(qpm[i,r], 0, force=true)
        end
    end

    @constraint(GTAPMod, eq_qgd[i in comm, r in reg; qgd0[i,r] != 0],
        qgd[i,r] - (qga[i,r]  * (pga[i,r] / pgd[i,r])^ESUBD[i,r]) ⟂ qgd[i,r])
    for i ∈ comm, r ∈ reg
        if qgd0[i,r] == 0
            fix(qgd[i,r], 0, force=true)
        end
    end

    @constraint(GTAPMod, eq_qgm[i in comm, r in reg; qgm0[i,r] != 0],
        qgm[i,r] - (qga[i,r]  * (pga[i,r] / pgm[i,r])^ESUBD[i,r]) ⟂ qgm[i,r])
    for i ∈ comm, r ∈ reg
        if qgm0[i,r] == 0
            fix(qgm[i,r], 0, force=true)
        end
    end

    @constraint(GTAPMod, eq_qid[i in comm, r in reg; qid0[i,r] != 0],
        qid[i,r] - (qia[i,r]  * (pia[i,r] / pid[i,r])^ESUBD[i,r]) ⟂ qid[i,r])
    for i ∈ comm, r ∈ reg
        if qid0[i,r] == 0
            fix(qid[i,r], 0, force=true)
        end
    end

    @constraint(GTAPMod, eq_qim[i in comm, r in reg; qim0[i,r] != 0],
        qim[i,r] - (qia[i,r]  * (pia[i,r] / pim[i,r])^ESUBD[i,r]) ⟂ qim[i,r])
    for i ∈ comm, r ∈ reg
        if qim0[i,r] == 0
            fix(qim[i,r], 0, force=true)
        end
    end

    #   Second level Armington -- sourcing by region
    #   Aggregate import demand
    @constraint(GTAPMod, eq_qms[i in comm, r in reg; qms0[i,r] != 0],
        qms[i,r] - (sum(qfm[i,a,r]*(qfm0[i,a,r]/qms0[i,r]) for a in acts) + qpm[i,r]*(qpm0[i,r]/qms0[i,r]) 
                        + qgm[i,r]*(qgm0[i,r]/qms0[i,r]) + qim[i,r]*(qim0[i,r]/qms0[i,r])) ⟂ qms[i,r])
    for i ∈ comm, r ∈ reg
        if qms0[i,r] == 0
            fix(qms[i,r], 0, force=true)
        end
    end

    #   Demand for imports by region d sourced from region s
    @constraint(GTAPMod, eq_qxs[i in comm, s in reg, d in reg; qxs0[i,s,d] != 0],
        qxs[i,s,d] - (qms[i,d] * (pms[i,d] / pmds[i,s,d])^ESUBM[i,d]) ⟂ qxs[i,s,d])
    for i ∈ comm, s ∈ reg, d ∈ reg
        if qxs0[i,s,d] == 0
            fix(qxs[i,s,d], 0, force=true)
        end
    end

    #   Aggregate price of imports
    @constraint(GTAPMod, eq_pms[i in comm, d in reg; qms0[i,d] != 0],
        pms[i,d]^(1 - ESUBM[i,d]) - (sum(shr_qxs[i,s,d]*pmds[i,s,d]^(1 - ESUBM[i,d]) for s in reg)) ⟂ pms[i,d])
    for i ∈ comm, d ∈ reg
        if qms0[i,d] == 0
            fix(pms[i,d], 1, force=true)
        end
    end

    #   Bilateral prices ---------------------------------------------------------------------------

    #   Export price at FOB
    @constraint(GTAPMod, eq_pfob[i in comm, s in reg, d in reg; qxs0[i,s,d] != 0],
        pfob[i,s,d] - ((1 + txs[i,s,d])*pds[i,s]*(pds0[i,s]/pfob0[i,s,d])) ⟂ pfob[i,s,d])
    for i ∈ comm, s ∈ reg, d ∈ reg
        if qxs0[i,s,d] == 0
            fix(pfob[i,s,d], 1, force=true)
        end
    end

    #   Import price at CIF
    @constraint(GTAPMod, eq_pcif[i in comm, s in reg, d in reg; qxs0[i,s,d] != 0],
        pcif[i,s,d] - (pfob[i,s,d]*(pfob0[i,s,d]/pcif0[i,s,d]) + ptrans[i,s,d]*tmarg[i,s,d]*(ptrans0[i,s,d]/pcif0[i,s,d])) ⟂ pcif[i,s,d])
    for i ∈ comm, s ∈ reg, d ∈ reg
        if qxs0[i,s,d] == 0
            fix(pcif[i,s,d], 1, force=true)
        end
    end

    #   Import price tariff inclusive
    @constraint(GTAPMod, eq_pmds[i in comm, s in reg, d in reg; qxs0[i,s,d] != 0],
        pmds[i,s,d] - ((1 + tms[i,s,d])*pcif[i,s,d]*(pcif0[i,s,d]/pmds0[i,s,d])) ⟂ pmds[i,s,d])
    for i ∈ comm, s ∈ reg, d ∈ reg
        if qxs0[i,s,d] == 0
            fix(pmds[i,s,d], 1, force=true)
        end
    end

    #   Global transport services ------------------------------------------------------------------

    @constraint(GTAPMod, eq_qtmfsd[m in comm, i in comm, s in reg, d in reg; qtmfsd0[m,i,s,d] !=0],
        qtmfsd[m,i,s,d] - tmarg[i,s,d]*qxs[i,s,d]*(qxs0[i,s,d]/qtmfsd0[m,i,s,d]) ⟂ qtmfsd[m,i,s,d])
    for m ∈ comm, i ∈ comm, s ∈ reg, d ∈ reg
        if qtmfsd0[m,i,s,d] == 0
            fix(qtmfsd[m,i,s,d], 0, force=true)
        end
    end

    @constraint(GTAPMod, eq_ptrans[i in comm, s in reg, d in reg; tmarg0[i,s,d] != 0],
        ptrans[i,s,d] - sum(shr_vtwr[m,i,s,d]*pt[m] for m in comm) ⟂ ptrans[i,s,d])
    for i ∈ comm, s ∈ reg, d ∈ reg
        if tmarg0[i,s,d] == 0
            fix(ptrans[i,s,d], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_qtm[m in comm; qtm0[m] != 0],
        qtm[m] - sum(qtmfsd[m,i,s,d]*(qtmfsd0[m,i,s,d]/qtm0[m]) for i in comm, s in reg, d in reg) ⟂ qtm[m])
    for m ∈ comm
        if qtm0[m] == 0
            fix(qtm[m], 0, force=true)
        end
    end

    @constraint(GTAPMod, eq_qst[m in comm, r in reg; qst0[m,r] != 0],
        qst[m,r] - (qtm[m] * (pt[m] / pds[m,r])^ESUBS[m]) ⟂ qst[m,r])
    for m ∈ comm, r ∈ reg
        if qst0[m,r] == 0
            fix(qst[m,r], 0, force=true)
        end
    end

    @constraint(GTAPMod, eq_pt[m in comm; qtm0[m] != 0],
        if ESUBS[m] == 1.0
            (pt[m] - (A_qtm[m]/pt0[m]) * prod((pds[m,r]*(pds0[m,r]/shr_qst[m,r]))^shr_qst[m,r] for r in reg))
        else
            (pt[m]^(1-ESUBS[m]) - sum(shr_qst[m,r]*pds[m,r]^(1-ESUBS[m]) for r in reg))  
        end
        ⟂ pt[m])
    for m ∈ comm
        if qtm0[m] == 0
            fix(pt[m], 1, force=true)
        end
    end

    #   Domestic goods market equilibrium ----------------------------------------------------------

    #   Total domestic demand for domestic goods
    @constraint(GTAPMod, eq_qds[i in comm, r in reg; qds0[i,r] != 0],
        qds[i,r] - (sum(qfd[i,a,r]*(qfd0[i,a,r]/qds0[i,r]) for a in acts) + qpd[i,r]*(qpd0[i,r]/qds0[i,r]) 
                    + qgd[i,r]*(qgd0[i,r]/qds0[i,r]) + qid[i,r]*(qid0[i,r])/qds0[i,r]) ⟂ qds[i,r])
    for i ∈ comm, r ∈ reg
        if qds0[i,r] == 0
            fix(qds[i,r], 0, force=true)
        end
    end

    #   Market clearing for domestically produced goods
    @constraint(GTAPMod, eq_qc[i in comm, r in reg; qc0[i,r] != 0],
        qc[i,r] - (qds[i,r]*(qds0[i,r]/qc0[i,r]) + sum(qxs[i,r,d]*(qxs0[i,r,d]/qc0[i,r]) for d in reg) + qst[i,r]*(qst0[i,r]/qc0[i,r])) ⟂ qc[i,r])
    for i ∈ comm, r ∈ reg
        if qc0[i,r] == 0
            fix(qc[i,r], 0, force=true)
        end
    end

    #   Factor market equilibrium ------------------------------------------------------------------

    #   Factor supply and equilibrium
    @constraint(GTAPMod, eq_pes[f in endw, a in acts, r in reg; qfe0[f,a,r] != 0],
        if mobFact[f]
            if ETRAE[f,r] != Inf
                qfe[f,a,r] - (qe[f,r] * (pes[f,a,r]/pe[f,r])^ETRAE[f,r])
            else
                pes[f,a,r] - pe[f,r]
            end
        else
            qfe[f,a,r] - (pes[f,a,r] / pfact[r])^ETAFF[f,a,r]
        end
        ⟂ pes[f,a,r])
    for f ∈ endw, a ∈ acts, r ∈ reg
        if qfe0[f,a,r] == 0
            fix(pes[f,a,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_pe[f in endw, r in reg; qe0[f,r] != 0],
        if mobFact[f]
            if ETRAE[f,r] != Inf
                pe[f,r]^(1 + ETRAE[f,r]) - sum(shr_qfes[f,a,r]*pes[f,a,r]^(1 + ETRAE[f,r]) for a in acts)
            else
                qe[f,r] - sum(shr_qfes[f,a,r]*qfe[f,a,r] for a in acts)
            end
        else
            pe[f,r] - sum(shr_qfes[f,a,r]*pes[f,a,r] for a in acts)
        end
        ⟂ pe[f,r])
    for f ∈ endw, r ∈ reg
        if qe0[f,r] == 0
            fix(pe[f,r], 1, force=true)
        end
    end

    @constraint(GTAPMod, eq_qe[f in endw, r in reg; qe0[f,r] != 0],
    if mobFact[f]
        #   By default, aggregate factor stocks are fixed
        qe[f,r] - chi_qe[f,r]
    else
        #   For sector specific factor, this is purely definitional
        pe[f,r]*qe[f,r] - sum(shr_qfes[f,a,r]*pes[f,a,r]*qfe[f,a,r] for a in acts)
    end
    ⟂ qe[f,r])
    for f ∈ endw, r ∈ reg
        if qe0[f,r] == 0
            fix(qe[f,r], 0, force=true)
        end
    end

    #   Factor returns at purchasers prices
    @constraint(GTAPMod, eq_pfe[f in endw, a in acts, r in reg; qfe0[f,a,r] != 0],
        pfe[f,a,r] - (1 + tfe[f,a,r])*peb[f,a,r]*(peb0[f,a,r]/pfe0[f,a,r]) ⟂ pfe[f,a,r])
    for f ∈ endw, a ∈ acts, r ∈ reg
        if qfe0[f,a,r] == 0
            fix(pfe[f,a,r], 1, force=true)
        end
    end

    #   Factor returns at basic prices
    @constraint(GTAPMod, eq_peb[f in endw, a in acts, r in reg; qfe0[f,a,r] != 0],
        pes[f,a,r] - (1 - tinc[f,a,r])*peb[f,a,r]*(peb0[f,a,r]/pes0[f,a,r]) ⟂ peb[f,a,r])
    for f ∈ endw, a ∈ acts, r ∈ reg
        if qfe0[f,a,r] == 0
            fix(peb[f,a,r], 1, force=true)
        end
    end

    #   Macro aggregates
    @constraint(GTAPMod, eq_gdpfc[r in reg],
        gdpfc[r] - y[r] ⟂ gdpfc[r])

    #   Capital account closure --------------------------------------------------------------------

    #   Non-normalized capital stock--change is equal to change in normalized capital stock
    #   !!!! In the GEMPACK code, isn't VES(ENDWC,r) / GROSSCAP(r) always equal to 1⟂
    @constraint(GTAPMod, eq_kb[r in reg],
        kb[r] - qe[capName,r] ⟂ kb[r])

    #   End-of-period capital stock
    @constraint(GTAPMod, eq_ke[r in reg],
        ke[r] - ((1 - depr[r])*kb[r]*(kb0[r]/ke0[r]) + qinv[r]*(qinv0[r]/ke0[r])) ⟂ ke[r])

    #   After-tax gross rate of return to capital
    @constraint(GTAPMod, eq_rental[r in reg],
        rental[r]*kb[r] - (sum(pes[capName,a,r]*qfe[capName,a,r]*(pes0[capName,a,r]*qfe0[capName,a,r]/(rental0[r]*kb0[r])) for a in acts)) ⟂ rental[r])

    #   Net rate of return to capital
    @constraint(GTAPMod, eq_rorc[r in reg],
        rorc[r] - ((rental0[r]/(rorc0[r]*pinv0[r]))*rental[r]/pinv[r] - depr[r]/rorc0[r]) ⟂ rorc[r])

    #   Expected rate of return to capital
    @constraint(GTAPMod, eq_rore[r in reg],
        rore[r] - (((rorc0[r]/rore0[r])*(kb0[r]/ke0[r])^RORFLEX[r])*rorc[r]*(kb[r]/ke[r])^RORFLEX[r]) ⟂ rore[r])

    #   Default capital account closure--determination of capital flows
    
    if savfClosure == :RoRFlex
        #   Equal rates of return (for all regions)
        @constraint(GTAPMod, eq_savf[r in reg],
            rorg - risk[r]*rore[r]*(rore0[r]/rorg0) ⟂ savf[r])
    elseif savfClosure == :capFix
        #   Fixed flows in real terms for all regions save residual
        @constraint(GTAPMod, eq_savf[r in reg; r != resName],
            savf[r] - pnum*savf0[r] ⟂ savf[r])
    elseif savfClosure == :capShrFix
        #   Fixed flows wrt GDP for all regions save residual
        @constraint(GTAPMod, eq_savf[r in reg; r != resName],
            savf[r] - savf0[r]*gdpfc[r] ⟂ savf[r])
    end

    #   Foreign saving must add to 0
    if savfClosure == :RoRFlex
        #   Determines the global rate of return
        @constraint(GTAPMod, eq_rorg,
            sum(savf[r] for r in reg) ⟂ rorg)
    else
        #   Determines foreign saving for the residual region
        @constraint(GTAPMod, eq_savfRes[r in reg; r == resName],
            sum(savf[s] for s in reg) ⟂ savf[r])
        #   rorg is the weighted average of all the rore's
        @constraint(GTAPMod, eq_rorg,
            rorg*rorg0*sum((pinv[r]*qinv[r]*pinv0[r]*qinv0[r] - depr[r]*kb[r]*kb0[r]) for r in reg) 
                - sum(rore[r]*rore0[r]*(pinv[r]*qinv[r]*pinv0[r]*qinv0[r] - depr[r]*kb[r]*kb0[r]) for r in reg) ⟂ rorg)
    end

    #   Savings = investment
    @constraint(GTAPMod, eq_yi[r in reg],
        yi[r] - (save[r]*(save0[r]/yi0[r]) + savf[r]/yi0[r] + depr[r]*pinv[r]*kb[r]*(pinv0[r]*kb0[r]/yi0[r])) + ifRes[r]*Walras/yi0[r] ⟂ yi[r])

    @constraint(GTAPMod, eq_chisave,
        chisave - sum(invwgt[r]*pinv[r]/pinvLag[r] for r in reg) / sum(savwgt[r]*psave[r]/psaveLag[r] for r in reg) ⟂ chisave)

    @constraint(GTAPMod, eq_psave[r in reg],
        psave[r]/psaveLag[r] - chisave*(pinv[r]/pinvLag[r]) ⟂ psave[r])

    @constraint(GTAPMod, eq_pfact[r in reg],
        pfact[r] - pfact0[r]*sqrt((
            (sum(peb0[f,a,r]*qfe0[f,a,r]*peb[f,a,r] for f in endw, a in acts))/
            (sum(peb0[f,a,r]*qfe0[f,a,r] for f in endw, a in acts)))
            *(
            (sum(peb0[f,a,r]*qfe0[f,a,r]*peb[f,a,r]*qfe[f,a,r] for f in endw, a in acts))/
            (sum(peb0[f,a,r]*qfe0[f,a,r]*qfe[f,a,r] for f in endw, a in acts)))) ⟂ pfact[r])

    @constraint(GTAPMod, eq_pfactw,
            pfactw - pfactw0*sqrt((
                (sum(peb0[f,a,r]*qfe0[f,a,r]*peb[f,a,r] for f in endw, a in acts, r in reg))/
                (sum(peb0[f,a,r]*qfe0[f,a,r] for f in endw, a in acts, r in reg))
                )*(
                (sum(peb0[f,a,r]*qfe0[f,a,r]*peb[f,a,r]*qfe[f,a,r] for f in endw, a in acts, r in reg))/
                (sum(peb0[f,a,r]*qfe0[f,a,r]*qfe[f,a,r] for f in endw, a in acts, r in reg))
                )) ⟂ pfactw)

        #   Numeraire definition
    @constraint(GTAPMod, eq_Walras,
        pnum - pfactw ⟂ Walras)

    return GTAPMod
end