#ifdef NO_PRE_ACC
void edge_pushing_s(short           tnum,
                    map<locint, map<locint,double> > *graph,
                    locint*         edge_index,
                    double*         edge_value,
                    unsigned int    edge_index_len,
                    unsigned int    edge_value_len,
                    unsigned int    max_index
){
#endif  // NO_PRE_ACC

#ifdef PRE_ACC
void edge_pushing_pre_s(short           tnum,
                    map<locint, map<locint,double> > *graph,
                    locint*         edge_index,
                    double*         edge_value,
                    unsigned int    edge_index_len,
                    unsigned int    edge_value_len,
                    unsigned int    max_index
){
#endif  // PRE_ACC
    unsigned char operation;
    derivative_info* info=new derivative_info();
    map<locint, double> *Adjoints= new map<locint, double>;
#ifdef PRE_ACC
    int edge_stmt_cnt = 0;
    EdgeLocalGraph *local_graph = new EdgeLocalGraph();
    derivative_info* dinfo=new derivative_info();
    dinfo->r=NULLLOC;dinfo->x=NULLLOC;dinfo->y=NULLLOC;
    dinfo->dx=0.0;dinfo->dy=0.0;
    dinfo->px=0.0;dinfo->py=0.0;dinfo->pxy=0.0;
    locint *dp= new locint[MAX_TEMP_ARRAY_SIZE];
    int dl=0;
    int r;
#endif

    unsigned int i,j;
    locint p;
    double w;
    double vx,vy,vr,vc;
    unsigned int tl;
    locint *tp = new locint[max_index];
    double *tw = new double[max_index];
    map<locint, double> *edges;

    init_rev_sweep(tnum);
    operation=get_op_r();
    while (operation != start_of_tape) { 
        info->opcode=operation;
        info->r=NULLLOC;info->x=NULLLOC;info->y=NULLLOC;
        info->dx=0.0;info->dy=0.0;
        info->px=0.0;info->py=0.0;info->pxy=0.0;
        switch (operation){
                /************************************************************/
                /*                                                  MARKERS */
                /*----------------------------------------------------------*/
            case end_of_op:                                    /* end_of_op */
                get_op_block_r();
                operation = get_op_r();
                /* Skip next operation, it's another end_of_op */
                break;
            case start_of_tape:                            /* start_of_tape */
            case end_of_tape:                                /* end_of_tape */
                break;
                /************************************************************/
                /*                                               COMPARISON */
                /*----------------------------------------------------------*/
            case eq_zero  :                                      /* eq_zero */
            case neq_zero :                                     /* neq_zero */
            case gt_zero  :                                      /* gt_zero */
            case lt_zero :                                       /* lt_zero */
            case ge_zero :                                       /* ge_zero */
            case le_zero :                                       /* le_zero */
                break;
                /************************************************************/
                /*                                              ASSIGNMENTS */
            case subscript:
            case ref_copyout:
            case ref_assign_a:
            case assign_a:     /* assign an adouble variable an    assign_a */
                /* adouble value. (=) */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                info->dx=1.0;
                edge_value_len-=2;
                break;
            
            case neq_a_a:
            case eq_a_a:
            case le_a_a:
            case ge_a_a:
            case lt_a_a:
            case gt_a_a:
            case ref_assign_d_zero:
            case ref_assign_d_one:
            case ref_assign_d:
            case assign_d:      /* assign an adouble variable a    assign_d */
            case assign_d_zero: /* assign an adouble a        assign_d_zero */
            case assign_d_one:  /* double value. (=)           assign_d_one */
                info->r=edge_index[--edge_index_len];
                edge_value_len-=1;
                break;

            case ref_assign_ind:
            case assign_ind:       /* assign an adouble variable an    assign_ind */
                /* independent double value (<<=) */
                edge_index_len-=1;
                edge_value_len-=1;
                break;
	
            case assign_dep:           /* assign a float variable a    assign_dep */
                /* dependent adouble value. (>>=) */
                (*Adjoints)[edge_index[--edge_index_len]]=1.0;
                edge_value_len-=1;
                break;

            /****************************************************************************/
            /*                                                   OPERATION + ASSIGNMENT */
            /*--------------------------------------------------------------------------*/
            case ref_eq_plus_d:
            case eq_plus_d:            /* Add a floating point to an    eq_plus_d */
                /* adouble. (+=) */
                break;

            case ref_eq_plus_a:
            case eq_plus_a:             /* Add an adouble to another    eq_plus_a */
                /* adouble. (+=) */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                info->y=edge_index[--edge_index_len];
                info->dx=1.0;info->dy=1.0;
                edge_value_len-=3;
                break;

            case ref_eq_min_d:
            case eq_min_d:       /* Subtract a floating point from an    eq_min_d */
                /* adouble. (-=) */
                break;

            case ref_eq_min_a:
            case eq_min_a:        /* Subtract an adouble from another    eq_min_a */
                /* adouble. (-=) */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                info->y=edge_index[--edge_index_len];
                info->dx=1.0;info->dy=-1.0;
                edge_value_len-=3;
                break;

            case ref_eq_mult_d:
            case eq_mult_d:              /* Multiply an adouble by a    eq_mult_d */
                /* flaoting point. (*=) */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                info->dx=edge_value[--edge_index_len];
                edge_value_len-=3;
                break;

            case ref_eq_mult_a:
            case eq_mult_a:       /* Multiply one adouble by another    eq_mult_a */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                info->y=edge_index[--edge_index_len];
                info->dx=edge_value[edge_index_len];
                info->dy=edge_value[edge_index_len+1];
                info->pxy=1.0;
                edge_value_len-=3;
                break;
            case eq_plus_prod:   /* increment a product of           eq_plus_prod */
                /* two adoubles (*) */
#ifdef NO_PRE_ACC
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                info->y=edge_index[--edge_index_len];
                info->dx=1.0;info->dy=1.0;
                compute_pushing(tl,tp,tw,info,graph);
                compute_adjoints(info,Adjoints);
                info->r=info->y;
                info->x=edge_index[--edge_index_len];
                info->y=edge_index[--edge_index_len];
                info->dx=edge_value[edge_index_len];
                info->dy=edge_value[edge_index_len+1];
                info->pxy=1.0;
#endif

#ifdef PRE_ACC
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                info->y=edge_index[--edge_index_len];
                info->dx=1.0;dinfo->dy=1.0;
                info->r=info->y;

                info->x=edge_index[--edge_index_len];
                info->y=edge_index[--edge_index_len];
                info->dx=edge_value[edge_index_len];
                info->dy=edge_value[edge_index_len+1];
                info->pxy=1.0;
#endif
                edge_value_len-=5;
                break;
            case eq_min_prod:   /* decrement a product of             eq_min_prod */
                /* two adoubles (*) */
#ifdef NO_PRE_ACC
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                info->y=edge_index[--edge_index_len];
                info->dx=1.0;info->dy=-1.0;
                compute_pushing(tl,tp,tw,info,graph);
                compute_adjoints(info,Adjoints);

                info->r=info->y;
                info->x=edge_index[--edge_index_len];
                info->y=edge_index[--edge_index_len];
                info->dx=edge_value[edge_index_len];
                info->dy=edge_value[edge_index_len+1];
                info->pxy=1.0;
                edge_value_len-=5;
#endif

#ifdef PRE_ACC
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                info->y=edge_index[--edge_index_len];
                info->dx=1.0;info->dy=-1.0;

                info->r=info->y;
                info->x=edge_index[--edge_index_len];
                info->y=edge_index[--edge_index_len];
                info->dx=edge_value[edge_index_len];
                info->dy=edge_value[edge_index_len+1];
                info->pxy=1.0;
                edge_value_len-=5;
#endif
                break;

                /*--------------------------------------------------------------------------*/
            case ref_incr_a:
            case ref_decr_a:
            case incr_a:                        /* Increment an adouble    incr_a */
            case decr_a:                        /* Increment an adouble    decr_a */
                break;

                /****************************************************************************/
                /*                                                        BINARY OPERATIONS */
                /*--------------------------------------------------------------------------*/
            case plus_a_a:                 /* : Add two adoubles. (+)    plus a_a */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                info->y=edge_index[--edge_index_len];
                info->dx=1.0;
                info->dy=1.0;
                edge_value_len-=3;
                break;

            case plus_d_a:             /* Add an adouble and a double    plus_d_a */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                info->dx=1.0;
                edge_value_len-=2;
                break;

            case min_a_a:              /* Subtraction of two adoubles    min_a_a */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                info->y=edge_index[--edge_index_len];
                info->dx=-1.0;
                info->dy=1.0;
                edge_value_len-=3;
                break;

            case min_d_a:                /* Subtract an adouble from a    min_d_a */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                info->dx=-1.0;
                edge_value_len-=2;
                break;

            case mult_a_a:               /* Multiply two adoubles (*)    mult_a_a */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                info->y=edge_index[--edge_index_len];
                info->dx=edge_value[edge_index_len];
                info->dy=edge_value[edge_index_len+1];
                info->pxy=1.0;
                edge_value_len-=3;
                break;

            case mult_d_a:         /* Multiply an adouble by a double    mult_d_a */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                info->dx=edge_value[--edge_index_len];
                edge_value_len-=3;
                break;

            case div_a_a:           /* Divide an adouble by an adouble    div_a_a */
                info->r=edge_index[--edge_index_len];
                info->y=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                vr=edge_value[--edge_value_len];
                vy=edge_value[--edge_value_len];
                vx=edge_value[--edge_value_len];
                info->dx=1.0/vy;
                info->dy=-vx/(vy*vy);
                info->py=2.0*vx/(vy*vy*vy);
                info->pxy=-1.0/(vy*vy);
                break;
            case div_d_a:             /* Division double - adouble (/)    div_d_a */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                edge_index_len-=1;
                vr=edge_value[--edge_value_len];
                vy=edge_value[--edge_value_len];
                vx=edge_value[--edge_value_len];
                info->dx=-vx/(vy*vy);
                info->px=2.0*vx/(vy*vy*vy);
                break;

                /****************************************************************************/
                /*                                                         SIGN  OPERATIONS */
                /*--------------------------------------------------------------------------*/
            case pos_sign_a:                                        /* pos_sign_a */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                info->dx=1.0;
                edge_value_len-=2;
                break;
		
            case neg_sign_a:                                        /* neg_sign_a */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                info->dx=-1.0;
                edge_value_len-=2;
                break;

                /****************************************************************************/
                /*                                                         UNARY OPERATIONS */
                /*--------------------------------------------------------------------------*/
            case exp_op:                          /* exponent operation    exp_op */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                vr=edge_value[--edge_value_len];
                info->dx=vr;
                info->px=vr;
                edge_value_len-=1;
                break;

            case sin_op:                              /* sine operation    sin_op */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                vr=edge_value[--edge_value_len];
                vx=edge_value[--edge_value_len];
                info->dx=cos(vx);
                info->px=-vr;
                break;

            case cos_op:                            /* cosine operation    cos_op */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                vr=edge_value[--edge_value_len];
                vx=edge_value[--edge_value_len];
                info->dx=-sin(vx);
                info->px=-vr;
                break;

            case atan_op:                                            /* atanh_op */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                vr=edge_value[--edge_value_len];
                vx=edge_value[--edge_value_len];
                info->dx=1.0/(1.0+vx*vx);
                info->px=-(2.0*vx)/((1.0+vx*vx)*(1.0+vx*vx));
                break;
            case asin_op:                                              /* asin_op */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                vr=edge_value[--edge_value_len];
                vx=edge_value[--edge_value_len];
                info->dx=1.0/sqrt(1.0-vx*vx);
                info->px=vx/((sqrt(1.0-vx*vx))*(1.0-vx*vx));
                break;
            case acos_op:                                              /* asin_op */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                vr=edge_value[--edge_value_len];
                vx=edge_value[--edge_value_len];
                info->dx=-1.0/sqrt(1.0-vx*vx);
                info->px=-vx/((sqrt(1.0-vx*vx))*(1.0-vx*vx));
                break;

            case asinh_op:                                            /* asinh_op */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                vr=edge_value[--edge_value_len];
                vx=edge_value[--edge_value_len];
                info->dx=1.0/sqrt(1.0+vx*vx);
                info->px=-vx*info->dx/(1.0+vx*vx);
                break;
            case acosh_op:                                            /* acosh_op */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                vr=edge_value[--edge_value_len];
                vx=edge_value[--edge_value_len];
                info->dx=1.0/sqrt(vx*vx-1.0);
                info->px=-vx*info->dx/(vx*vx-1.0);
                break;
            case atanh_op:                                            /* atanh_op */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                vr=edge_value[--edge_value_len];
                vx=edge_value[--edge_value_len];
                info->dx=1/(1-vx*vx);
                info->px=2.0*vx*info->dx*info->dx;
                break;
            case erf_op:                                              /* erf_op   */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                vr=edge_value[--edge_value_len];
                vx=edge_value[--edge_value_len];
                info->dx=2.0/sqrt(acos(-1.0))*exp(-vx*vx);
                info->px=-2.0*vx*info->dx;
                break;

            case log_op:                                                /* log_op */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                vr=edge_value[--edge_value_len];
                vx=edge_value[--edge_value_len];
                info->dx=1.0/vx;
                info->px=-1.0/(vx*vx);
                break;

            case pow_op:                                                /* pow_op */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                edge_index_len-=1;
                vr=edge_value[--edge_value_len];
                vx=edge_value[--edge_value_len];
                vc=edge_value[--edge_value_len];
                if (vx==0.0){
                    info->dx=0.0;info->px=0.0;
                }
                else{
                    info->dx=vc*(vr/vx);
                    info->px=(vc-1.0)*(info->dx/vx);
                }
                break;

            case sqrt_op:                                              /* sqrt_op */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                vr=edge_value[--edge_value_len];
                vx=edge_value[--edge_value_len];
                if (vx==0.0){
                    info->dx=0.0;info->px=0.0;
                }
                else{
                    info->dx=0.5*(vr/vx);
                    info->px=-0.5*(info->dx/vx);
                }
                break;

            case min_op:                                                /* min_op */
                info->r=edge_index[--edge_index_len];
                edge_index_len-=2;
                vr=edge_value[--edge_value_len];
                vx=edge_value[--edge_value_len];
                vy=edge_value[--edge_value_len];
                if (vx>vy){
                    info->x=edge_index[edge_value_len];
                    info->dx=1.0;
                }
                else if(vx<vy){
                    info->x=edge_index[edge_value_len+1];
                    info->dx=1.0;
                }
                else{  //tie here, what to do?
                    fprintf(DIAG_OUT,"WARNING: Tie in min_op\n");
                }
                break;
            case abs_val:                                              /* abs_val */
                info->r=edge_index[--edge_index_len];
                info->x=edge_index[--edge_index_len];
                vr=edge_value[--edge_value_len];
                vx=edge_value[--edge_value_len];
                if (vx>0.0){
                    info->dx=1.0;
                }
                else if(vx<0.0){
                    info->dx=-1.0;
                }
                else{
                    fprintf(DIAG_OUT,"WARNING: 0 in abs_val\n");
                }
                break;
                /*                                                             CONDITIONALS */
                /*--------------------------------------------------------------------------*/
            case ref_cond_assign:
            case cond_assign:                                      /* cond_assign */
                edge_value_len-=4;
                vc=edge_value[edge_value_len];
                if (vc>0.0){
                    info->r=edge_index[--edge_index_len];
                    edge_index_len-=1;
                    info->x=edge_index[--edge_index_len];
                    edge_index_len-=1;
                }
                else{
                    info->r=edge_index[--edge_index_len];
                    info->x=edge_index[--edge_index_len];
                    edge_index_len-=2;
                }
                info->dx=1.0;
                break;

            case ref_cond_assign_s:
            case cond_assign_s:                                  /* cond_assign_s */
                edge_value_len-=3;
                vc=edge_value[edge_value_len];
                if (vc>0.0){
                    info->r=edge_index[--edge_index_len];
                    info->x=edge_index[--edge_index_len];
                    info->dx=1.0;
                    edge_index_len-=1;
                }
                break;

            case subscript_ref:
	        case death_not:
	        case gen_quad:
	        case end_of_int:
	        case end_of_val:
	        case take_stock_op:
	        case ceil_op:
	        case floor_op:
	        case ext_diff:
	        case ext_diff_iArr:
		    break;
	        default:
                /* Die here, we screwed up */
                fprintf(DIAG_OUT,"EDGE_HESSIAN error: Fatal error in tape_doc for op %d\n",
                        operation);
                break;
 
        }//switch
#ifdef NO_PRE_ACC
        if (info->r!=NULLLOC){
//edge_check_info(info);
//pushing
            compute_pushing(tl,tp,tw,info,graph);
//creating
            compute_creating(info,Adjoints,graph);
//adjointing
            compute_adjoints(info,Adjoints);
//edge_check_graph(graph);
        }
#endif

#ifdef PRE_ACC
    switch (info->opcode){
        case assign_dep:
            dl=0;
            break;
        case ref_assign_ind:
        case assign_ind:
            break;
        case ref_assign_d:
        case ref_assign_d_one:
        case ref_assign_d_zero:
        case assign_d:
        case assign_d_one:
        case assign_d_zero:
            dp[dl++]=info->r;
            break;
        case ref_assign_a:
        case ref_eq_plus_a:
        case ref_eq_min_a:
        case ref_eq_mult_d:
        case ref_eq_mult_a:
        case ref_copyout:
        case ref_cond_assign:
        case ref_cond_assign_s:
        case assign_a:
        case eq_plus_a:
        case eq_min_a:
        case eq_mult_d:
        case eq_mult_a:
        case eq_plus_prod:
        case eq_min_prod:
        case cond_assign:
        case cond_assign_s:
//push previous result to Global Trace
//edge_check_graph(lGraph);
//            compute_global_pushing(tl,tp,tw,r,lAdjoints,graph);
//            compute_global_creating(r,lGraph,Adjoints,graph);
//            compute_global_adjoints(r,lAdjoints,Adjoints);
            ++edge_stmt_cnt;
            compute_global(tp, tw, local_graph, r, Adjoints, graph);
            for(i=0;i<dl;i++){
                dinfo->r=dp[i];
                compute_pushing(tl,tp,tw,dinfo,graph);
                compute_adjoints(dinfo,Adjoints);
            }
            dl=0;
            local_graph->reset();
            if ((info->opcode==eq_plus_prod)||(info->opcode==eq_min_prod)){
                dinfo->r=edge_index[edge_index_len+4];
                dinfo->x=edge_index[edge_index_len+3];
                dinfo->y=edge_index[edge_index_len+2];
                dinfo->dx=1.0;
                if (info->opcode==eq_plus_prod){
                    dinfo->dy=1.0;
                }
                else{
                    dinfo->dy=-1.0;
                }
                size_t ind_r = local_graph->AddLiveVar(dinfo->r);
                local_graph->adjoints[ind_r] = 1.0;
                r=dinfo->r;
                compute_local(tp, tw, dinfo, local_graph);
                dinfo->dx=0.0;dinfo->dy=0.0;
                dinfo->x=NULLLOC;dinfo->y=NULLLOC;
            } else {
                size_t ind_r = local_graph->AddLiveVar(info->r);
                local_graph->adjoints[ind_r] = 1.0;
                r=info->r;
            }
//edge_check_graph(graph);
//edge_check_adjoints(Adjoints,10);
            break;
        default:
            ; 
        }//switch
//edge_check_info(info);
//local_graph->Print();
        if (info->r!=NULLLOC){
// compute pushing, creating, adjoints
            compute_local(tp, tw, info, local_graph);
        }
//local_graph->Print();
#endif
        operation=get_op_r();
    }//while

#ifdef PRE_ACC
    compute_global(tp, tw, local_graph, r, Adjoints, graph);
    for(i=0;i<dl;i++){
        dinfo->r=dp[i];
        compute_pushing(tl,tp,tw,dinfo,graph);
        compute_adjoints(dinfo,Adjoints);
    }
    dl=0;
//edge_check_graph(graph);
//edge_check_adjoints(Adjoints,10);
#ifdef EDGE_DEBUG
    printf("edge_stmt_cnt = %d\n", edge_stmt_cnt);
#endif  // EDGE_DEBUG

    delete dinfo;
    delete local_graph;
    delete[] dp;
#endif
    end_sweep();
    delete info;
    delete Adjoints;
    delete[] tp;
    delete[] tw;
}


