# dyn.load("lib/Rfuncs.so")

linreg_path <- function(Y,X,val,idx,jdx,lambda_graph, 
						lambda_sparse,p,m,rho=1.0,
						eps_abs=1e-5,eps_rel=1e-4,
						max_num_iter=1e4,
						linearized_ADMM=FALSE, standard_ADMM=FALSE,
						nonadaptive_varying_rho=FALSE, constant_rho=FALSE,
						general_A=FALSE, graph_only=FALSE,
						graph_sparse=TRUE, reporting=FALSE)
{
	adaptive_varying_rho=TRUE
    if (nonadaptive_varying_rho && constant_rho) { stop("nonadaptive_varying_rho and constant_rho can't both be TRUE! \n") }
    if (constant_rho) { nonadaptive_varying_rho = FALSE; adaptive_varying_rho=FALSE }
	if (nonadaptive_varying_rho) { constant_rho = FALSE; adaptive_varying_rho=FALSE }
	special_ADMM = TRUE;
	if (standard_ADMM && linearized_ADMM) { stop("standard_ADMM and linearized_ADMM can't both be TRUE! \n") }
	if (standard_ADMM) { linearized_ADMM = FALSE; special_ADMM = FALSE }
	if (linearized_ADMM) { standard_ADMM = FALSE; special_ADMM = FALSE }	 
	
    out <- .Call('linreg',Y,t(X),val,
    as.integer(idx),as.integer(jdx),
	lambda_graph,lambda_sparse, 
	as.integer(p), as.integer(m),rho,
    eps_abs,eps_rel,as.integer(max_num_iter),
	 as.logical(special_ADMM), as.logical(linearized_ADMM), as.logical(standard_ADMM),
	 as.logical(adaptive_varying_rho),as.logical(nonadaptive_varying_rho),as.logical(constant_rho),
	 as.logical(general_A), as.logical(graph_only),
	 as.logical(graph_sparse),as.logical(reporting))
	 out$beta_path <- array(out$beta_path,c(p,length(lambda_sparse),length(lambda_graph)))
	 if (reporting) { out$fun_path <- array(out$fun_path,max_num_iter,c(length(lambda_sparse),length(lambda_graph))) }
	return (out)
}


filtering <- function(Y, val, idx, jdx, lambda, m, rho=1.0, eps=1e-3,max_iter=1e5,
                        variantADMM=TRUE,varying_rho=TRUE,general_A=FALSE,reporting=FALSE)
{
    out <- .Call('filtering',Y,val,as.integer(idx),as.integer(jdx),lambda,as.integer(m),
                rho,eps,as.integer(max_iter),as.logical(variantADMM),as.logical(varying_rho),
                as.logical(general_A),as.logical(reporting))
    out$beta_path <- matrix(out$beta_path,p,length(lambda))
    if (reporting) { out$fun_path <- matrix(out$fun_path, max_iter,length(lambda)) }
    return (out)
}



