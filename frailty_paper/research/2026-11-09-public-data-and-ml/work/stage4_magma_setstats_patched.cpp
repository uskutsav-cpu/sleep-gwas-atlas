/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#include "setstats.h"

using namespace SetStatsUtils;

void SetVariables::register_block(VariableIndex* index, int block, VariableType type) {
  CStringVariable& var_name = index->name; NumericID<int>& var_id = index->variable_id; var_id.init();
  for (int v = 0; v < index->data_total(); v++) {
    if (index->active_id(v)) var_id.set(v, variables.add_variable(index->data_id, var_name.get(v), block, v, type));
  }
}

void SetVariables::register_interactions(VariableIndex* index, int block, InteractionType type) {
  NumericID<int>& var_id = index->variable_id; var_id.init();
  NumericID<int> *main1 = index->typed_variable<NumericID<int> >("INTERACT_MAIN1"), *main2 = index->typed_variable<NumericID<int> >("INTERACT_MAIN2"); 
  if (!main1 || !main2) throw SetException("processing interactions", "INTERACT_MAIN variables are missing from VariableIndex");

  int total = index->data_used(true), offset = variables.extend_index(total);
  for (int v = 0; v < index->data_total(); v++) {
    if (index->active_id(v)) var_id.set(v, variables.set_interaction(offset++, main1->get(v), main2->get(v), block, v, type));
  }
}

int SetVariables::add_interactions(GeneSetData* ss) {
  int block = add_sets(ss);
  register_interactions(ss, block, it_InteractionSS);   
  return block;
}

int SetVariables::add_interactions(GeneCovarData<float>* sc) {
  int block = add_covar(sc);
  register_interactions(sc, block, it_InteractionSC);
  return block;
}


bool SetVariables::lock_variables() {
  variables.init(variables.no_vars(false).first);
  for (int s = 0; s < set_info.size(); s++) register_block(set_info[s], s, vt_GeneSet);
  for (int c = 0; c < covar_info.size(); c++) register_block(covar_info[c], c, vt_GeneCovar);  
  variables.build_index();
  
  return variables.get_duplicates().empty();
}

VariableInfo& SetVariables::get_variable(int var_id, bool check_data) {
  VariableInfo& info = variables.var_info(var_id); 
   
  if (!check_data || (info.is_valid() && info.get_column() >= 0)) return info;
  else return variables.var_info(-1);
}

int SetVariables::mean_sampsize() {
  if (!mean_N) {
    PositiveNumber<int>& nsamp = gene_info.nsamp;
    if (nsamp.initialised()) {
      int miss_code = nsamp.get_missing().second, count = 0; long long total = 0;
      BaseBuffer<int>& data = nsamp.get_buffer();
      for (int *curr = data.data(), *end = curr + data.size(); curr < end; curr++) {
        if (*curr != miss_code) {total += *curr; count++;}
      }
      mean_N = count > 0 ? round(total/double(count)) : -1;
    }  
  }
  return mean_N;
} 

double* SetVariables::get_outcome_data(int index) {           
  if (outcome_info.active_id(index) && outcome_info.is_valid(index)) {
    int col = outcome_info.storage_column.get(index);
    if (col >= 0) return outcome_info.data.get_buffer()[col];
  }
  return 0;
}

GeneMatrix& SetVariables::get_corrs(bool verbose) {
  if (!inv_corrs) {
    try {inv_corrs = new CompoundGeneMatrix(gene_info, gc_correction, pc_prune);}
    catch (const SetException& se) {_LOG.error("inverting gene correlation matrix") << se.get_msg() << endl; die();}
  }
  if (!inv_corrs->is_inverted()) inv_corrs->invert_matrix(verbose);
  return *inv_corrs;
}

void SetVariables::init_corrs() {
  _LOG << "Inverting gene-gene correlation matrix..." << endl;
  get_corrs(true);
}

void SetVariables::init_internal() {
  GeneCovarData<float>* var_info = new GeneCovarData<float>("COVAR_INTERNAL"); internal_block = add_covar(var_info); 
  NonNegativeNumber<short>* type = new NonNegativeNumber<short>; var_info->attach_variable(type, "INTERNAL_TYPE");

  vector<string> var_names = InternalVariable::names();
  int col = 0, no_core = InternalVariable::iv_COUNT, no_var = 2*no_core; var_info->init(no_var); 
  for (int log = false; log <= true; log++) {
    for (int i = 0; i < InternalVariable::iv_COUNT; i++) {
      int var_id = var_info->get_iid(string("#INTERNAL-") + (log ? "LOG" : "MAIN") + "##" + var_names[i]);
      var_info->storage_column.set(var_id, col++);
      type->set(var_id, i);
    }    
  }
  var_info->build_index(); var_info->link_data(gene_info); 
  set_internal(gene_info.nsnps, *var_info, InternalVariable::iv_GeneSize);
  set_internal(gene_info.nparam, *var_info, InternalVariable::iv_GeneDensity);
  set_internal(gene_info.nsamp, *var_info, InternalVariable::iv_SampleSize);
  set_internal(gene_info.mac, *var_info, InternalVariable::iv_InverseMac);

  BlockMultiVariable<NumericVariable<float> >& vars = var_info->data;
  float miss_code = vars.get_missing().second;
  for (int gene_id = 0; gene_id < gene_info.data_total(); gene_id++) {int count = 0;
    if (gene_info.active_id(gene_id)) {
      var_info->in_source.set(gene_id, true); 
      for (int i = 0; i < no_core; i++) {if (vars.get(gene_id, i) != miss_code) count += 2;}
    } 
    var_info->miss_count.set(gene_id, no_var - count);     
  }
}

///currently, ZSTAT is guaranteed to be complete
pair<int,int> SetVariables::process_missing(bool drop_genes, const string& impute_type, float impute_value) {
  vector<short> index_row(gene_info.data_total(), 0), index_any(gene_info.data_total(), 0); 

  for (int i = 0 ; i < covar_info.size(); i++) {GeneCovarData<float>& curr_info = *covar_info[i]; 
    BooleanVariable& observed = curr_info.in_source;
    NonNegativeNumber<int>& miss_count = curr_info.miss_count;
    BlockMultiVariable<NumericVariable<float> >& vars = curr_info.data;

    vector<int> invalid; float miss_code = vars.get_missing().second;

    for (int j = 0; j < curr_info.data_used(false); j++) {if (!curr_info.is_valid(j)) invalid.push_back(j);}
    for (int iid = 0; iid < gene_info.data_total(); iid++) {
      if (gene_info.active_id(iid)) {
        if (observed.get(iid)) {
          int miss = miss_count.get(iid);
          for (int j = 0; j < invalid.size(); j++) {if (vars.get(iid,invalid[j]) == miss_code) miss--;}
          if (miss > 0) index_any[iid] = true;  
        } else {index_row[iid] = true; index_any[iid] = true;}
      }
    }
  }
  int miss_genes = MathUtils::count(index_row), miss_partial = MathUtils::count(index_any) - miss_genes;

  if (impute_type == "drop") gene_info.filter_ids(index_any, false);
  else if (drop_genes) gene_info.filter_ids(index_row, false);
  gene_info.make_consecutive();    

  if (impute_type != "drop") {
    int nrow = gene_info.data_used(false); 
    bool use_mean = impute_type == "mean", use_median = impute_type == "median";

    for (int c = 0 ; c < covar_info.size(); c++) {
      GeneCovarData<float>& curr_info = *covar_info[c]; int ncol = curr_info.data_total();
      Buffer<float>& data = curr_info.data.get_buffer(); float miss_code = curr_info.data.get_missing().second; 
      curr_info.miss_count.get_data().assign_value(0);
      for (int v = 0; v < ncol; v++) {
        if (curr_info.is_valid(v)) {
          int nmiss = MathUtils::count_value(data[v], miss_code, nrow), nobs = nrow - nmiss;
          curr_info.values_observed.set(v, nobs);
          if (nobs > 0) {
            if (nmiss > 0) {
              float value = impute_value;
              if (use_mean) {double sum = 0;
                for (float *curr = data[v], *end = curr+nrow; curr < end; curr++) {if (*curr != miss_code) sum += *curr;}
                value = sum / nobs;
              } else if (use_median) {BaseBuffer<float> obs(nobs);
                for (float *read = data[v], *write = obs.data(), *end = read+nrow; read < end; read++) {if (*read != miss_code) *(write++) = *read;}              
                value = MathUtils::get_median(obs.data(), nobs, true);
              }
              for (float *curr = data[v], *end = curr+nrow; curr < end; curr++) {if (*curr == miss_code) *curr = value;}            
            }                
          } else Utils::fill_value(data[v], nrow, impute_value);
        }
      }      
    } 
  } else for (int c = 0; c < covar_info.size(); c++) {covar_info[c]->values_observed.get_data().assign_value(gene_info.data_used(false));}  
  
  return pair<int,int>(miss_genes,miss_partial);
}

void SetVariables::prep_sets(GeneSetData& set_block) {
  for (int v = 0; v < set_block.data_total(); v++) {
    int used = set_block.is_valid(v) ? set_block.data.count(v) : 0;
    set_block.used_genes.set(v, used);
  }
}

void SetVariables::prep_outcome(double trunc_low, double trunc_high) {
  outcome_info.storage_column.set(outcome_info.get_iid("ZSTAT"), 0);
  outcome_info.build_index(); outcome_info.link_data(gene_info); 

  BaseBuffer<double>& zstat = gene_info.zstat.get_buffer();
  outcome_info.data.get_buffer().assign_col(zstat.data(), 0);
  
  prep_covar(outcome_info, trunc_low, trunc_high, true);
}

void SelfContainedComputation::set_variable(int var_id) {
  error_msg.clear(); variable = &variable_data.get_variable(var_id, true);
  if (!variable || !variable->is_valid()) failure("trying to analyse variable that is unknown or unavailable");
  if (!variable->is_set()) failure(string("variable '") + variable->get_name() + "' is not a gene set");
  no_genes = variable->as_set()->get_ngenes();
}

void SelfContainedComputation::run() {
  set<int>& data = *variable->as_set()->get_data(); 
    
  CompoundGeneMatrix sigma(variable_data.gene_info, data, variable_data.gc_correction, variable_data.pc_prune);
  sigma.invert_matrix(false);
  
  buffer.resize(no_genes); sigma.product(buffer.data()); 
  double *read_i = buffer.data(), *read_z = variable_data.get_outcome_data(), iti = 0, itz = 0;
  for (set<int>::iterator it = data.begin(); it != data.end(); ++it, ++read_i) {iti += *read_i; itz += *read_i * read_z[*it];}
                                
  if (iti > 0) {
    try {
      mu = itz / iti;
      pval = converter.convert(mu * sqrt(iti));
    } catch (const StatException& se) {failure("unable to compute p-value (" + se.get_msg() + ")");}
  } else {
    mu = 0; pval = -1;
    failure("standard error is invalid");
  }
}



void RegressionComputation::failure(const string& msg, const string& sub_msg) {
  error_msg = msg; 
  if (!sub_msg.empty()) error_msg += ": " + sub_msg;
  throw SetException("using stat computation", msg);
}

void RegressionComputation::error(const string& msg, const string& sub_msg) {
  error_msg.clear(); 
  string full_msg = msg;
  if (!sub_msg.empty()) full_msg += "; " + sub_msg;
  throw SetException("using stat computation", full_msg);
}

bool RegressionComputation::check_init(bool has_init) {
  if (has_init && !initialised) error("trying to use an uninitialised computation module");
  if (!has_init && initialised) error("trying to modify an active computation module");  
  return (has_init == initialised);
}

double RegressionComputation::product(double* var, int len) {return MathUtils::sum(var, len);}
double RegressionComputation::product(float* var_A, double* var_B, int len) {return MathUtils::in_product(var_A, var_B, len);}
double RegressionComputation::product(double* var_A, double* var_B, int len) {return MathUtils::in_product(var_A, var_B, len);}
double RegressionComputation::product(set<int>& var_A, double* var_B) {double value = 0;
  for (set<int>::iterator it = var_A.begin(); it != var_A.end(); ++it) value += var_B[*it];
  return value;
}

void RegressionComputation::variable_sum(BaseBuffer<double>& target, VariableInfo& var, double scale) {
  double* write = target.data();
  if (var.is_set()) {
    set<int>& data = *var.as_set()->get_data();
    for (set<int>::iterator it = data.begin(); it != data.end(); ++it) write[*it] += scale;
  } else if (var.is_covar()) {
    float* data = var.as_covar()->get_data();
    for (double* end = write + no_rows; write < end; write++) *write += *(data++)*scale;
  } else error("trying to compute with undefined variable (sum)");
}
 
void RegressionComputation::init(int nvar) {
  if (check_init(false)) {initialised = true;
    collinear_core = collinear_mode = ct_None;
    no_rows = variable_data.gene_info.data_used(false);
  
    init_buffers();
    if (resid_core && variables.size() > 1) {
      init_model(nvar, true);
      residualize();
    }

    init_model(nvar, false);
    fit_model(mt_Intercept); 
    fit_model(mt_BaseModel);
  }
}

void RegressionComputation::init_buffers() {
  if (buffers.empty()) {for (int b = 0; b < bt_COUNT; b++) buffers.push_back(new Buffer<double>());}
  if (residuals.empty()) {for (int f = 0; f < mt_COUNT; f++) residuals.push_back(new BaseBuffer<double>());}
  if (computed.empty()) computed.assign(ct_COUNT, false);
  if (ssr.empty()) ssr.assign(mt_COUNT, -1);

  outcome.resize(no_rows); 
  zstat.assign(variable_data.get_outcome_data(), no_rows);
}

void RegressionComputation::init_model(int no_vars, bool prune_only) {
  GeneMatrix& inv_corrs = variable_data.get_corrs(); 
  eff_size = inv_corrs.size(false);
  
  init_params(no_vars);

  inv_corrs.product(predictors[0]);
  inv_corrs.product(zstat.data(), outcome.data());

  (*buffers[bt_XtX])(0,0) = product(predictors[0], no_rows);
  (*buffers[bt_XtZ])(0) = product(outcome.data(), no_rows);

  int core_size = param_count(prune_only ? pt_PruneCore : pt_Core);    
  for (int v = 1; v < core_size; v++) compute_core(v);
  
  invert_conditional(prune_only);
  set_degrees();
}

void RegressionComputation::init_params(int no_vars) {
  type_range.assign(pt_COUNT, pair<int,int>(0,0));  
  for (int i = 0; i < variable_types.size(); i++) {
    pair<int,int>& range = type_range[variable_types[i]];
    if (range.second > 0) {
      if (range.first + range.second != i) error("order of model parameters is invalid");  
      range.second++;
    } else range = pair<int,int>(i,1);
  }
  type_range[pt_Covariate] = range_sum(type_range[pt_InternalCovariate], range_sum(type_range[pt_HiddenCovariate], type_range[pt_UsedCovariate]));
  type_range[pt_PruneCore] = range_sum(type_range[pt_Intercept], type_range[pt_InternalCovariate]);
  type_range[pt_Core] = range_sum(type_range[pt_Intercept], type_range[pt_Covariate]);

  update_nvar(no_vars, true); 
  for (int i = 0; i < pt_COUNT; i++) {if (type_range[i].second < 0) error("order of model parameters is invalid");}

  buffers[bt_XtX]->set_size(no_param, no_param, true);
  buffers[bt_XtZ]->set_size(no_param, 1, true);
  buffers[bt_Coefficients]->set_size(no_param,1); 
  buffers[bt_StdErrors]->set_size(no_param,1);

  predictors.set_size(no_rows, no_param - 1);
}

void RegressionComputation::set_degrees() {
  df.resize(mt_COUNT);
  df[mt_Intercept] = eff_size - 1;
  df[mt_BaseModel] = eff_size - eff_core;
  df[mt_BaseModel] = eff_size - eff_core - (param_count(pt_Core) - param_count(pt_PruneCore));
  df[mt_FullModel] = df[mt_BaseModel] - param_count(pt_MainVariable);  
  converter.set_param("df", df[mt_FullModel]); 
} 

int RegressionComputation::update_nvar(int no_vars, bool is_init) {
  if (is_init || no_vars != param_count(pt_MainVariable)) {
    int core_size = param_count(pt_Core), diff = no_vars - param_count(pt_MainVariable);
    no_param = core_size + no_vars;
    variables.resize(no_param, 0);
    variable_types.resize(no_param, pt_MainVariable);
    type_range[pt_MainVariable] = pair<int,int>(core_size, no_vars);
    
    type_range[pt_Any] = pair<int,int>(0,no_param);
    type_range[pt_Variable] = range_sum(type_range[pt_UsedCovariate], type_range[pt_MainVariable]);  
    
    if (!is_init) {
      Buffer<double> swap_buffer(no_param, no_param, true);
      swap_buffer.matrix_sub(0, core_size) = buffers[bt_XtX]->matrix_sub(0, core_size);
      buffers[bt_XtX]->swap(swap_buffer);
      
      swap_buffer.set_size(no_param, 1, true);
      swap_buffer.matrix_sub(0, core_size, 0, 1) = buffers[bt_XtZ]->matrix_sub(0, core_size, 0, 1);
      buffers[bt_XtZ]->swap(swap_buffer);      

      buffers[bt_Coefficients]->set_size(no_param,1); 
      buffers[bt_StdErrors]->set_size(no_param,1);
      predictors.resize(no_param - 1);
 
      set_degrees();   
    }
    return diff;
  } else return 0;
}

pair<int,int> RegressionComputation::range_sum(pair<int,int> A, pair<int,int> B) {
  if (A.second < 0 || B.second < 0) return pair<int,int>(0,-1);
  if (A.first > B.first) return range_sum(B, A);

  if (A.second == 0) return B;
  else if (B.second == 0) return A;
  else if (A.first+A.second == B.first) return pair<int,int>(A.first, A.second+B.second);
  else return pair<int,int>(0,-1);
}

void RegressionComputation::invert_conditional(bool prune_only) {
  int core_size = param_count(pt_Core);
  bool non_pruned = !prune_only && core_size > param_count(pt_PruneCore);
  
  try {
    EigenInverter<double> eigen(*buffers[bt_XtX], param_count(pt_PruneCore));
    eigen.set_cutoff(PC_PRUNE_VARS, true);
    eigen.inverse(*buffers[non_pruned ? bt_XtXinvUse : bt_XtXinvStore]);
    eff_core = eigen.get_rank();
  } catch (const MathException& me) {error("unable to invert base model covariance matrix", me.get_msg());}

  if (non_pruned) {
    try {
      Buffer<double> cond(core_size, core_size);
      cond.matrix() = buffers[bt_XtX]->matrix_sub(0, core_size);
      collinear_core = inverter.invert(&cond, buffers[bt_XtXinvUse], buffers[bt_XtXinvStore]);
    } catch (const MathException& me) {error("could not invert design matrix of conditioned-on variables", me.get_msg());}
  }
}

void RegressionComputation::residualize() {
  BaseBuffer<double>& resid = resid_values(mt_PruneModel);
  double rsq = 1 - residual_variance(mt_PruneModel, false) / residual_variance(mt_Intercept, false);
  rsq_multiplier = 1 - max(min(rsq, 1.0), 0.0);
  zstat.assign(resid); resid.clear();

  int from = param_offset(pt_InternalCovariate), to = from + param_count(pt_InternalCovariate);
  variables.erase(variables.begin() + from, variables.begin() + to);
  variable_types.erase(variable_types.begin() + from, variable_types.begin() + to);
}

void RegressionComputation::run() {
  int offset = param_offset(pt_MainVariable);
  for (int i = 0; i < param_count(pt_MainVariable); i++) compute_core(offset+i);
  coefficients(); std_errors();
}

int RegressionComputation::convert_index(int& param_index, ParameterType type) {
  if (param_index >= 0) {
    pair<int,int>& range = type_range[type];
    if (param_index < range.second) param_index += range.first;
    else param_index = -1;
  }  
  if (param_index < 0) failure("trying to access unknown model parameter");     
  return param_index;
}


void RegressionComputation::set_variables(set<int>* var_ids) {
  error_msg.clear(); check_init(true); computed.assign(ct_COUNT, false);
  int index = param_offset(pt_MainVariable), total = param_count(pt_MainVariable);
  if (!var_ids || var_ids->size() != total) error("number of specified variables does not match model size");

  for (set<int>::iterator it = var_ids->begin(); it != var_ids->end(); ++it, ++index) {  
    variables[index] = &variable_data.get_variable(*it, true);
    if (!variables[index]->is_valid()) failure("trying to analyse variable that is unknown or unavailable");
  }
}

void RegressionComputation::set_variable(int var_id) {
  static set<int> wrapper; wrapper.clear();
  if (var_id >= 0) wrapper.insert(var_id);
  set_variables(&wrapper);
}

void RegressionComputation::set_conditional(set<int>& ids, ParameterType type) {check_init(false);
  if (type <= pt_Covariate || type >= pt_Variable) error("incorrect parameter type for conditional variable");
  for (set<int>::iterator it = ids.begin(); it != ids.end(); ++it) {
    VariableInfo& info = variable_data.get_variable(*it, true); 
    if (info.is_valid()) {
      variables.push_back(&info); 
      variable_types.push_back(type);
    }
    else error("conditioned variable is not available for analysis");
  }
}

void RegressionComputation::compute_core(int offset) {
  GeneMatrix& inv_corrs = variable_data.get_corrs(); Buffer<double>& xtx = *buffers[bt_XtX];
  VariableInfo& var = *variables[offset];
  if (var.is_set()) {
    set<int>& data = *var.as_set()->get_data();
    if (offset < no_param-1) inv_corrs.product(&data, predictors[offset]);
    (*buffers[bt_XtZ])(offset) = product(data, outcome.data());
    xtx(offset,offset) = inv_corrs.square(&data);
    for (int c = 0; c < offset; c++) xtx(c,offset) = xtx(offset,c) = product(data, predictors[c]);        
  } else if (var.is_covar()) {
    float* data = var.as_covar()->get_data(); 
    if (offset < no_param-1) inv_corrs.product(data, predictors[offset]);
    (*buffers[bt_XtZ])(offset) = product(data, outcome.data(), no_rows);
    xtx(offset,offset) = inv_corrs.square(data);    
    for (int c = 0; c < offset; c++) xtx(c,offset) = xtx(offset,c) = product(data, predictors[c], no_rows);            
  } else error("trying to compute with undefined variable (core)"); 
}

Buffer<double>& RegressionComputation::fit_coefficients(ModelType type, Buffer<double>& target) {
  if (type == mt_Intercept) {
    target.set_size(1,1);
    target(0) = xtz()(0) / xtx()(0);  
  } else if (type == mt_PruneModel || type == mt_BaseModel) {
    int size = param_count(type == mt_PruneModel ? pt_PruneCore : pt_Core);
    target.set_size(size, 1);
    target.matrix().noalias() = buffers[bt_XtXinvStore]->matrix() * xtz().matrix_vec(0,0,size);
  } else {
    if (!computed[ct_Coefficients]) {
      MathUtils::matrix_prod(xtx_inv(), xtz(), *buffers[bt_Coefficients]);
      computed[ct_Coefficients] = true;
    }  
    if (&target != buffers[bt_Coefficients]) target.assign(*buffers[bt_Coefficients]);
  }
  return target;
}

void RegressionComputation::fit_model(ModelType type) {
  BaseBuffer<double>& fit = *residuals[type];
  fit.assign(zstat.data(), no_rows);          

  Buffer<double>& coeff = fit_coefficients(type); double offset = coeff(0);
  for (double *curr = fit.data(), *end = curr+no_rows; curr < end; curr++) *curr -= offset;
  for (int j = 1; j < coeff.length(); j++) variable_sum(fit, *variables[j], -coeff(j));

  ssr[type] = MathUtils::in_product(fit.data(), outcome.data(), no_rows);
}

void RegressionComputation::fit_stderrors() {
  if (!computed[ct_StdErrors]) {
    Buffer<double> &inv = xtx_inv(), &se = *buffers[bt_StdErrors];    
    double resvar = residual_variance();

    for (int i = 0; i < no_param; i++) {
      double var = resvar * inv(i,i);
      se(i) = var > 0 ? sqrt(var) : 0;
    }  

    computed[ct_StdErrors] = true;
  }
}

Buffer<double>& RegressionComputation::xtx_inv() {
  if (!computed[ct_Inverse]) { try {
    collinear_mode = inverter.invert(buffers[bt_XtX], buffers[bt_XtXinvStore], buffers[bt_XtXinvUse]);
    if (collinear_core != ct_None) collinear_mode = (collinear_mode == ct_None) ? ct_Covariates : ct_All;
    computed[ct_Inverse] = true;
  } catch (const MathException& me) {failure("could not invert variable design matrix", me.get_msg());} }
  return *buffers[bt_XtXinvUse];
}

Buffer<double>& RegressionComputation::coefficients() {
  if (!computed[ct_Coefficients]) fit_coefficients(mt_FullModel);
  return *buffers[bt_Coefficients];  
}

Buffer<double>& RegressionComputation::std_errors() {
  if (!computed[ct_StdErrors]) fit_stderrors();
  return *buffers[bt_StdErrors];  
}

double RegressionComputation::residual_variance(ModelType type, bool scale) {
  resid_values(type);
  return scale ? ssr[type] / df[type] : ssr[type];
}

BaseBuffer<double>& RegressionComputation::resid_values(ModelType type) {
  if (type == mt_FullModel) {
    if (!computed[ct_Fitted]) {
      fit_model(mt_FullModel);
      computed[ct_Fitted] = true;
    }  
  } else if (residuals[type]->is_empty()) fit_model(type);
  return *residuals[type];
}

double RegressionComputation::scaling(int index) {
  VariableInfo& var = *variables[index];
  if (index > 0) {
    if (var.is_set()) return var.as_set()->compute_sd(no_rows);
    else if (var.is_covar()) return var.as_covar()->get_orig_sd();
  }
  return 1;
}

vector<VariableInfo*> RegressionComputation::get_variables(ParameterType type) {
  pair<int,int>& range = type_range[type];
  return vector<VariableInfo*>(variables.begin()+range.first, variables.begin()+range.first+range.second);
}

int RegressionComputation::get_ngenes(int index, ParameterType type) {
  if (convert_index(index, type) > 0) {
    VariableInfo& var = *variables[index];
    if (var.is_set()) return var.as_set()->get_ngenes();
    else if (var.is_covar()) return var.as_covar()->get_obs_genes();
    else return 0;
  } else return no_rows;
} 

double RegressionComputation::get_coeff(int index, ParameterType type) {
  if (convert_index(index, type) > 0 && variables[index]->is_covar()) {
    double sd = scaling(index);
    return sd > 0 ? coefficients()(index)/sd : 0;    
  } else return coefficients()(index);
}

double RegressionComputation::get_coeff_std(int index, ParameterType type) {
  if (convert_index(index, type) > 0 && variables[index]->is_set()) {
    double sd = scaling(index);
    return  sd > 0 ? coefficients()(index)*sd : 0;    
  } else return coefficients()(index);
}

double RegressionComputation::get_se(int index, ParameterType type) {
  if (convert_index(index, type) > 0 && variables[index]->is_covar()) {
    double sd = scaling(index);
    return  sd > 0 ? std_errors()(index)/sd : 0;    
  } else return std_errors()(index);
}

double RegressionComputation::get_pval(int index, short direction, ParameterType type) {
  convert_index(index, type);
  double se = std_errors()(index);  
  if (se <= 0) failure("standard error is zero or negative");
  double stat = coefficients()(index)/se, pval; 
  
  try {
    if (direction != 0) {
      if (direction < 0) stat *= -1;
      pval = converter.convert(stat);
    } else pval = converter.convert(abs(stat))*2;
  } catch (const StatException& se) {failure("unable to compute p-value (" + se.get_msg() + ")");}
  return pval;
}

bool RegressionComputation::get_collinear(int index, ParameterType type) {
  if (collinear_mode == ct_All) return true;
  else if (collinear_mode == ct_None) return false;
  else {
    convert_index(index, type);
    if (collinear_mode == ct_Covariates && has_type(index, pt_Core, pt_Any)) return true;
    else if (collinear_mode == ct_Variables && has_type(index, pt_MainVariable, pt_Any)) return true;
    else return false;
  } 
}

bool RegressionComputation::has_type(int index, ParameterType type, ParameterType select_type) {
	pair<int,int>& range = type_range[type]; convert_index(index, select_type);
  return (index >= range.first && index < range.first+range.second);
}

double RegressionComputation::get_rsq(bool base_rsq, bool rescale) {
  double rsq = 1 - residual_variance(base_rsq ? mt_BaseModel : mt_FullModel, false) / residual_variance(mt_Intercept, false);
  return rescale ? rsq * rsq_multiplier : rsq;
}

void RegressionComputation::load_fitted(Buffer<double>& target, ModelType type) {
  BaseBuffer<double>& resid = resid_values(type);
  target.set_size(no_rows, 2);
  target.assign_col(resid.data(), 1);
  MathUtils::vec_diff(zstat.data(), resid.data(), target[0], no_rows);
}
