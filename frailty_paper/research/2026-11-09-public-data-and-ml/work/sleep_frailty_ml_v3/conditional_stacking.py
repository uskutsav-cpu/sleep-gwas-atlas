"""Small trainable combination of a DNA expert and permitted AVI scores.

Quantile references, imputation and coefficients fit only the supplied training
rows. Prediction for a row is independent of other rows in its input batch.
"""
import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import QuantileTransformer
from sklearn.linear_model import LogisticRegression
def inputs(base,avi,include_dna):
    a=np.asarray(avi,dtype='float64');b=np.asarray(base,dtype='float64');assert len(a)==len(b);x=np.c_[b,a[:,0]] if include_dna else a[:,:1];return x,a[:,2:3]
def fit_models(base,avi,y,study,training,include_dna):
    x,flag=inputs(base,avi,include_dna);models={}
    for s in range(2):
        ix=training[study[training]==s];assert len(np.unique(y[ix]))==2;imputer=SimpleImputer(strategy='median',keep_empty_features=True);z=imputer.fit_transform(x[ix]);qt=QuantileTransformer(n_quantiles=min(256,len(ix)),output_distribution='normal',subsample=None,random_state=42);q=qt.fit_transform(z);classifier=LogisticRegression(C=1.,max_iter=2000,random_state=42);classifier.fit(np.c_[q,flag[ix]],y[ix]);models[s]={'imputer':imputer,'quantile':qt,'classifier':classifier,'training_n':len(ix),'training_calls':int(y[ix].sum()),'include_DNA_expert':include_dna}
    return models
def predict_models(models,base,avi,study,include_dna):
    x,flag=inputs(base,avi,include_dna);pred=np.empty(len(base),dtype='float64')
    for s in range(2):
        ix=np.flatnonzero(study==s)
        if not len(ix):continue
        m=models[s];assert m['include_DNA_expert']==include_dna;z=m['quantile'].transform(m['imputer'].transform(x[ix]));pred[ix]=m['classifier'].predict_proba(np.c_[z,flag[ix]])[:,1]
    assert np.isfinite(pred).all();return pred
