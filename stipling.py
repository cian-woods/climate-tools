import numpy as np

def stipling(Nb,xx=None,yy=None,x=None,y=None,thresh=0.8):
	# Stipling for > thresh of models with same bias sign
	# Input array with models on first dimension, also axes of data
	n           = 1.*Nb.shape[0]
	thresh      = 2*thresh - 1
	Nb_sign     = np.sign(Nb)
	Nb_sign_zer = np.abs(np.abs(Nb_sign)-1)
	Nb_pos      = np.abs(np.ma.masked_where(Nb_sign<=0,Nb_sign).sum(axis=0)).data
	Nb_neg      = np.abs(np.ma.masked_where(Nb_sign>=0,Nb_sign).sum(axis=0)).data
	Nb_zer      = np.abs(np.ma.masked_where(Nb_sign_zer==0,Nb_sign_zer).sum(axis=0)).data
	Nb_diff     = np.abs(Nb_pos-Nb_neg)/n
	# If axes are 1d not 2d
        if (xx,yy) == (None,None): xx,yy = np.meshgrid(x,y)
	stipx,stipy = np.ma.masked_where(Nb_diff<thresh,xx),np.ma.masked_where(Nb_diff<thresh,yy)
	return stipx,stipy
