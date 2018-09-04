from UnPickle import *
from KDE import *
from scipy import stats

import glob,sys
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.colors as colors

def select_markers(n):
        kernel = ['o','+','v','s','^','x']
        marks  = []
        for i in range(n):
                c = i%len(kernel)
                marks.append(kernel[c])
        return marks

def stagger_levels(lev0,dlev,n,N):
	# n = number of levs
	# N = number of points
	levs   = [lev0 + i*dlev for i in range(n)]
	points = []
	for i in range(N):
		c = i%n
		points.append(levs[c])
	return points
	

# Input parameters
N          = int(sys.argv[1])
Gradient   = str(sys.argv[2])=='True'
Models_not = ['CMCC-CM','GFDL-ESM2M']
Models     = [g[33:] for g in glob.glob('/mnt/climstorage/cian/historical/*') if g[33:] not in Models_not]
# Sea-ice from Johan model (March 1980-2018 time series)
ci_y = range(1980,2018+1,1)
ci = [0.4293079022311042, 0.4591951242091214, 0.4645827117891332,\
 0.36947709681221536, 0.3807608434449665, 0.3551705207349827,\
 0.40107250460590915, 0.4319571162032606, 0.42006771226291006,\
 0.41116043062494445, 0.3883145164095926, 0.37796831613172, 0.3501907727510132,\
 0.3862038046350324, 0.4070822611339657, 0.3503428942059941,\
 0.3778855535402496, 0.41153637168938234, 0.4698939765739258,\
 0.4273113137753377, 0.3654033529557372, 0.4183206170351418,\
 0.41079619013495167, 0.45060716711721316, 0.41373399409613326,\
 0.3798740628205584, 0.3537437721915546, 0.32327573603759285,\
 0.3510563606059502, 0.3520895069904999, 0.3799395464228133,\
 0.3470848312343721, 0.2745806226379127, 0.3667244449140885,\
 0.29552216115608254, 0.32954975494906014, 0.27082991519288346,\
 0.3168441468066702, 0.34950260619227635]
m3e = 1e2*np.array([10*np.polyfit(ci_y[ci_y.index(i):ci_y.index(i)+N],ci[ci_y.index(i):ci_y.index(i)+N],1)[0] for i in range(1980,2018-N+1+1,1)])

if Gradient:
	grad,reg_end = '.grad','.70-80N_340-15E'
else:
	grad,reg_end = '',''
year_end = 2018 - N + 1
# CMIP data
M2h,M1_end,M2_end,Model_end = [],[],[],[]
M2r = []
for Model in Models:
	try:
		y,m1,m2 = unpick('/mnt/climstorage/cian/trendfiles/%s/%s.rcp85%s.1850-%s.NDJF.NDJF.%syears.75-83N_20-80E_65-75N_60-100E%s.p' % (N,Model,grad,2018-N+1,N,reg_end))
		print Model, len([i for i in m2 if i != -999]),len(y)
		m1_end  = m1[y.index(year_end)]
		m2_end  = m2[y.index(year_end)]
		m2      = [i/100. for i in m2 if i!=-999]
		M2h     = M2h + m2
		if (m1_end != -999) and (m2_end != -999):
			M1_end.append(m1_end)
			M2_end.append(m2_end/100)
			Model_end.append(Model)
	except:
		pass
        try:
		y,m1,m2 = unpick('/mnt/climstorage/cian/trendfiles/%s/%s.rcp85%s.2019-%s.NDJF.NDJF.%syears.75-83N_20-80E_65-75N_60-100E%s.p' % (N,Model,grad,2100-N+1,N,reg_end))
		print Model, len([i for i in m2 if i != -999]),len(y)
		m2      = [i/100. for i in m2 if i!=-999]
		M2r     = M2r + m2
        except:
                pass

M2h,M2r     = np.array(M2h),np.array(M2r)
M2h_m,M2r_m = M2h.mean(),M2r.mean()
M2h_s,M2r_s = M2h.std(),M2r.std()
print M2h_m,M2r_m
print M2h_s,M2r_s
print Model_end
print len(Model_end)
print len(M2h),len(M2r)

# ERAInt data
ye,m1e,m2e = unpick('/mnt/climstorage/cian/trendfiles/%s/%s.historical%s.1980-%s.NDJF.NDJF.%syears.75-83N_20-80E_65-75N_60-100E%s.p' % (N,'ERAInt',grad,2018-N+1,N,reg_end))
m2e = [i/100. for i in m2e]
print m2e

print stats.percentileofscore(M2h,m2e)

# PDF edges
emin,emax = -3.0,3.0
edges     = np.arange(emin,emax+0.10,0.10)
x_grid    = np.arange(emin,emax+0.01,0.01)

# Best fit kernel bandwith
grid   = GridSearchCV(KernelDensity(),{'bandwidth': np.linspace(0.1, 0.4, 30)},cv=20) # 20-fold cross-validation
grid.fit(M2h[:, None])
kde_h = grid.best_estimator_
grid.fit(M2r[:, None])
kde_r = grid.best_estimator_
print kde_h.bandwidth
print kde_r.bandwidth

# Make PDF anf get percentiles
pdf_h = np.exp(kde_h.score_samples(x_grid[:, None]))
cdf_h = np.cumsum(pdf_h)
x50_h   = x_grid[np.argmin((cdf_h-50)**2)]
x95_h   = x_grid[np.argmin((cdf_h-95)**2)]
x99_h   = x_grid[np.argmin((cdf_h-99)**2)]
pdf_r = np.exp(kde_r.score_samples(x_grid[:, None]))
cdf_r = np.cumsum(pdf_r)
x50_r   = x_grid[np.argmin((cdf_r-50)**2)]
x95_r   = x_grid[np.argmin((cdf_r-95)**2)]
x99_r   = x_grid[np.argmin((cdf_r-99)**2)]


scattermap = pl.get_cmap('nipy_spectral')
norm       = colors.Normalize(vmin=0, vmax=1)
markers    = select_markers(len(Model_end))
ixs        = np.linspace(0.03,0.97,len(Model_end))

# Plot
fig, ax1 = pl.subplots(num=1,figsize=(8,8))
#ax2     = ax1.twinx()
ax1.plot(x_grid, pdf_h, 'grey', linestyle='solid' , linewidth=3.2)#, alpha=0.5)
ax1.plot(x_grid, pdf_r, 'grey', linestyle='dashed', linewidth=3.2)#, alpha=0.5)
ax1.hist(M2h,edges, fc='gray', histtype='stepfilled', alpha=0.3, normed=True,align='mid')
ymin1,ymax1 = ax1.get_ylim()
y_points    = stagger_levels(lev0=ymax1,dlev=0.012*ymax1,n=1,N=len(M2_end))
#Model_end   = [i[1] for i in sorted(zip(M2_end,Model_end))]
#M2_end      = sorted(M2_end)
for i in range(len(M2_end)):
	#ax1.plot(M2_end[i],y_points[i],marker=markers[i],markerfacecolor="None",markeredgecolor=scattermap(norm(ixs[i])),markersize=6.5,mew=2,label=Model_end[i],linestyle="None")
	ax1.plot(M2_end[i],y_points[i],marker='|',markerfacecolor='grey',markeredgecolor='grey',markersize=10,mew=1,label=Model_end[i],linestyle="None")
ax1.plot(np.mean(M2_end),np.mean(y_points),marker='|',markeredgecolor="grey",markerfacecolor='grey',markersize=20,mew=4,label='Multi-model mean',linestyle="None")
ax1.plot(np.mean(m2e),np.mean(y_points),marker='*',markeredgecolor="None",markerfacecolor='k',markersize=12,mew=2,label='ERA-Interim',linestyle="None")
pl.xlim(edges[0],edges[-1])
xmin1,xmax1 = ax1.get_xlim()
ymin1,ymax1 = ax1.get_ylim()
ymax1 = 0.9
ax1.plot([x95_h,x95_h],[ymin1,1.1*ymax1],'k',linewidth=1.1,alpha=0.3)
ax1.plot([x99_h,x99_h],[ymin1,1.1*ymax1],'k',linewidth=1.1,alpha=0.3)
#ax1.plot([x50_h,x50_h],[ymin1,1.1*ymax1],'b--',linewidth=1.1,alpha=0.3)
#ax1.plot([x50_r,x50_r],[ymin1,1.1*ymax1],'r--',linewidth=1.1,alpha=0.3)
#ax1.plot([M2h_m,M2h_m],[ymin1,pdf_h[np.argmin((x_grid-M2h_m)**2)]],'b--',linewidth=1.1,alpha=0.3)
#ax1.plot([M2r_m,M2r_m],[ymin1,pdf_r[np.argmin((x_grid-M2r_m)**2)]],'r--',linewidth=1.1,alpha=0.3)
ax1.plot(M2h_m,pdf_h[np.argmin((x_grid-M2h_m)**2)],marker='o',markeredgecolor='grey',markerfacecolor='grey',markersize=6.5)
ax1.plot(M2r_m,pdf_r[np.argmin((x_grid-M2r_m)**2)],marker='o',markeredgecolor='grey',markerfacecolor='grey',markersize=6.5)
ax1.set_ylim(ymin1,ymax1)
fs = 14
if grad == '':
	ax1.set_xlabel('Sea-level pressure trend [hPa decade$^{-1}$]',fontsize=fs)
else:
	ax1.set_xlabel('Sea-level pressure gradient trend [hPa decade$^{-1}$]',fontsize=fs)
ax1.set_ylabel('Frequency',fontsize=fs)
pl.title('%s year windows' % (N))
ax1.set_aspect(abs((xmax1-xmin1)/(ymax1-ymin1))*1)
ax1.set_xticks([-3,-2,-1,0,1,2,3])
ax1.set_xticklabels(['-3','-2','-1','0','1','2','3'],fontsize=fs)
ax1.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
ax1.set_yticklabels(['0.0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9'],fontsize=fs)
pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/dists/%syears%s.pdf' % (N,grad),format='pdf')
# New figure for legend
#fig_leg = pl.figure(num=2,figsize=(3,6))
#ax_leg = fig_leg.add_subplot(111)
#ax_leg.legend(*ax1.get_legend_handles_labels(), loc='center',frameon=True,prop={'size':8},ncol=1,columnspacing=1,handletextpad=0.2,numpoints=1)
#ax_leg.axis('off')
#pl.savefig('/mnt/climstorage/cian/scripts/figs/BSO/dists/%syears.legend.pdf' % (N),format='pdf')

pl.show()
