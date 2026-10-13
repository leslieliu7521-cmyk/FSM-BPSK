"""Actual real passband -> sampled PN wipeoff -> resample -> LPF -> estimator.
Random known 31-chip code, not the unavailable original mseq/Gold implementation.
"""
from pathlib import Path
source=Path(__file__).with_name('known_code_validate.py').read_text().split('errors=[]')[0]
# Reuse identical estimator and defaults, no equivalent-baseband trial loop.
exec(source)
from scipy.signal import resample_poly
rawfs=24.8e6; fc=1e5; rawN=round(rawfs*T)
tr=(np.arange(rawN)-(rawN-1)/2)/rawfs
errors=[]
for k in range(args.trials):
 fd=rng.uniform(-45000,45000); rate=rng.uniform(-900,900); phase=rng.uniform(-np.pi,np.pi)
 pnchips=rng.choice([-1.,1.],31); pn=pnchips[(np.arange(rawN)//8)%31]
 bits=rng.choice([-1.,1.],int(np.ceil(rawN/248)))
 data=bits[np.arange(rawN)//248]
 clean=pn*data*np.cos(2*np.pi*(fc+fd)*tr+np.pi*rate*tr**2+phase)
 variance=np.mean(clean*clean)*rawfs/(2*rb*10**(-7/10))
 raw=clean+rng.normal(size=rawN)*np.sqrt(variance)
 # Time coordinates match raw midpoint; PN excludes data.
 base=resample_poly(2*raw*pn*np.exp(-2j*np.pi*fc*tr),1,62)
 y=lfilter(h,[1],base)[2*delay:]
 # Recenter local estimator times after trimming and known FIR group delay.
 tq=(np.arange(len(y))+delay)*62/rawfs-(rawN-1)/(2*rawfs)
 center=(tq[0]+tq[-1])/2
 oldt=t;oldw=w;oldbank=bank
 t=tq-center;w=np.hanning(len(y));bank=np.exp(-2j*np.pi*rates[:,None]*t[None,:]**2)
 fh,ah=estimate(y)
 errors.append([fh-(fd+rate*center),ah-rate])
 t=oldt;w=oldw;bank=oldbank
 print('raw trial',k+1,flush=True)
e=np.array(errors); ok=(abs(e[:,0])<=20)&(abs(e[:,1])<=200)
r={'trials':len(e),'duration_s':T,'eb_n0_db':-7,'carrier_success':float((abs(e[:,0])<=20).mean()),'rate_success':float((abs(e[:,1])<=200).mean()),'joint_success':float(ok.mean()),'frequency_rmse_hz':float(np.sqrt(np.mean(e[:,0]**2))),'rate_rmse_hz_per_s':float(np.sqrt(np.mean(e[:,1]**2))),'errors':e.tolist(),'pipeline':'real passband, sampled known PN wipeoff, downconversion, 62x resampling, lowpass, square chirp estimator','limitations':'exact code timing, random 31-chip code not original Gold; unknown data; linear Doppler only'}
Path(__file__).with_name('known_code_raw_results.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
