"""True frequency-domain FSM: frame FFT cross-products, temporal cyclic search,
then exact polynomial phase shifts and complex spectral averaging. No x**2.
Prototype checks, not MATLAB execution nor a 1000-trial performance claim.
"""
import numpy as np
from scipy.signal import resample_poly
from scipy.optimize import minimize
from pathlib import Path
import json
FS=800000.; L=128; H=64; F_RANGE=(50000.,150000.); A_RANGE=(-1000.,1000.)
def setup(x):
 n=len(x); starts=np.arange(0,n-L+1,H); ix=starts[:,None]+np.arange(L)
 w=np.hanning(L); frames=(x[ix]*w).T
 tc=(ix-(n-1)/2)/FS
 return frames,tc.T,(starts-(n-1)/2)/FS, (n-1)/2/FS

def coarse_search(frames,tc,ts):
 m=frames.shape[1]; df=FS/L; X=np.fft.fft(frames,axis=0)
 fidx=np.arange(-8,9); nt=2**int(np.ceil(np.log2(2*m)))
 delta=np.fft.fftfreq(nt,H/FS)
 rates=np.arange(-1000,1000.1,100.)
 best=(-np.inf,None,None)
 candidate_list=[]
 for k in range(int(np.floor(F_RANGE[0]/df)),int(np.ceil(F_RANGE[1]/df))+1):
  f0=k*df
  p=np.mean(X[(fidx+k)%L]*np.conj(X[(fidx-k)%L]),axis=0)
  p=p*np.exp(-4j*np.pi*f0*ts)
  values=abs(np.fft.fft(p[None,:]*np.exp(-2j*np.pi*rates[:,None]*ts[None,:]**2),nt,axis=1))**2
  freq=f0+delta/2
  valid=(delta>=-df)&(delta<df)&(freq>=F_RANGE[0])&(freq<=F_RANGE[1])
  ids=np.flatnonzero(valid); ir,j=np.unravel_index(np.argmax(values[:,valid]),values[:,valid].shape)
  v=values[ir,ids[j]]
  candidate_list.append((v,freq[ids[j]],rates[ir]))
  if v>best[0]:best=(v,freq[ids[j]],rates[ir])
 return best[1],best[2],sorted(candidate_list,reverse=True)[:5]

def spectrum(frames,tc,f,a):
 phase=2*np.pi*f*tc+np.pi*a*tc**2
 xp=np.fft.fft(frames*np.exp(-1j*phase),axis=0)
 xm=np.fft.fft(frames*np.exp(1j*phase),axis=0)
 return np.mean(xp*np.conj(xm),axis=1)/(FS*np.sum(np.hanning(L)**2))

def exact_score(frames,tc,f,a):
 s=spectrum(frames,tc,f,a)
 return abs(np.mean(s[np.arange(-8,9)%L]))**2

def estimate(x):
 frames,tc,ts,ref=setup(x); f0,a0,candidates=coarse_search(frames,tc,ts); T=len(x)/FS
 # Confirm several fixed-prior branch candidates with the EXACT FSM score.
 exact=[(exact_score(frames,tc,f,a),f,a) for _,f,a in candidates]
 _,f0,a0=max(exact)
 scale=exact_score(frames,tc,f0,a0)
 def fun(p):
  f=f0+p[0]/T; a=a0+p[1]/T**2
  if not F_RANGE[0]<=f<=F_RANGE[1] or not A_RANGE[0]<=a<=A_RANGE[1]:return 1e6
  return -exact_score(frames,tc,f,a)/max(scale,1e-30)
 opt=minimize(fun,[0,0],method='Nelder-Mead',options={'xatol':1e-5,'fatol':1e-8,'maxiter':180})
 return f0,a0,f0+opt.x[0]/T,a0+opt.x[1]/T**2,bool(opt.success)

def gold_code():
 def seq(taps):
  s=[1]*5
  for i in range(31):s.append(sum(s[i+k] for k in taps)%2)
  assert s[31:36]==s[:5]
  return np.array(s[:31])
 return 2*np.logical_xor(seq([0,2]),np.roll(seq([0,2,3,4]),1)).astype(float)-1

if __name__=='__main__':
 rng=np.random.default_rng(20261014); rawfs=24.8e6; rb=1e5; T=.1; N=round(rawfs*T)
 t=(np.arange(N)-(N-1)/2)/rawfs; pn=gold_code()[(np.arange(N)//8)%31]
 cases=[(50000.,-1000.,None),(149999.,999.,None),(100173.7,413.9,None),
        (76543.2,-873.1,-7),(123456.7,731.4,-7),(143210.3,-531.2,-7)]
 results=[]
 for f,a,eb in cases:
  data=np.repeat(rng.choice([-1.,1.],round(rb*T)),248)
  clean=data*pn*np.cos(2*np.pi*f*t+np.pi*a*t*t+rng.uniform(-np.pi,np.pi))
  x=clean.copy()
  if eb is not None:x+=rng.normal(size=N)*np.sqrt(np.mean(clean*clean)*rawfs/(2*rb*10**(eb/10)))
  # Known, exact PN. Real low-pass resampling: preserve bands through 400 kHz.
  y=resample_poly(x*pn,1,31)
  fc,ac,ff,af,ok=estimate(y)
  # resampled sample midpoint differs from raw by <= 1/fs_out; correct truth.
  ref_shift=(len(y)-1)/(2*FS)-(N-1)/(2*rawfs)
  row={'f_true_hz':f+a*ref_shift,'rate_true_hz_s':a,'eb_n0_db':eb,'coarse_f_hz':fc,'coarse_rate_hz_s':ac,'fine_f_hz':ff,'fine_rate_hz_s':af,'fine_f_error_hz':ff-(f+a*ref_shift),'fine_rate_error_hz_s':af-a,'converged':ok}
  results.append(row);print(json.dumps(row),flush=True)
 Path(__file__).with_name('fsm_frequency_checks.json').write_text(json.dumps(results,indent=2))
