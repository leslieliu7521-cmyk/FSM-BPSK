"""Equivalent complex-baseband validation AFTER ideal synchronized known-code wipeoff.
Not a validation of unknown code/bit timing acquisition or original mseq.
Complex noise variance = fs/(Rb*EbN0), amplitude=1, bits unknown to estimator.
"""
import numpy as np
from scipy.signal import firwin, lfilter
from scipy.optimize import minimize
from pathlib import Path
import json, argparse
parser=argparse.ArgumentParser(); parser.add_argument('--trials',type=int,default=100); parser.add_argument('--duration',type=float,default=.1); args=parser.parse_args()
rng=np.random.default_rng(20261014); fs=400000.; rb=100000.; T=args.duration
N=round(fs*T); t=(np.arange(N)-(N-1)/2)/fs
rates=np.arange(-1000,1000+1,max(50,1/T**2))
if rates[-1]!=1000: rates=np.append(rates,1000)
bank=np.exp(-2j*np.pi*rates[:,None]*t[None,:]**2)
h=firwin(129,90000,fs=fs); delay=(len(h)-1)//2
nfft=2**int(np.ceil(np.log2(2*N))); f=np.fft.fftfreq(nfft,1/fs)/2
sel=np.flatnonzero(abs(f)<=50000); w=np.hanning(N)
def estimate(y):
 q=y*y*w
 score=abs(np.fft.fft(q[None,:]*bank,nfft,axis=1)[:,sel])**2
 ir,j=np.unravel_index(score.argmax(),score.shape); f0=f[sel[j]]; a0=rates[ir]
 def fun(p):
  fh=f0+p[0]/T; ah=a0+p[1]/T**2
  if abs(fh)>50000 or abs(ah)>1000:return 1e99
  return -abs(np.sum(q*np.exp(-1j*(4*np.pi*fh*t+2*np.pi*ah*t*t))))**2
 res=minimize(fun,[0,0],method='Nelder-Mead',options={'xatol':1e-6,'maxiter':250})
 return f0+res.x[0]/T,a0+res.x[1]/T**2
errors=[]
for k in range(args.trials):
 fd=rng.uniform(-45000,45000); rate=rng.uniform(-900,900); phase=rng.uniform(-np.pi,np.pi)
 # Guard symbols avoid artificial zero-extension boundary transients.
 extN=N+2*delay; te=(np.arange(extN)-delay-(N-1)/2)/fs
 bits=rng.choice([-1.,1.],int(np.ceil(extN/4)))
 clean=np.repeat(bits,4)[:extN]*np.exp(1j*(2*np.pi*fd*te+np.pi*rate*te*te+phase))
 variance=fs/(rb*10**(-7/10))
 noise=(rng.normal(size=extN)+1j*rng.normal(size=extN))*np.sqrt(variance/2)
 y=lfilter(h,[1],clean+noise)[2*delay:2*delay+N]
 fh,ah=estimate(y); errors.append([fh-fd,ah-rate])
 if (k+1)%50==0: print(f'{T}s completed {k+1}',flush=True)
e=np.array(errors); ok=(abs(e[:,0])<=20)&(abs(e[:,1])<=200); n=len(e);p=float(ok.mean());z=1.96
c=(p+z*z/(2*n))/(1+z*z/n);d=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
r={'trials':n,'duration_s':T,'eb_n0_db':-7,'carrier_tolerance_hz':20,'rate_tolerance_hz_per_s':200,'fd_test_range_hz':[-45000,45000],'rate_test_range_hz_per_s':[-900,900],'carrier_success':float((abs(e[:,0])<=20).mean()),'rate_success':float((abs(e[:,1])<=200).mean()),'joint_success':p,'joint_wilson_95':[c-d,c+d],'frequency_rmse_hz':float(np.sqrt(np.mean(e[:,0]**2))),'rate_rmse_hz_per_s':float(np.sqrt(np.mean(e[:,1]**2))),'assumptions':['ideal known-code and bit-boundary synchronization','complex baseband equivalence, not original real-passband pipeline','unknown independent BPSK bits','AWGN; constant amplitude; linear Doppler','no code-rate Doppler, multipath, oscillator drift, or jerk']}
Path(__file__).with_name(f'known_code_{T:g}s_results.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
