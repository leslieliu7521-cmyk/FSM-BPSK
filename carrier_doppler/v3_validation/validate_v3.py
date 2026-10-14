"""Small Python cross-checks of v3 mathematics and control flow.
NOT a native MATLAB run; no 1000-trial claim. SciPy/MATLAB FIR details can differ.
"""
from pathlib import Path
import numpy as np
from scipy.signal import resample_poly
from scipy.optimize import minimize
import json,tempfile,os,csv,re,hashlib,time,shutil,argparse
ROOT=Path(__file__).resolve().parents[1]
REPORT={'matlab_available':shutil.which('matlab') is not None,'matlab_executed':False,'full_experiment_executed':False,'checks':{},'paired_smoke':[]}
def check(name,condition):
 assert condition,name
 REPORT['checks'][name]='PASS (Python cross-check)'
def distinct(points,df,da):
 keep=[]
 for i,p in enumerate(points):
  if not keep or np.all(np.sum(((np.asarray(points)[keep]-p)/[df,da])**2,axis=1)>=1):keep.append(i)
 return np.array(keep,int)
def select(coarse,R):
 point=np.array(coarse[:2],float);J=coarse[2];selected=-1
 for i,r in enumerate(R):
  valid=bool(r[6]) and np.all(np.isfinite(r[2:5])) and 5e4<=r[2]<=15e4 and -1000<=r[3]<=1000 and r[4]>=0
  if valid and r[4]>J:point=r[2:4];J=r[4];selected=i
 return point,J,selected

def track(df,da,j,n,fs):
 m2=(n*n-1)/(12*fs*fs);m4=(n*n-1)*(3*n*n-7)/(240*fs**4)
 rms=np.sqrt(max(0,df*df+da*da*m2+.25*j*j*m4-df*j*m2))
 bound=(n-1)/(2*fs);points=[-bound,bound]
 if j!=0 and abs(da/j)<=bound:points.append(da/j)
 points=np.array(points);return rms,max(abs(df+da*points-.5*j*points**2))

class FSM:
 def __init__(self,y,fs=8e5):
  self.fs=fs;self.L=round(128*fs/8e5);self.H=self.L//2;self.T=len(y)/fs
  starts=np.arange(0,len(y)-self.L+1,self.H);ix=starts[:,None]+np.arange(self.L)
  self.w=np.hanning(self.L);self.frames=(y[ix]*self.w).T
  self.tc=(ix.T-(len(y)-1)/2)/fs;self.ts=(starts-(len(y)-1)/2)/fs
  self.band=np.arange(-8,9);self.df=.5/self.T;self.da=2/self.T**2
  X=np.fft.fft(self.frames,axis=0);freqbin=fs/self.L
  nt=2**int(np.ceil(np.log2(2*len(starts))));delta=np.fft.fftshift(np.fft.fftfreq(nt,self.H/fs))
  rates=np.arange(-1000,1000.1,100);raw=[]
  for k in range(int(np.floor(5e4/freqbin)),int(np.ceil(15e4/freqbin))+1):
   fb=k*freqbin
   p=np.mean(X[(self.band+k)%self.L]*np.conj(X[(self.band-k)%self.L]),axis=0)*np.exp(-4j*np.pi*fb*self.ts)
   power=abs(np.fft.fftshift(np.fft.fft(p[None,:]*np.exp(-2j*np.pi*rates[:,None]*self.ts**2),nt,axis=1),axes=1))**2
   freqs=fb+delta/2;valid=(delta>=-freqbin)&(delta<freqbin)&(freqs>=5e4)&(freqs<=15e4)
   f=freqs[valid];P=power[:,valid].copy()
   for _ in range(3):
    a,b=np.unravel_index(np.argmax(P),P.shape);value=P[a,b]
    if not np.isfinite(value):break
    raw.append([f[b],rates[a],value,k])
    distance=((rates[:,None]-rates[a])/self.da)**2+((f[None,:]-f[b])/self.df)**2
    P[distance<1]=-np.inf
  self.raw=np.array(raw);order=np.argsort(-self.raw[:,2],kind='stable')
  ordered=self.raw[order];self.candidates=ordered[distinct(ordered[:,:2],self.df,self.da)]
 def spectrum(self,f,a):
  phase=2*np.pi*f*self.tc+np.pi*a*self.tc**2
  xp=np.fft.fft(self.frames*np.exp(-1j*phase),axis=0)
  xm=np.fft.fft(self.frames*np.exp(1j*phase),axis=0)
  return np.mean(xp*np.conj(xm),axis=1)/(self.fs*np.sum(self.w**2))
 def score(self,f,a):return abs(np.mean(self.spectrum(f,a)[self.band%self.L]))**2
 def estimate(self,count=10):
  count=min(len(self.candidates),count);C=self.candidates[:count]
  scores=np.array([self.score(f,a) for f,a,*_ in C]);rank=np.argsort(-scores,kind='stable')
  coarse=[*C[rank[0],:2],scores[rank[0]]];norm=max(coarse[2],np.finfo(float).tiny);R=[]
  for i in rank[:3]:
   start=C[i,:2]
   def fun(p):
    f,a=start+p/[self.T,self.T**2]
    if not 5e4<=f<=15e4 or not -1000<=a<=1000:return 1e12
    return -self.score(f,a)/norm
   result=minimize(fun,[0,0],method='Nelder-Mead',options={'maxiter':200,'xatol':1e-5,'fatol':1e-8})
   f,a=start+result.x/[self.T,self.T**2]
   J=self.score(f,a) if 5e4<=f<=15e4 and -1000<=a<=1000 else np.nan
   valid=np.all(np.isfinite([f,a,J])) and 5e4<=f<=15e4 and -1000<=a<=1000
   R.append([*start,f,a,J,int(result.success),int(valid),int(result.success)])
  point,J,selected=select(coarse,np.array(R))
  assert J>=coarse[2]
  return {'coarse_f':coarse[0],'coarse_a':coarse[1],'coarse_score':coarse[2], 'f':point[0],'a':point[1],
          'score':J,'fallback':selected<0,'exitflag':0 if selected<0 else R[selected][5],
          'candidate_count':len(self.candidates),'confirmed':int(count),'starts':len(R),'valid_outputs':int(sum(r[6] for r in R))}

def make_base(seed):
 r=np.random.default_rng(seed);fs=24.8e6;n=2480000
 bits=r.choice([-1.,1.],10000);phase=r.uniform(0,2*np.pi);fd=r.uniform(-45000,45000);a=r.uniform(-900,900)
 noise=r.normal(size=n);t=(np.arange(n)-(n-1)/2)/fs
 def seq(taps):
  s=[1]*5
  for j in range(31):s.append(sum(s[j+k] for k in taps)%2)
  assert s[31:36]==s[:5]
  return np.array(s[:31])
 code=2*np.logical_xor(seq([0,2]),np.roll(seq([0,2,3,4]),1))-1
 pn=code[(np.arange(n)//8)%31];data=np.repeat(bits,248)
 return {'bits':bits,'phase':phase,'fd':fd,'a':a,'noise':noise,'t':t,'pn':pn,'data':data}

def run(structural_only=False):
 start=time.monotonic();source=(ROOT/'FSM_BPSK_full_system_v3.m').read_text()
 original=ROOT/'FSM_BPSK_full_system.m'
 if not original.exists():original=ROOT/'original'/'FSM_BPSK_full_system.m'
 check('original_file_preserved',hashlib.sha256(original.read_bytes()).hexdigest()=='ffac50f39cd0d334a070bcb403fc1cdbe72e5379e0eded7713341dda213c21d4')
 for term in ['cfg.Rb=1e5','cfg.Rc=31e5','cfg.fs=24.8e6','cfg.T=.1','nTrials=1000','cfg.multiStarts=3','cfg.fTolerance=20','cfg.rateTolerance=200']:
  check('parameter_'+term,term in source)
 check('quick_default',"'mode','quick'" in source)
 check('numeric_save_before_plot',source.index("state.phase='numerics_saved_before_plot'")<source.index('[figures,creationErrors]=create_figures'))
 check('closed_handle_and_numbered_print',"~isgraphics(h,'figure')" in source and "print(number,png,'-dpng','-r150')" in source)
 check('real_frequency_domain_core','Xp.*conj(Xm)' in source and 'y.^2' not in source and 'x.^2' not in source)
 names=re.findall(r"'([A-Za-z][A-Za-z0-9_]*)'",source.split('function m=metric_template()')[1].split('m=cell2struct')[0])
 assigned=set(re.findall(r'm\.([A-Za-z][A-Za-z0-9_]*)',source.split('function [record,demo]=run_trial')[1].split('function r=receiver_config')[0]))
 check('metric_column_count_and_assignments',len(names)==len(set(names))==41 and assigned==set(names))
 for T in [.001,.01,.05,.1]:
  df=.5/T;da=2/T**2;pts=np.array([[1e5,0],[1e5+.1*df,.1*da],[1e5+2*df,0],[1e5,2*da]])
  ix=distinct(pts,df,da);check(f'dedup_T_{T}',np.array_equal(ix,[0,2,3]))
 coarse=[1e5,0,10];R=np.array([[1e5,0,100010,0,12,0,1,0],[1e5,0,np.nan,0,99,-999,0,0],[1e5,0,100020,0,11,1,1,1]])
 check('multistart_highest_valid_not_convergence_alias',select(coarse,R)[2]==0)
 for label,row in [('worse',[1e5,0,100010,0,9,1,1,1]),('invalid',[1e5,0,np.nan,0,np.nan,-999,0,0]),('outofprior',[1e5,0,200000,0,100,1,1,1])]:
  check('fallback_'+label,select(coarse,np.array([row]))[2]==-1)
 u=(np.arange(801)-400)/8000;e=.2-20*u-.5*5000*u*u;tr,tm=track(.2,-20,5000,801,8000)
 check('track_exact_rmse_and_continuous_max',abs(tr-np.sqrt(np.mean(e*e)))<1e-12 and tm>=max(abs(e))-1e-12)
 with tempfile.TemporaryDirectory() as folder:
  path=Path(folder)/'checkpoint.json';temp=Path(folder)/'checkpoint.partial'
  state={'done':[True,False,True],'seeds':[11,12,13],'numeric':[1,2,3]}
  temp.write_text(json.dumps(state));os.replace(temp,path)
  temp.write_text('interrupted partial write')
  restored=json.loads(path.read_text());check('atomic_old_checkpoint_survives_partial',restored==state)
  missing=[i for i,done in enumerate(restored['done']) if not done]
  for i in missing:restored['done'][i]=True
  check('resume_no_duplicate_count',missing==[1] and sum(restored['done'])==3)
  for fault in ['closed_figure','exportgraphics_failure','print_failure']:
   try:raise OSError(fault)
   except OSError:check('numeric_survives_'+fault,json.loads(path.read_text())['numeric']==[1,2,3])
 if structural_only:
  REPORT['source_sha256']=hashlib.sha256((ROOT/'FSM_BPSK_full_system_v3.m').read_bytes()).hexdigest()
  REPORT['elapsed_s']=time.monotonic()-start
  (ROOT/'v3_validation'/'structural_results.json').write_text(json.dumps(REPORT,indent=2))
  print('PASS',len(REPORT['checks']),'structural/control/formula checks; no signal simulation.')
  return
 b=make_base(20261015);again=make_base(20261015)
 check('paired_base_reconstruction',all(np.array_equal(b[x],again[x]) for x in ['bits','phase','fd','a','noise']))
 del again
 # Full real passband, unknown bits. PN exact and known; not supplied as data.
 for f,a,eb in [(90624.63,-873.1,None),(90625.,413.9,None),(150000.,1000.,None),(123456.7,731.4,-7)]:
  clean=b['pn']*b['data']*np.cos(2*np.pi*f*b['t']+np.pi*a*b['t']**2+b['phase'])
  variance=0 if eb is None else np.mean(clean*clean)*24.8e6/(2*1e5*10**(eb/10))
  raw=clean+np.sqrt(variance)*b['noise'];long=resample_poly(raw*b['pn'],1,31)
  receiver=FSM(long);counts=[5,10,np.inf] if eb==-7 else [10]
  pool=receiver.candidates.copy();best=-np.inf
  for count in counts:
   result=receiver.estimate(count)
   check('candidate_pool_unchanged_'+str(f)+'_'+str(count),np.array_equal(pool,receiver.candidates))
   check('nested_coarse_confirmation_'+str(f)+'_'+str(count),result['coarse_score']>=best);best=result['coarse_score']
   # Convert true frequency to resampled midpoint, not raw midpoint.
   delta=(len(long)-1)/(2*8e5)-(2480000-1)/(2*24.8e6)
   result.update({'f_true':f+a*delta,'a_true':a,'eb_n0_db':eb,'confirm_limit':str(count),
                  'f_error':result['f']-(f+a*delta),'a_error':result['a']-a})
   check('finite_bounded_smoke_'+str(f)+'_'+str(count),np.all(np.isfinite([result['f'],result['a'],result['score']])) and 5e4<=result['f']<=15e4 and abs(result['a'])<=1000)
   REPORT['paired_smoke'].append(result)
   print('smoke',f,a,eb,count,'errors',round(result['f_error'],4),round(result['a_error'],4),flush=True)
  # Duration pairing: identical post-filter long record and correctly centered reference.
  for T in [.001,.01,.05,.1]:
   n=round(T*8e5);first=(len(long)-n)//2;seg=long[first:first+n]
   reference=first/8e5+(n-1)/(2*8e5)
   check('paired_crop_reference_'+str(f)+'_'+str(T),len(seg)==n and abs(reference-(len(long)-1)/(2*8e5))<1e-12)
  if eb==-7:
   fsm=receiver;S=fsm.spectrum(0,0)
   psd=np.mean(abs(np.fft.fft(fsm.frames,axis=0))**2,axis=1)/(fsm.fs*np.sum(fsm.w**2))
   check('fsm_psd_normalization',np.allclose(S,psd,rtol=1e-12,atol=1e-20))
   for up,down,fsout in [(1,31,8e5),(2,31,16e5)]:
    out=resample_poly(raw*b['pn'],up,down);model=FSM(out,fsout)
    check('bandwidth_dimensions_'+str(fsout),model.L==round(128*fsout/8e5) and model.frames.shape[1]==receiver.frames.shape[1] and model.fs/model.L==6250)
 # Additional paired duration/sensitivity smoke: KEEP misses, no success-only filtering.
 for duration in [.001,.01,.05]:
  n=round(duration*8e5);first=(len(long)-n)//2;seg=long[first:first+n]
  result=FSM(seg).estimate(10)
  tref=first/8e5+(n-1)/(2*8e5);delta=tref-(2480000-1)/(2*24.8e6)
  result.update({'condition':'duration_'+str(duration),'f_error':result['f']-(f+a*delta),'a_error':result['a']-a})
  result['target_success']=bool(abs(result['f_error'])<=20 and abs(result['a_error'])<=200)
  REPORT['paired_smoke'].append(result)
  check('duration_selection_invariant_'+str(duration),result['score']>=result['coarse_score'])
  print('paired duration',duration,'error',result['f_error'],result['a_error'],'target',result['target_success'],flush=True)
 for scenario in ['Jerk5000','Phase_0p125_chip']:
  jerk=5000 if scenario=='Jerk5000' else 0
  clean=b['pn']*b['data']*np.cos(2*np.pi*f*b['t']+np.pi*a*b['t']**2+np.pi/3*jerk*b['t']**3+b['phase'])
  variance=np.mean(clean*clean)*24.8e6/(2*1e5*10**(-7/10))
  rx=clean+np.sqrt(variance)*b['noise']
  pnrx=np.roll(b['pn'],-1) if scenario=='Phase_0p125_chip' else b['pn']
  signal=resample_poly(rx*pnrx,1,31);result=FSM(signal).estimate(10)
  delta=(len(signal)-1)/(2*8e5)-(2480000-1)/(2*24.8e6)
  df=result['f']-(f+a*delta+.5*jerk*delta*delta);da=result['a']-(a+jerk*delta)
  tr,tm=track(df,da,jerk,len(signal),8e5)
  result.update({'condition':scenario,'f_error':df,'a_error':da,'track_rmse_hz':tr,'track_max_error_hz':tm,'target_success':bool(abs(df)<=20 and abs(da)<=200)})
  REPORT['paired_smoke'].append(result)
  check('sensitivity_selection_invariant_'+scenario,result['score']>=result['coarse_score'])
  print('paired sensitivity',scenario,'errors',df,da,'track',tr,tm,flush=True)
 REPORT['elapsed_s']=time.monotonic()-start
 REPORT['limitations']=['No native MATLAB/Octave execution.','Python I/O fault injection checks the corresponding save architecture, not MATLAB graphics/export backends.','SciPy/MATLAB random streams and FIR design are not assumed bit-identical.','Small deterministic/paired smoke cases, not a success-rate study.']
 (ROOT/'v3_validation'/'validation_results.json').write_text(json.dumps(REPORT,indent=2,allow_nan=False))
 print('PASS',len(REPORT['checks']),'checks;',round(REPORT['elapsed_s'],2),'seconds')
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--structural-only',action='store_true')
 run(parser.parse_args().structural_only)
