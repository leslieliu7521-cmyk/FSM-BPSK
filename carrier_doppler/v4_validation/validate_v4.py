"""v4 structural audit and Python paired diagnostic cross-check; NEVER MATLAB results.
Run: PYTHONPATH=/tmp/fsm_v4_parser python .../validate_v4.py [--structural-only]
The optional tree-sitter parser is only syntax parsing, not MATLAB name resolution.
"""
from pathlib import Path
import sys, importlib.util, hashlib, json, re, shutil, time, argparse
import numpy as np
from scipy.signal import resample_poly
from scipy.optimize import minimize
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('v3checks',ROOT/'v3_validation/validate_v3.py')
v3=importlib.util.module_from_spec(spec);spec.loader.exec_module(v3)
SOURCE=(ROOT/'FSM_BPSK_full_system_v4.m').read_text()
REPORT={'native_matlab_available':shutil.which('matlab'),'native_octave_available':shutil.which('octave'),
        'native_matlab_executed':False,'full_experiment_executed':False,'checks':{},'paired_diagnostics':[],
        'limitations':['Python FFT/random/FIR implementation is not MATLAB execution.',
        'Small deterministic cases do not establish success rates or general causes.']}
def check(name,yes):
 if not yes:raise AssertionError(name)
 REPORT['checks'][name]='PASS'
def function_text(source,name):
 starts=list(re.finditer(r'^function\b.*',source,re.M))
 for i,start in enumerate(starts):
  if re.search(r'\b'+name+r'\(',start.group()):
   text=source[start.start():starts[i+1].start() if i+1<len(starts) else len(source)]
   ends=list(re.finditer(r'^end\s*$',text,re.M))
   return text[:ends[-1].end()].rstrip()+'\n'
 raise KeyError(name)
def audit():
 old=(ROOT/'FSM_BPSK_full_system_v3.m').read_bytes()
 check('v3_preserved',hashlib.sha256(old).hexdigest()=='dfbfb1e4fa03a0a77c56418e53a014a80872eeba79deb6bd6a81d3e3ec676114')
 # Compare the scientifically critical functions byte for byte.
 for f in ['make_base','make_gold','mseq','fsm_spectrum','fsm_score','branch_candidates','distinct_indices','choose_final','refine_objective','track_error','metric_template','error_stats','wilson','atomic_state_save','atomic_bundle_save','export_figures']:
  check('preserved_'+f,function_text(SOURCE,f)==function_text(old.decode(),f))
 for term in ['cfg.Rb=1e5','cfg.Rc=31e5','cfg.fs=24.8e6','cfg.fc=1e5','cfg.T=.1','nTrials=1000','cfg.fRange=[50e3 150e3]','cfg.rateRange=[-1000 1000]','cfg.fTolerance=20','cfg.rateTolerance=200']:
  check('parameter_'+term,term in SOURCE)
 check('diagnostic_default',"opts=struct('mode','diagnostic'" in SOURCE)
 check('native_trials_30_100_configurable','nTrials=opts.diagnosticTrials' in SOURCE)
 check('noise_free_first','eb=[eb(isinf(eb)),eb(isfinite(eb))]' in SOURCE)
 check('full_retains_original_conditions','conditions=make_conditions(cfg,nTrials,opts.mode)' in SOURCE)
 check('diagnostic_modes_do_not_feedback_truth',SOURCE.index('search=fsm_receiver(y,receiverCfg)')<SOURCE.index('diagnosis=evaluate_stages('))
 receiver=function_text(SOURCE,'fsm_receiver')
 for forbidden in ['trueF','trueA','b.bits','scanF','rxRatePPM']:
  check('no_search_truth_'+forbidden,forbidden not in receiver)
 white=function_text(SOURCE,'receiver_config')
 for forbidden in ["'scanF'","'testCase'","'bits'","'rxRatePPM'","'fdTest'"]:
  check('whitelist_excludes_'+forbidden,forbidden not in white)
 check('direct_fixed_prior',"cfg.directRange=[90000 91250]" in SOURCE and 'c.fRange=c.directRange' in receiver)
 check('all_generated_not_global','全局最优' in function_text(SOURCE,'direct_exact_grid') or 'global optimality' in SOURCE)
 check('numerics_before_figures',SOURCE.index("atomic_state_save(state,fullfile(outputDir,'results.mat'))",SOURCE.index('state.summary=make_summary'))<SOURCE.index('create_diagnostic_figures(state)'))
 check('pair_signature_assertion','isequaln(a.receiverSignature,b.receiverSignature)' in SOURCE)
 check('phase_quantization_documented','PhaseResolutionChip' in SOURCE and 'ContinuousFractionalDelay' in SOURCE)
 check('no_adc_claim','ADCClockMismatch' in SOURCE and 'LocalPNIndexRateMismatch_NotADCClock' in SOURCE)
 check('separate_frequency_stage',"nearG=abs(G(:,1)-trueF)<=d.FrequencyNeighborhood_Hz" in SOURCE)
 check('dedup_identity_mapping','dedupRepresentativeGeneratedIndex' in SOURCE)
 check('all_branch_surfaces','info.ApproxPower=power(:,valid)' in SOURCE)
 check('rate_profile_no_rejection','RateProfileAtSelected' in SOURCE and 'RateProfileRelativeContrast' in SOURCE)
 script=SOURCE.split('\nfunction ',1)[0]
 script=re.sub(r"'(?:''|[^'\n])*'","''",script);script=re.sub(r'%[^\n]*','',script)
 variables=set(re.findall(r'\b([A-Za-z]\w*)\s*=(?!=)',script))
 for outputs in re.findall(r'\[([^]\n]+)\]\s*=(?!=)',script):variables.update(re.findall(r'\b[A-Za-z]\w*\b',outputs))
 functions=re.findall(r'^function\s+(?:(?:\[[^]\n]*\]|\w+)\s*=\s*)?(\w+)\s*\(',SOURCE,re.M)
 check('no_duplicate_local_function_names',len(functions)==len(set(functions)))
 check('no_script_function_name_collision',not variables.intersection(functions))
 metrics=function_text(SOURCE,'metric_template')
 names=re.findall(r"'([A-Za-z][A-Za-z0-9_]*)'",metrics.split('m=cell2struct')[0])
 assigned=set(re.findall(r'm\.([A-Za-z][A-Za-z0-9_]*)',function_text(SOURCE,'run_trial')))
 check('metric_columns_41',len(names)==len(set(names))==41 and set(names)==assigned)
 try:
  from tree_sitter import Language,Parser
  import tree_sitter_matlab
  parser=Parser(Language(tree_sitter_matlab.language()));tree=parser.parse(SOURCE.encode())
  check('tree_sitter_syntax_only',not tree.root_node.has_error)
  REPORT['parser']='tree-sitter-matlab; not MATLAB semantic validation'
 except ImportError:REPORT['parser']='not installed; parser check skipped'
 for T in [.001,.01,.05,.1]:
  n=round(T*8e5);first=(80000-n)//2;ref=first/8e5+(n-1)/(2*8e5)
  check('center_reference_'+str(T),abs(ref-79999/(2*8e5))<1e-12)
  df=max(20,.5/T);da=max(200,2/T**2)
  check('frequency_independent_of_rate_'+str(T),abs(100)<df if T<=.001 else df<=50)
 for phase,samples in [(0,0),(.0625,0),(.125,1),(.1875,1),(.5,4)]:
  idx=np.arange(24800);ci=np.floor(idx/8+phase).astype(int)
  check('pn_phase_quantization_'+str(phase),np.array_equal(ci,np.floor((idx+samples)/8).astype(int)))
 for count in [3,30,100]:
  p=1.;z=1.96;lo=(p+z*z/(2*count)-z*np.sqrt(z*z/(4*count*count)))/(1+z*z/count)
  check('wilson_all_success_not_certainty_'+str(count),0<lo<1)

class DiagnosticFSM(v3.FSM):
 def __init__(self,y,branch_peaks=3,dedup_factor=1.):
  # First build frame data with the unmodified mathematical v3 cross-check.
  super().__init__(y)
  self.df*=dedup_factor;self.da*=dedup_factor
  if branch_peaks==3 and dedup_factor==1:return
  X=np.fft.fft(self.frames,axis=0);freqbin=self.fs/self.L
  nt=2**int(np.ceil(np.log2(2*len(self.ts))));delta=np.fft.fftshift(np.fft.fftfreq(nt,self.H/self.fs))
  rates=np.arange(-1000,1000.1,100);raw=[]
  for k in range(8,25):
   fb=k*freqbin
   p=np.mean(X[(self.band+k)%self.L]*np.conj(X[(self.band-k)%self.L]),axis=0)*np.exp(-4j*np.pi*fb*self.ts)
   P=abs(np.fft.fftshift(np.fft.fft(p[None,:]*np.exp(-2j*np.pi*rates[:,None]*self.ts**2),nt,axis=1),axes=1))**2
   freqs=fb+delta/2;valid=(delta>=-freqbin)&(delta<freqbin)&(freqs>=5e4)&(freqs<=15e4)
   f=freqs[valid];P=P[:,valid].copy()
   for j in range(branch_peaks):
    ai,fi=np.unravel_index(np.argmax(P),P.shape);value=P[ai,fi]
    if not np.isfinite(value):break
    raw.append([f[fi],rates[ai],value,k])
    dist=((rates[:,None]-rates[ai])/self.da)**2+((f[None,:]-f[fi])/self.df)**2;P[dist<1]=-np.inf
  self.raw=np.array(raw);ordered=self.raw[np.argsort(-self.raw[:,2],kind='stable')]
  self.candidates=ordered[v3.distinct(ordered[:,:2],self.df,self.da)]
 def estimate(self,count=10):
  lo,hi=getattr(self,'search_range',(5e4,15e4))
  count=min(len(self.candidates),count);C=self.candidates[:count]
  scores=np.array([self.score(f,a) for f,a,*_ in C]);rank=np.argsort(-scores,kind='stable')
  coarse=[*C[rank[0],:2],scores[rank[0]]];norm=max(coarse[2],np.finfo(float).tiny);R=[]
  for i in rank[:3]:
   start=C[i,:2]
   def fun(p):
    f,a=start+p/[self.T,self.T**2]
    if not lo<=f<=hi or not -1000<=a<=1000:return 1e12
    return -self.score(f,a)/norm
   out=minimize(fun,[0,0],method='Nelder-Mead',options={'maxiter':200,'xatol':1e-5,'fatol':1e-8})
   f,a=start+out.x/[self.T,self.T**2]
   valid=np.all(np.isfinite([f,a])) and lo<=f<=hi and -1000<=a<=1000
   J=self.score(f,a) if valid else np.nan
   valid=bool(valid and np.isfinite(J) and J>=0)
   R.append([*start,f,a,J,int(out.success),int(valid),int(out.success)])
  point,J,selected=v3.select(coarse,np.array(R))
  assert J>=coarse[2] and lo<=point[0]<=hi
  return {'coarse_f':coarse[0],'coarse_a':coarse[1],'coarse_score':coarse[2],
          'f':point[0],'a':point[1],'score':J,'fallback':selected<0,
          'refinements':R,'candidate_count':len(self.candidates)}
 def diagnose(self,result,true_f,true_a,count):
  near=lambda x:abs(x[:,0]-true_f)<=max(20,.5/self.T)
  G=self.raw;C=self.candidates;conf=C[:min(len(C),count)]
  gs=np.array([self.score(f,a) for f,a,*_ in G]);nearG=near(G)
  out={'generated_frequency_present':bool(any(nearG)),'dedup_frequency_present':bool(any(near(C))),
       'confirmed_frequency_present':bool(any(near(conf))),'near_frequency_scale_hz':max(20,.5/self.T),
       'final_frequency_near':bool(abs(result['f']-true_f)<=max(20,.5/self.T)),
       'truth_point_score':self.score(true_f,true_a),'selected_score':result['score'],
       'best_near_generated_exact_score':float(max(gs[nearG],default=np.nan)),
       'best_other_generated_exact_score':float(max(gs[~nearG],default=np.nan))}
  stage='None'
  for name,flag in [('Generation',out['generated_frequency_present']),('Deduplication',out['dedup_frequency_present']),('ConfirmationCutoff',out['confirmed_frequency_present'])]:
   if not flag:stage=name;break
  if stage=='None' and abs(result['coarse_f']-true_f)>out['near_frequency_scale_hz']:stage='ExactMetricRanking'
  elif stage=='None' and not out['final_frequency_near']:stage='RefinementSelection'
  out['first_missing_frequency_stage']=stage
  return out

def smoke():
 b=v3.make_base(20261015);fs=24.8e6
 # Reproduce the historical 1ms miss, then a fixed diagnostic frequency near a branch boundary.
 cases=[('Historical',123456.7,731.4),('FixedOffGrid',90643.7,731.4)]
 for case,f,a in cases:
  clean=b['pn']*b['data']*np.cos(2*np.pi*f*b['t']+np.pi*a*b['t']**2+b['phase'])
  for eb in [None,-7]:
   var=0 if eb is None else np.mean(clean*clean)*fs/(2*1e5*10**(eb/10))
   raw=clean+np.sqrt(var)*b['noise']
   for scene,phase,ppm in [('IdealSync',0,0),('Phase_0p125_chip',.125,0),('Rate_5_ppm',0,5)]:
    ci=np.floor(np.arange(len(raw))*(1+ppm*1e-6)/8+phase).astype(int)
    code=b['pn'][np.arange(31)*8];pn=code[ci%31]
    long=resample_poly(raw*pn,1,31)
    times=[.001,.01]+([.05,.1] if case=='Historical' and eb==-7 and scene=='IdealSync' else [])
    for T in times:
     n=round(T*8e5);first=(len(long)-n)//2;y=long[first:first+n]
     tref=first/8e5+(n-1)/(2*8e5);delta=tref-(len(raw)-1)/(2*fs);truth=f+a*delta
     baseline=DiagnosticFSM(y)
     for method in ['V3','ConfirmAll','Branch9','DedupHalf','DirectExact']:
      if T>.01 and method!='V3':continue
      if method=='DirectExact':
       # Deliberately fixed limited reference even when Historical lies outside this scope.
       model=DiagnosticFSM(y);ff=np.unique(np.r_[np.arange(90000,91250+.0001,.25/T),91250]);aa=np.arange(-1000,1000.1,500)
       points=np.array([[fp,ap,model.score(fp,ap),0] for fp in ff for ap in aa]);model.raw=points.copy()
       model.candidates=points[np.argsort(-points[:,2],kind='stable')];count=len(points)
       model.search_range=(90000,91250)
      else:
       model=baseline if method in ['V3','ConfirmAll'] else DiagnosticFSM(y,branch_peaks=9 if method=='Branch9' else 3,dedup_factor=.5 if method=='DedupHalf' else 1.)
       count=len(model.candidates) if method=='ConfirmAll' else 10
      clock=time.monotonic();result=model.estimate(count);elapsed=time.monotonic()-clock
      diag=model.diagnose(result,truth,a,count)
      if method=='DirectExact' and not 90000<=truth<=91250:diag['first_missing_frequency_stage']='OutsideFixedReferenceScope'
      check('objective_nondecrease_'+str(len(REPORT['paired_diagnostics'])),result['score']>=result['coarse_score'])
      tr,tm=v3.track(result['f']-truth,result['a']-a,0,n,8e5)
      row={'case':case,'eb_n0_db':eb,'scene':scene,'duration_s':T,'method':method,'f_error_hz':result['f']-truth,
           'rate_error_hz_per_s':result['a']-a,'track_rmse_hz':tr,'track_max_hz':tm,'coarse_f_error_hz':result['coarse_f']-truth,
           'fallback':result['fallback'],'refinements':result['refinements'],'selection_elapsed_s':elapsed,'diagnosis':diag}
      REPORT['paired_diagnostics'].append(row)
      print(case,eb,scene,T,method,round(row['f_error_hz'],3),diag['first_missing_frequency_stage'],flush=True)

def main():
 args=argparse.ArgumentParser();args.add_argument('--structural-only',action='store_true');o=args.parse_args()
 start=time.monotonic();audit()
 if not o.structural_only:smoke()
 REPORT['elapsed_s']=time.monotonic()-start;REPORT['source_sha256']=hashlib.sha256(SOURCE.encode()).hexdigest()
 name='structural_results.json' if o.structural_only else 'python_diagnostic_results.json'
 # NaN denotes unavailable diagnostics and is converted to JSON null explicitly.
 def clean(x):
  if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
  if isinstance(x,list):return [clean(v) for v in x]
  if isinstance(x,float) and not np.isfinite(x):return None
  return x
 (ROOT/'v4_validation'/name).write_text(json.dumps(clean(REPORT),ensure_ascii=False,indent=2)+'\n')
 print('PASS',len(REPORT['checks']),'checks;',len(REPORT['paired_diagnostics']),'Python diagnostic rows; no MATLAB run.')
if __name__=='__main__':main()
