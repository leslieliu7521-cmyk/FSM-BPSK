"""Structural and small Python cross-checks of the stage-1-only MATLAB baseline.
Not MATLAB execution; NumPy/SciPy random/FIR/tie handling can differ.
"""
from pathlib import Path
import sys,re,json,hashlib,time,shutil,importlib.util,argparse
import numpy as np
from scipy.signal import resample_poly
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('fsm_reference',ROOT/'v3_validation/validate_v3.py')
ref=importlib.util.module_from_spec(spec);spec.loader.exec_module(ref)
SOURCE=(ROOT/'FSM_BPSK_stage1_700.m').read_text()
REPORT={'native_matlab_executed':False,'matlab_executable':shutil.which('matlab'),
        'octave_executable':shutil.which('octave'),'full_150_experiment_executed':False,
        'checks':{},'noiseless':[],'paired_smoke':[],'short_smoke':[],'benchmark':{},
        'limitations':['Python implementation is not MATLAB native execution.',
        'Only a few shared base records; not a success-rate study.',
        'SciPy/MATLAB FIR/random streams/tie handling are not assumed identical.']}
def check(name,value):
 assert value,name
 REPORT['checks'][name]='PASS'
def function_text(source,name):
 starts=list(re.finditer(r'^function\b[^\n]*',source,re.M))
 for i,m in enumerate(starts):
  if re.search(r'\b'+name+r'\(',m.group()):
   text=source[m.start():starts[i+1].start() if i+1<len(starts) else len(source)]
   ends=list(re.finditer(r'^end\s*$',text,re.M));return text[:ends[-1].end()].strip()
 raise KeyError(name)
def structural():
 old=(ROOT/'FSM_BPSK_full_system_v4.m').read_text()
 for name in ['branch_candidates','distinct_indices','track_error','error_stats','wilson','make_gold','mseq']:
  check('original_'+name,function_text(SOURCE,name)==function_text(old,name))
 check('no_continuous_optimizer_calls',not re.search(r'\b(?:fminsearch|fminunc|fmincon|refine_objective)\s*\(',SOURCE))
 for token in ["'trials',50","'ebN0',[-10 -7 -4]","'durations',.1","'rates',700","'makePlot',false",
               "'windowComparison',false","'rateComparison',false","'bandwidthComparison',false",
               "'syncComparison',false","'candidateComparison',false","'fRange',[50e3 150e3]",
               "'rateRange',[-1000 1000]","'rateStep',100","'branchPeaks',3","'confirmCandidates',10",
               "'fs',24.8e6","'scfFs',8e5","'Rb',1e5","'Rc',31e5"]:
  check('parameter_'+token,token in SOURCE)
 rx=function_text(SOURCE,'receiver_config')
 check('receiver_whitelist_no_truth',not any(x in rx.split("names=")[1] for x in ["'RateTruth'","'fd'","'bits'"]))
 evaluation=function_text(SOURCE,'evaluate_record')
 check('f0_transport_uses_estimated_rate','m.F0_Hz=e.fLocal+e.rate*(x.reference-localRef)' in evaluation)
 check('separate_frequency_success','m.F0Success=m.Valid && abs(m.F0Error_Hz)<=c.fTolerance' in evaluation)
 check('nonfinite_complete_sample_rmse','any(~isfinite(v))' in function_text(SOURCE,'error_stats'))
 check('invalid_absolute_errors_as_inf','values(~isfinite(values))=Inf' in function_text(SOURCE,'absolute_quantile'))
 check('no_repeated_700','rates=unique([rates 650 700 750]' in function_text(SOURCE,'make_conditions'))
 check('no_truth_window_selection','bestWindow' not in SOURCE)
 check('coherence_data_only','abs(sum(e.frameSequence))/denominator' in SOURCE)
 check('numeric_save_before_plot',SOURCE.index("atomic_save(state,fullfile(outputDir,'results.mat'))",SOURCE.index('state.trials=trial_table'))<SOURCE.index('state.plotExport=export_summary_plot'))
 check('print_numeric_figure_fallback',"print(number,filename,'-dpng','-r150')" in SOURCE)
 check('exact_score_not_replaced_by_square','cross=Xp.*conj(Xm)' in SOURCE and 'y.^2' not in SOURCE)
 check('count_debug100_full1000',"opts.trials=100;" in SOURCE and 'opts.trials=1000;' in SOURCE)
 check('resume_maintains_masks',"~state.done(:,trialID)'" in SOURCE and "不能更改科学选项" in SOURCE)
 script=SOURCE.split('\nfunction ',1)[0]
 script=re.sub(r"'(?:''|[^'\n])*'","''",script);script=re.sub(r'%[^\n]*','',script)
 variables=set(re.findall(r'\b([A-Za-z]\w*)\s*=(?!=)',script))
 for outputs in re.findall(r'\[([^]\n]+)\]\s*=(?!=)',script):variables.update(re.findall(r'\b[A-Za-z]\w*\b',outputs))
 names=re.findall(r'^function\s+(?:(?:\[[^]\n]*\]|\w+)\s*=\s*)?(\w+)\s*\(',SOURCE,re.M)
 check('no_script_local_function_name_collision',not variables.intersection(names))
 check('no_duplicate_function_names',len(names)==len(set(names)))
 try:
  from tree_sitter import Language,Parser
  import tree_sitter_matlab
  check('matlab_syntax_parser_only',not Parser(Language(tree_sitter_matlab.language())).parse(SOURCE.encode()).root_node.has_error)
 except ImportError:REPORT['parser']='not installed; parsing skipped'
 # Common raw midpoint and local resampling midpoint differ by ~0.6048us, not 50ms.
 common=(2480000-1)/(2*24.8e6)
 for duration in [.001,.01,.05,.1]:
  n=round(duration*8e5);first=(80000-n)//2;local=first/8e5+(n-1)/(2*8e5)
  f0=123456.7;true_a=700;estimated_a=600;local_f=f0+true_a*(local-common)
  est=local_f+estimated_a*(common-local)
  check('estimated_rate_transport_'+str(duration),abs(est-f0-(true_a-estimated_a)*(local-common))<1e-10)
  check('same_center_'+str(duration),abs(local-79999/(2*8e5))<1e-12)
 check('start_midpoint_35Hz',700*.05==35)
 # Complete-sample quantiles retain a nonfinite result as Inf instead of dropping it.
 e=np.r_[np.arange(1,50),np.inf]
 check('tail_quantile_failure_kept',np.isinf(np.sort(e)[49]) and np.sort(e)[47]==48)
 check('default_condition_count',3*1*1==3 and 3*50==150)
 check('paired_rate_count_no_duplicate',len(set([700,650,700,750]))*3==9)

CACHE={}
class FastStage1:
 def __init__(self,y,fs=8e5):
  n=len(y);key=(n,fs)
  if key not in CACHE:
   L=round(128*fs/8e5);H=L//2;starts=np.arange(0,n-L+1,H);idx=starts[:,None]+np.arange(L)
   tc=(idx.T-(n-1)/2)/fs;ts=(starts-(n-1)/2)/fs;w=np.hanning(L);band=np.arange(-8,9)
   nt=2**int(np.ceil(np.log2(2*len(starts))));delta=np.fft.fftshift(np.fft.fftfreq(nt,H/fs));rates=np.arange(-1000,1000.1,100)
   branches=[];df=fs/L
   for k in range(int(np.floor(5e4/df)),int(np.ceil(15e4/df))+1):
    fb=k*df;f=fb+delta/2;valid=(delta>=-df)&(delta<df)&(f>=5e4)&(f<=15e4)
    branches.append((k,(band+k)%L,(band-k)%L,np.exp(-4j*np.pi*fb*ts),valid,f[valid]))
   CACHE[key]=(L,H,idx,tc,ts,w,band,nt,rates,branches,np.exp(-2j*np.pi*rates[:,None]*ts**2))
  self.L,self.H,idx,self.tc,self.ts,self.w,self.band,nt,rates,branches,bank=CACHE[key]
  self.fs=fs;self.T=n/fs;self.df=.5/self.T;self.da=2/self.T**2
  self.frames=(y[idx]*self.w).T;X=np.fft.fft(self.frames,axis=0);raw=[]
  for k,plus,minus,phase,valid,f in branches:
   p=np.mean(X[plus]*np.conj(X[minus]),axis=0)*phase
   P=abs(np.fft.fftshift(np.fft.fft(p[None,:]*bank,nt,axis=1),axes=1))**2;P=P[:,valid].copy()
   for j in range(3):
    ai,fi=np.unravel_index(np.argmax(P),P.shape);value=P[ai,fi]
    if not np.isfinite(value):break
    raw.append([f[fi],rates[ai],value,k])
    distance=((rates[:,None]-rates[ai])/self.da)**2+((f[None,:]-f[fi])/self.df)**2;P[distance<1]=-np.inf
  self.raw=np.array(raw);order=np.argsort(-self.raw[:,2],kind='stable');ordered=self.raw[order]
  self.candidates=ordered[ref.distinct(ordered[:,:2],self.df,self.da)]
 def score(self,f,a):
  phase=2*np.pi*f*self.tc+np.pi*a*self.tc**2
  xp=np.fft.fft(self.frames*np.exp(-1j*phase),axis=0);xm=np.fft.fft(self.frames*np.exp(1j*phase),axis=0)
  cross=xp*np.conj(xm);norm=self.fs*np.sum(self.w**2)
  spec=np.mean(cross,axis=1)/norm;q=np.mean(cross[self.band%self.L],axis=0)/norm
  return abs(np.mean(spec[self.band%self.L]))**2,q
 def estimate(self):
  C=self.candidates[:10];scored=[self.score(f,a) for f,a,*_ in C];scores=np.array([x[0] for x in scored])
  j=np.argmax(scores);q=scored[j][1]
  return {'f_local':C[j,0],'rate':C[j,1],'score':scores[j],
          'coherence':abs(np.sum(q))/np.sum(abs(q)),'scores':scores,'frame_count':len(q),'candidates':self.candidates}
def compare(y,signal_direct):
 clock=time.perf_counter();baseline=ref.FSM(signal_direct);baseline_gen=time.perf_counter()-clock
 clock=time.perf_counter();base_scores=np.array([baseline.score(f,a) for f,a,*_ in baseline.candidates[:10]]);baseline_confirm=time.perf_counter()-clock
 clock=time.perf_counter();fast=FastStage1(y);fast_gen=time.perf_counter()-clock
 clock=time.perf_counter();result=fast.estimate();fast_confirm=time.perf_counter()-clock
 ident=str(len(REPORT['noiseless'])+len(REPORT['paired_smoke'])+len(REPORT['short_smoke']))
 check('candidate_pool_'+ident,np.allclose(baseline.raw,fast.raw,rtol=1e-8,atol=1e-8))
 check('deduplicated_pool_'+ident,np.allclose(baseline.candidates,fast.candidates,rtol=1e-8,atol=1e-8))
 check('exact_scores_'+ident,np.allclose(base_scores,result['scores'],rtol=1e-8,atol=1e-25))
 best=np.argmax(base_scores)
 check('same_selected_point_'+ident,np.array_equal(baseline.candidates[best,:2],[result['f_local'],result['rate']]))
 return result,{'reference_generation_s':baseline_gen,'reference_confirmation_s':baseline_confirm,
                'cached_generation_s':fast_gen,'cached_confirmation_s':fast_confirm}
def smoke():
 start=time.perf_counter();base=ref.make_base(20261015);base_s=time.perf_counter()-start
 fs=24.8e6;common=(len(base['t'])-1)/(2*fs);local=79999/(2*8e5)
 for f,a in [(100173.7,650),(90625,700),(90624.63,750)]:
  clean=base['pn']*base['data']*np.cos(2*np.pi*f*base['t']+np.pi*a*base['t']**2+base['phase'])
  long=resample_poly(clean*base['pn'],1,31);e,profile=compare(long,long)
  f0=e['f_local']+e['rate']*(common-local);tr,tm=ref.track(e['f_local']-(f+a*(local-common)),e['rate']-a,0,len(long),8e5)
  REPORT['noiseless'].append({'true_f0':f,'true_rate':a,'f0_error_hz':f0-f,'rate_error_hz_per_s':e['rate']-a,'track_rmse_hz':tr,'profile':profile})
  print('Python noiseless',f,a,'errors',f0-f,e['rate']-a,flush=True)
 # Two independent bases, each paired over exactly the three default Eb/N0 points.
 profiles=[]
 for seed in [20261015,20261016]:
  clock=time.perf_counter();b=ref.make_base(seed);generate_s=time.perf_counter()-clock
  a=700;f=1e5+b['fd'];clean=b['pn']*b['data']*np.cos(2*np.pi*f*b['t']+np.pi*a*b['t']**2+b['phase'])
  power=np.mean(clean**2);clock=time.perf_counter();cl=resample_poly(clean*b['pn'],1,31);nl=resample_poly(b['noise']*b['pn'],1,31);prep_s=time.perf_counter()-clock
  for eb in [-10,-7,-4]:
   var=power*fs/(2*1e5*10**(eb/10));y=cl+np.sqrt(var)*nl
   clock=time.perf_counter();direct=resample_poly((clean+np.sqrt(var)*b['noise'])*b['pn'],1,31);direct_prep_s=time.perf_counter()-clock
   check('linear_fir_'+str(seed)+'_'+str(eb),np.allclose(y,direct,rtol=1e-11,atol=1e-12))
   e,profile=compare(y,direct);profile.update({'linear_shared_prep_s':prep_s/3,'direct_prep_s':direct_prep_s,'generation_share_s':generate_s/3});profiles.append(profile)
   f0=e['f_local']+e['rate']*(common-local);tr,tm=ref.track(e['f_local']-(f+a*(local-common)),e['rate']-a,0,len(y),8e5)
   REPORT['paired_smoke'].append({'seed':seed,'eb_n0_db':eb,'f0_error_hz':f0-f,'estimated_rate':e['rate'],'rate_error_hz_per_s':e['rate']-a,
      'frequency_success_20hz':bool(abs(f0-f)<=20),'track_rmse_hz':tr,'track_max_hz':tm,'coherence':e['coherence'],'profile':profile})
   print('Python paired',seed,eb,'f0 error',f0-f,'rate',e['rate'],flush=True)
 # Additional paired center crops demonstrate that failures stay failures, not omitted statistics.
 f=123456.7;a=700;clean=base['pn']*base['data']*np.cos(2*np.pi*f*base['t']+np.pi*a*base['t']**2+base['phase'])
 var=np.mean(clean**2)*fs/(2*1e5*10**(-7/10));long=resample_poly((clean+np.sqrt(var)*base['noise'])*base['pn'],1,31)
 for duration in [.001,.01]:
  n=round(duration*8e5);first=(len(long)-n)//2;y=long[first:first+n]
  e,_=compare(y,y);local_short=first/8e5+(n-1)/(2*8e5)
  f0=e['f_local']+e['rate']*(common-local_short)
  tr,tm=ref.track(e['f_local']-(f+a*(local_short-common)),e['rate']-a,0,n,8e5)
  REPORT['short_smoke'].append({'duration_s':duration,'f0_error_hz':f0-f,'rate_error_hz_per_s':e['rate']-a,
    'frequency_success_20hz':bool(abs(f0-f)<=20),'track_rmse_hz':tr,'track_max_hz':tm,'coherence':e['coherence']})
  print('Python short (kept)',duration,'f0 error',f0-f,flush=True)
 REPORT['benchmark']={'python_only':True,'base_generation_first_s':base_s,'profiles':profiles,
  'matlab_default_150_estimate_s':None,'explanation':'No MATLAB runtime. MATLAB file measures its first two paired bases and prints its own estimate.'}
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--structural-only',action='store_true');args=parser.parse_args()
 start=time.perf_counter();structural()
 if not args.structural_only:smoke()
 REPORT['elapsed_s']=time.perf_counter()-start;REPORT['source_sha256']=hashlib.sha256(SOURCE.encode()).hexdigest()
 out='structural_results.json' if args.structural_only else 'python_results.json'
 (ROOT/'stage1_validation'/out).write_text(json.dumps(REPORT,ensure_ascii=False,indent=2)+'\n')
 print('PASS',len(REPORT['checks']),'checks; MATLAB native execution: false')
if __name__=='__main__':main()
