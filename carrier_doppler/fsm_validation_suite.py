"""Full-passband Python cross-check of FSM core; MATLAB remains unexecuted."""
from fsm_frequency_validate import estimate,gold_code,FS
from scipy.signal import resample_poly
import numpy as np
import json,csv,argparse,time
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--trials',type=int,default=100);args=p.parse_args()
rng=np.random.default_rng(20261015);rawfs=24.8e6;rb=1e5;T=.1;N=round(rawfs*T)
t=(np.arange(N)-(N-1)/2)/rawfs;pn=gold_code()[(np.arange(N)//8)%31]
ref_shift=(round(FS*T)-1)/(2*FS)-(N-1)/(2*rawfs)
def run(f,a,eb):
 data=np.repeat(rng.choice([-1.,1.],round(rb*T)),248)
 s=data*pn*np.cos(2*np.pi*f*t+np.pi*a*t*t+rng.uniform(-np.pi,np.pi))
 x=s.copy()
 if eb is not None:x+=rng.normal(size=N)*np.sqrt(np.mean(s*s)*rawfs/(2*rb*10**(eb/10)))
 fc,ac,ff,af,ok=estimate(resample_poly(x*pn,1,31))
 true=f+a*ref_shift
 return [true,a,fc,ac,ff,af,fc-true,ac-a,ff-true,af-a,int(ok)]
start=time.monotonic();scan=[]
for f in [50000,50037.7,100173.7,149963.3,150000]:
 for a in [-1000,-873.1,0,731.4,1000]:scan.append(run(f,a,None))
print('Completed noiseless 25-point boundary/off-grid scan',flush=True)
rows=[]
for k in range(args.trials):
 rows.append(run(rng.uniform(55000,145000),rng.uniform(-900,900),-7))
 if (k+1)%10==0: print('FSM -7 dB trial',k+1,'elapsed_s',round(time.monotonic()-start,1),flush=True)
r=np.array(rows);ss=np.array(scan)
cok=(abs(r[:,6])<=20)&(abs(r[:,7])<=200)
fok=(abs(r[:,8])<=20)&(abs(r[:,9])<=200)&(r[:,10]>0)
prob=float(fok.mean());n=len(rows);z=1.96;den=1+z*z/n;c=(prob+z*z/(2*n))/den;h=z*np.sqrt(prob*(1-prob)/n+z*z/(4*n*n))/den
summary={'noise_trials':n,'eb_n0_db':-7,'duration_s':T,'algorithm':'Frequency-domain FSM shifted FFT cross-products, complex time averaging and 17-bin frequency smoothing; polynomial phase rate compensation. Not signal squaring.','preprocessing':'known sampled PN, exact synchronization, real 31x resampling to 800 kHz','coarse_joint_success':float(cok.mean()),'fine_joint_success':prob,'fine_wilson_95':[max(0,c-h),min(1,c+h)],'coarse_f_rmse_hz':float(np.sqrt(np.mean(r[:,6]**2))),'coarse_rate_rmse_hz_s':float(np.sqrt(np.mean(r[:,7]**2))),'fine_f_rmse_hz':float(np.sqrt(np.mean(r[:,8]**2))),'fine_rate_rmse_hz_s':float(np.sqrt(np.mean(r[:,9]**2))),'noiseless_scan_points':len(scan),'noiseless_max_f_error_hz':float(max(abs(ss[:,8]))),'noiseless_max_rate_error_hz_s':float(max(abs(ss[:,9]))),'matlab_executed':False,'frequency_tolerance_hz':20,'rate_tolerance_hz_s':200,'exact_coarse_confirmation_candidates':5,'fine_f_bias_hz':float(np.mean(r[:,8])),'fine_rate_bias_hz_s':float(np.mean(r[:,9])),'fine_f_std_hz':float(np.std(r[:,8],ddof=1)),'fine_rate_std_hz_s':float(np.std(r[:,9],ddof=1)),'strict_frequency_improvement':float(np.mean(abs(r[:,8])<abs(r[:,6]))),'strict_rate_improvement':float(np.mean(abs(r[:,9])<abs(r[:,7]))),'elapsed_s':time.monotonic()-start}
root=Path(__file__).parent
(root/'fsm_validation_summary.json').write_text(json.dumps(summary,indent=2))
header=['true_f_hz','true_rate_hz_s','coarse_f_hz','coarse_rate_hz_s','fine_f_hz','fine_rate_hz_s','coarse_f_error','coarse_rate_error','fine_f_error','fine_rate_error','converged']
for name,records in [('fsm_validation_trials.csv',rows),('fsm_noiseless_scan.csv',scan)]:
 with (root/name).open('w',newline='') as f:
  writer=csv.writer(f,lineterminator="\n");writer.writerow(header);writer.writerows(records)
print(json.dumps(summary,indent=2))
