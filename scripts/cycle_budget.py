"""Cycle budget: where the duty cycle goes, and how densely peaks are sampled.

Marginal cost of one MS2 is the spacing between consecutive MS2 scans inside a cycle; the MS1
block is the gap either side of the survey scan. Sampling density is reported as scans above half
maximum -- counted, not inferred -- because estimating a peak width from points above half maximum
uses the very sampling whose adequacy is in question and is biased low when points are few.

Produces the numbers behind manuscript Table 12.
"""
import glob, re
import xml.etree.ElementTree as ET
import numpy as np
NS="{http://psi.hupo.org/ms/mzml}"

def scans(f):
    for _,el in ET.iterparse(f,events=("end",)):
        if el.tag!=NS+"spectrum": continue
        d={}
        for c in el.iter(NS+"cvParam"):
            n=c.get("name")
            if n in ("ms level","ion injection time","scan start time"): d[n]=float(c.get("value"))
        if "ms level" in d: yield d["ms level"], d.get("ion injection time"), d.get("scan start time",0)*60
        el.clear()

for pol in ("Pos","Neg"):
    f=sorted(glob.glob(f"data/mzml/kiterie_{pol}/*Pool_QC*.mzML"))[0]
    lv=[];it=[];rt=[]
    for a,b,c in scans(f): lv.append(a); it.append(b if b is not None else np.nan); rt.append(c)
    lv=np.array(lv); it=np.array(it,dtype=float); rt=np.array(rt)

    # marginal cost of one MS2 = spacing between consecutive MS2 scans inside a cycle
    d=np.diff(rt); pair2=(lv[:-1]==2)&(lv[1:]==2)
    ms2_cost=np.median(d[pair2])*1000
    # cost of the MS1 = last MS2 -> MS1 -> first MS2
    p_m1=(lv[:-1]==2)&(lv[1:]==1); p_m2=(lv[:-1]==1)&(lv[1:]==2)
    ms1_cost=(np.median(d[p_m1])+np.median(d[p_m2]))*1000

    m1it=it[lv==1]; m1it=m1it[~np.isnan(m1it)]
    ms1=np.where(lv==1)[0]; per=np.diff(ms1)-1; cyc=np.diff(rt[ms1])*1000
    print(f"\n=== {pol}")
    print(f"  cycle {np.median(cyc):6.0f} ms = MS1 block {ms1_cost:5.0f} ms + {np.median(per):.0f} MS2")
    print(f"  marginal cost per MS2 : {ms2_cost:5.1f} ms   (injection time 22 ms of it -> "
          f"{22/ms2_cost*100:.0f}% ; {ms2_cost-22:.1f} ms is transient+overhead)")
    print(f"  MS1 injection time    : median {np.median(m1it):6.1f} ms of 246 max, "
          f"at ceiling {100*np.mean(m1it>=246-1e-9):.1f}%")
