import pandas as pd
import numpy as np
d1 = pd.read_csv('temp/LEDD_CATDenovoPD_vs_LEDD_CATControl.csv')
d2 = pd.read_csv('temp/LEDD_CATPD_vs_LEDD_CATControl.csv')
df = pd.merge(d1[['Probe', 'logFC', 't', 'adj.P.Val']], 
              d2[['Probe', 'logFC', 't', 'adj.P.Val']], on='Probe')
df['tt'] = df['t_x'] * df['t_y']
df.head()
d1sigProbe = d1[d1['adj.P.Val'] < 0.05].Probe
d2sigProbe = d2[d2['adj.P.Val'] < 0.05].Probe
dfsigProbe = np.intersect1d(d1sigProbe, d2sigProbe)
df[df.Probe.isin(dfsigProbe)]
len(dfsigProbe)