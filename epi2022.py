import numpy as np
import pandas as pd
import pylab as pl
import scipy.integrate as integrate
import scipy.optimize as optimize
from scipy.ndimage.filters import uniform_filter1d as ufilt1d
import sys
import getpass as gp

state = "Utah"
#state = "Idaho"
#state = "Florida"
#state = "Missouri"

# choose a state at the command line?
if len(sys.argv)>1:
    state = ' '.join(sys.argv[1:])

# get population try US Census Bureau...

urlpop = "https://www2.census.gov/programs-surveys/popest/datasets/2010-2019/national/totals/nst-est2019-alldata.csv"
df = pd.read_csv(urlpop)
df.index = df["NAME"]
# set up a series:
pop = pd.Series(df["POPESTIMATE2019"].values, index=df["NAME"]).astype(int)
print(f'{state} population (2019): {pop[state]}')

# get case counts...

urlcovid = "https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv"
df = pd.read_csv(urlcovid)

print('cols',df.columns)
df = df.loc[df["state"]==state].copy() # just keep this state's data

df['dcases'] = (df['cases']-df['cases'].shift(1,fill_value=0)) # daily cases


# smooth out some of the reporting/download irregularities...
df['dcases7dayavg'] = ufilt1d(df.dcases, size=7,mode='nearest')

# create a time coord in days since first case
df.index = np.arange(len(df.index))
df['dday'] = df['date'].to_numpy().astype(np.datetime64)
df['dday'] = (df.dday-df.dday[0]).astype('timedelta64[D]').astype(float)

# turn to model:

def R0(b,g): return b/g

def dsir(t,Y, N, b0, g0):
    b, g = b0, g0
    S,I,R = Y
    dS, dI, dR = -b*I*S/N, b*I*S/N - g*I, g*I
    return dS, dI, dR

def sir(t,tinit,N,b0,g0):
    Iinit = 1.0; Sinit = N-Iinit; Rinit = 0.0
    S, I, R = np.ones(len(t))*Sinit, np.zeros(len(t)), np.zeros(len(t))
    tmsk = t>=tinit # start time is different from t[0], maybe
    sol = integrate.solve_ivp(dsir,[tinit,t[-1]],[Sinit,Iinit,Rinit],                              t_eval=t[tmsk],args=(N,b0,g0))
    S[tmsk],I[tmsk],R[tmsk] = sol.y
    return S,I,R

def dcasesmodel(t,tinit,N,b0,g0):
    S,I,R = sir(t,tinit,N,b0,g0)
    dS = -(S[1:]-S[:-1]) # change in suseptible from previous day
    return np.append(np.array([0]),dS)

def chi2ish(p,dday,dcases,tepi_start,tfit_start,tfit_end,N):
    b0,g0 = p
    model = dcasesmodel(dday.to_numpy(),tepi_start,N,b0,g0)
    tmsk = (dday>=tfit_start)&(dday<tfit_end)
    return np.sum((model[tmsk]-dcases[tmsk])**2)

def chi2isht(pt,dday,dcases,tfitstart,tfitend,N):
    b0,g0,tepi_start = pt
    if b0<0.05 or g0<0.05 or tepi_start < tfitstart-60: return 1e99
    model = dcasesmodel(dday.to_numpy(),tepi_start,N,b0,g0)
    tmsk = (dday>=tfitstart)&(dday<tfitend)
    return np.sum((model[tmsk]-dcases[tmsk])**2)

g0 = 1/3.; b0 = 1.1*g0
N = pop[state]
df['dmod'] = dcasesmodel(df.dday.to_numpy(),0,N,b0,g0)

# wild/alpha
bgguess = [b0,g0]
t0wild,tfbeg,tfend = 0,0,420 
args=(df.dday,df.dcases,t0wild,tfbeg,tfend,N)
sol = optimize.minimize(chi2ish,bgguess,args=args,method='Nelder-Mead')        
bgwild = sol.x
print('b,g fit:',sol.x)
df['dmod'] = dcasesmodel(df.dday.to_numpy(),0,N,*bgwild)
Rtot = np.sum(df.dmod)
print(Rtot)

# delta
bgguess = bgwild
t0delta,tfbeg,tfend = 320,320,df.index[-1] 
t0delta,tfbeg,tfend = 320,320,630
Ndelta = 0.5*N-Rtot
fitguess = np.append(bgguess,[t0delta])
args=(df.dday,df.dcases,tfbeg,tfend,Ndelta,)
# fitting for start time, too!
sol = optimize.minimize(chi2isht,fitguess,args=args,method='Nelder-Mead')        
bgdelta = sol.x[:2]
t0delta = sol.x[2] 
print('b,g fit:',sol.x,R0(*bgdelta))
df['dmod'] += dcasesmodel(df.dday.to_numpy(),t0delta,Ndelta,*bgdelta)

# add in omicron, do not fit, just use delta and match
bgom = np.array([1.25*bgdelta[0],1.1*bgdelta[1]])
t0om = 0
print('t index, date',df.index[-1],df.date[df.index[-1]])
# find index at 2022-01-01:
tmark = df.index[df.date=='2022-01-01'][0]
#tmark = min(690,df.index[-1])
Nom = 0.5*N
dcasemark = df['dcases7dayavg'][tmark]
dmodmark = df['dmod'][tmark]

dcom = dcasesmodel(df.dday.to_numpy(),t0om,Nom,*bgom)
tx = np.argmax(dcom>dcasemark-dmodmark)
if tx<=0: tx=tmark
t0om = t0om + (tmark-tx)


savethedate = df.index[-1]
xdays = 200
datex = pd.date_range(start=df['date'].iloc[0], freq="1d",periods=len(df.index)+xdays)
df = df.reindex(np.arange(len(datex)))
df['date'] = datex
df['dday'] = df['date'].to_numpy().astype(np.datetime64)
df['dday'] = (df.dday-df.dday[0]).astype('timedelta64[D]').astype(float)

df.dmod = dcasesmodel(df.dday.to_numpy(),t0wild,N,*bgwild)
df.dmod += dcasesmodel(df.dday.to_numpy(),t0delta,Ndelta,*bgdelta)

#print(savethedate)
#bgom = 2*np.array([bgdelta[0],bgdelta[1]])

dcom = dcasesmodel(df.dday.to_numpy(),t0om,Nom,*bgom)
print('Omicron peak:',df.date[np.argmax(dcom)])
df['dmod'] += dcom

# shift to case load on Jan 4


print(bgom, R0(*bgom))


# show the world...
#ax = df.plot('date','dcases',color='#aaaaff', label='daily cases')
#ax.plot(df['date'],df['dcases7dayavg'],'-k',label='7-day running average')
#ax.plot(df.date,df.dmod,'k:',label='model')

ax = df.plot('date',['dcases','dcases7dayavg','dmod'],color=['#aaaaff','#0000ff','#663333'],
             label=['daily cases','7-day running average','model'])

lines = ax.get_lines()
zo=[0,2,1]
for i, line in enumerate(lines, -len(lines)):
        line.set_zorder(zo[i])
        

ax.set_ylabel('case count')
ax.figure.autofmt_xdate()
ax.set_title(state+' ')
ax.figure.savefig('tmp.png')

# comment this out!?
if gp.getuser() == 'bromley':
    import os
    os.system('convert tmp.png ~/www/tmp.jpg')

print('done')
quit()


