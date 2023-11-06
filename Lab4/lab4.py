import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp


def funk(X, t):  # PDL desines puses funkcija
  h1=X[0];h2=X[1];v1=X[2];v2=X[3]; 
  m1=0.2; m2=0.4; k1=0.02; k2=0.005; k = 0.015       
  rez=np.zeros(4,dtype=float)  # k = 1, p = 2

  if 1 > t:   # iki atsiskyrimo
    rez[0]=v1;rez[1]=v2;
    rez[2]=-9.81+(-k*v1**2*np.sign(v1))/(m2+m1);
    rez[3]=rez[2];
  else:       # po atsiskyrimo
    rez[0]=v1;rez[2]=-9.81-k1*v1**2*np.sign(v1)/m1
    rez[1]=v2; rez[3]=-9.81-k2*v2**2*np.sign(v2)/m2;
    
  return rez

ttt=15  # skaiciavimo laikas  (s)
step = 0.01
t=np.arange(0, 15, 0.1)
dt=t[1]-t[0]; #laiko momentai
h0=0  # pradinis aukstis (m)
N=len(t); rez=np.zeros([4,N],dtype=float) # rezultatu masyvas
rez[:,0]=np.array([h0,h0,80,80]);           # pradines salygos

if 0:   # sprendimas Eulerio metodu
  for i in range (N-1):  
   rez[:,i+1]=rez[:,i]+funk(rez[:,i],t[i])*dt
   if rez[0,i+1] <= 0 : rez[[0,2],i+1]=0  
   if rez[1,i+1] <= 0 : rez[[1,3],i+1]=0
else: 
  for i in range (N-1) :  # sprendimas IV RK metodu
    fz=rez[:,i]+funk(rez[:,i],t[i])*dt/2
    fzz=rez[:,i]+funk(fz,t[i]+dt/2)*dt/2
    fzzz=rez[:,i]+funk(fzz,t[i]+dt/2)*dt
    rez[:,i+1]=rez[:,i]+dt/6*(funk(rez[:,i],t[i])+2*funk(fz,t[i]+dt/2)+2*funk(fzz,t[i]+dt/2)+funk(fzzz,t[i]+dt))
    if rez[0,i+1] <= 0 : rez[[0,2],i+1]=0  
    if rez[1,i+1] <= 0 : rez[[1,3],i+1]=0



# rezultatu pavaizdavimas: 
fig1=plt.figure(1); ax1=fig1.add_subplot(1,1,1); ax1.set_xlabel('t, s');ax1.set_ylabel('h, m');ax1.grid();plt.title('Aukštis')
ax1.plot(t,rez[0,:],'b-');ax1.plot(t,rez[1,:],'r-'); plt.legend(['M1','M2']);
maxy1 = np.max(rez[0])
result1 = np.where(rez[0] == maxy1)
maxt1 = t[result1[0]][0]
ax1.plot(maxt1, maxy1, marker='x', color='green')

maxy = np.max(rez[1])
result = np.where(rez[1] == maxy)
maxt = t[result[0]][0]
ax1.plot(maxt, maxy, marker='x', color='green')
plt.show()

print(f"M1 pasiekė aukščiausia tašką {maxy1}m. - {maxt1}s.")
print(f"M2 pasiekė aukščiausia tašką {maxy}m. - {maxt}s.")

fig2=plt.figure(1); ax2=fig2.add_subplot(1,1,1); ax2.set_xlabel('t, s');ax2.set_ylabel('v, m/s');ax2.grid();plt.title('Greitis')
ax2.plot(t,rez[2,:],'b-');ax2.plot(t,rez[3,:],'r-');plt.legend(['M1','M2']);plt.show()








# # patikrinimas
# t = np.array([0,15], dtype="float")
# function_rez = solve_ivp(funk, t, [0,0,80,80])
# H1 = function_rez.y
# T = function_rez.t
# fig2=plt.figure(1); ax2=fig2.add_subplot(1,1,1); ax2.set_xlabel('t');ax2.set_ylabel('v');ax2.grid();plt.title('Aukštis')
# ax2.plot(T,H1[0,:],'b-');ax2.plot(T,H1[1,:],'r-');plt.legend(['M1','M2']);plt.show();

# fig2=plt.figure(1); ax2=fig2.add_subplot(1,1,1); ax2.set_xlabel('t');ax2.set_ylabel('v');ax2.grid();plt.title('Greitis')
# ax2.plot(T,H1[2,:],'b-');ax2.plot(T,H1[3,:],'r-');plt.legend(['M1','M2']);plt.show();






# print(function_rez.y)