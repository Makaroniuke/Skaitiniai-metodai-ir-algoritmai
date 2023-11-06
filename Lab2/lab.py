import numpy as np


# QR skaidos, Gauso - Zeidelio metodai

# Pirma lygtis
# 2𝑥1 + 𝑥2 + 𝑥3 + 𝑥4 = 6
# 𝑥1 + 3𝑥2 + 𝑥3 − 3 𝑥4 = −4
# 𝑥1 + 𝑥2 + 5𝑥3 + 𝑥4 = 4
# 2𝑥1 + 3𝑥2 − 3𝑥3 − 2 𝑥4 = 0

# Antra lygtis
# 3𝑥1 + 𝑥2 − 𝑥3 + 5𝑥4 = 8
# −3𝑥1 + 4𝑥2 − 8𝑥3 − 𝑥4 = 10
# 𝑥1 − 3𝑥2 + 7𝑥3 + 6𝑥4 = 11
# 5𝑥2 − 9𝑥3 + 4𝑥4 = 1

A24 = np.matrix([[2,1,1,1],[1,3,1,-3],[1,1,5,1],[2,3,-3,-2]], dtype=float);
B24 = np.matrix([[6],[-4],[4],[0]], dtype=float);
n24=(np.shape(A24))[0] 
nb24=(np.shape(B24))[1]



A27 = np.matrix([[3,1,-1,5],[-3,4,-8,-1],[1,-3,7,6],[0,5,-9,4]], dtype=float);
B27 = np.matrix([[8],[10],[11],[1]], dtype=float);
n27=(np.shape(A27))[0] 
nb27=(np.shape(B27))[1]



# ==== QR skaidos metodas ====


def qr_skaida(A, b, n, nb):
  Ap = A;
  Q=np.identity(n)

  for i in range(0,n-1):
      z=A[i:n,i] # A matricos stulpelis
      zp=np.zeros(np.shape(z)); # atspindetas vektorius
      zp[0]=np.linalg.norm(z)
      omega=z-zp; 
      omega=omega/np.linalg.norm(omega)
      Qi=np.identity(n-i)-2*omega*omega.transpose()
      A[i:n,:]=Qi.dot(A[i:n,:])
      Q[:,i:n]=Q[:,i:n].dot(Qi) # kaupia visas q sandauga
      print("Iteracija: ", i+1)
      print(A)

  # atgalinis etapas:
  b1=Q.transpose().dot(b);
  x=np.zeros(shape=(n,nb));
  for i in range (n-1,-1,-1):  
      x[i,:]=(b1[i,:]-A[i,i+1:n]*x[i+1:n,:])/A[i,i];
  print(x)

  liekana=Ap.dot(x)-b1;
  print(liekana);
  print(np.linalg.norm(liekana)/ np.linalg.norm(x))


qr_skaida(A24, B24, n24, nb24);
# qr_skaida(A27, B27,n27, nb27);


def check_with_numpy(a, b):
    q, r = np.linalg.qr(a)
    y = np.dot(q.T, b)
    return np.linalg.solve(r, y)


# x_check = check_with_numpy(A24,B24);
# print("Patikrintos reikšmės:")
# print(x_check)

# x_check = check_with_numpy(A27,B27);
# print("Patikrintos reikšmės:")
# print(x_check)


# ==== Gauso - Zeidelio metodas ====

def gauso_zeidelio(A, b, n):
  Ap = A;
  alpha = np.array([1, 1, 1, -10]);
  # alpha = np.array([10, 10, 100, 100]);
  Atld=np.diag(1./np.diag(A)).dot(A)-np.diag(alpha)
  btld=np.diag(1./np.diag(A)).dot(b)
  nitmax=1000; 
  eps=1e-12
  x=np.zeros(shape=(n,1));
  x1=np.zeros(shape=(n,1));

  for it in range (0,nitmax):
    for i in range (0,n) : 
      # A padauginta is x1 vektoriaus
      x1[i]=(btld[i]-Atld[i,:].dot(x1))/alpha[i];
    prec=(np.linalg.norm(x1-x)/(np.linalg.norm(x)+np.linalg.norm(x1)))
    if prec < eps : break
    x[:]=x1[:]

  print(x)
  liekana=Ap.dot(x)-b;
  print(liekana);
  print(np.linalg.norm(liekana)/ np.linalg.norm(x))


# gauso_zeidelio(A24, B24, n24)
# gauso_zeidelio(A27, B27, n27)


# X=np.linalg.solve(A24,B24)
# print(X, "patikrintas")

# X=np.linalg.solve(A27,B27);
# print(X, "patikrintas")