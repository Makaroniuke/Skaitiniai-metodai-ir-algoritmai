import matplotlib.pyplot as plt
import numpy as np
import math
import sympy as sym

# x0 = 0   # pradinis artinys

STEP = 0.001
EPS = 1e-10
maxI = 100 # max allowed iterations



coeff_f = [0.1, -0.05, -1.95, 1.75, 5.18, -2.14]
# 0.10𝑥^5 − 0.05𝑥^4 − 1.95𝑥^3 + 1.75𝑥^2 + 5.18𝑥 − 2.14 
# def f(x):
#   return (0.10 * pow(x, 5)) - (0.05 * pow(x, 4)) - (1.95 * pow(x 3)) + (1.75 * pow(x, 2)) + (5.18 * x) - 2.14;
 

# ((x+1)^2*(x-3)^2)/(x^3 + 2) + (x-2)^3*cos(x);   0<=x<=15
def g(x):
  return (((x+1)**2*(x-3)**2)/(pow(x, 3) + 2)) + (x-2)**3*math.cos(x); 


def v(x):
  return 100*math.exp(-0.25/x)+(9.8*x/0.05)*(math.exp(-0.25/x) - 1) - 22



x, fx, df, gx, dg, vx, dv = sym.symbols(('x', 'fx', 'df', 'gx', 'dg', 'vx', 'dv'))
fx= 0.10*x**5 - 0.05*x**4 - 1.95*x**3 + 1.75*x**2 + 5.18*x - 2.14 
df = fx.diff(x)
gx = ((x+1)**2*(x-3)**2)/(x**3 + 2) + (x-2)**3*sym.cos(x)
dg = gx.diff(x)
vx =100*sym.exp(-0.25/x)+((9.8*x)/0.05)*(sym.exp(-0.25/x) - 1) - 22
dv = vx.diff(x)


# suskaičiuoja visus galimus taškus funcijoje duotame intervale
def create_function_points(function, intervals):
    points = []
    # arrange - values are generated within the interval, with spacing between values given by step
    for number in np.arange(intervals[0], intervals[1], STEP):
        points.append(function(number))
    return points


# atspausdina duota funkciją
def show_graph(intervals, func):
    # surandami taškai
    points = create_function_points(func, intervals)
    # add a horizontal line across the Axes
    plt.axhline(y=0, color='r', linestyle=':')
    plt.plot(np.arange(intervals[0], intervals[1], STEP), points)
    plt.plot(intervals[0], 0, marker='x', color='red', markersize=12)
    plt.plot(intervals[1], 0, marker='x', color='red', markersize=12)
    plt.grid(True)
    plt.show()

def grubus(coef):
  largest_number = 0
  for co in coef:
    if abs(co) > largest_number:
        largest_number = abs(co)
  result = 1+(largest_number/coef[0])
  return [-result, result]

def tikslesnis_teigiamas(coef):
  smallest_number = 0
  index = 0
  for co in coef:
    if co < 0:
      if abs(co) > smallest_number:
        smallest_number = co
  for c in range(len(coef)):
    if coef[c] < 0:
      index = c + 1;
      
  b = abs(smallest_number) 
  k = len(coef) - 1 -(len(coef) - index) 
  result = 1 + (b/coef[0])**(1/k)
  return result;

def tikslesnis_neigiamas(coef):
  
  for c in range(len(coef)):
    if c % 2 == 0:
      coef[c] *= -1;
      
  if coef[0] < 0:
    for co in range(len(coef)):
      coef[co] *= -1;

  smallest_number = 0
  for co in coef:
    if co < 0:
      if abs(co) > smallest_number:
        smallest_number = abs(co)
  for c in range(len(coef)):
    if coef[c] < 0:
      index = c + 1;      
  
  b = smallest_number
  k = len(coef) - 1 -(len(coef) - index)
  result = 1 + (b/coef[0])**(1/k)
  return result;


coef = [10, 0, -1, 0, -2, 0]
print('rezai')
print(tikslesnis_neigiamas(coef))
print(grubus(coef))

# suranda visus šaknų intervalus 
def scan(interval, func):
  x1 = interval[0];
  x2 = x1 + STEP;
  root_intervals = [];
  while x1 < interval[1]:
    if (np.sign(func(x1)) != np.sign(func(x2))):
      root_intervals.append([x1, x2]) 
    x1 = x1 + STEP
    x2 = x2 + STEP
  return root_intervals;

# skenavimo metodas su mažėjančiu žingsniu
def scan_method(root_interval, func):
  scan_step = STEP;
  iterations = 0;
  x1 = root_interval[0];
  x2 = root_interval[1];
  precision = 1;
  while precision > EPS and maxI > iterations:
    iterations += 1;
    if (np.sign(func(x1)) != np.sign(func(x2))):
        scan_step = (scan_step / 2);
        x1 -= scan_step;
        x1 = x1 + scan_step;
        x2 = x1 + scan_step;
    elif func(x1) is 0:
      return x1
    else:
      x1 = x1 + scan_step;
      x2 = x2 + scan_step;
    precision = abs(func(x1))
  return x1, precision, iterations

def simple_iteration(interval, func, i):
  x = interval[0] ;
  precision = 1;
  iterations = 0;
  while precision > EPS and maxI > iterations:
    iterations += 1;
    x_next = (func(x)/alpha[i]) + x;
    
    precision = abs(x - x_next)
    x = x_next
  return x, iterations, precision, alpha[i]

def newton(interval, fx, df, func):
  precision = 1;
  iterations = 0;
  x0 = interval[0];
  while precision > EPS and maxI > iterations:
    iterations += 1;  
    x1 = x0 - fx.subs(x,x0).evalf()/df.subs(x,x0).evalf()
    precision = abs(x1-x0)
    x0 = x1
  return x1, iterations, precision;



# paleidžiama f(x) funkcija ir skaičiavimai
def run_function(f, fx, df, intervals):
  show_graph(intervals, f)

  root_interval = scan(intervals, f);

  print("Skenavimo mažėjančiu žingsniu metodas")
  print(("| {0:^42} | {1:^20} | {2:^10} | {3:^25} |").format("Intervalas", "Šaknis", "Iteracijos", "Tikslumas"))
  for it in root_interval:
    x, precisions, iteration = scan_method(it, f)
    print(("| {0:<20}- {1:^20} | {2:^20} | {3:^10} | {4:^25} |").format(it[0],it[1], x, iteration, precisions))
  print("\n")

  print("Niutono liestinių metodas")
  print(("| {0:^42} | {1:^20} | {2:^10} | {3:^25} |").format("Intervalas", "Šaknis", "Iteracijos", "Tikslumas"))
  for it in root_interval:
    x, iteration, pr = newton(it, fx, df, f)
    print(("| {0:<20}- {1:^20} | {2:^20} | {3:^10} | {4:^25} |").format(it[0], it[1], x, iteration, pr or 0))
  print("\n")

  print("Parastųjų iteracijų metodas")
  print(("| {0:^42} | {1:^20} | {2:^10} | {3:^25} | {4:^7} |").format("Intervalas", "Šaknis", "Iteracijos", "Tikslumas", "Alpha"))
  for it in range(len(root_interval)):
   x, iteration, precisions, alpha = simple_iteration(root_interval[it], f, it)
   print(("| {0:<20}- {1:^20} | {2:^20} | {3:^10} | {4:^25} | {5:^7} |").format(root_interval[it][0], root_interval[it][1], x, iteration, precisions, alpha))
   

 
def run_v_function():
  intervals = [0.3, 1];
  show_graph(intervals, v)

  root_interval = scan(intervals, v);

  print("Niutono liestinių metodas")
  print(("| {0:^42} | {1:^20} | {2:^10} | {3:^25} |").format("Intervalas", "Šaknis", "Iteracijos", "Tikslumas"))
  for it in root_interval:
    x, iteration, pr = newton(it, vx, dv, v)
    print(("| {0:<20}- {1:^20} | {2:^20} | {3:^10} | {4:^25} |").format(it[0], it[1], x, iteration, pr))

def fx_intervals():
  print("Grubus įvertis") 
  print(grubus(coeff_f));
  print("\n")
  print("Tikslesnis teigiamas įvertis")
  print(tikslesnis_teigiamas(coeff_f))
  print("Tikslesnis neigiamas įvertis")
  print(-tikslesnis_neigiamas(coeff_f))


# fx_intervals();


alpha = [-100, 10, -10, 10, -10]  # daugiklio reiksme
# run_function(f, fx, df, [-5.4,5])

# alpha = [-30, 1, -20, 200, -600, 1000]
# run_function(g, gx, dg, [0,15])


# run_v_function()
















