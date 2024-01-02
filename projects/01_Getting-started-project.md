---
jupytext:
  formats: notebooks//ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.4
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('fivethirtyeight')
```

# Computational Mechanics Project #01 - Heat Transfer in Forensic Science

We can use our current skillset for a macabre application. We can predict the time of death based upon the current temperature and change in temperature of a corpse. 

Forensic scientists use Newton's law of cooling to determine the time elapsed since the loss of life, 

$\frac{dT}{dt} = -K(T-T_a)$,

where $T$ is the current temperature, $T_a$ is the ambient temperature, $t$ is the elapsed time in hours, and $K$ is an empirical constant. 

Suppose the temperature of the corpse is 85$^o$F at 11:00 am. Then, 45
min later the temperature is 80$^{o}$F. 

Assume ambient temperature is a constant 65$^{o}$F.

1. Use Python to calculate $K$ using a finite difference approximation, $\frac{dT}{dt} \approx \frac{T(t+\Delta t)-T(t)}{\Delta t}$.

+++

## 1.
First, we must rearrange Newton's Law of Cooling to solve for K:

$K = -\frac{dT}{dt}\frac{1}{(T-T_a)}$

Then we can replace $\frac{dT}{dt}$ using $\frac{dT}{dt} \approx \frac{T(t+\Delta t)-T(t)}{\Delta t}$ and then the equation for K becomes:

$K = -\frac{T(t+\Delta t)-T(t)}{\Delta t}\frac{1}{(T-T_a)}$

```{code-cell} ipython3
T_1 = 85 # degrees Farenheit
T_2 = 80 # degrees Farenheit
T_a = 65 # degrees Farenheit

def timetomins(t, ampm):
    '''Takes time input and turns it into absolute minutes from 00:00 am
    
    Arguments
    ---------
    t = time input integer in the form of hhmm (hour, hour, minute, minute)
    ampm = am or pm
    
    Returns
    ---------
    min = the total number of minutes that make up the time input since 00:00 am '''
    
    mins = int(str(t)[0:2])*60 + int(str(t)[2:])
    if ampm == 'pm':
        mins += 720
    else:
        mins = mins
    
    return mins

t1 = timetomins(1100, 'am')
t2 = timetomins(1145, 'pm') 
dt = t2 - t1

dTdt = (T_2 - T_1)/dt

K = -(dTdt*(1/(T_2 - T_a)))
print('The value of K =',K)
```

2. Change your work from problem 1 to create a function that accepts the temperature at two times, ambient temperature, and the time elapsed to return $K$.

```{code-cell} ipython3
def newtoncoolingk(T1, T2, T_a, dt):
    '''Takes two temperarture inputs with two corresponding time values in absolute 
    minutes and returns a K value following Newton's Law of Cooling.
    
    Arguments
    ---------
    T1 = first temperature
    T2 = second temperature
    T_a = ambient temperature
    dt = time elapsed in minutes
    
    Returns
    ---------
    K = empirical constant in Newton's Law of Cooling equation '''
    
    K = ((T2 - T1)/(dt))*(1/(T2 - T_a))
    
    return K

print('The value of K =', newtoncoolingk(85, 80, 65, timetomins(1100, 'am') - timetomins(1145, 'pm')))
```

3. A first-order thermal system has the following analytical solution, 

    $T(t) =T_a+(T(0)-T_a)e^{-Kt}$

    where $T(0)$ is the temperature of the corpse at t=0 hours i.e. at the time of discovery and $T_a$ is a constant ambient temperature. 

    a. Show that an Euler integration converges to the analytical solution as the time step is decreased. Use the constant $K$ derived above and the initial temperature, T(0) = 85$^o$F. 

    b. What is the final temperature as t$\rightarrow\infty$?
    
    c. At what time was the corpse 98.6$^{o}$F? i.e. what was the time of death?

+++

## Part a
Creating a plot with both the analytical and numerical models along an axis of time (minutes) we can see that the two curves converge and idenctical with a timestep of 20 minutes as show below.

```{code-cell} ipython3
T1 = 85 # degrees Farenheit
T2 = 80 # degrees Farenheit
T_a = 65 # degrees Farenheit

t = np.arange(0, 1440, 20)
dt = t[1] - t[0]

N = len(t)

temp_analytical = T_a + (T1 - T_a)*np.exp(-K*t)


temp_numerical = np.zeros(N)
temp_numerical[0] = 85.0
for i in range(1,N):
    temp_numerical[i] = temp_numerical[i-1] - K*(temp_numerical[i-1]-T_a)*dt
    

plt.plot(t, temp_analytical,'-', label='Analytical');
plt.plot(t, temp_numerical,'o-', label='Numerical');
plt.legend(loc='best');
plt.xlabel('Time (min)');
plt.ylabel('Temperature ($^oF$)');
```

## Part b
Running the same code again, this time increasing the time range to a large value such as 100,000 minutes (333 hours), we can clearly see that temperature values approach a limit of 65 $^oF$ as time appraoches infinity.

```{code-cell} ipython3
T1 = 85 # degrees Farenheit
T2 = 80 # degrees Farenheit
T_a = 65 # degrees Farenheit

t = np.arange(0, 20000, 20)
dt = t[1] - t[0]

N = len(t)

temp_analytical = T_a + (T1 - T_a)*np.exp(-K*t)


temp_numerical = np.zeros(N)
temp_numerical[0] = 85.0
for i in range(1,N):
    temp_numerical[i] = temp_numerical[i-1] - K*(temp_numerical[i-1]-T_a)*dt
    

plt.plot(t, temp_analytical,'-', label='Analytical');
plt.plot(t, temp_numerical,'o-', label='Numerical');
plt.legend(loc='best');
plt.xlabel('Time (min)');
plt.ylabel('Temperature ($^oF$)');
```

## Part c
To find the time of death, first we must rearrange the analytical solution to solve for time, t.

$T(t) =T_a+(T(0)-T_a)e^{-Kt}$

$T(t)-T_a =(T(0)-T_a)e^{-Kt}$

$\frac {T(t)-T_a}{T(0)-T_a} =e^{-Kt}$

$\ln{\frac{T(t)-T_a}{T(0)-T_a}} =-Kt$

$t =-\frac{1}{K}\ln{\frac{T(t)-T_a}{T(0)-T_a}}$

```{code-cell} ipython3
T_d = 98.6 # Temp at time of death
T1 = 85 # Initial temperature
T_a = 65 # Ambient Temperature


t_an = -(1/K)*np.log((T_d - T_a)/(T1 - T_a))
t_d_min = timetomins(1100, 'am') + t_an # time of death minutes
t_d = 12 - abs(t_d_min/60)

time = str(t_d)
mins = int(time[2:4])*60/100


print('The death occured {:.2f} minutes before 11:00 am. In other words the death occurred at {:.0f}:0{:.0f} pm the previous day.'.format(abs(t_an), int(t_d), mins))
```

```{code-cell} ipython3

```
