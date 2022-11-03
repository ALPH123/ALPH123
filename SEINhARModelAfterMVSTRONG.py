import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import odeint
#import plotly.graph_objects as go
#import plotly.io as pio
#import requests
#import kaleido
from lmfit import minimize, Parameters, Parameter, report_fit
#import os
#import sys
import plotly.graph_objects as go
import plotly.io as pio
#import requests
#from IPython.display import HTML
from ipywidgets.widgets import interact, IntSlider, FloatSlider, Layout, ToggleButton, ToggleButtons
plt.rcParams['font.sans-serif'] = ['Times New Roman']
plt.rcParams['axes.unicode_minus'] = False
style = {'description_width': '100px'}
slider_layout = Layout(width='99%')
#import panel as pn
pio.renderers.default = "notebook"
#pn.extension()
#plt.style.use('ggplot')
#plt.style.available()
#plt.style.use()

C0, C1, C2, C3= 0, 0.96, 0.96, 0.95,
def ode_model(y,t,beta, gamma, sigma, alpha, epsilon, lamb, eta, tau,omega,lota):
    S, E, I, Cd, A, R = y
    N=S+E+I+Cd+A+R
    dSdt0 = (eta*N-beta * S*(E+A)/ (N-Cd) - lamb*C0*S -eta*S-lota*S)*0.10
    dSdt1 = (-beta * S * (E + A) / (N-Cd) + lamb*C0 * S - lamb*C1 * S -eta*S-lota*S)*0.10
    dSdt2 = (-beta * S * (E + A) / (N-Cd) + lamb*C0 * S + lamb*C1 * S - lamb*C2 * S -eta*S-lota*S)*0.10
    dSdt3 = (-beta * S * (E + A) / (N-Cd) + lamb*C0 * S + lamb*C1 * S + lamb*C2 * S - lamb*C3 * S -eta*S-lota*S)*0.10
    dSdt = dSdt0+dSdt1+dSdt2+dSdt3
    dEdt = (beta * S*(E+A)/(N-Cd) - (sigma+eta)*E)*0.10#+(1-C1)*S
    dIdt = (sigma*alpha*E - (gamma+eta) * I)*0.10
    dCddt = (lota*S+tau * I - (omega + eta) * Cd - (gamma + eta) * Cd)*0.10
    dAdt= (sigma*epsilon*E - (gamma+eta) * A)#*0.10
    dRdt = (gamma * (I+A+Cd)-eta*R)*0.10
    return np.array([dSdt, dEdt, dIdt, dCddt,dAdt, dRdt])






def ode_solver(t, initial_conditions, params):
    initE, initI, initCd, initA, initR, initN = initial_conditions
    beta, gamma, sigma, alpha, epsilon,lamb,eta,tau,omega,lota = params['beta'].value,params['gamma'].value, params['sigma'].value, params['alpha'].value, params['epsilon'].value, params['lamb'].value, params['eta'].value, params['tau'].value,params['omega'].value,params['lota'].value
    initS = initN - (initE + initI + initR + initA+initCd)
    res = odeint(ode_model, [initS, initE, initI, initCd, initA, initR], t, args=(beta, gamma, sigma, alpha, epsilon,lamb,eta,tau,omega,lota))
    return res


initN = 2000000  # 1380000000
# S0 = 966000000
initE = 100000
initI = 1  # 47
initR = 1
initA = 0
initCd = 0
sigma = 0.18  # 1/5.2
gamma = 0.06 # 1/2.9
eta = 0.00
alpha = 0.08
tau = 0.008
epsilon = 0.06
lamb = 0.97
omega = 0.12
lota = 0.03
R0 = 4
beta = R0 * gamma
days = 200  # 180

params = Parameters()
params.add('beta', value=beta, min=0, max=10)
params.add('sigma', value=sigma, min=0, max=10)
params.add('gamma', value=gamma, min=0, max=10)
params.add('eta', value=eta, min=0, max=10)
params.add('alpha', value=alpha, min=0, max=10)
params.add('epsilon', value=epsilon, min=0, max=10)
params.add('lamb', value=lamb, min=0, max=10)
params.add('tau', value=tau, min=0, max=10)
params.add('omega', value=omega, min=0, max=10)
params.add('lota', value=lota, min=0, max=10)

# Simulation


def main(initE, initI, initCd, initA, initR, initN, beta, gamma, sigma, alpha, epsilon, lamb, eta, tau, days,
         param_fitting):
    initial_conditions = [initE, initI, initCd, initA, initR, initN]
    params['beta'].value, params['gamma'].value, params['sigma'].value, params['alpha'].value, params['epsilon'].value, \
    params['lamb'].value, params['eta'].value, params['tau'].value, params['omega'].value, params['lota'].value = [beta,
                                                                                                                   gamma,
                                                                                                                   sigma,
                                                                                                                   alpha,
                                                                                                                   epsilon,
                                                                                                                   lamb,
                                                                                                                   eta,
                                                                                                                   tau,
                                                                                                                   omega,                                                                                                              lota]
    tspan = np.arange(0, days, 1)
    sol = ode_solver(tspan, initial_conditions, params)
    S, E, I, Cd, A, R = sol[:, 0], sol[:, 1], sol[:, 2], sol[:, 3], sol[:, 4], sol[:, 5]

    # Create traces
    fig = go.Figure()
    if not param_fitting:
        fig.add_trace(
            go.Scatter(x=tspan, y=S, mode='lines+markers', name='Susceptible', line=dict(width=4), fill='tozeroy'))
        fig.add_trace(
            go.Scatter(x=tspan, y=E, mode='lines+markers', name='Exposed', line=dict(width=4), fill='tozeroy'))
        # fig.add_vline(x=2.5, line_width=3, line_dash="dash", line_color="green",row=np.max(E),col=0)
        # fig.add_hline(y=np.max(E), line_width=3, line_dash="dash", line_color="green")



params = Parameters()
params.add('beta', value=beta, min=0, max=10)
params.add('sigma', value=sigma, min=0, max=10)
params.add('gamma', value=gamma, min=0, max=10)
params.add('eta', value=eta, min=0, max=10)
params.add('alpha', value=alpha, min=0, max=10)
params.add('epsilon', value=epsilon, min=0, max=10)
params.add('lamb', value=lamb, min=0, max=10)
params.add('tau', value=tau, min=0, max=10)
params.add('omega', value=omega, min=0, max=10)
params.add('lota', value=lota, min=0, max=10)

# Simulation


def main(initE, initI, initCd, initA, initR, initN, beta, gamma, sigma, alpha, epsilon, lamb, eta, tau, days,
         param_fitting):
    initial_conditions = [initE, initI, initCd, initA, initR, initN]
    params['beta'].value, params['gamma'].value, params['sigma'].value, params['alpha'].value, params['epsilon'].value, \
    params['lamb'].value, params['eta'].value, params['tau'].value, params['omega'].value, params['lota'].value = [beta,
                                                                                                                   gamma,
                                                                                                                   sigma,
                                                                                                                   alpha,
                                                                                                                   epsilon,
                                                                                                                   lamb,
                                                                                                                   eta,
                                                                                                                   tau,
                                                                                                                   omega,
                                                                                                                   lota]
    tspan = np.arange(0, days, 1)
    sol = ode_solver(tspan, initial_conditions, params)
    S, E, I, Cd, A, R = sol[:, 0], sol[:, 1], sol[:, 2], sol[:, 3], sol[:, 4], sol[:, 5]

    # Create traces
    fig = go.Figure()
    if not param_fitting:
        fig.add_trace(
            go.Scatter(x=tspan, y=S, mode='lines+markers', name='Susceptible', line=dict(width=4), fill='tozeroy'))
        fig.add_trace(
            go.Scatter(x=tspan, y=E, mode='lines+markers', name='Exposed', line=dict(width=4), fill='tozeroy'))
        # fig.add_vline(x=2.5, line_width=3, line_dash="dash", line_color="green",row=np.max(E),col=0)
        # fig.add_hline(y=np.max(E), line_width=3, line_dash="dash", line_color="green")
        fig.add_trace(
            go.Scatter(x=tspan, y=I, mode='lines+markers', name='Infected', line=dict(width=4), fill='tozeroy'))
        fig.add_trace(go.Scatter(x=tspan, y=Cd, mode='lines+markers', name='Notified and hospitalized infection',
                                 line=dict(width=4), fill='tozeroy'))
        fig.add_trace(
            go.Scatter(x=tspan, y=A, mode='lines+markers', name='Asymptomatic', line=dict(width=4), fill='tozeroy'))
        fig.add_trace(
            go.Scatter(x=tspan, y=R, mode='lines+markers', name='Recovered', line=dict(width=4), fill='tozeroy'))
    # fig.show()
    if param_fitting:
        # fig.add_trace(go.Scatter(x=tspan, y=df_covid_history.infected, mode='lines+markers', \
        # name='Infections Observed', line=dict(dash='dash')))
        fig.add_trace(go.Scatter(x=tspan, y=S, mode='lines+markers', \
                                 name='Susceptible Observed', line=dict(dash='dash')))
        fig.add_trace(go.Scatter(x=tspan, y=E, mode='lines+markers', \
                                 name='Exposed Observed', line=dict(dash='dash')))
        fig.add_trace(go.Scatter(x=tspan, y=I, mode='lines+markers', \
                                 name='Infections Observed', line=dict(dash='dash')))
        fig.add_trace(go.Scatter(x=tspan, y=Cd, mode='lines+markers', \
                                 name='Notified and hospitalized infection Observed', line=dict(dash='dash')))
        fig.add_trace(go.Scatter(x=tspan, y=A, mode='lines+markers', \
                                 name='Asymptomatic Observed', line=dict(dash='dash')))
        fig.add_trace(go.Scatter(x=tspan, y=R, mode='lines+markers', \
                                 name='Recovered Observed', line=dict(dash='dash')))

    if days <= 30:
        step = 1
    elif days <= 90:
        step = 7
    else:
        step = 30

        # Edit the layout
    fig.update_layout(title='Simulation of SEINhAR Model',
                      xaxis_title='Day',
                      yaxis_title='Counts',
                      title_x=0.4, font=dict(size=30, family='Times New Roman'),
                      width=1300, height=600
                      )
    fig.update_xaxes(tickangle=0, tickformat=None, tickmode='array', tickvals=np.arange(0, days + 1, step))
    if not os.path.exists("F:/Project_bioinformatic_COVID19/Second_paper/NewFigure"):
        os.mkdir("F:/Project_bioinformatic_COVID19/Second_paper/NewFigure")
    pio.write_image(fig, "F:/Project_bioinformatic_COVID19/Second_paper/NewFigure/seird_simulation11.pdf")
    fig.show()




interact(main,
         initE=IntSlider(min=0, max=10000, step=1, value=initE, description='initE', style=style, layout=slider_layout),
         initI=IntSlider(min=0, max=10000, step=10, value=initI, description='initI', style=style, layout=slider_layout),
         initCd=IntSlider(min=0, max=10000, step=10, value=initCd, description='initCd', style=style,layout=slider_layout),
         initA=IntSlider(min=0, max=10000, step=10, value=initA, description='initA', style=style,layout=slider_layout),
         initR=IntSlider(min=0, max=10000, step=10, value=initR, description='initR', style=style, layout=slider_layout),

         initN=IntSlider(min=0, max=2000000, step=100000, value=initN, description='initN', style=style, layout=slider_layout),
         beta=FloatSlider(min=0, max=4, step=0.01, value=beta, description='Infection rate', style=style, layout=slider_layout),
         sigma=FloatSlider(min=0, max=4, step=0.01, value=sigma, description='Incubation rate', style=style, layout=slider_layout),
         gamma=FloatSlider(min=0, max=4, step=0.01, value=gamma, description='Recovery rate', style=style, layout=slider_layout),
         alpha=FloatSlider(min=0, max=1, step=0.001, value=alpha, description='Mortality rate', style=style, layout=slider_layout),
         epsilon=FloatSlider(min=0, max=1, step=0.001, value=epsilon, description='Mortality rate', style=style, layout=slider_layout),
         lamb=FloatSlider(min=0, max=1, step=0.001, value=lamb, description='Mortality rate', style=style, layout=slider_layout),
         eta=FloatSlider(min=0, max=1, step=0.001, value=eta, description='Mortality rate', style=style, layout=slider_layout),
         tau=FloatSlider(min=0, max=1, step=0.001, value=tau, description='Mortality rate', style=style, layout=slider_layout),
         omega=FloatSlider(min=0, max=1, step=0.001, value=omega, description='Mortality rate', style=style, layout=slider_layout),
         lota=FloatSlider(min=0, max=1, step=0.001, value=lota, description='Mortality rate', style=style, layout=slider_layout),

         days=IntSlider(min=0, max=600, step=7, value=days, description='Days', style=style, layout=slider_layout),
         param_fitting=ToggleButton(value=False, description='Fitting Mode', disabled=False, button_style='', \
             tooltip='Click to show fewer plots', icon='check-circle')
        )


#interactive_plot

initial_conditions = [initE, initI, initCd, initA, initR, initN]
"""sigma = 0.09
gamma = 0.09
eta = 0.00
alpha = 0.085
tau = 0.19
epsilon = 0.04
lamb = 0.097
R0 = 3
beta=R0*gamma
omega=0.04
lota=0.0173
R0 = 4"""
params['beta'].value = beta
params['sigma'].value = sigma
params['gamma'].value = gamma
params['eta'].value = eta
params['alpha'].value = alpha
params['tau'].value = tau
params['epsilon'].value = epsilon
params['lamb'].value = lamb
params['omega'].value = omega
params['lota'].value = lota
days = 200
tspan = np.arange(0, days, 1)
sol = ode_solver(tspan, initial_conditions, params)

t=tspan


ordered_time = [time_list for (score,time_list) in sorted(zip(sol[:,1],t))]
best_time = ordered_time[-1]
max_coords = '('+ np.str(best_time)+', ' + str("%.4f" % (np.max(sol[:,1])))+')'
#max_point = plt.plot(best_time, np.max(sol[:,1]), 'bo', label="(Opt. Time, Max Score)")
#plt.text(best_time, np.max(sol[:,1]), max_coords)


ordered_time1 = [time_list for (score1,time_list) in sorted(zip(sol[:,2],t))]
best_time1 = ordered_time1[-1]
max_coords1 = '('+ np.str(best_time1)+', ' + str("%.4f" % (np.max(sol[:,2])))+')'

ordered_time2 = [time_list for (score2,time_list) in sorted(zip(sol[:,3],t))]
best_time2 = ordered_time2[-1]
max_coords2 = '('+ np.str(best_time2)+', ' + str("%.4f" % (np.max(sol[:,3])))+')'


ordered_time3 = [time_list for (score3,time_list) in sorted(zip(sol[:,4],t))]
best_time3 = ordered_time3[-1]
max_coords3 = '('+ np.str(best_time3)+', ' + str("%.4f" % (np.max(sol[:,4])))+')'



fig = plt.figure()#facecolor='w')
ax = fig.add_subplot(111, axisbelow=True)
ax.fill_between(t, sol[:,0], step="pre", alpha=0.4)
ax.plot(t, sol[:,0], 'b', alpha=0.5, lw=3, label='Susceptible')
ax.fill_between(t, sol[:,1], step="pre", alpha=0.4)
ax.plot(t, sol[:,1], 'r', alpha=0.5, lw=3, label='Exposed')
ax.fill_between(t, sol[:,2], step="pre", alpha=0.4)
ax.plot(t, sol[:,2], 'y', alpha=0.5, lw=3, label='Infected')
ax.fill_between(t, sol[:,3], step="pre", alpha=0.4)
ax.plot(t,sol[:,3],'#00FF00', alpha=0.6,lw=3,label="Notified and hospitalized infection")
ax.fill_between(t, sol[:,4], step="pre", alpha=0.4)
ax.plot(t, sol[:,4], '#FF6103', alpha=0.5,lw=3, label='Asymptomatic')
ax.fill_between(t, sol[:,5], step="pre", alpha=0.4)
ax.plot(t, sol[:,5], 'g', alpha=0.5, lw=3, label='Recovered')
ax.set_xlabel('Days',fontsize=16,fontweight='bold')
ax.set_ylabel('Counts',fontsize=16,fontweight='bold')
ax.plot(best_time, np.max(sol[:,1]), 'bo')
ax.text(best_time, np.max(sol[:,1]), max_coords,ha="right")
plt.vlines(best_time,0,np.max(sol[:,1]),colors="r",linestyles="dashed")
plt.hlines(np.max(sol[:,1]),0,best_time,colors="r",linestyles="dashed")

ax.plot(best_time1, np.max(sol[:,2]), 'bo')
ax.text(best_time1, np.max(sol[:,2]), max_coords1,ha="left")
plt.vlines(best_time1,0,np.max(sol[:,2]),colors="y",linestyles="dashed")
plt.hlines(np.max(sol[:,2]),0,best_time1,colors="y",linestyles="dashed")

ax.plot(best_time2, np.max(sol[:,3]), 'bo')
ax.text(best_time2, np.max(sol[:,3]), max_coords2)
plt.vlines(best_time2,0,np.max(sol[:,3]),colors='#00FF00',linestyles="dashed")
plt.hlines(np.max(sol[:,3]),0,best_time2,colors='#00FF00',linestyles="dashed")

ax.plot(best_time3, np.max(sol[:,4]), 'bo')
ax.text(best_time3, np.max(sol[:,4]), max_coords3)
plt.vlines(best_time3,0,np.max(sol[:,4]),colors='#FF6103',linestyles="dashed")
plt.hlines(np.max(sol[:,4]),0,best_time3,colors='#FF6103',linestyles="dashed")


ax.yaxis.set_tick_params(length=2,labelsize=12)
ax.xaxis.set_tick_params(length=2,labelsize=12)

legend = ax.legend(fontsize=2,prop={'weight':'bold'})
legend.get_frame().set_linewidth(0.0)
for spine in ('top', 'right','bottom', 'left'):
    ax.spines[spine].set_visible(True)

plt.tick_params(axis="both",labelsize=15)
plt.title('Simulation of SEINhAR Model',fontweight='bold',fontsize=12)
plt.show()
