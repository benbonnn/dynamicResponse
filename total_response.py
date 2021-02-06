# --------------------------------------------------------------
# PURPOSE : Compute total dynamic response of a damped system
#
# INPUT   : m = Mass [Kg]
#           c =  damping [N/(m/s)]
#           k =  stiffness [N/m]
#           t  -  time vector
#
# OUPUT   : u(t) = Total response of the system
#
# CALL    :  total_response(m,c,k,p0)
#--------------------------------------------------------------
# Last modified : Benjamin Bondsman     2021-02-06
#--------------------------------------------------------------
import math
import numpy as np
import matplotlib.pyplot as plt
# ---------------Total response of a damped system ------------
def total_response(m,c,k,p0):
    # NATURAL CIRCULAR FREQUENCY [rad/s]
    w_n = math.sqrt(k/m)
    # CRITICAL DAMPING COEFFICIENT
    Ccr = 2*m*w_n
    # RELATIVE CRITICAL VISCOUS DAMPING
    zeta = c/Ccr
    # DAMPED CIRCULAR FREQUENCY  [rad/s]
    w_d = w_n * math.sqrt(1-zeta**2)
    # NATURAL TIME PERIOD [s]
    T_n = (2*math.pi)
    # Def time vector
    t = np.arange (0,0.5*T_n,0.01)
    #
    # Solution of a system u = u_h + u_p
    # Solve Particular solution u_p
    u_p = p0/k
    # Solve homogeneous solution u_h
    A = u_p
    B = (w_n*zeta*A)/w_d
    # Initialize vectors
    response = []
    time = []
    # Homogeneous solution for all time steps
    for t in np.arange (0,0.5*T_n,0.01):
        u_tot = math.exp(-zeta*w_n*t)*((A*math.cos(w_n*t))+(B*math.sin(w_d*t)))+u_p
        response.append(u_tot)
        time.append(t)
    # Plot the results
    plt.figure()
    plt.plot(time,response,color='teal',marker=11)
    plt.xlabel('Time [-]')
    plt.ylabel('Displacement [mm]')
    plt.grid(True)
    plt.title('Total Response')
    plt.show()
#--------------------------------------------------------------
# ------------------ Call the function ------------------------
total_response(m,c,k,p0)